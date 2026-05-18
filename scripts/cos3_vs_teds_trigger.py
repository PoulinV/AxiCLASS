#!/usr/bin/env python3
"""
Study: can (1-cos)^3 work as a DM-triggered EDE potential?

Central question: if we accept tuning theta_i close to pi (hilltop), does
the coupling term (not the bare slope) control the EDE release for the
(1-cos)^3 potential?  We compare against the tEDS benchmark.

Produces four figures saved to analysis/cos3_vs_teds/:
  fig1_fede_comparison.pdf   -- f_EDE(z) and theta(z) side-by-side
  fig2_force_budget.pdf      -- coupling vs bare slope vs z, varying theta_i
  fig3_phase_diagram.pdf     -- (delta, g) phase diagram: coupling fraction at peak
  fig4_coupling_fraction.pdf -- coupling fraction vs z for best cases
"""

from __future__ import annotations

import math
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import get_cmap

REPO = Path(__file__).resolve().parents[1]
CLASS = REPO / "class"
OUT = REPO / "analysis" / "cos3_vs_teds"

# Cosmological parameters (matching existing scripts)
H0 = 72.81
OMEGA_B = 0.02251
OMEGA_DM_TODAY = 0.1280   # present-day IDM-EDE (h^2); CLASS shoots omega_ini internally
Z_EQ = 3400.0

# (1-cos)^3 parameters
N_AXION = 3
F_AXION = 1.0     # Mpl
M_AXION = 20000.0  # in H0 units

# tEDS benchmark parameters
F_TEDS = 0.05
P_TEDS = 8.0
EV4_TO_CLASS = (1.56e29 / 2.435e27) ** 2  # 1 eV^4 in CLASS units; no /3 since background.c divides rho_scf=(K+V)/3 internally
V0_TEDS = 0.12 * EV4_TO_CLASS  # Lin+ 2212.08098 Fig 1 caption states V0 = 0.12 eV^4

# Common ini footer
COMMON_FOOTER = """N_ur = 2.0328
N_ncdm = 1
deg_ncdm = 1
m_ncdm = 0.06
T_ncdm = 0.71611
scf_evolve_as_fluid = no
scf_evolve_like_axionCAMB = no
do_shooting = yes
do_shooting_scf = no
scf_has_perturbations = no
attractor_ic_scf = no
modes = s
gauge = synchronous
lensing = no
"""

# No omega_cdm (all DM is IDM-EDE, matching Lin+ 2212.08098).
# No Omega_Lambda: CLASS derives it from the flat-universe condition.
# omega_idm_ede is the TODAY value; CLASS internal shooting finds omega_ini_idm_ede.
COMMON_HEADER = """output = tCl
background_evolver = 0
evolver = 0
write background = yes
write parameters = no
omega_b = {omega_b}
omega_cdm = 0
omega_idm_ede = {omega_dm}
H0 = {H0}
tau_reio = 0.068
A_s = 2.191e-9
n_s = 0.9860
coupling_type = 1
idm_ede_mass_form = poly
g_idm_ede = {g}
""".format


def axion_ini(root: Path, theta_i: float, g: float, omega_dm: float = OMEGA_DM_TODAY) -> str:
    header = COMMON_HEADER(omega_b=OMEGA_B, omega_dm=omega_dm, H0=H0, g=g)
    return (
        f"root = {root}\n"
        + header
        + f"scf_potential = axion\n"
        + f"n_axion = {N_AXION}\n"
        + f"scf_parameters = {theta_i},0.0\n"
        + f"f_axion = {F_AXION}\n"
        + f"m_axion = {M_AXION}\n"
        + COMMON_FOOTER
    )


def teds_ini(root: Path, theta_i: float, g: float, omega_dm: float = OMEGA_DM_TODAY) -> str:
    phi_i = theta_i * F_TEDS
    header = COMMON_HEADER(omega_b=OMEGA_B, omega_dm=omega_dm, H0=H0, g=g)
    return (
        f"root = {root}\n"
        + header
        + f"scf_potential = teds\n"
        + f"scf_parameters = {phi_i},0.0\n"
        + f"V0_teds = {V0_TEDS}\n"
        + f"f_teds = {F_TEDS}\n"
        + f"p_teds = {P_TEDS}\n"
        + COMMON_FOOTER
    )


def load_bg(path: Path) -> Optional[dict[str, np.ndarray]]:
    try:
        header = ""
        with path.open() as fh:
            for line in fh:
                if line.startswith("#") and "1:z" in line:
                    header = line[1:]
                    break
        if not header:
            return None
        names = [e.split(":", 1)[1] for e in header.split() if ":" in e]
        data = np.loadtxt(path, comments="#")
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return {n: data[:, i] for i, n in enumerate(names)}
    except Exception as exc:
        print(f"  load error: {exc}", file=sys.stderr)
        return None


def coupling_force_arr(bg: dict, g: float) -> np.ndarray:
    """3 rho_IDM d(ln m)/d(phi) — same units as V'_scf."""
    phi = bg["phi_scf"]
    rho = bg["(.)rho_idm_ede"]
    dlnm = 2.0 * g * phi / (1.0 + g * phi**2)
    return 3.0 * rho * dlnm


def _run_class_once(ini_text: str, root: Path, timeout: int = 90) -> bool:
    ini_path = Path(str(root) + ".ini")
    ini_path.write_text(ini_text)
    for suf in ("00_background.dat", "00_cl.dat", "00_parameters.ini", "00_unused_parameters"):
        p = Path(str(root) + suf)
        if p.exists():
            p.unlink()
    res = subprocess.run(
        [str(CLASS), str(ini_path)],
        cwd=REPO, capture_output=True, text=True, timeout=timeout,
    )
    if res.returncode != 0:
        snippet = (res.stderr or res.stdout or "")[-400:].strip()
        print(f"  CLASS failed: {snippet}", file=sys.stderr)
    return res.returncode == 0


def run_class_normalized(
    ini_builder,
    root: Path,
    g: float = 0.0,
    phi_i: float = 0.0,
    omega_dm: float = OMEGA_DM_TODAY,
    timeout: int = 120,
) -> bool:
    """Run CLASS with internal shooting for omega_idm_ede (omega_dm is the today target)."""
    return _run_class_once(ini_builder(root), root, timeout)


@dataclass
class Result:
    label: str
    potential: str
    theta_i: float
    g: float
    ok: bool
    log10_z_peak: float = 0.0
    f_ede_peak: float = 0.0
    coupling_frac_peak: float = 0.0
    bg: dict = field(default_factory=dict, repr=False)


def analyse(bg: dict, g: float) -> dict:
    z = bg["z"]
    f_ede = bg["(.)Omega_scf"]
    v_prime = bg["V'_scf"]
    coup = coupling_force_arr(bg, g)

    peak = int(np.argmax(f_ede))
    bare_abs = np.abs(v_prime)
    coup_abs = np.abs(coup)
    total = bare_abs + coup_abs + 1e-300
    frac = coup_abs / total

    return {
        "log10_z_peak": math.log10(float(z[peak])) if z[peak] > 0 else 0.0,
        "f_ede_peak": float(f_ede[peak]),
        "coupling_frac_peak": float(frac[peak]),
        "bare_abs": bare_abs,
        "coup_abs": coup_abs,
        "coup_frac": frac,
    }


# ---------------------------------------------------------------------------
# Scan cases
# ---------------------------------------------------------------------------

# (1-cos)^3: theta_i = pi - delta for various delta and g values
DELTA_VALUES = [2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001]
G_VALUES = [0.3, 0.56, 0.68, 1.0, 2.0]

# tEDS benchmarks
TEDS_CASES = [
    (8.0, 0.75),
    (6.0, 0.68),
    (4.0, 0.56),
]

# Representative (1-cos)^3 cases for the comparison figures (chosen after scan)
# We pick a few coupling-dominated cases that bracket g and delta values
COS3_SHOWCASE = [
    (math.pi - 0.05,  0.68),
    (math.pi - 0.1,   0.68),
    (math.pi - 0.5,   0.68),
    (math.pi - 0.05,  1.0),
    (math.pi - 0.5,   1.0),
]


def run_axion_scan() -> list[Result]:
    results = []
    total = len(DELTA_VALUES) * len(G_VALUES)
    idx = 0
    for delta in DELTA_VALUES:
        theta_i = math.pi - delta
        for g in G_VALUES:
            idx += 1
            tag = f"axion_d{delta:.4g}_g{g:.4g}".replace(".", "p")
            root = OUT / tag
            phi_i = theta_i * F_AXION
            print(f"[{idx}/{total}] (1-cos)^3 delta={delta:.4g} g={g:.4g} ...", end=" ", flush=True)
            ok = run_class_normalized(
                lambda r, _ti=theta_i, _g=g: axion_ini(r, _ti, _g),
                root,
            )
            r = Result(
                label=rf"$\delta={delta:.3g},\,g={g:.2f}$",
                potential="axion",
                theta_i=theta_i, g=g, ok=ok,
            )
            if ok:
                bg = load_bg(Path(str(root) + "00_background.dat"))
                if bg is not None:
                    info = analyse(bg, g)
                    r.log10_z_peak = info["log10_z_peak"]
                    r.f_ede_peak = info["f_ede_peak"]
                    r.coupling_frac_peak = info["coupling_frac_peak"]
                    r.bg = bg
                    print(f"ok  f_EDE={r.f_ede_peak:.3f}  z_peak=10^{r.log10_z_peak:.2f}  "
                          f"coup_frac={r.coupling_frac_peak:.2f}")
                else:
                    r.ok = False
                    print("load failed")
            else:
                print("failed")
            results.append(r)
    return results


def run_teds_benchmarks() -> list[Result]:
    results = []
    for theta_i, g in TEDS_CASES:
        tag = f"teds_th{theta_i:.4g}_g{g:.4g}".replace(".", "p")
        root = OUT / tag
        phi_i = theta_i * F_TEDS
        print(f"tEDS theta_i={theta_i:.4g} g={g:.4g} ...", end=" ", flush=True)
        ok = run_class_normalized(
            lambda r, _ti=theta_i, _g=g: teds_ini(r, _ti, _g),
            root,
        )
        r = Result(
            label=rf"tEDS $\theta_i={theta_i:.4g},\,g={g:.2f}$",
            potential="teds", theta_i=theta_i, g=g, ok=ok,
        )
        if ok:
            bg = load_bg(Path(str(root) + "00_background.dat"))
            if bg is not None:
                info = analyse(bg, g)
                r.log10_z_peak = info["log10_z_peak"]
                r.f_ede_peak = info["f_ede_peak"]
                r.coupling_frac_peak = info["coupling_frac_peak"]
                r.bg = bg
                print(f"ok  f_EDE={r.f_ede_peak:.3f}  z_peak=10^{r.log10_z_peak:.2f}  "
                      f"coup_frac={r.coupling_frac_peak:.2f}")
            else:
                r.ok = False
                print("load failed")
        else:
            print("failed")
        results.append(r)
    return results


def run_cos3_showcase() -> list[Result]:
    """A smaller set of (1-cos)^3 cases for the comparison figure."""
    results = []
    for theta_i, g in COS3_SHOWCASE:
        delta = math.pi - theta_i
        phi_i = theta_i * F_AXION
        tag = f"axion_show_d{delta:.4g}_g{g:.4g}".replace(".", "p")
        root = OUT / tag
        print(f"showcase (1-cos)^3 delta={delta:.4g} g={g:.4g} ...", end=" ", flush=True)
        ok = run_class_normalized(
            lambda r, _ti=theta_i, _g=g: axion_ini(r, _ti, _g),
            root,
        )
        r = Result(
            label=rf"$(1\!-\!\cos)^3$, $\delta\pi={delta:.3g}$, $g={g:.2f}$",
            potential="axion", theta_i=theta_i, g=g, ok=ok,
        )
        if ok:
            bg = load_bg(Path(str(root) + "00_background.dat"))
            if bg is not None:
                info = analyse(bg, g)
                r.log10_z_peak = info["log10_z_peak"]
                r.f_ede_peak = info["f_ede_peak"]
                r.coupling_frac_peak = info["coupling_frac_peak"]
                r.bg = bg
                print(f"ok  f_EDE={r.f_ede_peak:.3f}  z_peak=10^{r.log10_z_peak:.2f}  "
                      f"coup_frac={r.coupling_frac_peak:.2f}")
            else:
                r.ok = False
                print("load failed")
        else:
            print("failed")
        results.append(r)
    return results


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

STYLE = dict(lw=1.8, solid_capstyle="round")


def fig1_fede_comparison(cos3_show: list[Result], teds: list[Result]) -> None:
    """
    Side-by-side: f_EDE(z) and theta(z) for tEDS (left) and (1-cos)^3 (right).
    """
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 5.5), sharex=True)
    fig.subplots_adjust(hspace=0.08, wspace=0.28)

    colors_teds = plt.cm.Blues(np.linspace(0.45, 0.85, len(teds)))
    colors_cos3 = plt.cm.Oranges(np.linspace(0.45, 0.85, len([r for r in cos3_show if r.ok])))

    ax_fede_t, ax_fede_c = axes[0]
    ax_th_t, ax_th_c = axes[1]

    ok_cos3 = [r for r in cos3_show if r.ok]

    for color, r in zip(colors_teds, teds):
        if not r.ok:
            continue
        z = r.bg["z"]
        log10z = np.log10(np.maximum(z, 1e-3))
        ax_fede_t.plot(log10z, r.bg["(.)Omega_scf"], color=color, label=r.label, **STYLE)
        ax_th_t.plot(log10z, r.bg["phi_scf"] / F_TEDS, color=color, **STYLE)

    for color, r in zip(colors_cos3, ok_cos3):
        z = r.bg["z"]
        log10z = np.log10(np.maximum(z, 1e-3))
        ax_fede_c.plot(log10z, r.bg["(.)Omega_scf"], color=color, label=r.label, **STYLE)
        ax_th_c.plot(log10z, r.bg["phi_scf"] / F_AXION, color=color, **STYLE)

    for ax in axes.flat:
        ax.axvline(math.log10(Z_EQ), color="red", ls="--", lw=0.9, alpha=0.7, label=r"$z_{\rm eq}$")
        ax.set_xlim(1.5, 6.5)
        ax.tick_params(direction="in", top=True, right=True)

    ax_fede_t.set_ylabel(r"$f_{\rm EDE}(z)$")
    ax_fede_c.set_ylabel(r"$f_{\rm EDE}(z)$")
    ax_th_t.set_ylabel(r"$\theta(z) = \phi/f$")
    ax_th_c.set_ylabel(r"$\theta(z) = \phi/f$")
    ax_th_t.set_xlabel(r"$\log_{10} z$")
    ax_th_c.set_xlabel(r"$\log_{10} z$")

    ax_fede_t.set_title("tEDS", fontsize=10)
    ax_fede_c.set_title(r"$(1-\cos)^3$", fontsize=10)

    for ax in [ax_fede_t, ax_fede_c]:
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")

    fig.savefig(OUT / "fig1_fede_comparison.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig1_fede_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved fig1_fede_comparison")


def fig2_force_budget(cos3_show: list[Result]) -> None:
    """
    For fixed g=0.68, vary delta = pi - theta_i.
    Show |V'| and |coupling| vs log10(z) with coupling fraction.
    Three panels stacked.
    """
    target_g = 0.68
    cases = [r for r in cos3_show if r.ok and abs(r.g - target_g) < 0.01]
    if not cases:
        # fall back to any g value
        cases = [r for r in cos3_show if r.ok]
    if not cases:
        print("fig2: no ok cases, skipping")
        return

    n = len(cases)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, n))

    fig, axes = plt.subplots(3, 1, figsize=(5.5, 7.0), sharex=True)
    fig.subplots_adjust(hspace=0.08)

    for color, r in zip(colors, sorted(cases, key=lambda x: x.theta_i)):
        z = r.bg["z"]
        log10z = np.log10(np.maximum(z, 1e-3))
        bare = np.abs(r.bg["V'_scf"])
        coup = np.abs(coupling_force_arr(r.bg, r.g))
        frac = coup / (bare + coup + 1e-300)
        delta = math.pi - r.theta_i
        label = rf"$\delta\pi = {delta:.3g}$"

        axes[0].semilogy(log10z, bare, color=color, ls="-", label=label, **STYLE)
        axes[1].semilogy(log10z, coup, color=color, ls="-", **STYLE)
        axes[2].plot(log10z, frac, color=color, ls="-", **STYLE)

    for ax in axes:
        ax.axvline(math.log10(Z_EQ), color="red", ls="--", lw=0.9, alpha=0.7)
        ax.set_xlim(1.5, 6.5)
        ax.tick_params(direction="in", top=True, right=True)

    axes[2].axhline(0.5, color="black", ls=":", lw=1.0)
    axes[0].set_ylabel(r"$|V_{,\phi}|$ [CLASS units]")
    axes[1].set_ylabel(r"$|3\rho_{\rm IDM}\partial_\phi \ln m|$")
    axes[2].set_ylabel(r"$f_{\rm coup}$")
    axes[2].set_xlabel(r"$\log_{10} z$")
    axes[0].legend(frameon=False, fontsize=8, title=rf"$(1-\cos)^3$, $g={target_g:.2f}$",
                   title_fontsize=8)
    axes[2].set_ylim(-0.05, 1.05)

    fig.savefig(OUT / "fig2_force_budget.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig2_force_budget.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved fig2_force_budget")


def fig3_phase_diagram(scan: list[Result]) -> None:
    """
    Phase diagram in (log10 delta, g) space, colored by coupling fraction at peak.
    """
    ok = [r for r in scan if r.ok and r.potential == "axion"]
    if not ok:
        print("fig3: no ok axion scan results, skipping")
        return

    deltas = np.array([math.pi - r.theta_i for r in ok])
    gs = np.array([r.g for r in ok])
    fracs = np.array([r.coupling_frac_peak for r in ok])
    zpeaks = np.array([r.log10_z_peak for r in ok])

    log_deltas = np.log10(deltas)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
    fig.subplots_adjust(wspace=0.32)

    # Left: coupling fraction
    sc = axes[0].scatter(log_deltas, gs, c=fracs, cmap="RdYlGn",
                         vmin=0, vmax=1, s=90, edgecolors="gray", lw=0.4, zorder=3)
    cb = fig.colorbar(sc, ax=axes[0], label=r"$f_{\rm coup}$ at $f_{\rm EDE}$ peak")
    axes[0].axhline(0.68, color="blue", ls=":", lw=1, alpha=0.6, label=r"$g=0.68$ (tEDS benchmark)")
    axes[0].set_xlabel(r"$\log_{10}(\pi - \theta_i)$")
    axes[0].set_ylabel(r"$g$")
    axes[0].set_title(r"Coupling fraction at release: $(1-\cos)^3$")
    axes[0].legend(frameon=False, fontsize=8)

    # Mark coupling-dominated region (frac > 0.5)
    dominated = fracs > 0.5
    if dominated.any():
        axes[0].scatter(log_deltas[dominated], gs[dominated],
                        marker="*", s=120, color="green", zorder=4, label=r"$f_{\rm coup}>0.5$")
        axes[0].legend(frameon=False, fontsize=8)

    # Right: z_peak
    sc2 = axes[1].scatter(log_deltas, gs, c=zpeaks, cmap="plasma",
                          vmin=3.0, vmax=4.5, s=90, edgecolors="gray", lw=0.4, zorder=3)
    cb2 = fig.colorbar(sc2, ax=axes[1], label=r"$\log_{10} z_{\rm peak}$")
    axes[1].axhline(math.log10(Z_EQ), color="red", ls="--", lw=1.0, alpha=0.7)
    axes[1].axhline(0.68, color="blue", ls=":", lw=1, alpha=0.6)
    axes[1].set_xlabel(r"$\log_{10}(\pi - \theta_i)$")
    axes[1].set_ylabel(r"$g$")
    axes[1].set_title(r"EDE release redshift: $(1-\cos)^3$")

    # Draw z_eq contour line
    log10_zeq = math.log10(Z_EQ)
    axes[1].contour(
        np.unique(log_deltas), np.unique(gs),
        zpeaks.reshape(len(DELTA_VALUES), len(G_VALUES)).T,
        levels=[log10_zeq], colors=["red"], linewidths=[1.5],
    )

    for ax in axes:
        ax.tick_params(direction="in", top=True, right=True)

    fig.savefig(OUT / "fig3_phase_diagram.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig3_phase_diagram.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved fig3_phase_diagram")


def fig4_coupling_fraction_vs_z(cos3_show: list[Result], teds: list[Result]) -> None:
    """
    Coupling fraction f_coup(z) for both potentials.
    One panel for tEDS, one for (1-cos)^3 best cases.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0), sharey=True)
    fig.subplots_adjust(wspace=0.12)

    colors_t = plt.cm.Blues(np.linspace(0.45, 0.85, len(teds)))
    ok_cos3 = [r for r in cos3_show if r.ok]
    colors_c = plt.cm.Oranges(np.linspace(0.45, 0.85, len(ok_cos3)))

    for color, r in zip(colors_t, teds):
        if not r.ok:
            continue
        z = r.bg["z"]
        log10z = np.log10(np.maximum(z, 1e-3))
        bare = np.abs(r.bg["V'_scf"])
        coup = np.abs(coupling_force_arr(r.bg, r.g))
        frac = coup / (bare + coup + 1e-300)
        axes[0].plot(log10z, frac, color=color, label=r.label, **STYLE)

    for color, r in zip(colors_c, ok_cos3):
        z = r.bg["z"]
        log10z = np.log10(np.maximum(z, 1e-3))
        bare = np.abs(r.bg["V'_scf"])
        coup = np.abs(coupling_force_arr(r.bg, r.g))
        frac = coup / (bare + coup + 1e-300)
        axes[1].plot(log10z, frac, color=color, label=r.label, **STYLE)

    for ax in axes:
        ax.axvline(math.log10(Z_EQ), color="red", ls="--", lw=0.9, alpha=0.7,
                   label=r"$z_{\rm eq}$")
        ax.axhline(0.5, color="black", ls=":", lw=1.0)
        ax.set_xlim(1.5, 6.5)
        ax.set_ylim(-0.05, 1.05)
        ax.tick_params(direction="in", top=True, right=True)
        ax.set_xlabel(r"$\log_{10} z$")
        ax.legend(frameon=False, fontsize=7.5)

    axes[0].set_ylabel(r"$f_{\rm coup}(z) = \frac{|\mathrm{coupling}|}{|\mathrm{bare}|+|\mathrm{coupling}|}$")
    axes[0].set_title("tEDS")
    axes[1].set_title(r"$(1-\cos)^3$")

    fig.savefig(OUT / "fig4_coupling_fraction.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig4_coupling_fraction.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved fig4_coupling_fraction")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Running tEDS benchmarks")
    print("=" * 60)
    teds = run_teds_benchmarks()

    print()
    print("=" * 60)
    print("Running (1-cos)^3 showcase cases")
    print("=" * 60)
    cos3_show = run_cos3_showcase()

    print()
    print("=" * 60)
    print("Running (1-cos)^3 scan (delta x g grid)")
    print("=" * 60)
    scan = run_axion_scan()

    print()
    print("=" * 60)
    print("Making figures")
    print("=" * 60)
    fig1_fede_comparison(cos3_show, teds)
    fig2_force_budget(cos3_show)
    fig3_phase_diagram(scan)
    fig4_coupling_fraction_vs_z(cos3_show, teds)

    print()
    print("Summary of scan results:")
    print(f"{'delta':>10}  {'g':>5}  {'ok':>4}  {'f_EDE':>7}  {'log10z':>7}  {'f_coup':>7}")
    for r in scan:
        delta = math.pi - r.theta_i
        ok_str = "ok" if r.ok else "FAIL"
        if r.ok:
            print(f"{delta:10.4g}  {r.g:5.2f}  {ok_str:>4}  "
                  f"{r.f_ede_peak:7.4f}  {r.log10_z_peak:7.3f}  {r.coupling_frac_peak:7.3f}")
        else:
            print(f"{delta:10.4g}  {r.g:5.2f}  {ok_str:>4}")

    print()
    print(f"Output directory: {OUT}")


if __name__ == "__main__":
    main()
