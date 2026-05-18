#!/usr/bin/env python3
"""
Make final figures for cos3 vs tEDS study from pre-computed CLASS data.

This script reads the background .dat files from analysis/cos3_vs_teds/
and produces four publication-quality figures.

Key physical insight being illustrated:
  - tEDS: coupling triggers the ONSET of rolling (coup_frac high at trigger z),
    then the field rolls off the plateau to phi~0 where the coupling dies.
    coup_frac at peak is ~0 (bare slope dominates at peak itself).
  - (1-cos)^3: phi stays near pi*f throughout the release, so the coupling
    remains active at and after the peak. But whether the RELEASE is actually
    driven by the coupling (vs the bare slope) is tested by comparing z_peak
    at fixed m with and without coupling.

The physically correct diagnostic for "coupling-triggered release" is:
  Does removing the coupling (g->0) significantly shift z_peak?
We approximate this by showing z_peak vs delta for each g value, and identifying
the region where z_peak ~ z_eq (the coincidence-problem condition).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "analysis" / "cos3_vs_teds"
OUT  = REPO / "analysis" / "cos3_vs_teds"
FIGS = REPO / "report" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

Z_EQ = 3400.0
F_AXION = 1.0
F_TEDS = 0.05

DELTA_VALUES = [2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001]
G_VALUES = [0.3, 0.56, 0.68, 1.0, 2.0]

TEDS_CASES = [(8.0, 0.75), (6.0, 0.68), (4.0, 0.56)]

# Colours consistent across figures
G_COLORS = {0.3: "#e6194B", 0.56: "#f58231", 0.68: "#3cb44b", 1.0: "#4363d8", 2.0: "#911eb4"}
TEDS_COLORS = ["#0077bb", "#33bbee", "#009988"]


def load_bg(path: Path) -> dict[str, np.ndarray] | None:
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
    except Exception:
        return None


def coupling_force(bg: dict, g: float) -> np.ndarray:
    phi = bg["phi_scf"]
    rho = bg["(.)rho_idm_ede"]
    dlnm = 2.0 * g * phi / (1.0 + g * phi**2)
    return 3.0 * rho * dlnm


def trigger_idx(bg: dict, frac: float = 0.05) -> int:
    """Index where field first moves by frac * |phi_i| from initial value."""
    phi = bg["phi_scf"]
    phi0 = phi[0]
    threshold = frac * abs(phi0) if abs(phi0) > 1e-10 else frac
    candidates = np.where(np.abs(phi - phi0) > threshold)[0]
    return int(candidates[0]) if candidates.size else int(np.argmax(bg["(.)Omega_scf"]))


def build_scan_table() -> list[dict]:
    rows = []
    for delta in DELTA_VALUES:
        theta_i = math.pi - delta
        tag_base = f"axion_d{delta:.4g}".replace(".", "p")
        for g in G_VALUES:
            tag = f"{tag_base}_g{g:.4g}".replace(".", "p")
            bg_path = DATA / f"{tag}00_background.dat"
            if not bg_path.exists():
                # try alternate naming
                bg_path = list(DATA.glob(f"{tag}*background.dat"))
                bg_path = bg_path[0] if bg_path else None
            if bg_path is None or not bg_path.exists():
                continue
            bg = load_bg(bg_path)
            if bg is None:
                continue
            z = bg["z"]
            f_ede = bg["(.)Omega_scf"]
            peak = int(np.argmax(f_ede))
            trig = trigger_idx(bg)
            bare = np.abs(bg["V'_scf"])
            coup = np.abs(coupling_force(bg, g))
            frac_arr = coup / (bare + coup + 1e-300)
            rows.append({
                "delta": delta,
                "theta_i": theta_i,
                "g": g,
                "z_peak": float(z[peak]),
                "f_ede_peak": float(f_ede[peak]),
                "coup_frac_peak": float(frac_arr[peak]),
                "coup_frac_trigger": float(frac_arr[trig]),
                "z_trigger": float(z[trig]),
                "bg": bg,
            })
    return rows


def build_teds_table() -> list[dict]:
    rows = []
    for theta_i, g in TEDS_CASES:
        tag = f"teds_th{theta_i:.4g}_g{g:.4g}".replace(".", "p")
        bg_path = DATA / f"{tag}00_background.dat"
        if not bg_path.exists():
            bg_path = list(DATA.glob(f"{tag}*background.dat"))
            bg_path = bg_path[0] if bg_path else None
        if bg_path is None or not bg_path.exists():
            continue
        bg = load_bg(bg_path)
        if bg is None:
            continue
        z = bg["z"]
        f_ede = bg["(.)Omega_scf"]
        peak = int(np.argmax(f_ede))
        trig = trigger_idx(bg, frac=0.05)
        bare = np.abs(bg["V'_scf"])
        phi = bg["phi_scf"]
        rho = bg["(.)rho_idm_ede"]
        dlnm = 2.0 * g * phi / (1.0 + g * phi**2)
        coup = np.abs(3.0 * rho * dlnm)
        frac_arr = coup / (bare + coup + 1e-300)
        rows.append({
            "theta_i": theta_i,
            "g": g,
            "z_peak": float(z[peak]),
            "f_ede_peak": float(f_ede[peak]),
            "coup_frac_peak": float(frac_arr[peak]),
            "coup_frac_trigger": float(frac_arr[trig]),
            "z_trigger": float(z[trig]),
            "bg": bg,
        })
    return rows


STYLE = dict(lw=1.9, solid_capstyle="round")


# ---------------------------------------------------------------------------
# Fig 1: z_peak vs delta for each g — the main "does it work?" figure
# ---------------------------------------------------------------------------

def fig_zpeak_vs_delta(scan: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    fig.subplots_adjust(wspace=0.28)

    ax_z, ax_f = axes

    for g in G_VALUES:
        rows = sorted([r for r in scan if r["g"] == g], key=lambda r: r["delta"])
        deltas = [r["delta"] for r in rows]
        zpeaks = [math.log10(r["z_peak"]) if r["z_peak"] > 0 else 0 for r in rows]
        fpeaks = [r["f_ede_peak"] for r in rows]
        color = G_COLORS[g]
        ax_z.semilogx(deltas, zpeaks, color=color, marker="o", ms=5, label=rf"$g={g:.2f}$", **STYLE)
        ax_f.semilogx(deltas, fpeaks, color=color, marker="o", ms=5, label=rf"$g={g:.2f}$", **STYLE)

    ax_z.axhline(math.log10(Z_EQ), color="red", ls="--", lw=1.2, label=r"$z_{\rm eq}$")
    ax_f.axhline(0.10, color="black", ls=":", lw=1.2, label=r"$f_{\rm EDE}=0.1$")

    ax_z.set_xlabel(r"$\delta \equiv \pi - \theta_i$")
    ax_z.set_ylabel(r"$\log_{10} z_{\rm peak}$")
    ax_z.set_title(r"Release redshift vs initial field offset")
    ax_z.legend(frameon=False, fontsize=8.5)
    ax_z.tick_params(direction="in", top=True, right=True)

    ax_f.set_xlabel(r"$\delta \equiv \pi - \theta_i$")
    ax_f.set_ylabel(r"$f_{\rm EDE,peak}$")
    ax_f.set_title(r"Peak EDE fraction vs initial field offset")
    ax_f.legend(frameon=False, fontsize=8.5)
    ax_f.tick_params(direction="in", top=True, right=True)

    for f in ("fig1_zpeak_vs_delta.pdf", "fig1_zpeak_vs_delta.png"):
        fig.savefig(OUT / f, bbox_inches="tight", dpi=180)
        fig.savefig(FIGS / f, bbox_inches="tight", dpi=180)
    plt.close(fig)
    print("saved fig1_zpeak_vs_delta")


# ---------------------------------------------------------------------------
# Fig 2: f_EDE(z) histories — tEDS and best (1-cos)^3 side by side
# ---------------------------------------------------------------------------

def fig_histories(scan: list[dict], teds: list[dict]) -> None:
    # Best (1-cos)^3 cases: g=0.68, delta in {0.5, 0.2, 0.05, 0.01}
    cos3_cases = [r for r in scan if abs(r["g"] - 0.68) < 0.01 and r["delta"] in {0.5, 0.2, 0.05, 0.01}]
    cos3_cases.sort(key=lambda r: r["delta"], reverse=True)

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.0), sharex=True)
    fig.subplots_adjust(hspace=0.08, wspace=0.28)

    ax_ft, ax_fc = axes[0]
    ax_tt, ax_tc = axes[1]

    teds_colors = plt.cm.Blues(np.linspace(0.45, 0.88, len(teds)))
    cos3_colors = plt.cm.Oranges(np.linspace(0.40, 0.90, len(cos3_cases)))

    for color, r in zip(teds_colors, teds):
        bg = r["bg"]
        z = bg["z"]
        lz = np.log10(np.maximum(z, 1e-3))
        ax_ft.plot(lz, bg["(.)Omega_scf"], color=color,
                   label=rf"$\theta_i={r['theta_i']:.0f},\,g={r['g']:.2f}$", **STYLE)
        ax_tt.plot(lz, bg["phi_scf"] / F_TEDS, color=color, **STYLE)

    for color, r in zip(cos3_colors, cos3_cases):
        bg = r["bg"]
        z = bg["z"]
        lz = np.log10(np.maximum(z, 1e-3))
        ax_fc.plot(lz, bg["(.)Omega_scf"], color=color,
                   label=rf"$\delta={r['delta']:.3g},\,g={r['g']:.2f}$", **STYLE)
        ax_tc.plot(lz, bg["phi_scf"] / F_AXION, color=color, **STYLE)

    for ax in axes.flat:
        ax.axvline(math.log10(Z_EQ), color="red", ls="--", lw=0.9, alpha=0.7)
        ax.set_xlim(1.8, 6.5)
        ax.tick_params(direction="in", top=True, right=True)

    ax_ft.set_ylabel(r"$f_{\rm EDE}(z)$")
    ax_fc.set_ylabel(r"$f_{\rm EDE}(z)$")
    ax_tt.set_ylabel(r"$\theta = \phi/f$")
    ax_tc.set_ylabel(r"$\theta = \phi/f$")
    ax_tt.set_xlabel(r"$\log_{10} z$")
    ax_tc.set_xlabel(r"$\log_{10} z$")
    ax_ft.set_title("tEDS", fontsize=10)
    ax_fc.set_title(r"$(1-\cos)^3$", fontsize=10)
    ax_ft.legend(frameon=False, fontsize=8)
    ax_fc.legend(frameon=False, fontsize=8)

    for f in ("fig2_histories.pdf", "fig2_histories.png"):
        fig.savefig(OUT / f, bbox_inches="tight", dpi=180)
        fig.savefig(FIGS / f, bbox_inches="tight", dpi=180)
    plt.close(fig)
    print("saved fig2_histories")


# ---------------------------------------------------------------------------
# Fig 3: Coupling fraction at TRIGGER point (onset of rolling), not at peak
# ---------------------------------------------------------------------------

def fig_trigger_coupling(scan: list[dict], teds: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2), sharey=True)
    fig.subplots_adjust(wspace=0.12)

    # Left: tEDS — show f_coup(z) with trigger marked
    teds_colors = plt.cm.Blues(np.linspace(0.45, 0.88, len(teds)))
    for color, r in zip(teds_colors, teds):
        bg = r["bg"]
        z = bg["z"]
        lz = np.log10(np.maximum(z, 1e-3))
        phi = bg["phi_scf"]
        bare = np.abs(bg["V'_scf"])
        coup_arr = np.abs(coupling_force(bg, r["g"]))
        frac = coup_arr / (bare + coup_arr + 1e-300)
        trig = trigger_idx(bg, 0.05)
        fct = r["coup_frac_trigger"]
        axes[0].plot(lz, frac, color=color,
                     label=(rf"$\theta_i={r['theta_i']:.0f},\,g={r['g']:.2f}$"
                            rf"  $f_\mathrm{{coup}}^\mathrm{{trig}}={fct:.2f}$"),
                     **STYLE)
        axes[0].axvline(math.log10(z[trig]), color=color, ls=":", lw=1.2, alpha=0.7)

    # Right: (1-cos)^3 best cases (g=0.68, various delta)
    cos3_cases = [r for r in scan if abs(r["g"] - 0.68) < 0.01 and r["delta"] in {0.5, 0.2, 0.05, 0.01}]
    cos3_cases.sort(key=lambda r: r["delta"], reverse=True)
    cos3_colors = plt.cm.Oranges(np.linspace(0.40, 0.90, len(cos3_cases)))
    for color, r in zip(cos3_colors, cos3_cases):
        bg = r["bg"]
        z = bg["z"]
        lz = np.log10(np.maximum(z, 1e-3))
        bare = np.abs(bg["V'_scf"])
        coup_arr = np.abs(coupling_force(bg, r["g"]))
        frac = coup_arr / (bare + coup_arr + 1e-300)
        trig = trigger_idx(bg, 0.05)
        fct = r["coup_frac_trigger"]
        axes[1].plot(lz, frac, color=color,
                     label=(rf"$\delta={r['delta']:.3g},\,g={r['g']:.2f}$"
                            rf"  $f_\mathrm{{coup}}^\mathrm{{trig}}={fct:.2f}$"),
                     **STYLE)
        axes[1].axvline(math.log10(z[trig]), color=color, ls=":", lw=1.2, alpha=0.7)

    for ax in axes:
        ax.axvline(math.log10(Z_EQ), color="red", ls="--", lw=0.9, alpha=0.7, label=r"$z_{\rm eq}$")
        ax.axhline(0.5, color="black", ls=":", lw=1.0)
        ax.set_xlim(1.8, 6.5)
        ax.set_ylim(-0.05, 1.05)
        ax.tick_params(direction="in", top=True, right=True)
        ax.set_xlabel(r"$\log_{10} z$")
        ax.legend(frameon=False, fontsize=7.5)

    axes[0].set_ylabel(r"$f_{\rm coup}(z)$")
    axes[0].set_title("tEDS  [dotted = trigger point]")
    axes[1].set_title(r"$(1-\cos)^3$  [dotted = trigger point]")

    for f in ("fig3_trigger_coupling.pdf", "fig3_trigger_coupling.png"):
        fig.savefig(OUT / f, bbox_inches="tight", dpi=180)
        fig.savefig(FIGS / f, bbox_inches="tight", dpi=180)
    plt.close(fig)
    print("saved fig3_trigger_coupling")


# ---------------------------------------------------------------------------
# Fig 4: Force budget for (1-cos)^3 — bare vs coupling for various delta
# ---------------------------------------------------------------------------

def fig_force_budget(scan: list[dict]) -> None:
    # Fix g=0.68; show delta = {2, 0.5, 0.1, 0.01}
    target_g = 0.68
    target_deltas = [2.0, 0.5, 0.1, 0.01]
    cases = [r for r in scan if abs(r["g"] - target_g) < 0.01 and r["delta"] in target_deltas]
    cases.sort(key=lambda r: r["delta"], reverse=True)
    if not cases:
        print("fig4: no cases found, skipping")
        return

    colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(cases)))

    fig, axes = plt.subplots(3, 1, figsize=(5.5, 7.0), sharex=True)
    fig.subplots_adjust(hspace=0.07)

    for color, r in zip(colors, cases):
        bg = r["bg"]
        z = bg["z"]
        lz = np.log10(np.maximum(z, 1e-3))
        bare = np.abs(bg["V'_scf"])
        coup_arr = np.abs(coupling_force(bg, r["g"]))
        frac = coup_arr / (bare + coup_arr + 1e-300)
        lbl = rf"$\delta={r['delta']:.3g}$"
        axes[0].semilogy(lz, bare, color=color, lw=1.8, label=lbl)
        axes[1].semilogy(lz, coup_arr, color=color, lw=1.8)
        axes[2].plot(lz, frac, color=color, lw=1.8)

    for ax in axes:
        ax.axvline(math.log10(Z_EQ), color="red", ls="--", lw=0.9, alpha=0.7)
        ax.set_xlim(1.8, 6.5)
        ax.tick_params(direction="in", top=True, right=True)

    axes[2].axhline(0.5, color="black", ls=":", lw=1.0)
    axes[0].set_ylabel(r"$|V_{,\phi}|$ [CLASS]")
    axes[1].set_ylabel(r"$|3\rho_{\rm IDM}\partial_\phi\!\ln m|$")
    axes[2].set_ylabel(r"$f_{\rm coup}(z)$")
    axes[2].set_xlabel(r"$\log_{10} z$")
    axes[2].set_ylim(-0.05, 1.05)
    axes[0].legend(frameon=False, fontsize=9,
                   title=rf"$(1-\cos)^3$, $g={target_g:.2f}$", title_fontsize=8.5)

    for f in ("fig4_force_budget.pdf", "fig4_force_budget.png"):
        fig.savefig(OUT / f, bbox_inches="tight", dpi=180)
        fig.savefig(FIGS / f, bbox_inches="tight", dpi=180)
    plt.close(fig)
    print("saved fig4_force_budget")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("Loading scan data...")
    scan = build_scan_table()
    teds = build_teds_table()
    print(f"  loaded {len(scan)} axion runs, {len(teds)} tEDS runs")

    print("\nTrigger coupling fractions:")
    print(f"  tEDS:")
    for r in teds:
        print(f"    theta_i={r['theta_i']:.1f} g={r['g']:.2f}  "
              f"f_coup_trigger={r['coup_frac_trigger']:.3f}  "
              f"f_coup_peak={r['coup_frac_peak']:.3f}  "
              f"z_trigger=10^{math.log10(r['z_trigger']):.2f}  "
              f"z_peak=10^{math.log10(r['z_peak']):.2f}")

    print(f"  (1-cos)^3 g=0.68 representative cases:")
    for r in sorted([r for r in scan if abs(r["g"]-0.68)<0.01], key=lambda x: x["delta"]):
        print(f"    delta={r['delta']:.4g}  "
              f"f_coup_trigger={r['coup_frac_trigger']:.3f}  "
              f"f_coup_peak={r['coup_frac_peak']:.3f}  "
              f"z_peak=10^{math.log10(r['z_peak']):.2f}  "
              f"f_EDE={r['f_ede_peak']:.3f}")

    print("\nMaking figures...")
    fig_zpeak_vs_delta(scan)
    fig_histories(scan, teds)
    fig_trigger_coupling(scan, teds)
    fig_force_budget(scan)
    print(f"\nFigures saved to {FIGS} and {OUT}")


if __name__ == "__main__":
    main()
