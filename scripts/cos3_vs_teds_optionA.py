#!/usr/bin/env python3
"""
Option A: omega_cdm + subdominant IDM-EDE.

Key changes vs the original setup:
  - Add omega_cdm (standard CDM, dominant DM component)
  - Set omega_ini_idm_ede = F_IDM * omega_DM * (1+g*phi_i^2)
    where F_IDM is the IDM-EDE fraction of total DM today
  - Do NOT fix Omega_Lambda -- CLASS derives it from the flat-universe condition
    (otherwise Omega_Lambda=0.685 overrides omega_IDM regardless of omega_ini)
  - Scale g up by 1/F_IDM relative to the original runs to keep the same trigger strength

Why this matters:
  In the original setup (Omega_Lambda fixed, omega_cdm=0), the IDM-EDE today was
  determined entirely by the flat-universe residual (~0.194), not by omega_ini.
  Removing Omega_Lambda as input and specifying omega_cdm lets omega_IDM_today track
  the formula omega_ini/(1+g*phi_i^2) correctly.

  With F_IDM = 0.10 (IDM-EDE = 10% of DM):
    omega_cdm = 0.90 * omega_DM = 0.90 * 0.128 = 0.1152
    omega_ini_idm_ede = 0.10 * 0.128 * (1+g*phi_i^2) = 0.0128 * (1+g*phi_i^2)
    g must be ~10x larger to maintain the same trigger strength
"""

from __future__ import annotations

import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

REPO = Path(__file__).resolve().parents[1]
CLASS = REPO / "class"
OUT = REPO / "analysis" / "cos3_vs_teds_optA"
FIGS = REPO / "report" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

H0 = 72.81
OMEGA_B = 0.02251
OMEGA_DM_TOTAL = 0.1280     # total DM in h^2 units (CDM + IDM)
F_IDM = 0.10                # fraction of total DM that is IDM-EDE
OMEGA_CDM = (1.0 - F_IDM) * OMEGA_DM_TOTAL   # 0.1152
OMEGA_IDM_TODAY = F_IDM * OMEGA_DM_TOTAL      # 0.0128

# tEDS parameters
F_TEDS = 0.05
P_TEDS = 8.0
EV4_TO_CLASS = (1.56e29 / 2.435e27) ** 2 / 3.0
V0_TEDS = 0.12 * EV4_TO_CLASS

# axion parameters
N_AXION = 3
F_AXION = 1.0
M_AXION = 20000.0

# Rescaled g values: original was ~0.68 for full-DM IDM.
# With F_IDM=0.10, need g ~ 0.68/0.10 = 6.8 for same trigger strength.
TEDS_CASES = [
    (8.0, 7.5),
    (6.0, 6.8),
    (4.0, 5.6),
]
COS3_CASES = [
    (math.pi - 0.1,  6.8),
    (math.pi - 0.5,  6.8),
    (math.pi - 0.5,  10.0),
    (math.pi - 1.0,  6.8),
]

# LCDM reference (no IDM-EDE at all)
LCDM_BG = REPO / "analysis" / "axion_cubed_trigger_report" / "lcdm_reference00_background.dat"

COMMON_TAIL = f"""N_ur = 2.0328
N_ncdm = 1
deg_ncdm = 1
m_ncdm = 0.06
T_ncdm = 0.71611
scf_evolve_as_fluid = no
scf_evolve_like_axionCAMB = no
do_shooting = no
do_shooting_scf = no
scf_has_perturbations = no
attractor_ic_scf = no
modes = s
gauge = synchronous
lensing = no
"""


def axion_ini(root: Path, theta_i: float, g: float) -> str:
    phi_i = theta_i * F_AXION
    omega_ini = OMEGA_IDM_TODAY * (1.0 + g * phi_i**2)
    return (
        f"root = {root}\n"
        f"output = tCl\nbackground_evolver = 0\nevolver = 0\n"
        f"write background = yes\nwrite parameters = no\n"
        f"omega_b = {OMEGA_B}\n"
        f"omega_cdm = {OMEGA_CDM}\n"
        f"omega_ini_idm_ede = {omega_ini}\n"
        # No Omega_Lambda: CLASS derives it from flatness
        f"H0 = {H0}\n"
        f"tau_reio = 0.068\n"
        f"A_s = 2.191e-9\nn_s = 0.9860\n"
        f"coupling_type = 1\nidm_ede_mass_form = poly\ng_idm_ede = {g}\n"
        f"scf_potential = axion\nn_axion = {N_AXION}\n"
        f"scf_parameters = {theta_i},0.0\n"
        f"f_axion = {F_AXION}\nm_axion = {M_AXION}\n"
        + COMMON_TAIL
    )


def teds_ini(root: Path, theta_i: float, g: float) -> str:
    phi_i = theta_i * F_TEDS
    omega_ini = OMEGA_IDM_TODAY * (1.0 + g * phi_i**2)
    return (
        f"root = {root}\n"
        f"output = tCl\nbackground_evolver = 0\nevolver = 0\n"
        f"write background = yes\nwrite parameters = no\n"
        f"omega_b = {OMEGA_B}\n"
        f"omega_cdm = {OMEGA_CDM}\n"
        f"omega_ini_idm_ede = {omega_ini}\n"
        # No Omega_Lambda: CLASS derives it from flatness
        f"H0 = {H0}\n"
        f"tau_reio = 0.068\n"
        f"A_s = 2.191e-9\nn_s = 0.9860\n"
        f"coupling_type = 1\nidm_ede_mass_form = poly\ng_idm_ede = {g}\n"
        f"scf_potential = teds\n"
        f"scf_parameters = {phi_i},0.0\n"
        f"V0_teds = {V0_TEDS}\nf_teds = {F_TEDS}\np_teds = {P_TEDS}\n"
        + COMMON_TAIL
    )


def load_bg(path: Path) -> Optional[dict[str, np.ndarray]]:
    try:
        header = ""
        with path.open() as f:
            for line in f:
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


def run_class(ini_text: str, root: Path, timeout: int = 90) -> bool:
    ini = Path(str(root) + ".ini")
    ini.write_text(ini_text)
    for suf in ("00_background.dat", "00_cl.dat", "00_parameters.ini", "00_unused_parameters"):
        p = Path(str(root) + suf)
        if p.exists():
            p.unlink()
    res = subprocess.run(
        [str(CLASS), str(ini)],
        cwd=REPO, capture_output=True, text=True, timeout=timeout,
    )
    if res.returncode != 0:
        print(f"  CLASS failed: {(res.stderr or res.stdout or '')[-300:].strip()}", file=sys.stderr)
    return res.returncode == 0


def summarise(bg: dict, label: str) -> None:
    z = bg["z"]
    f_ede = bg["(.)Omega_scf"]
    rho_idm = bg["(.)rho_idm_ede"]
    rho_cdm = bg.get("(.)rho_cdm", np.zeros_like(z))
    rho_tot = bg["(.)rho_tot"]
    peak = int(np.argmax(f_ede))
    omega_idm_today = float(rho_idm[-1] / rho_tot[-1])
    omega_cdm_today = float(rho_cdm[-1] / rho_tot[-1]) if rho_cdm.any() else 0.0
    print(f"  {label:50s}  f_EDE={f_ede[peak]:.4f}  z_peak={z[peak]:.0f}  "
          f"Omega_IDM_today={omega_idm_today:.4f}  Omega_CDM_today={omega_cdm_today:.4f}")


STYLE = dict(lw=1.8, solid_capstyle="round")
Z_EQ = 3400.0


def make_fig(teds_bgs, cos3_bgs, lcdm):
    """Four-panel: Omega_IDM/Omega_tot and f_EDE for tEDS (left) and (1-cos)^3 (right)."""
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 6.0), sharex=True)
    fig.subplots_adjust(hspace=0.07, wspace=0.28)
    ax_ot, ax_oc = axes[0]
    ax_ft, ax_fc = axes[1]

    z_ref = lcdm["z"]
    omega_cdm_ref = lcdm["(.)rho_cdm"] / lcdm["(.)rho_tot"]

    colors_t = plt.cm.Blues(np.linspace(0.45, 0.88, len(teds_bgs)))
    colors_c = plt.cm.Oranges(np.linspace(0.40, 0.90, len(cos3_bgs)))

    for color, (label, bg) in zip(colors_t, teds_bgs):
        if bg is None:
            continue
        z = bg["z"]
        x = 1 + z
        omega_idm = bg["(.)rho_idm_ede"] / bg["(.)rho_tot"]
        ax_ot.loglog(x, omega_idm, color=color, lw=1.8, label=label)
        ax_ft.loglog(x, np.maximum(bg["(.)Omega_scf"], 1e-7), color=color, lw=1.8)

    for color, (label, bg) in zip(colors_c, cos3_bgs):
        if bg is None:
            continue
        z = bg["z"]
        x = 1 + z
        omega_idm = bg["(.)rho_idm_ede"] / bg["(.)rho_tot"]
        ax_oc.loglog(x, omega_idm, color=color, lw=1.8, label=label)
        ax_fc.loglog(x, np.maximum(bg["(.)Omega_scf"], 1e-7), color=color, lw=1.8)

    # LCDM reference
    for ax in [ax_ot, ax_oc]:
        ax.loglog(1 + z_ref, omega_cdm_ref, color="black",
                  lw=1.4, ls="--", label=r"$\Lambda$CDM $\Omega_{\rm cdm}$", zorder=2)
        ax.axvline(1 + Z_EQ, color="gray", ls=":", lw=1.0, alpha=0.8)

    for ax in axes.flat:
        ax.axvline(1 + Z_EQ, color="gray", ls=":", lw=1.0, alpha=0.8)
        ax.set_xlim(1.5, 2e6)
        ax.tick_params(direction="in", top=True, right=True, which="both")
        ax.invert_xaxis()

    ax_ot.set_ylabel(r"$\Omega_{\rm IDM\text{-}EDE}/\Omega_{\rm tot}$")
    ax_oc.set_ylabel(r"$\Omega_{\rm IDM\text{-}EDE}/\Omega_{\rm tot}$")
    ax_ft.set_ylabel(r"$f_{\rm EDE}(z)$")
    ax_fc.set_ylabel(r"$f_{\rm EDE}(z)$")
    ax_ft.set_xlabel(r"$1+z$")
    ax_fc.set_xlabel(r"$1+z$")
    ax_ot.set_title(r"tEDS  ($F_{\rm IDM}=10\%$ of DM)", fontsize=10)
    ax_oc.set_title(r"$(1-\cos)^3$  ($F_{\rm IDM}=10\%$ of DM)", fontsize=10)
    ax_ot.legend(frameon=False, fontsize=7.5)
    ax_oc.legend(frameon=False, fontsize=7.5)
    ax_ft.legend(frameon=False, fontsize=7.5)
    ax_fc.legend(frameon=False, fontsize=7.5)

    note = (rf"$\omega_{{\rm cdm}}={OMEGA_CDM:.4f}$, $\omega_{{\rm IDM,today}}\approx{OMEGA_IDM_TODAY:.4f}$,"
            r" $\Omega_\Lambda$ derived from flatness.")
    fig.text(0.5, 0.005, note, ha="center", fontsize=8, color="gray")

    for fname in ("fig_omega_optA.pdf", "fig_omega_optA.png"):
        fig.savefig(OUT / fname, bbox_inches="tight")
        fig.savefig(FIGS / fname, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("saved fig_omega_optA")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Option A: omega_cdm={OMEGA_CDM:.4f}, omega_IDM_today={OMEGA_IDM_TODAY:.4f}, F_IDM={F_IDM:.0%}")
    print(f"  g values rescaled by 1/F_IDM = {1/F_IDM:.0f}x relative to original\n")

    # Run tEDS cases
    teds_bgs = []
    for theta_i, g in TEDS_CASES:
        tag = f"teds_th{theta_i:.4g}_g{g:.4g}".replace(".", "p")
        root = OUT / tag
        ini = teds_ini(root, theta_i, g)
        label = rf"tEDS $\theta_i={theta_i:.0f},\,g={g:.1f}$"
        print(f"tEDS theta_i={theta_i} g={g} ...", end=" ", flush=True)
        ok = run_class(ini, root)
        if ok:
            bg = load_bg(Path(str(root) + "00_background.dat"))
            if bg is not None:
                summarise(bg, label)
                teds_bgs.append((label, bg))
            else:
                print("load failed"); teds_bgs.append((label, None))
        else:
            print("CLASS failed"); teds_bgs.append((label, None))

    print()

    # Run (1-cos)^3 cases
    cos3_bgs = []
    for theta_i, g in COS3_CASES:
        delta = math.pi - theta_i
        tag = f"axion_d{delta:.4g}_g{g:.4g}".replace(".", "p")
        root = OUT / tag
        ini = axion_ini(root, theta_i, g)
        label = rf"$(1-\cos)^3$ $\delta={delta:.2g},\,g={g:.1f}$"
        print(f"(1-cos)^3 delta={delta:.3g} g={g} ...", end=" ", flush=True)
        ok = run_class(ini, root)
        if ok:
            bg = load_bg(Path(str(root) + "00_background.dat"))
            if bg is not None:
                summarise(bg, label)
                cos3_bgs.append((label, bg))
            else:
                print("load failed"); cos3_bgs.append((label, None))
        else:
            print("CLASS failed"); cos3_bgs.append((label, None))

    print()
    lcdm = load_bg(LCDM_BG)
    make_fig(teds_bgs, cos3_bgs, lcdm)
    print(f"\nFigures saved to {FIGS}")


if __name__ == "__main__":
    main()
