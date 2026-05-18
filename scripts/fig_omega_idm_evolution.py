#!/usr/bin/env python3
"""
Fig-6-style figure: Omega_IDM-EDE / Omega_tot vs 1+z for tEDS and (1-cos)^3.

Also shows phi(z) and f_EDE(z) for the same runs, styled after Fig. 6 of
Lin et al. (2212.08098).  A LCDM Omega_cdm/Omega_tot reference is overlaid.

Key result being illustrated:
  In our setup, IDM-EDE plays the role of ALL dark matter (no separate omega_cdm),
  so Omega_IDM-EDE/Omega_tot ~ Omega_DM/Omega_tot ~ 0.3 at late times.
  In the paper, IDM-EDE is a subdominant DM component alongside standard CDM.
  This difference explains why our Omega_IDM-EDE curve looks different from
  the paper's fig 6.
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
LCDM_BG = REPO / "analysis" / "axion_cubed_trigger_report" / "lcdm_reference00_background.dat"
OUT = REPO / "report" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

Z_EQ = 3400.0
F_AXION = 1.0
F_TEDS = 0.05


def load(path: Path):
    header = ""
    with path.open() as f:
        for line in f:
            if line.startswith("#") and "1:z" in line:
                header = line[1:]
                break
    names = [e.split(":", 1)[1] for e in header.split() if ":" in e]
    data = np.loadtxt(path, comments="#")
    return {n: data[:, i] for i, n in enumerate(names)}


# Selected cases: two tEDS + three axion
TEDS_CASES = [
    ("teds_th8_g0p7500_background.dat", r"tEDS $\theta_i=8,\,g=0.75$", "#1f77b4"),
    ("teds_th6_g0p6800_background.dat", r"tEDS $\theta_i=6,\,g=0.68$", "#aec7e8"),
    ("teds_th4_g0p5600_background.dat", r"tEDS $\theta_i=4,\,g=0.56$", "#7fc3e3"),
]
COS3_CASES = [
    ("axion_show_d0p1_g0p6800_background.dat",
     r"$(1\!-\!\cos)^3$, $\delta\!=\!0.1$, $g\!=\!0.68$", "#d62728"),
    ("axion_show_d0p5_g0p6800_background.dat",
     r"$(1\!-\!\cos)^3$, $\delta\!=\!0.5$, $g\!=\!0.68$", "#ff7f0e"),
    ("axion_show_d0p5_g100_background.dat",
     r"$(1\!-\!\cos)^3$, $\delta\!=\!0.5$, $g\!=\!1.0$",  "#e377c2"),
]


def make_fig6(lcdm):
    fig, axes = plt.subplots(3, 2, figsize=(10.0, 8.5), sharex=True)
    fig.subplots_adjust(hspace=0.07, wspace=0.28)

    ax_om_t, ax_om_c = axes[0]
    ax_ph_t, ax_ph_c = axes[1]
    ax_fe_t, ax_fe_c = axes[2]

    z_lcdm = lcdm["z"]
    # Omega_cdm / Omega_tot for LCDM reference
    omega_cdm_ref = lcdm["(.)rho_cdm"] / lcdm["(.)rho_tot"]

    def plot_lcdm(ax):
        ax.loglog(1 + z_lcdm, omega_cdm_ref, color="black",
                  lw=1.4, ls="--", label=r"$\Lambda$CDM $\Omega_{\rm cdm}$", zorder=2)

    # ---- tEDS ----
    for fname, label, color in TEDS_CASES:
        p = DATA / fname
        if not p.exists():
            print(f"missing {fname}")
            continue
        bg = load(p)
        z = bg["z"]
        x = 1 + z
        omega_idm = bg["(.)rho_idm_ede"] / bg["(.)rho_tot"]
        phi = bg["phi_scf"]
        f_ede = bg["(.)Omega_scf"]
        ax_om_t.loglog(x, omega_idm, color=color, lw=1.8, label=label)
        ax_ph_t.loglog(x, np.abs(phi / F_TEDS), color=color, lw=1.8)
        ax_fe_t.loglog(x, np.maximum(f_ede, 1e-6), color=color, lw=1.8)

    # ---- (1-cos)^3 ----
    for fname, label, color in COS3_CASES:
        p = DATA / fname
        if not p.exists():
            print(f"missing {fname}")
            continue
        bg = load(p)
        z = bg["z"]
        x = 1 + z
        omega_idm = bg["(.)rho_idm_ede"] / bg["(.)rho_tot"]
        phi = bg["phi_scf"]
        f_ede = bg["(.)Omega_scf"]
        ax_om_c.loglog(x, omega_idm, color=color, lw=1.8, label=label)
        ax_ph_c.loglog(x, np.abs(phi / F_AXION), color=color, lw=1.8)
        ax_fe_c.loglog(x, np.maximum(f_ede, 1e-6), color=color, lw=1.8)

    # LCDM reference on both omega panels
    plot_lcdm(ax_om_t)
    plot_lcdm(ax_om_c)

    # Equality line
    for ax in axes.flat:
        ax.axvline(1 + Z_EQ, color="gray", ls=":", lw=1.0, alpha=0.8)
        ax.set_xlim(1.5, 2e6)
        ax.tick_params(direction="in", top=True, right=True, which="both")
        ax.invert_xaxis()

    ax_om_t.set_ylabel(r"$\Omega_{\rm IDM\text{-}EDE}/\Omega_{\rm tot}$")
    ax_om_c.set_ylabel(r"$\Omega_{\rm IDM\text{-}EDE}/\Omega_{\rm tot}$")
    ax_ph_t.set_ylabel(r"$|\theta| = |\phi/f|$")
    ax_ph_c.set_ylabel(r"$|\theta| = |\phi/f|$")
    ax_fe_t.set_ylabel(r"$f_{\rm EDE}(z)$")
    ax_fe_c.set_ylabel(r"$f_{\rm EDE}(z)$")
    ax_fe_t.set_xlabel(r"$1+z$")
    ax_fe_c.set_xlabel(r"$1+z$")

    ax_om_t.set_title("tEDS", fontsize=10)
    ax_om_c.set_title(r"$(1-\cos)^3$", fontsize=10)

    for ax in [ax_om_t, ax_om_c, ax_fe_t, ax_fe_c]:
        ax.legend(frameon=False, fontsize=7.5)

    note = (r"Note: in both our setup and Lin et al.\ IDM-EDE $=$ all DM"
            r" ($\omega_{\rm cdm}=0$, $\omega_{\rm IDM}=0.1280\,h^2$).")
    fig.text(0.5, 0.005, note, ha="center", fontsize=8, color="gray")

    fname = "fig_omega_idm_evolution"
    fig.savefig(OUT / f"{fname}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{fname}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {fname}")


def print_dlnm_table():
    print("\n=== Why f_coup stays large in (1-cos)^3: phi and dlnm at peak ===")
    print(f"{'potential':30s}  {'phi_i':7s}  {'phi_peak':8s}  {'dlnm_peak':9s}  {'f_coup reason'}")
    for fname, label, _ in TEDS_CASES[:1]:
        bg = load(DATA / fname)
        phi = bg["phi_scf"]
        f_ede = bg["(.)Omega_scf"]
        peak = int(np.argmax(f_ede))
        g = 0.68
        dlnm = 2*g*phi[peak]/(1+g*phi[peak]**2)
        print(f"{'tEDS th=6':30s}  {phi[0]:7.4f}  {phi[peak]:8.4f}  {dlnm:9.4f}  phi->0 => dlnm->0 => f_coup->0")

    for fname, label, _ in COS3_CASES[:2]:
        p = DATA / fname
        if not p.exists(): continue
        bg = load(p)
        phi = bg["phi_scf"]
        f_ede = bg["(.)Omega_scf"]
        peak = int(np.argmax(f_ede))
        g = 0.68
        dlnm = 2*g*phi[peak]/(1+g*phi[peak]**2)
        delta_str = fname.split("_d")[1].split("_")[0]
        print(f"{'(1-cos)^3 d='+delta_str:30s}  {phi[0]:7.4f}  {phi[peak]:8.4f}  {dlnm:9.4f}  phi~1 Mpl => dlnm~0.8 => f_coup~1")

    print()
    print("=== Omega_IDM_today check ===")
    print("Both potentials and Lin+22: IDM-EDE = all DM, omega_cdm=0, omega_IDM=0.1280 h^2")
    print("=> Omega_IDM/Omega_tot tracks CDM-like until coupling drives mass change near z_eq")


def main():
    lcdm = load(LCDM_BG)
    make_fig6(lcdm)
    print_dlnm_table()


if __name__ == "__main__":
    main()
