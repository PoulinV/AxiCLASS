#!/usr/bin/env python3
"""Background-level tEDS trigger check for arXiv:2212.08098.

This reproduces the qualitative Fig. 2 test: for the plateau power-law
potential with p=8 and f=0.05 Mpl, an order-unity quadratic DM-mass coupling
triggers the field near equality for several initial plateau positions.
"""

from __future__ import annotations

import csv
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt


REPO = Path(__file__).resolve().parents[1]
OUT_ROOT = REPO / "analysis" / "teds_trigger_check"
CLASS = REPO / "class"

H0 = 72.81
OMEGA_IDM_EDE = 0.1280
F_TEDS = 0.05
P_TEDS = 8.0

# CLASS stores scalar-field densities as (K+V)/3. In this branch the
# historical AxiCLASS convention has been to pass the eV^4/Mpl^2 potential
# normalization with the same factor of 1/3 used for CLASS densities.
EV4_TO_CLASS = (1.56e29 / 2.435e27) ** 2  # 1 eV^4 in CLASS units (m_Pl^2/Mpc^2); no /3 since background.c already divides rho_scf = (K+V)/3
V0_TEDS = 0.12 * EV4_TO_CLASS  # Lin+ 2212.08098 Fig 1 caption states V0 = 0.12 eV^4


@dataclass(frozen=True)
class Case:
    theta_i: float
    g: float


CASES = (
    Case(theta_i=8.0, g=0.75),
    Case(theta_i=6.0, g=0.68),
    Case(theta_i=4.0, g=0.56),
)


def build_ini(root: Path, case: Case, coupled: bool = True) -> str:
    g = case.g if coupled else 0.0
    # omega_cdm = 0: all DM is IDM-EDE (Lin+ 2212.08098).
    # omega_idm_ede is the TODAY target; CLASS internal shooting finds omega_ini_idm_ede.
    return f"""root = {root}
output = tCl
background_evolver = 0
evolver = 0
write background = yes
write parameters = no
omega_b = 0.02251
omega_cdm = 0
omega_idm_ede = {OMEGA_IDM_EDE}
H0 = {H0}
tau_reio = 0.068
A_s = 2.191e-9
n_s = 0.9860
coupling_type = 1
idm_ede_mass_form = poly
g_idm_ede = {g}
N_ur = 2.0328
N_ncdm = 1
deg_ncdm = 1
m_ncdm = 0.06
T_ncdm = 0.71611
scf_potential = teds
scf_parameters = {case.theta_i * F_TEDS},0.0
V0_teds = {V0_TEDS}
f_teds = {F_TEDS}
p_teds = {P_TEDS}
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


def load_background(path: Path) -> dict[str, list[float]]:
    header = ""
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") and "1:z" in line:
                header = line[1:]
                break
    if not header:
        raise RuntimeError(f"could not find column header in {path}")

    names: list[str] = []
    for entry in header.split():
        if ":" in entry:
            names.append(entry.split(":", 1)[1])

    columns = {name: [] for name in names}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            values = [float(value) for value in line.split()]
            for name, value in zip(names, values):
                columns[name].append(value)
    return columns


def run_case(case: Case, coupled: bool = True) -> tuple[dict[str, float], dict[str, list[float]]]:
    tag = f"theta{case.theta_i:g}_g{case.g:g}_{'coupled' if coupled else 'free'}".replace(".", "p")
    root = OUT_ROOT / tag
    g = case.g if coupled else 0.0
    phi_i = case.theta_i * F_TEDS

    ini_path = OUT_ROOT / f"{tag}.ini"
    ini_path.write_text(build_ini(root, case, coupled))
    for suffix in ("00_background.dat", "00_cl.dat", "00_parameters.ini", "00_unused_parameters"):
        p = Path(f"{root}{suffix}")
        if p.exists():
            p.unlink()
    subprocess.run([str(CLASS), str(ini_path)], cwd=REPO,
                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                   text=True, check=True, timeout=120)

    bg = load_background(Path(f"{root}00_background.dat"))
    omega = bg["(.)Omega_scf"]
    z = bg["z"]
    peak = max(range(len(omega)), key=omega.__getitem__)
    phi = bg["phi_scf"]
    mass_shift = [(case.g if coupled else 0.0) * value * value for value in phi]
    bare = [abs(value) for value in bg["V'_scf"]]
    if "dVadd_contribution" in bg:
        coupling = [abs(value) for value in bg["dVadd_contribution"]]
    else:
        coupling = [0.0 for _ in bare]
    ratio = [c / max(b, 1e-300) for c, b in zip(coupling, bare)]
    early_ratio = [value for value, zi in zip(ratio, z) if zi > 100.0]

    return {
        "theta_i": case.theta_i,
        "g": case.g if coupled else 0.0,
        "log10_z_peak": math.log10(z[peak]),
        "f_ede_peak": omega[peak],
        "force_ratio_at_peak": ratio[peak],
        "max_force_ratio": max(early_ratio),
        "dm_mass_shift_max": max(mass_shift),
        "phi_final": phi[-1],
    }, bg


def plot_figure2(histories: list[tuple[Case, dict[str, list[float]]]], summary: list[dict[str, float]]) -> None:
    colors = ("#ff7f0e", "#2ca02c", "#d62728")
    fig, axes = plt.subplots(3, 1, figsize=(5.2, 6.4), sharex=True)

    peak_log10z = summary[0]["log10_z_peak"]

    for color, (case, bg) in zip(colors, histories):
        z = np.asarray(bg["z"])
        phi = np.asarray(bg["phi_scf"])
        theta = phi / F_TEDS
        log10z = np.log10(np.maximum(z, 1e-300))

        mass_shift = case.g * phi**2

        label = rf"$\theta_i={case.theta_i:g},\,g={case.g:.2f}$"
        axes[0].plot(log10z, theta, color=color, lw=1.8, label=label)
        axes[1].plot(log10z, mass_shift, color=color, lw=1.8)
        axes[2].plot(log10z, bg["(.)Omega_scf"], color=color, lw=1.8)

    for ax in axes:
        ax.axvline(peak_log10z, color="red", ls="--", lw=0.9, alpha=0.75)
        ax.set_xlim(2.0, 6.0)
        ax.tick_params(direction="in", top=True, right=True)

    axes[0].axhline(1.0, color="black", lw=0.8)
    axes[0].set_ylim(-0.4, 8.2)
    axes[0].set_ylabel(r"$\theta=\phi/f$")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")

    axes[1].set_ylim(-0.005, 0.125)
    axes[1].set_ylabel(r"$\Delta m_{\rm DM}/m_{\rm DM}$")

    axes[2].set_ylim(-0.005, 0.12)
    axes[2].set_ylabel(r"$f_{\rm EDE}$")
    axes[2].set_xlabel(r"$\log_{10} z$")

    fig.tight_layout(h_pad=0.15)
    fig.savefig(OUT_ROOT / "figure2_reproduction.png", dpi=220)
    fig.savefig(OUT_ROOT / "figure2_reproduction.pdf")


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    runs = [run_case(case, coupled=True) for case in CASES]
    rows = [row for row, _ in runs]
    histories = [(case, bg) for case, (_, bg) in zip(CASES, runs)]
    plot_figure2(histories, rows)

    csv_path = OUT_ROOT / "summary.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"V0_teds = {V0_TEDS:.6g} CLASS units")
    print("theta_i,g,log10_z_peak,f_ede_peak,force_ratio_at_peak,max_force_ratio,dm_mass_shift_max,phi_final")
    for row in rows:
        print(",".join(f"{row[key]:.6g}" for key in row))
    print(f"summary = {csv_path}")
    print(f"figure = {OUT_ROOT / 'figure2_reproduction.png'}")


if __name__ == "__main__":
    main()
