#!/usr/bin/env python3
"""Check whether a hilltop (1-cos)^3 potential can be DM-triggered.

The script runs background-only CLASS jobs for

    V(phi) = (m H0)^2 f^2 [1 - cos(phi/f)]^3

with initial theta close to pi and a type-1 quadratic IDM-EDE mass coupling.
It writes a compact PDF report with f_EDE, Omega_idm_ede/Omega_tot, and a
reference LCDM Omega_cdm/Omega_tot curve.
"""

from __future__ import annotations

import csv
import math
import subprocess
import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import numpy as np

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


REPO = Path(__file__).resolve().parents[1]
CLASS = REPO / "class"
OUT = REPO / "analysis" / "axion_cubed_trigger_report"

H0 = 72.81
OMEGA_B = 0.02251
OMEGA_DM_TODAY = 0.1280   # present-day target; CLASS shoots omega_ini internally
Z_EQ = 3400.0
THETA_I = math.pi - 1.0e-3
N_AXION = 3

# Approximate target for f_EDE ~ 0.1 at equality in this branch's
# convention: rho_scf = V/3, V(pi) = 8 (m H0)^2 f^2.
M_TIMES_F_TARGETS = (1.0e4, 2.0e4, 3.0e4)


@dataclass(frozen=True)
class Case:
    m_axion: float
    f_axion: float
    q_initial: float
    mf_target: float

    @property
    def phi_i(self) -> float:
        return THETA_I * self.f_axion

    @property
    def g(self) -> float:
        return self.q_initial / (self.phi_i * self.phi_i)

    @property
    def omega_ini_idm_ede(self) -> float:
        return OMEGA_DM_TODAY * (1.0 + self.q_initial)


def tag_number(value: float) -> str:
    return f"{value:.6g}".replace("-", "m").replace(".", "p")


def class_ini(root: Path, case: Case) -> str:
    # omega_cdm = 0: all DM is IDM-EDE (Lin+ 2212.08098).
    # omega_idm_ede is the TODAY target; CLASS internal shooting finds omega_ini_idm_ede.
    return f"""root = {root}
output = tCl
background_evolver = 0
evolver = 0
write background = yes
write parameters = no
omega_b = {OMEGA_B}
omega_cdm = 0
omega_idm_ede = {OMEGA_DM_TODAY}
H0 = {H0}
tau_reio = 0.068
A_s = 2.191e-9
n_s = 0.9860
coupling_type = 1
idm_ede_mass_form = poly
g_idm_ede = {case.g}
N_ur = 2.0328
N_ncdm = 1
deg_ncdm = 1
m_ncdm = 0.06
T_ncdm = 0.71611
scf_potential = axion
n_axion = {N_AXION}
scf_parameters = {THETA_I},0.0
f_axion = {case.f_axion}
m_axion = {case.m_axion}
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


def lcdm_ini(root: Path) -> str:
    return f"""root = {root}
output = tCl
background_evolver = 0
evolver = 0
write background = yes
write parameters = no
omega_b = {OMEGA_B}
omega_cdm = {OMEGA_DM_TODAY}
Omega_Lambda = {OMEGA_LAMBDA}
H0 = {H0}
tau_reio = 0.068
A_s = 2.191e-9
n_s = 0.9860
N_ur = 2.0328
N_ncdm = 1
deg_ncdm = 1
m_ncdm = 0.06
T_ncdm = 0.71611
modes = s
gauge = synchronous
lensing = no
"""


def load_background(path: Path) -> dict[str, np.ndarray]:
    header = ""
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") and "1:z" in line:
                header = line[1:]
                break
    if not header:
        raise RuntimeError(f"missing background header in {path}")

    names = [entry.split(":", 1)[1] for entry in header.split() if ":" in entry]
    data = np.loadtxt(path)
    return {name: data[:, index] for index, name in enumerate(names)}


def cleanup(root: Path) -> None:
    for suffix in ("00_background.dat", "00_cl.dat", "00_parameters.ini", "00_unused_parameters"):
        path = Path(f"{root}{suffix}")
        if path.exists():
            path.unlink()


def run_ini(ini: Path, root: Path, timeout: int = 35) -> tuple[bool, str]:
    cleanup(root)
    run = subprocess.run(
        [str(CLASS), str(ini)],
        cwd=REPO,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if run.returncode != 0:
        return False, (run.stderr or run.stdout or "").strip()[-600:]
    return True, ""


def run_case(case: Case) -> dict[str, object]:
    tag = f"m{tag_number(case.m_axion)}_f{tag_number(case.f_axion)}_q{tag_number(case.q_initial)}"
    root = OUT / tag
    ini_path = OUT / f"{tag}.ini"
    ini_path.write_text(class_ini(root, case))
    ok, error = run_ini(ini_path, root)
    row: dict[str, object] = {
        "ok": ok,
        "m_axion": case.m_axion,
        "f_axion": case.f_axion,
        "q_initial": case.q_initial,
        "mf_target": case.mf_target,
        "g_idm_ede": case.g,
        "omega_ini_idm_ede": case.omega_ini_idm_ede,
        "background": str(Path(f"{root}00_background.dat")),
        "error": error,
    }
    if not ok:
        return row

    bg = load_background(Path(f"{root}00_background.dat"))
    z = bg["z"]
    f_ede = bg["(.)Omega_scf"]
    omega_idm = bg["(.)rho_idm_ede"] / bg["(.)rho_tot"]
    phi = bg["phi_scf"]
    peak = int(np.argmax(f_ede))
    eq = int(np.argmin(np.abs(np.log1p(z) - np.log1p(Z_EQ))))
    release_candidates = np.where(np.abs(phi - phi[0]) > 0.1 * max(abs(phi[0]), 1e-300))[0]
    release = int(release_candidates[0]) if release_candidates.size else peak

    row.update(
        {
            "log10_z_peak": math.log10(z[peak]),
            "f_ede_peak": f_ede[peak],
            "f_ede_eq": f_ede[eq],
            "log10_z_release": math.log10(z[release]),
            "omega_idm_frac_peak": omega_idm[peak],
            "omega_idm_frac_eq": omega_idm[eq],
            "phi_final": phi[-1],
        }
    )
    return row


def score(row: dict[str, object]) -> float:
    if not row["ok"]:
        return float("inf")
    return (
        abs(math.log(float(row["f_ede_peak"]) / 0.1))
        + 0.75 * abs(float(row["log10_z_peak"]) - math.log10(Z_EQ))
    )


def run_lcdm_reference() -> dict[str, np.ndarray]:
    root = OUT / "lcdm_reference"
    ini = OUT / "lcdm_reference.ini"
    ini.write_text(lcdm_ini(root))
    ok, error = run_ini(ini, root)
    if not ok:
        raise RuntimeError(error)
    return load_background(Path(f"{root}00_background.dat"))


def make_cases() -> list[Case]:
    cases: list[Case] = []
    for f_axion, m_values, q_values in (
        (1.0, (7000.0, 10000.0, 20000.0, 30000.0), (1.0, 3.0, 10.0, 30.0)),
        (0.3333333333, (30000.0,), (0.1, 0.15, 0.2, 0.3, 1.0)),
    ):
        for m_axion in m_values:
            for q_initial in q_values:
                cases.append(Case(m_axion, f_axion, q_initial, m_axion * f_axion))
    return cases


def write_summary(rows: list[dict[str, object]]) -> Path:
    path = OUT / "summary.csv"
    fields = [
        "ok",
        "m_axion",
        "f_axion",
        "q_initial",
        "mf_target",
        "g_idm_ede",
        "omega_ini_idm_ede",
        "log10_z_peak",
        "f_ede_peak",
        "f_ede_eq",
        "log10_z_release",
        "omega_idm_frac_peak",
        "omega_idm_frac_eq",
        "phi_final",
        "background",
        "error",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})
    return path


def read_summary(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            row["ok"] = row["ok"] == "True"
            for key in (
                "m_axion",
                "f_axion",
                "q_initial",
                "mf_target",
                "g_idm_ede",
                "omega_ini_idm_ede",
                "log10_z_peak",
                "f_ede_peak",
                "f_ede_eq",
                "log10_z_release",
                "omega_idm_frac_peak",
                "omega_idm_frac_eq",
                "phi_final",
            ):
                if row.get(key):
                    row[key] = float(row[key])
            rows.append(row)
    return rows


def plot_report(selected: list[dict[str, object]], lcdm: dict[str, np.ndarray]) -> Path:
    report = OUT / "axion_cubed_trigger_report.pdf"
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(selected)))
    z_lcdm = lcdm["z"]
    cdm_frac = lcdm["(.)rho_cdm"] / lcdm["(.)rho_tot"]

    with PdfPages(report) as pdf:
        fig, axes = plt.subplots(3, 1, figsize=(7.4, 8.6), sharex=True)
        for color, row in zip(colors, selected):
            bg = load_background(Path(str(row["background"])))
            z = bg["z"]
            label = (
                rf"$m/H_0={float(row['m_axion']):g}$, "
                rf"$f={float(row['f_axion']):.3g}$, "
                rf"$q_i={float(row['q_initial']):g}$"
            )
            axes[0].semilogx(1.0 + z, bg["(.)Omega_scf"], lw=1.9, color=color, label=label)
            axes[1].semilogx(1.0 + z, bg["(.)rho_idm_ede"] / bg["(.)rho_tot"], lw=1.9, color=color)
            axes[2].semilogx(1.0 + z, bg["phi_scf"] / float(row["f_axion"]), lw=1.6, color=color)

        axes[1].semilogx(1.0 + z_lcdm, cdm_frac, color="black", lw=2.0, ls="--", label=r"LCDM $\Omega_{\rm cdm}/\Omega_{\rm tot}$")

        for ax in axes:
            ax.axvline(1.0 + Z_EQ, color="red", ls="--", lw=1.0, alpha=0.7)
            ax.set_xlim(1.0e6, 1.0)
            ax.grid(alpha=0.25)
            ax.tick_params(direction="in", top=True, right=True)

        axes[0].set_ylabel(r"$f_{\rm EDE}=\Omega_\phi$")
        axes[1].set_ylabel(r"$\Omega_{\rm IDM-EDE}/\Omega_{\rm tot}$")
        axes[2].set_ylabel(r"$\theta=\phi/f$")
        axes[2].set_xlabel(r"$1+z$")
        axes[0].legend(frameon=False, fontsize=8)
        axes[1].legend(frameon=False, fontsize=8)
        fig.tight_layout()
        pdf.savefig(fig)
        fig.savefig(OUT / "selected_histories.png", dpi=180)
        plt.close(fig)

        fig = plt.figure(figsize=(7.4, 5.0))
        ax = fig.add_subplot(111)
        ok_rows = [row for row in selected if row["ok"]]
        for row in ok_rows:
            ax.scatter(float(row["log10_z_peak"]), float(row["f_ede_peak"]), s=60)
            ax.annotate(
                f"m={float(row['m_axion']):g}, q={float(row['q_initial']):g}",
                (float(row["log10_z_peak"]), float(row["f_ede_peak"])),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=8,
            )
        ax.axvline(math.log10(Z_EQ), color="red", ls="--", lw=1)
        ax.axhline(0.1, color="black", ls=":", lw=1)
        ax.set_xlabel(r"$\log_{10} z_{\rm peak}$")
        ax.set_ylabel(r"$f_{\rm EDE,peak}$")
        ax.grid(alpha=0.25)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--from-summary",
        action="store_true",
        help="Only rebuild plots from analysis/axion_cubed_trigger_report/summary.csv.",
    )
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    if args.from_summary:
        summary = OUT / "summary.csv"
        rows = read_summary(summary)
    else:
        rows = []
        cases = make_cases()
        for index, case in enumerate(cases, start=1):
            row = run_case(case)
            rows.append(row)
            status = "ok" if row["ok"] else "failed"
            peak = row.get("f_ede_peak", "")
            zpeak = row.get("log10_z_peak", "")
            print(
                f"[{index}/{len(cases)}] m={case.m_axion:g} f={case.f_axion:g} "
                f"q={case.q_initial:g} mf={case.mf_target:g}: {status} f_peak={peak} log10z={zpeak}",
                flush=True,
            )
        summary = write_summary(rows)

    ok_rows = [row for row in rows if row["ok"]]
    ok_rows.sort(key=score)
    selected = ok_rows[:5]
    lcdm = run_lcdm_reference()
    report = plot_report(selected, lcdm)

    print(f"summary = {summary}")
    print(f"report = {report}")
    print("selected:")
    for row in selected:
        print(
            f"  m={float(row['m_axion']):g}, f={float(row['f_axion']):g}, "
            f"q={float(row['q_initial']):g}, mf={float(row['mf_target']):g}, "
            f"g={float(row['g_idm_ede']):.4e}, "
            f"log10_z_peak={float(row['log10_z_peak']):.4f}, "
            f"f_peak={float(row['f_ede_peak']):.4f}, "
            f"Omega_idm/Om_tot(eq)={float(row['omega_idm_frac_eq']):.4f}"
        )


if __name__ == "__main__":
    main()
