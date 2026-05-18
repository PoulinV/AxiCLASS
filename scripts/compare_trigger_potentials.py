#!/usr/bin/env python3
"""Compare DM-triggered release for axion-like and tEDS potentials.

The comparison is intentionally background-only.  It scans a small grid for
axion-like EDE with n=2 and n=3 and compares it to the plateau power-law tEDS
potential using the same quadratic DM mass coupling.
"""

from __future__ import annotations

import csv
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
CLASS = REPO / "class"
OUT_ROOT = REPO / "analysis" / "potential_trigger_compare"

H0 = 72.81
OMEGA_INI_IDM_EDE = 0.1280
OMEGA_LAMBDA = 0.685
ZEQ = 3400.0

EV4_TO_CLASS = (1.56e29 / 2.435e27) ** 2 / 3.0
V0_TEDS = 0.12 * EV4_TO_CLASS
F_TEDS = 0.05
P_TEDS = 8.0


@dataclass(frozen=True)
class Case:
    model: str
    n: float
    theta_i: float
    f: float
    m: float
    g: float
    V0: float = 0.0
    p: float = 0.0


def common_ini(root: Path, g: float) -> str:
    return f"""root = {root}
output = tCl
background_evolver = 0
evolver = 0
write background = yes
write parameters = no
omega_b = 0.02251
omega_ini_idm_ede = {OMEGA_INI_IDM_EDE}
Omega_Lambda = {OMEGA_LAMBDA}
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


def build_ini(root: Path, case: Case) -> str:
    if case.model == "teds":
        scf = f"""scf_potential = teds
scf_parameters = {case.theta_i * case.f},0.0
V0_teds = {case.V0}
f_teds = {case.f}
p_teds = {case.p}
"""
    else:
        scf = f"""scf_potential = axion
n_axion = {case.n}
scf_parameters = {case.theta_i},0.0
f_axion = {case.f}
m_axion = {case.m}
"""
    return common_ini(root, case.g) + scf


def load_background(path: Path) -> dict[str, list[float]]:
    header = ""
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") and "1:z" in line:
                header = line[1:]
                break
    if not header:
        raise RuntimeError(f"missing background header in {path}")

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


def cleanup(root: Path) -> None:
    for suffix in ("00_background.dat", "00_cl.dat", "00_parameters.ini", "00_unused_parameters"):
        path = Path(f"{root}{suffix}")
        if path.exists():
            path.unlink()


def run_case(case: Case, timeout: int = 25) -> dict[str, float | str]:
    tag = (
        f"{case.model}_n{case.n:g}_th{case.theta_i:.8g}_f{case.f:g}_"
        f"m{case.m:g}_g{case.g:g}"
    ).replace(".", "p").replace("-", "m")
    root = OUT_ROOT / tag
    ini = OUT_ROOT / f"{tag}.ini"
    ini.write_text(build_ini(root, case))
    cleanup(root)

    run = subprocess.run(
        [str(CLASS), str(ini)],
        cwd=REPO,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if run.returncode != 0:
        return {
            "ok": "no",
            "model": case.model,
            "n": case.n,
            "theta_i": case.theta_i,
            "f": case.f,
            "m": case.m,
            "g": case.g,
            "error": (run.stderr or run.stdout or "").strip().replace("\n", " ")[-240:],
        }

    bg = load_background(Path(f"{root}00_background.dat"))
    z = bg["z"]
    omega = bg["(.)Omega_scf"]
    phi = bg["phi_scf"]
    bare = [abs(value) for value in bg["V'_scf"]]
    coupling = [abs(value) for value in bg.get("dVadd_contribution", [0.0] * len(z))]
    ratio = [c / max(b, 1e-300) for c, b in zip(coupling, bare)]

    peak = max(range(len(omega)), key=omega.__getitem__)
    eq = min(range(len(z)), key=lambda i: abs(math.log1p(z[i]) - math.log1p(ZEQ)))
    phi0 = phi[0]
    scale = max(abs(phi0), 1e-300)
    release = next((i for i, value in enumerate(phi) if abs(value - phi0) > 0.1 * scale), peak)

    mrel = [1.0 + case.g * value * value for value in phi]
    score = abs(math.log((z[peak] + 1e-12) / ZEQ))
    score += 1.5 * abs(math.log((omega[peak] + 1e-12) / 0.1))
    if ratio[release] < 1.0:
        score += 3.0
    if z[release] < 0.3 * ZEQ or z[release] > 3.0 * ZEQ:
        score += 1.0

    return {
        "ok": "yes",
        "model": case.model,
        "n": case.n,
        "theta_i": case.theta_i,
        "f": case.f,
        "m": case.m,
        "g": case.g,
        "log10_z_peak": math.log10(z[peak]),
        "f_ede_peak": omega[peak],
        "log10_z_release": math.log10(z[release]),
        "force_ratio_release": ratio[release],
        "force_ratio_eq": ratio[eq],
        "dm_mass_shift_max": max(abs(value / mrel[0] - 1.0) for value in mrel),
        "phi_final": phi[-1],
        "score": score,
        "background": str(Path(f"{root}00_background.dat")),
        "error": "",
    }


def cases() -> list[Case]:
    output: list[Case] = []
    for n in (2.0, 3.0):
        for delta in (3e-3, 1e-3):
            for f in (0.3, 1.0):
                for m in (100.0,):
                    for g in (0.0, 0.5, 1.0):
                        output.append(Case("axion", n, math.pi - delta, f, m, g))

    for theta_i, g in ((8.0, 0.75), (6.0, 0.68), (4.0, 0.56)):
        output.append(Case("teds", 3.0, theta_i, F_TEDS, 0.0, g, V0_TEDS, P_TEDS))
        output.append(Case("teds", 3.0, theta_i, F_TEDS, 0.0, 0.0, V0_TEDS, P_TEDS))
    return output


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "ok", "model", "n", "theta_i", "f", "m", "g", "log10_z_peak", "f_ede_peak",
        "log10_z_release", "force_ratio_release", "force_ratio_eq", "dm_mass_shift_max",
        "phi_final", "score", "background", "error",
    ]
    csv_path = OUT_ROOT / "summary.csv"
    rows = []
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        all_cases = cases()
        for index, case in enumerate(all_cases, start=1):
            row = run_case(case)
            rows.append(row)
            writer.writerow(row)
            handle.flush()
            print(f"[{index}/{len(all_cases)}] {case.model} n={case.n:g} theta={case.theta_i:.8g} f={case.f:g} m={case.m:g} g={case.g:g}: {row['ok']}")

    ok_rows = [row for row in rows if row["ok"] == "yes"]
    for model, n in (("axion", 2.0), ("axion", 3.0), ("teds", 3.0)):
        subset = [row for row in ok_rows if row["model"] == model and float(row["n"]) == n]
        subset.sort(key=lambda row: float(row["score"]))
        print(f"\n{model} n={n:g}: {len(subset)} successful runs")
        for row in subset[:5]:
            print(
                "  "
                f"theta={float(row['theta_i']):.8g} f={float(row['f']):g} "
                f"m={float(row['m']):g} g={float(row['g']):g} "
                f"log10_z_peak={float(row['log10_z_peak']):.3f} "
                f"f_peak={float(row['f_ede_peak']):.3f} "
                f"log10_z_rel={float(row['log10_z_release']):.3f} "
                f"R_rel={float(row['force_ratio_release']):.3g} "
                f"score={float(row['score']):.3f}"
            )

    print(f"\nsummary = {csv_path}")


if __name__ == "__main__":
    main()
