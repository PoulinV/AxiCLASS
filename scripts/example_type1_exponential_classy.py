#!/usr/bin/env python3

"""Minimal classy example for the type-1 exponential DM/EDE coupling.

Run from the repository root after building the local wrapper:

    make class libclass.a
    /usr/local/anaconda3/bin/python python/setup.py build_ext --inplace
    /usr/local/anaconda3/bin/python scripts/example_type1_exponential_classy.py

The model choice is

    coupling_type = 1
    idm_ede_mass_form = exp
    c_idm_ede = 0.1

The script uses classy for the actual CLASS call, then reads the generated
background file to print a compact diagnostic summary.
"""

from __future__ import annotations

import ctypes
import math
import os
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "output" / "example_type1_exp_classy_"


def preload_openmp_runtime() -> None:
    """Preload libgomp when using the local macOS gcc/OpenMP build."""
    candidates = [
        os.environ.get("CLASSY_GOMP_LIB"),
        "/usr/local/Cellar/gcc@12/12.3.0/lib/gcc/12/libgomp.1.dylib",
        "/opt/homebrew/opt/gcc/lib/gcc/current/libgomp.1.dylib",
        "/usr/local/opt/gcc/lib/gcc/current/libgomp.1.dylib",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            ctypes.CDLL(candidate, mode=ctypes.RTLD_GLOBAL)
            return


def read_background(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            rows.append([float(x) for x in line.split()])
    return rows


def nearest_row(rows: list[list[float]], z_target: float) -> list[float]:
    target = math.log1p(z_target)
    return min(rows, key=lambda row: abs(math.log1p(row[0]) - target))


def phi_over_hubble(row: list[float]) -> float:
    z = row[0]
    a = 1.0 / (1.0 + z)
    H = row[3]
    return row[24] / (a * H)


def main() -> None:
    preload_openmp_runtime()
    sys.path.insert(0, str(REPO))

    from classy import Class

    params = {
        "root": str(ROOT),
        "output": "tCl",
        "write background": "yes",
        "write parameters": "no",
        "omega_b": 0.02251,
        "omega_ini_idm_ede": 0.12,
        "H0": 72.81,
        "tau_reio": 0.068,
        "A_s": 2.191e-9,
        "n_s": 0.9860,
        "coupling_type": 1,
        "idm_ede_mass_form": "exp",
        "c_idm_ede": 0.1,
        "N_ur": 2.0328,
        "N_ncdm": 1,
        "deg_ncdm": 1,
        "m_ncdm": 0.06,
        "T_ncdm": 0.71611,
        "scf_potential": "axion",
        "n_axion": 3,
        "scf_parameters": "3.141592,0.0",
        "f_axion": 0.3,
        "m_axion": 300.0,
        "scf_evolve_as_fluid": "no",
        "scf_evolve_like_axionCAMB": "no",
        "do_shooting": "yes",
        "do_shooting_scf": "no",
        "scf_has_perturbations": "no",
        "attractor_ic_scf": "no",
        "modes": "s",
        "gauge": "synchronous",
        "lensing": "no",
    }

    cosmo = Class()
    cosmo.set(**params)
    cosmo.compute()
    cosmo.struct_cleanup()
    cosmo.empty()

    bg_files = sorted(ROOT.parent.glob(ROOT.name + "*background.dat"), key=lambda path: path.stat().st_mtime)
    if not bg_files:
        raise RuntimeError("CLASS finished but no background file was produced")

    rows = read_background(bg_files[-1])
    eq = nearest_row(rows, 3400.0)
    peak = max(rows, key=lambda row: row[18])

    print(f"background = {bg_files[-1]}")
    print(f"z_peak(fEDE) = {peak[0]:.6g}")
    print(f"fEDE_peak = {peak[18]:.6g}")
    print(f"fEDE(z_eq ~= 3400) = {eq[18]:.6g}")
    print(f"phi_prime/(aH) at z_eq = {phi_over_hubble(eq):.6g}")
    print(f"Veff_phi/(H^2 phi) at z_eq = {eq[29]:.6g}")


if __name__ == "__main__":
    main()
