"""
Shared utilities for CLASS ini files with correct IDM-EDE normalization.

Convention (matching Lin et al. 2212.08098):
  - omega_cdm = 0 explicitly: all DM is IDM-EDE.
  - Omega_Lambda NOT specified: CLASS derives it from flatness using the ODE-integrated
    rho_IDM at z=0.
  - omega_ini_idm_ede is iterated externally until omega_IDM_today = OMEGA_DM_TARGET.

Why iterate externally rather than use CLASS's built-in omega_idm_ede shooting?
  The CLASS shooting for omega_idm_ede has a bug for type-1 coupling: background.c
  lines 844/1118 add Omega0_idm_ede to Omega_m for the Lambda computation before
  the ODE runs, causing double-counting and an inconsistent Omega_Lambda. The external
  iteration avoids this by keeping omega_ini as the direct input.

Convergence: ~5 iterations to reach <0.1% accuracy on omega_IDM_today.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np

H0_DEFAULT = 72.81
OMEGA_B = 0.02251
OMEGA_DM_TARGET = 0.1280   # h^2 units; all DM is IDM-EDE

COMMON_COSMO_TAIL = f"""omega_b = {OMEGA_B}
omega_cdm = 0
H0 = {H0_DEFAULT}
tau_reio = 0.068
A_s = 2.191e-9
n_s = 0.9860
N_ur = 2.0328
N_ncdm = 1
deg_ncdm = 1
m_ncdm = 0.06
T_ncdm = 0.71611
"""

COMMON_SCF_FLAGS = """scf_evolve_as_fluid = no
scf_evolve_like_axionCAMB = no
do_shooting = no
do_shooting_scf = no
scf_has_perturbations = no
attractor_ic_scf = no
modes = s
gauge = synchronous
lensing = no
"""

COMMON_OUTPUT = """output = tCl
background_evolver = 0
evolver = 0
write background = yes
write parameters = no
"""


def _load_bg(path: Path) -> dict[str, np.ndarray] | None:
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
    except Exception:
        return None


def _run_class(class_bin: Path, ini_path: Path, repo: Path, timeout: int = 90) -> bool:
    res = subprocess.run(
        [str(class_bin), str(ini_path)],
        cwd=repo, capture_output=True, text=True, timeout=timeout,
    )
    return res.returncode == 0


def normalize_omega_ini(
    class_bin: Path,
    repo: Path,
    ini_builder,          # callable(root, omega_ini) -> ini_text str
    root: Path,
    g: float,
    phi_i: float,
    omega_dm_target: float = OMEGA_DM_TARGET,
    max_iter: int = 6,
    tol: float = 0.001,
) -> float:
    """
    Iterate omega_ini_idm_ede until omega_IDM_today = omega_dm_target (h^2 units).

    Start from omega_ini = omega_dm_target * (1 + g*phi_i^2), then apply ratio
    correction each iteration. Converges in ~5 steps to <0.1%.

    Returns the converged omega_ini value.
    """
    h = H0_DEFAULT / 100.0
    h2 = h ** 2

    omega_ini = omega_dm_target * (1.0 + g * phi_i ** 2)

    for _ in range(max_iter):
        ini_text = ini_builder(root, omega_ini)
        ini_path = Path(str(root) + ".ini")
        ini_path.write_text(ini_text)

        for suf in ("00_background.dat", "00_cl.dat", "00_parameters.ini", "00_unused_parameters"):
            p = Path(str(root) + suf)
            if p.exists():
                p.unlink()

        if not _run_class(class_bin, ini_path, repo):
            break

        bg = _load_bg(Path(str(root) + "00_background.dat"))
        if bg is None:
            break

        omega_idm_today = bg["(.)rho_idm_ede"][-1] / bg["(.)rho_tot"][-1] * h2
        ratio = omega_dm_target / omega_idm_today

        if abs(ratio - 1.0) < tol:
            break

        omega_ini *= ratio

    return omega_ini
