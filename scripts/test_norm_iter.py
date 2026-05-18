#!/usr/bin/env python3
"""Test iterative omega_ini normalization with omega_cdm=0 explicitly."""
import subprocess, numpy as np
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CLASS = REPO / 'class'
H0 = 72.81; h = H0/100; h2 = h**2; OMEGA_DM = 0.1280

BASE_INI = """output = tCl
background_evolver = 0
evolver = 0
write background = yes
write parameters = no
omega_b = 0.02251
omega_cdm = 0
omega_ini_idm_ede = {omega_ini}
H0 = 72.81
tau_reio = 0.068
A_s = 2.191e-9
n_s = 0.9860
coupling_type = 1
idm_ede_mass_form = poly
g_idm_ede = 0.68
N_ur = 2.0328
N_ncdm = 1
deg_ncdm = 1
m_ncdm = 0.06
T_ncdm = 0.71611
scf_potential = teds
scf_parameters = 0.3,0.0
V0_teds = 164.17659980857525
f_teds = 0.05
p_teds = 8.0
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

def run(omega_ini, tag):
    ini = REPO / 'analysis' / f'{tag}.ini'
    ini.parent.mkdir(parents=True, exist_ok=True)
    ini.write_text(BASE_INI.format(omega_ini=omega_ini) + f'root = {REPO}/analysis/{tag}\n')
    for s in ('00_background.dat','00_cl.dat','00_parameters.ini','00_unused_parameters'):
        p = REPO/'analysis'/f'{tag}{s}'
        if p.exists(): p.unlink()
    r = subprocess.run([str(CLASS), str(ini)], cwd=REPO,
                       capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        print('FAILED:', r.stderr[-300:])
        return None
    bg_path = REPO / 'analysis' / f'{tag}00_background.dat'
    header = ''
    with open(bg_path) as f:
        for line in f:
            if line.startswith('#') and '1:z' in line: header = line[1:]; break
    names = [e.split(':',1)[1] for e in header.split() if ':' in e]
    data = np.loadtxt(bg_path, comments='#')
    bg = {n: data[:,i] for i,n in enumerate(names)}
    z = bg['z']; fe = bg['(.)Omega_scf']; pk = int(np.argmax(fe))
    tot = bg['(.)rho_tot'][-1]
    w_idm = bg['(.)rho_idm_ede'][-1] / tot * h2
    OmLam = bg['(.)rho_lambda'][-1] / tot
    OmCDM = bg['(.)rho_cdm'][-1] / tot if '(.)rho_cdm' in bg else 0.0
    return {'w_idm': w_idm, 'OmLam': OmLam, 'OmCDM': OmCDM,
            'z_peak': z[pk], 'fede': fe[pk]}

omega_ini = OMEGA_DM * (1 + 0.68 * 0.3**2)
for i in range(4):
    r = run(omega_ini, f'norm_test_iter{i}')
    if r is None: break
    print(f'iter{i}: omega_ini={omega_ini:.6f}  '
          f'w_IDM={r["w_idm"]:.5f} (target {OMEGA_DM})  '
          f'z_peak={r["z_peak"]:.0f}  '
          f'OmLam={r["OmLam"]:.4f}  OmCDM={r["OmCDM"]:.4f}')
    if abs(r['w_idm'] / OMEGA_DM - 1) < 0.001:
        print('CONVERGED')
        break
    omega_ini *= OMEGA_DM / r['w_idm']
