# Report Handoff Notes

This file is a compact state dump for continuing the report work with another LLM.

## Current Goal

Build a PRD-style report in `report.tex` about interacting EDE-DM in this AxiCLASS fork, with:

- an introduction following the model logic of `2112.09128v2`
- a section on the `(1-cos)^3` potential
- a section on the tEDS potential
- plots showing the tEDS scan reproduction of Fig. 2 from `2212.08098v1`
- plots showing the `(1-cos)^3` result for the discussed `m, f, g` region
- plots comparing `Omega_idm_ede/Omega_tot` between tEDS and `(1-cos)^3`
- an appendix that explicitly documents the debugging history and code changes

## Main Report File

The report is now drafted in:

- [`report.tex`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/report.tex)

It is self-contained now. I removed the missing `\input{...}` dependency files and wrote the whole paper directly into `report.tex`.

## Content Already Added To `report.tex`

### Introduction

The intro is aligned with the logic of `2112.09128v2`:

- EDE defined through `f_EDE(z) = rho_phi / rho_tot`
- baseline axion-like EDE potential written as
  `V(phi) = m^2 f^2 [1 - cos(phi/f)]^3`
- type-1 coupled DM mass written as
  `m_IDM(phi) = m0 (1 + g phi^2)`
- coupled KG equation written with the force term
  `V_phi + 3 rho_IDM d ln m / d phi`

### `(1-cos)^3` Section

I included:

- model definition
- expected hilltop behavior
- explanation that coupling modifies the effective slope
- the diagnostic caveat that raw `rho_phi` at very high `z` can be misleading because the field can be position-frozen but still have a nonzero terminal `phi'`

### tEDS Section

I included:

- the tEDS potential
  `V(phi) = V0 theta^6 / (1 + theta^(6p))^(1/p)`
- the interpretation as a broad plateau plus steep transition
- the role of tEDS as a more flexible testbed for coupling-triggered release

### Figures Embedded In Report

The report now references these figures directly:

- tEDS Fig. 2-style reproduction:
  - [`analysis/teds_trigger_check/figure2_reproduction.png`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/analysis/teds_trigger_check/figure2_reproduction.png)
- `(1-cos)^3` natural-parameter scan:
  - [`analysis/axion_natural_trigger_report/selected_histories.png`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/analysis/axion_natural_trigger_report/selected_histories.png)
- direct `Omega_idm_ede` comparison between tEDS and `(1-cos)^3`:
  - [`analysis/report_assets/omega_idm_compare_teds_vs_cos3.png`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/analysis/report_assets/omega_idm_compare_teds_vs_cos3.png)
- earlier debugging scan plots:
  - [`analysis/type1_trigger_scan/plots/scan_scatter.png`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/analysis/type1_trigger_scan/plots/scan_scatter.png)
  - [`analysis/type1_trigger_scan/plots/top_histories.png`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/analysis/type1_trigger_scan/plots/top_histories.png)

## Code Changes Made

These are the key code changes already in the repo.

### 1) Coupled slow-roll IC for type-1 runs

I added a new option:

- `coupled_scf_slowroll_ic = yes/no`

Default behavior is now effectively `yes` for type-1 coupling, with explicit override support.

Relevant locations:

- [`include/background.h`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/include/background.h:220)
- [`source/input.c`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/source/input.c:3822)
- [`source/input.c`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/source/input.c:4914)
- [`source/input.c`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/source/input.c:7707)
- [`source/background.c`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/source/background.c:2966)

The IC uses the overdamped terminal-velocity branch:

```c
phi_prime_ini = -0.5 * a * (V_phi + 3 rho_idm dlnm_dphi) / H
```

This was added because starting from `phi_prime = 0` produced an artificial initial transient in `rho_scf`.

### 2) Type-1 DM continuity and KG force

The coupled background equations now use:

```text
d rho_idm / d ln a = -3 rho_idm + (phi' / aH) (d ln m / d phi) rho_idm
d phi' / d ln a   = -2 phi' - (a/H) [V_phi + 3 rho_idm d ln m / d phi]
```

Relevant lines:

- [`source/background.c`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/source/background.c:3382)
- [`source/background.c`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/source/background.c:3415)

### 3) Polynomial coupling derivatives

For `m(phi) = m0 (1 + g phi^2)`, the code uses:

```text
d ln m / d phi = 2 g phi / (1 + g phi^2)
d^2 ln m / d phi^2 = 2 g (1 - g phi^2) / (1 + g phi^2)^2
```

Relevant line:

- [`source/background.c`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/source/background.c:4186)

### 4) Type-1 mass-form parser

The parser accepts:

- `idm_ede_mass_form = poly`
- `idm_ede_mass_form = exp`
- `idm_ede_mass_form = phantom` / `andriot` / `am`

Relevant lines:

- [`source/input.c`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/source/input.c:3821)
- [`source/input.c`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/source/input.c:3831)
- [`source/input.c`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/source/input.c:3837)

## Debugging History To Preserve

The earlier `type1_trigger_scan` work is important because it established the main model-debugging story.

Main points from [`analysis/type1_trigger_scan/type1_trigger_scan.tex`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/analysis/type1_trigger_scan/type1_trigger_scan.tex):

- corrected the type-1 DM continuity equation
- corrected the type-1 KG force
- corrected the polynomial second derivative
- fixed perturbation-level coupling sourcing to use the actual background derivatives
- added the exponential coupling option
- found that polynomial coupling saturates and tends to fail as a clean equality trigger
- found that exponential coupling improves trigger control but still does not give a robust equality-era `f_EDE ~ 0.1` solution in the scanned region

The corresponding scan plots are already in:

- [`analysis/type1_trigger_scan/plots/scan_scatter.png`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/analysis/type1_trigger_scan/plots/scan_scatter.png)
- [`analysis/type1_trigger_scan/plots/top_histories.png`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/analysis/type1_trigger_scan/plots/top_histories.png)

## Relevant Numerical Results

### tEDS

The tEDS reproduction file already exists and should be used as the report figure source:

- [`analysis/teds_trigger_check/figure2_reproduction.png`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/analysis/teds_trigger_check/figure2_reproduction.png)

The tEDS benchmark used in the direct comparison figure is:

- `theta_i = 6`
- `g = 0.68`
- file: `analysis/teds_trigger_check/theta6_g0p68_coupled00_background.dat`

### `(1-cos)^3`

The natural-grid scan that produced the report figure is in:

- `analysis/axion_natural_trigger_report/`

Representative cases that look best for the report:

- `m/H0 = 1, 3, 10, 30, 100`
- `f = 1`
- `g = 1`

These are the histories shown in `selected_histories.png`.

### Direct `Omega_idm_ede` comparison

The comparison figure was generated from:

- tEDS `theta_i=6, g=0.68`
- `(1-cos)^3` `m/H0=30, f=1, g=1`
- `(1-cos)^3` `m/H0=30, f=1, g=0.3`

The image lives at:

- [`analysis/report_assets/omega_idm_compare_teds_vs_cos3.png`](/Users/vpoulin/Dropbox/Labo/ProgrammeCMB/AxiCLASS_EDE_DM_coupling/analysis/report_assets/omega_idm_compare_teds_vs_cos3.png)

## Current Status / Known Gaps

- `report.tex` is written, but I did not successfully compile it here because this environment lacks `revtex4-2.cls`.
- The appendix line references are now corrected, but if you make further code edits they may shift again.
- The report currently references the generated plots directly from `analysis/...`; those files are already present in the repo.

## Suggested Next Step

If the next LLM continues the report work, the most useful actions are:

1. Compile `report.tex` in your normal TeX environment.
2. Decide whether to add a short results subsection summarizing the numerical outcomes of the scans.
3. Optionally move the generated figures into a cleaner `figures/` or `paper_figures/` directory if the report is going to be finalized.

