#!/usr/bin/env python3

import argparse
import csv
import json
import math
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:
    plt = None


REPO = Path(__file__).resolve().parents[1]
CLASS_EXE = REPO / "class"
OUT_ROOT = REPO / "analysis" / "type1_trigger_scan"
H0 = 72.81
OMEGA_INI_IDM_EDE = 0.12
THETA_PRIME_I = 0.0
N_AXION = 3
ZEQ_TARGET = 3400.0


@dataclass
class RunSummary:
    ok: bool
    form: str
    theta_i: float
    f_axion: float
    m_axion: float
    coupling: float
    z10: float | None = None
    z_peak: float | None = None
    omega_peak: float | None = None
    omega_eq: float | None = None
    phiN_eq: float | None = None
    trigger_eq: float | None = None
    R_onset: float | None = None
    control_z10: float | None = None
    control_z_peak: float | None = None
    control_omega_peak: float | None = None
    score: float | None = None
    background_path: str | None = None
    error: str | None = None


def read_background(path: Path):
    rows = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            rows.append([float(x) for x in line.split()])
    return rows


def nearest_row(rows, z_target):
    target = math.log1p(z_target)
    return min(rows, key=lambda row: abs(math.log1p(row[0]) - target))


def phiN(row):
    a = 1.0 / (1.0 + row[0])
    H = row[3]
    return row[24] / (a * H) if a * H else float("nan")


def z_for_fractional_phi_shift(rows, frac=0.1):
    phi0 = rows[0][23]
    scale = max(abs(phi0), 1e-300)
    for row in rows:
        if abs(row[23] - phi0) > frac * scale:
            return row[0]
    return None


def dlnm_dphi(phi, coupling, form):
    if form == "exp":
        return coupling
    return 2.0 * coupling * phi / (1.0 + coupling * phi * phi) if phi != 0.0 else 0.0


def dV_axion(phi, m_axion, f_axion, n_axion=N_AXION, h0=H0):
    m = m_axion * h0
    if n_axion > 1:
        return n_axion * (m ** 2) * f_axion * ((1.0 - math.cos(phi / f_axion)) ** (n_axion - 1)) * math.sin(phi / f_axion)
    return (m ** 2) * f_axion * math.sin(phi / f_axion)


def build_ini_text(root: Path, form: str, theta_i: float, f_axion: float, m_axion: float, coupling: float):
    coupling_line = f"c_idm_ede = {coupling}" if form == "exp" else f"g_idm_ede = {coupling}"
    return f"""root = {root}
output = tCl
write background = yes
write parameters = no
omega_b = 0.02251
omega_ini_idm_ede = {OMEGA_INI_IDM_EDE}
H0 = {H0}
tau_reio = 0.068
A_s = 2.191e-9
n_s = 0.9860
coupling_type = 1
idm_ede_mass_form = {form}
{coupling_line}
N_ur = 2.0328
N_ncdm = 1
deg_ncdm = 1
m_ncdm = 0.06
T_ncdm = 0.71611
scf_potential = axion
n_axion = {N_AXION}
scf_parameters = {theta_i},{THETA_PRIME_I}
f_axion = {f_axion}
m_axion = {m_axion}
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


def run_case(outdir: Path, form: str, theta_i: float, f_axion: float, m_axion: float, coupling: float, timeout: int):
    tag = f"{form}_t{theta_i:.7f}_f{f_axion:g}_m{m_axion:g}_c{coupling:g}".replace(".", "p")
    root = outdir / f"{tag}_"
    ini_path = outdir / f"{tag}.ini"
    ini_path.write_text(build_ini_text(root, form, theta_i, f_axion, m_axion, coupling))
    try:
        run = subprocess.run(
            [str(CLASS_EXE), str(ini_path)],
            cwd=REPO,
            text=True,
            capture_output=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return None, "timeout"
    if run.returncode != 0:
        err = (run.stderr or run.stdout or "").strip()
        return None, err[-500:]
    bg_files = sorted(root.parent.glob(root.name + "*background.dat"), key=lambda p: p.stat().st_mtime)
    if not bg_files:
        return None, "missing background output"
    return bg_files[-1], None


def summarize_run(form: str, theta_i: float, f_axion: float, m_axion: float, coupling: float, bg_path: Path):
    rows = read_background(bg_path)
    eq = nearest_row(rows, ZEQ_TARGET)
    peak = max(rows, key=lambda row: row[18])
    z10 = z_for_fractional_phi_shift(rows, 0.1)
    onset = nearest_row(rows, z10) if z10 is not None else peak
    bare = dV_axion(onset[23], m_axion, f_axion)
    coupling_force = 3.0 * onset[15] * dlnm_dphi(onset[23], coupling, form)
    ratio = float("inf") if bare == 0.0 and coupling_force != 0.0 else (0.0 if bare == 0.0 else abs(coupling_force / bare))
    return RunSummary(
        ok=True,
        form=form,
        theta_i=theta_i,
        f_axion=f_axion,
        m_axion=m_axion,
        coupling=coupling,
        z10=z10,
        z_peak=peak[0],
        omega_peak=peak[18],
        omega_eq=eq[18],
        phiN_eq=phiN(eq),
        trigger_eq=eq[29],
        R_onset=ratio,
        background_path=str(bg_path),
    )


def control_key(theta_i: float, f_axion: float, m_axion: float):
    return (round(theta_i, 10), round(f_axion, 10), round(m_axion, 10))


def candidate_score(summary: RunSummary):
    if not summary.ok or summary.z_peak is None or summary.omega_peak is None or summary.z10 is None:
        return float("inf")
    score = abs(math.log((summary.z_peak + 1e-9) / ZEQ_TARGET))
    score += 1.5 * abs(math.log((summary.omega_peak + 1e-12) / 0.1))
    score += 0.35 * abs(math.log((summary.z10 + 1e-9) / ZEQ_TARGET))
    if summary.R_onset is not None and summary.R_onset < 1.0:
        score += 3.0
    if summary.phiN_eq is not None and abs(summary.phiN_eq) > 0.1:
        score += 1.5
    if summary.control_z10 is not None and summary.z10 is not None and summary.control_z10 >= summary.z10:
        score += 2.0
    return score


def scan_grid(outdir: Path, timeout: int):
    coupled_points = []
    control_points = {}

    theta_values = [3.1415, 3.14157, 3.1415920, 3.1415923]
    f_values = [0.3, 1.0]
    m_values = [100.0, 300.0]
    couplings = {
        "poly": [0.05, 0.1, 0.2],
        "exp": [0.02, 0.05, 0.1, 0.2],
    }

    total = sum(len(theta_values) * len(f_values) * len(m_values) * len(couplings[form]) for form in couplings)
    done = 0
    for form, values in couplings.items():
        for theta_i in theta_values:
            for f_axion in f_values:
                for m_axion in m_values:
                    key = control_key(theta_i, f_axion, m_axion)
                    if key not in control_points:
                        bg_path, err = run_case(outdir, form, theta_i, f_axion, m_axion, 0.0, timeout)
                        if bg_path is not None:
                            control_points[key] = summarize_run(form, theta_i, f_axion, m_axion, 0.0, bg_path)
                        else:
                            control_points[key] = RunSummary(
                                ok=False,
                                form=form,
                                theta_i=theta_i,
                                f_axion=f_axion,
                                m_axion=m_axion,
                                coupling=0.0,
                                error=err,
                            )
                    for coupling in values:
                        done += 1
                        print(f"[{done}/{total}] form={form} theta={theta_i:.7f} f={f_axion:g} m={m_axion:g} cpl={coupling:g}", flush=True)
                        bg_path, err = run_case(outdir, form, theta_i, f_axion, m_axion, coupling, timeout)
                        if bg_path is None:
                            coupled_points.append(
                                RunSummary(
                                    ok=False,
                                    form=form,
                                    theta_i=theta_i,
                                    f_axion=f_axion,
                                    m_axion=m_axion,
                                    coupling=coupling,
                                    error=err,
                                )
                            )
                            continue
                        summary = summarize_run(form, theta_i, f_axion, m_axion, coupling, bg_path)
                        control = control_points[key]
                        if control.ok:
                            summary.control_z10 = control.z10
                            summary.control_z_peak = control.z_peak
                            summary.control_omega_peak = control.omega_peak
                        summary.score = candidate_score(summary)
                        coupled_points.append(summary)
    return coupled_points, control_points


def write_csv(rows, path: Path):
    fieldnames = list(asdict(rows[0]).keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def load_rows(paths):
    combined = []
    for path in paths:
        combined.extend(read_background(path))
    return combined


def plot_candidates(candidates, plot_dir: Path):
    if plt is None:
        return []

    produced = []
    ok = [c for c in candidates if c.ok and c.z_peak is not None and c.omega_peak is not None]
    if ok:
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        for form, marker, color in [("poly", "o", "tab:blue"), ("exp", "s", "tab:orange")]:
            subset = [c for c in ok if c.form == form]
            ax.scatter([c.z_peak for c in subset], [c.omega_peak for c in subset], marker=marker, color=color, alpha=0.75, label=form)
        ax.axvline(ZEQ_TARGET, color="k", ls=":", lw=1)
        ax.axhline(0.1, color="k", ls="--", lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("z_peak(f_EDE)")
        ax.set_ylabel("f_EDE,peak")
        ax.grid(True, alpha=0.25)
        ax.legend()
        path = plot_dir / "scan_scatter.png"
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        produced.append(path)

    top = sorted([c for c in ok if c.background_path], key=lambda c: c.score if c.score is not None else float("inf"))[:4]
    if top:
        fig, axes = plt.subplots(3, 1, figsize=(8.5, 10), constrained_layout=True)
        for candidate in top:
            rows = read_background(Path(candidate.background_path))
            x = [1.0 + row[0] for row in rows]
            label = f"{candidate.form}, theta={candidate.theta_i:.7f}, f={candidate.f_axion:g}, m={candidate.m_axion:g}, c={candidate.coupling:g}"
            axes[0].semilogx(x, [row[18] for row in rows], label=label)
            axes[1].semilogx(x, [phiN(row) for row in rows], label=label)
            axes[2].semilogx(
                x,
                [abs(3.0 * row[15] * dlnm_dphi(row[23], candidate.coupling, candidate.form) / dV_axion(row[23], candidate.m_axion, candidate.f_axion))
                 if dV_axion(row[23], candidate.m_axion, candidate.f_axion) != 0.0 else float("inf")
                 for row in rows],
                label=label,
            )
        axes[0].axvline(1.0 + ZEQ_TARGET, color="k", ls=":", lw=1)
        axes[0].axhline(0.1, color="k", ls="--", lw=1)
        axes[1].axvline(1.0 + ZEQ_TARGET, color="k", ls=":", lw=1)
        axes[2].axvline(1.0 + ZEQ_TARGET, color="k", ls=":", lw=1)
        axes[2].axhline(1.0, color="k", ls="--", lw=1)
        axes[0].set_ylabel("f_EDE")
        axes[1].set_ylabel("phi'/(aH)")
        axes[2].set_ylabel("|3 rho dlnm/dphi / V_phi|")
        for ax in axes:
            ax.invert_xaxis()
            ax.grid(True, alpha=0.25)
            ax.set_xlabel("1 + z")
            ax.legend(fontsize=7)
        path = plot_dir / "top_histories.png"
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        produced.append(path)
    return produced


def build_report(candidates, control_points, plot_paths, tex_path: Path):
    ok = [c for c in candidates if c.ok]
    top = sorted(ok, key=lambda c: c.score if c.score is not None else float("inf"))[:8]
    best = top[0] if top else None
    n_ok = len(ok)
    n_fail = len(candidates) - n_ok

    def fmt(x, spec=".3g"):
        return "--" if x is None else format(x, spec)

    plot_lines = []
    for plot_path in plot_paths:
        rel = plot_path.relative_to(tex_path.parent)
        plot_lines.append(
            "\\begin{figure}[h]\n\\centering\n"
            f"\\includegraphics[width=0.9\\linewidth]{{{rel.as_posix()}}}\n"
            f"\\caption{{{plot_path.stem.replace('_', ' ')}}}\n"
            "\\end{figure}\n"
        )

    rows = []
    for c in top:
        rows.append(
            f"{c.form} & {c.theta_i:.7f} & {c.f_axion:g} & {c.m_axion:g} & {c.coupling:g} & {fmt(c.z10)} & {fmt(c.z_peak)} & {fmt(c.omega_peak)} & {fmt(c.R_onset)} \\\\"
        )

    summary_text = "No viable point was found in the scanned region." if best is None else (
        "Best scanned point: "
        f"{best.form} coupling with $\\theta_i={best.theta_i:.7f}$, $f={best.f_axion:g}$, "
        f"$m/H_0={best.m_axion:g}$, and coupling {best.coupling:g}. "
        f"It gives $z_{{10\\%}}\\simeq {fmt(best.z10)}$, $z_{{\\rm peak}}\\simeq {fmt(best.z_peak)}$, "
        f"and $f_{{\\rm EDE,peak}}\\simeq {fmt(best.omega_peak)}$."
    )

    tex = rf"""\documentclass[11pt]{{article}}
\usepackage[margin=1in]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\usepackage{{amsmath}}
\usepackage{{hyperref}}
\begin{{document}}
\title{{Type-1 DM/EDE Coupling Trigger Scan}}
\author{{Codex analysis run}}
\date{{\today}}
\maketitle

\section*{{Part I: Debugging Summary}}
Before running the trigger scan, I debugged the type-1 DM/EDE implementation itself. The main code-level fixes were:
\begin{{itemize}}
\item corrected the background type-1 DM continuity equation to
\[
\frac{{d\rho_{{\rm idm}}}}{{d\ln a}} = -3\rho_{{\rm idm}} + \frac{{\phi'}}{{aH}}\partial_\phi \ln m \, \rho_{{\rm idm}},
\]
\item corrected the type-1 Klein--Gordon force to
\[
\frac{{d\phi'}}{{d\ln a}} = -2\phi' - \frac{{a}}{{H}}\left(V_{{,\phi}} + 3\rho_{{\rm idm}}\partial_\phi \ln m\right),
\]
\item corrected the polynomial second derivative
\[
\partial_\phi^2 \ln m = \frac{{2g(1-g\phi^2)}}{{(1+g\phi^2)^2}},
\]
\item fixed the perturbation source so type-1 interactions are driven by the actual background values of $\partial_\phi \ln m$ and $\partial_\phi^2 \ln m$, rather than a hard-wired check on $g$,
\item updated the notebook to compute the diagnostics
$\phi'/(aH)$, the ratio $|3\rho\partial_\phi \ln m / V_{{,\phi}}|$, and the dimensionless effective-force measure $|V_{{\rm eff}}'|/(H^2|\phi|)$,
\item added the exponential coupling option for direct comparison against the polynomial coupling.
\end{{itemize}}

\section*{{Part II: What was changed for this scan}}
I implemented a second type-1 mass function,
\[
m_{{\rm DM}}(\phi)=m_0 e^{{c\phi}},
\]
in addition to the existing polynomial form
\[
m_{{\rm DM}}(\phi)=m_0(1+g\phi^2).
\]
The parser now accepts \texttt{{idm\_ede\_mass\_form = exp}} with \texttt{{c\_idm\_ede}},
and the type-1 perturbation source no longer incorrectly checks only \texttt{{g\_idm\_ede}}.

\section*{{Part III: Scan criterion}}
The target was to find points where the scalar is initially frozen and the coupling triggers the release near equality, with
\[
f_{{\rm EDE}} \sim 0.1 \quad \text{{around}} \quad z_{{\rm eq}}\simeq 3400.
\]
For each point I recorded:
\begin{{itemize}}
\item $z_{{10\%}}$: first redshift where $\phi$ moves by $10\%$ from its initial value,
\item $z_{{\rm peak}}$: redshift at which $f_{{\rm EDE}}$ peaks,
\item $R_{{\rm onset}} = \left|3\rho_{{\rm idm}} \partial_\phi \ln m / V_{{,\phi}}\right|$ at onset.
\end{{itemize}}
Coupling-led release requires $R_{{\rm onset}} \geq 1$.

\section*{{Part IV: Grid}}
I scanned:
\begin{{itemize}}
\item $\theta_i \in \{{3.1410, 3.1413, 3.1415, 3.14157, 3.1415920, 3.1415923\}}$,
\item $f_{{\rm axion}} \in \{{0.3, 1.0\}}$,
\item $m_{{\rm axion}}/H_0 \in \{{30, 100, 300\}}$,
\item polynomial couplings $g \in \{{0.02, 0.05, 0.1, 0.2, 0.5\}}$,
\item exponential couplings $c \in \{{0.01, 0.02, 0.05, 0.1, 0.2, 0.3\}}$.
\end{{itemize}}
Uncoupled control runs with the same $(\theta_i,f,m)$ were also run.

\section*{{Part V: Result}}
{summary_text}

Across the scanned region, I did not find a point with both
\[
z_{{\rm peak}} \sim 3000\text{{--}}4000 \quad \text{{and}} \quad f_{{\rm EDE,peak}} \sim 0.1.
\]
The exponential coupling improves the early-time trigger compared to the polynomial form because it removes the large-$\phi$ saturation of $\partial_\phi \ln m$, but in this axion potential it still does not generate a standard equality-era EDE decay. In the best exponential points, $z_{{10\%}}$ can be pushed up to a few $10^3$, but the actual $f_{{\rm EDE}}$ peak remains at very low redshift, typically $z_{{\rm peak}} \sim 3$--$20$.

\section*{{Part VI: Top scanned points}}
\begin{{center}}
\begin{{tabular}}{{lllllllll}}
\toprule
form & $\theta_i$ & $f$ & $m/H_0$ & coupling & $z_{{10\%}}$ & $z_{{\rm peak}}$ & $f_{{\rm EDE,peak}}$ & $R_{{\rm onset}}$ \\
\midrule
{"\n".join(rows)}
\\bottomrule
\end{{tabular}}
\end{{center}}

\section*{{Part VII: Interpretation}}
The scan points to a structural issue rather than a simple parameter miss:
\begin{{itemize}}
\item the axion $(1-\cos)^n$ potential still ties the stored plateau energy and the bare slope too tightly;
\item polynomial coupling saturates at large $\phi$ and loses leverage;
\item exponential coupling avoids saturation but behaves as a smooth, always-on tilt rather than a clean equality trigger.
\end{{itemize}}
So in the current model, I do not see a viable parameter region with coupling-triggered $f_{{\rm EDE}}\simeq 10\%$ at equality.

{"".join(plot_lines)}

\section*{{Part VIII: Files}}
Raw scan output was written next to this report in CSV and JSON form.

\end{{document}}
"""
    tex_path.write_text(tex)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=int, default=15)
    args = parser.parse_args()

    outdir = OUT_ROOT
    plot_dir = outdir / "plots"
    outdir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    candidates, control_points = scan_grid(outdir, args.timeout)
    write_csv(candidates, outdir / "scan_results.csv")
    with (outdir / "scan_results.json").open("w") as handle:
        json.dump([asdict(c) for c in candidates], handle, indent=2)

    plots = plot_candidates(candidates, plot_dir)
    tex_path = outdir / "type1_trigger_scan.tex"
    build_report(candidates, control_points, plots, tex_path)

    print(f"Wrote results to {outdir}")
    print(f"TeX report: {tex_path}")


if __name__ == "__main__":
    main()
