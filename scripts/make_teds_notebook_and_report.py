#!/usr/bin/env python3
"""Create the tEDS trigger notebook and PDF report.

The notebook is a runnable analysis notebook.  The PDF report is generated from
the current CLASS background outputs in analysis/teds_trigger_check and
analysis/potential_trigger_compare.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import nbformat as nbf
import numpy as np


REPO = Path(__file__).resolve().parents[1]
TEDS_DIR = REPO / "analysis" / "teds_trigger_check"
COMPARE_DIR = REPO / "analysis" / "potential_trigger_compare"
NOTEBOOK = REPO / "notebooks" / "tEDS_trigger_reproduction.ipynb"
REPORT = TEDS_DIR / "tEDS_trigger_report.pdf"
FIG_DIR = TEDS_DIR / "plots"
ZEQ = 3400.0


def read_summary(path: Path) -> list[dict[str, float | str]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    output: list[dict[str, float | str]] = []
    for row in rows:
        clean: dict[str, float | str] = {}
        for key, value in row.items():
            if value is None or value == "":
                clean[key] = value or ""
                continue
            try:
                clean[key] = float(value)
            except ValueError:
                clean[key] = value
        output.append(clean)
    return output


def load_background(path: Path) -> dict[str, np.ndarray]:
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
    data = np.loadtxt(path)
    return {name: data[:, index] for index, name in enumerate(names)}


def teds_backgrounds() -> list[tuple[dict[str, float | str], dict[str, np.ndarray]]]:
    rows = read_summary(TEDS_DIR / "summary.csv")
    out = []
    for row in rows:
        theta = f"{row['theta_i']:g}".replace(".", "p")
        g = f"{row['g']:g}".replace(".", "p")
        path = TEDS_DIR / f"theta{theta}_g{g}_coupled00_background.dat"
        out.append((row, load_background(path)))
    return out


def plot_potential(ax: plt.Axes) -> None:
    theta = np.logspace(-2, 1.2, 800)
    for p, color in [(2, "#6c757d"), (4, "#2a9d8f"), (8, "#264653")]:
        v = theta**6 / (1.0 + theta ** (6 * p)) ** (1.0 / p)
        ax.loglog(theta, v, color=color, lw=2, label=f"p={p}")
    ax.axvspan(4, 8, color="#e9c46a", alpha=0.18, label="paper IC range")
    ax.set_xlabel(r"$\theta=\phi/f$")
    ax.set_ylabel(r"$V/V_0$")
    ax.set_title("tEDS plateau power-law potential")
    ax.legend(frameon=False)
    ax.grid(True, which="both", alpha=0.25)


def make_figures() -> dict[str, Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    data = teds_backgrounds()
    colors = {8.0: "#264653", 6.0: "#2a9d8f", 4.0: "#e76f51"}
    paths: dict[str, Path] = {}

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    plot_potential(ax)
    paths["potential"] = FIG_DIR / "teds_potential.pdf"
    fig.tight_layout()
    fig.savefig(paths["potential"])
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.6), sharex=True)
    for row, bg in data:
        theta_i = float(row["theta_i"])
        color = colors[theta_i]
        z = bg["z"]
        x = 1.0 + z
        label = rf"$\theta_i={theta_i:g}$, $g={float(row['g']):.2f}$"
        axes[0, 0].loglog(x, bg["(.)Omega_scf"], color=color, lw=2.2, label=label)
        axes[0, 1].semilogx(x, bg["phi_scf"] / 0.05, color=color, lw=2.2)
        bare = np.abs(bg["V'_scf"])
        coupling = np.abs(bg.get("dVadd_contribution", np.zeros_like(bare)))
        axes[1, 0].loglog(x, coupling / np.maximum(bare, 1e-300), color=color, lw=2.2)
        mrel = 1.0 + float(row["g"]) * bg["phi_scf"] ** 2
        axes[1, 1].semilogx(x, mrel / mrel[0] - 1.0, color=color, lw=2.2)

    for ax in axes.flat:
        ax.axvline(1 + ZEQ, color="0.35", ls="--", lw=1)
        ax.invert_xaxis()
        ax.grid(True, which="both", alpha=0.25)
    axes[0, 0].set_ylabel(r"$f_{\rm EDE}=\Omega_\phi$")
    axes[0, 1].set_ylabel(r"$\theta=\phi/f$")
    axes[1, 0].set_ylabel(r"$|3\rho_{\rm DM}\partial_\phi\ln m|/|V'|$")
    axes[1, 1].set_ylabel(r"$m_{\rm DM}/m_{\rm DM,ini}-1$")
    axes[1, 0].set_xlabel(r"$1+z$")
    axes[1, 1].set_xlabel(r"$1+z$")
    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[0, 0].set_title("EDE fraction peaks near equality")
    axes[0, 1].set_title("Field release from the plateau")
    axes[1, 0].set_title("Coupling dominates the release force")
    axes[1, 1].set_title("Small DM mass variation")
    paths["histories"] = FIG_DIR / "teds_histories.pdf"
    fig.tight_layout()
    fig.savefig(paths["histories"])
    plt.close(fig)

    compare_rows = [row for row in read_summary(COMPARE_DIR / "summary.csv") if row.get("ok") == "yes"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.9))
    markers = {("axion", 2.0): ("o", "#457b9d", "axion n=2"), ("axion", 3.0): ("s", "#1d3557", "axion n=3"), ("teds", 3.0): ("^", "#e76f51", "tEDS")}
    for key, (marker, color, label) in markers.items():
        model, n = key
        subset = [row for row in compare_rows if row["model"] == model and float(row["n"]) == n]
        if not subset:
            continue
        axes[0].scatter(
            [row["log10_z_peak"] for row in subset],
            [row["f_ede_peak"] for row in subset],
            marker=marker,
            color=color,
            s=54,
            label=label,
        )
        axes[1].scatter(
            [row["dm_mass_shift_max"] for row in subset],
            [row["f_ede_peak"] for row in subset],
            marker=marker,
            color=color,
            s=54,
            label=label,
        )
    axes[0].axvline(math.log10(ZEQ), color="0.35", ls="--", lw=1)
    axes[0].axhspan(0.08, 0.14, color="#2a9d8f", alpha=0.12)
    axes[0].set_xlabel(r"$\log_{10} z_{\rm peak}$")
    axes[0].set_ylabel(r"$f_{\rm EDE,peak}$")
    axes[0].set_title("Trigger target: equality and 10% EDE")
    axes[1].set_xlabel("max fractional DM mass shift")
    axes[1].set_ylabel(r"$f_{\rm EDE,peak}$")
    axes[1].set_title("Energy release vs DM mass change")
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend(frameon=False, fontsize=9)
    paths["comparison"] = FIG_DIR / "potential_comparison.pdf"
    fig.tight_layout()
    fig.savefig(paths["comparison"])
    plt.close(fig)

    examples = [
        ("standard EDE, no coupling", "axion", 3.0, 0.3, 100.0, 0.0, "#6c757d", "-"),
        ("standard EDE + coupling", "axion", 3.0, 0.3, 100.0, 1.0, "#1d3557", "--"),
        ("standard EDE + coupling", "axion", 3.0, 1.0, 100.0, 0.5, "#457b9d", "-."),
        ("tEDS + coupling", "teds", 3.0, 0.05, 0.0, 0.68, "#e76f51", "-"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for label, model, n, f, m, g, color, ls in examples:
        subset = [
            row for row in compare_rows
            if row["model"] == model and float(row["n"]) == n and float(row["f"]) == f
            and float(row["m"]) == m and float(row["g"]) == g and row.get("background")
        ]
        if not subset:
            continue
        bg = load_background(Path(str(subset[0]["background"])))
        ax.loglog(1.0 + bg["z"], bg["(.)Omega_scf"], color=color, ls=ls, lw=2.2, label=label)
    ax.axvline(1 + ZEQ, color="0.35", ls=":", lw=1.3, label="equality")
    ax.axhspan(0.08, 0.14, color="#2a9d8f", alpha=0.12, label="10% target band")
    ax.invert_xaxis()
    ax.set_xlabel(r"$1+z$")
    ax.set_ylabel(r"$f_{\rm EDE}=\Omega_\phi$")
    ax.set_title("Standard EDE vs coupling-triggered tEDS")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    paths["standard_comparison"] = FIG_DIR / "standard_vs_teds_histories.pdf"
    fig.tight_layout()
    fig.savefig(paths["standard_comparison"])
    plt.close(fig)

    return paths


def write_report() -> None:
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    figures = make_figures()
    summary = read_summary(TEDS_DIR / "summary.csv")
    compare_rows = read_summary(COMPARE_DIR / "summary.csv")

    with PdfPages(REPORT) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.08, 0.93, "tEDS Coupling-Triggered EDE Release", fontsize=20, weight="bold")
        fig.text(0.08, 0.895, "Background reproduction of the plateau power-law trigger in arXiv:2212.08098", fontsize=11)
        body = [
            r"tEDS potential: $V=V_0\theta^6/(1+\theta^{6p})^{1/p}$, $\theta=\phi/f$.",
            r"Standard axion-like EDE potential: $V=m^2f^2[1-\cos(\phi/f)]^n$.",
            r"DM mass coupling: $m_{\rm DM}(\phi)=m_0(1+g\phi^2)$, so $\partial_\phi\ln m=2g\phi/(1+g\phi^2)$.",
            r"Background KG force term in this implementation: $V' + 3\rho_{\rm DM}\partial_\phi\ln m$.",
            r"Fiducial tEDS values: $p=8$, $f=0.05\,M_{\rm Pl}$, $V_0=0.12\,{\rm eV}^4$ converted to CLASS units.",
        ]
        fig.text(0.08, 0.84, "\n".join(body), fontsize=11, va="top", linespacing=1.5)

        y = 0.68
        fig.text(0.08, y, "Numerical summary", fontsize=14, weight="bold")
        y -= 0.035
        header = "theta_i    g       log10 z_peak    f_EDE_peak    max DM mass shift"
        fig.text(0.08, y, header, family="monospace", fontsize=10)
        for row in summary:
            y -= 0.027
            line = f"{row['theta_i']:>6.1f}  {row['g']:>5.2f}      {row['log10_z_peak']:>8.3f}      {row['f_ede_peak']:>8.3f}        {row['dm_mass_shift_max']:>8.3f}"
            fig.text(0.08, y, line, family="monospace", fontsize=10)
        fig.text(0.08, 0.36, "Conclusion", fontsize=14, weight="bold")
        conclusion = (
            "Compared with the axion-like hilltop potential, the tEDS power-law plateau is structurally better "
            "for coupling-triggered release.  The plateau stores enough energy while keeping the bare slope tiny; "
            "the coupling then supplies the release force near equality.  In the axion-like case, tuning the "
            "hilltop slope, stored energy, and release redshift remains entangled."
        )
        fig.text(0.08, 0.33, conclusion, fontsize=11, va="top", wrap=True)
        pdf.savefig(fig)
        plt.close(fig)

        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.08, 0.93, "Model Details and Standard-EDE Comparison", fontsize=18, weight="bold")
        details = (
            "Standard EDE without an extra coupling releases when Hubble friction becomes small compared with the "
            "bare curvature/slope of the axion-like potential.  The same potential with a DM coupling changes the "
            "force balance, but not the underlying link between plateau energy, hilltop tuning, and bare slope.  "
            "That is why the axion-like cases split into two bad regimes: equality-era timing with too small "
            "f_EDE, or order-10% f_EDE with a peak far below equality.\n\n"
            "The tEDS potential was designed to break that link.  For large theta it has a high, flat plateau, "
            "while near the minimum it behaves as theta^6 and redshifts with the same averaged w=1/2 as an n=3 "
            "axion-like potential.  On the plateau, the bare slope is extremely suppressed, so the DM coupling "
            "can be the physical trigger rather than a small perturbation to ordinary Hubble release."
        )
        fig.text(0.08, 0.87, details, fontsize=10.5, va="top", wrap=True, linespacing=1.45)
        y = 0.58
        fig.text(0.08, y, "Representative outcomes", fontsize=13, weight="bold")
        y -= 0.035
        fig.text(0.08, y, "case                              log10 z_peak   f_EDE_peak   interpretation", family="monospace", fontsize=8.7)

        representative = [
            ("axion n=3, f=0.3, g=0", lambda r: r.get("ok") == "yes" and r["model"] == "axion" and float(r["n"]) == 3.0 and float(r["f"]) == 0.3 and float(r["g"]) == 0.0),
            ("axion n=3, f=0.3, g=1", lambda r: r.get("ok") == "yes" and r["model"] == "axion" and float(r["n"]) == 3.0 and float(r["f"]) == 0.3 and float(r["g"]) == 1.0),
            ("axion n=3, f=1.0, g=0.5", lambda r: r.get("ok") == "yes" and r["model"] == "axion" and float(r["n"]) == 3.0 and float(r["f"]) == 1.0 and float(r["g"]) == 0.5),
            ("tEDS theta=6, g=0.68", lambda r: r.get("ok") == "yes" and r["model"] == "teds" and float(r["theta_i"]) == 6.0 and float(r["g"]) == 0.68),
        ]
        interpretations = [
            "too late and too large",
            "near equality but too small",
            "right amplitude but too late",
            "near equality and ~10%",
        ]
        for (name, pred), interp in zip(representative, interpretations):
            match = next((row for row in compare_rows if pred(row)), None)
            y -= 0.027
            if match:
                line = f"{name:<33} {match['log10_z_peak']:>8.3f}      {match['f_ede_peak']:>8.3f}   {interp}"
            else:
                line = f"{name:<33} {'failed':>8}      {'--':>8}   no viable background for same ICs"
            fig.text(0.08, y, line, family="monospace", fontsize=8.7)

        failed_teds = [row for row in compare_rows if row.get("model") == "teds" and float(row.get("g", -1)) == 0.0 and row.get("ok") != "yes"]
        if failed_teds:
            y -= 0.055
            fig.text(0.08, y, "Uncoupled tEDS controls", fontsize=12, weight="bold")
            y -= 0.03
            fig.text(0.08, y, "With g=0 the tEDS plateau cases in this setup do not yield the paper-like release; they fail the background viability checks before producing a controlled equality-era peak.", fontsize=10, wrap=True)
        pdf.savefig(fig)
        plt.close(fig)

        for key in ("potential", "histories", "comparison"):
            fig = plt.figure(figsize=(8.27, 11.69))
            image = plt.imread(str(figures[key]).replace(".pdf", ".png")) if False else None
            plt.close(fig)

        # Recreate the figures as native PDF pages instead of embedding rasterized PDFs.
        fig, ax = plt.subplots(figsize=(8.0, 5.2))
        plot_potential(ax)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        data = teds_backgrounds()
        colors = {8.0: "#264653", 6.0: "#2a9d8f", 4.0: "#e76f51"}
        fig, axes = plt.subplots(2, 2, figsize=(8.0, 7.2), sharex=True)
        for row, bg in data:
            theta_i = float(row["theta_i"])
            color = colors[theta_i]
            z = bg["z"]
            x = 1.0 + z
            axes[0, 0].loglog(x, bg["(.)Omega_scf"], color=color, lw=2.0, label=rf"$\theta_i={theta_i:g}$")
            axes[0, 1].semilogx(x, bg["phi_scf"] / 0.05, color=color, lw=2.0)
            bare = np.abs(bg["V'_scf"])
            coupling = np.abs(bg.get("dVadd_contribution", np.zeros_like(bare)))
            axes[1, 0].loglog(x, coupling / np.maximum(bare, 1e-300), color=color, lw=2.0)
            mrel = 1.0 + float(row["g"]) * bg["phi_scf"] ** 2
            axes[1, 1].semilogx(x, mrel / mrel[0] - 1.0, color=color, lw=2.0)
        for ax in axes.flat:
            ax.axvline(1 + ZEQ, color="0.35", ls="--", lw=1)
            ax.invert_xaxis()
            ax.grid(True, which="both", alpha=0.25)
        axes[0, 0].set_ylabel(r"$\Omega_\phi$")
        axes[0, 1].set_ylabel(r"$\phi/f$")
        axes[1, 0].set_ylabel(r"coupling force / bare force")
        axes[1, 1].set_ylabel(r"$m_{\rm DM}/m_{\rm DM,ini}-1$")
        axes[1, 0].set_xlabel(r"$1+z$")
        axes[1, 1].set_xlabel(r"$1+z$")
        axes[0, 0].legend(frameon=False)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        compare_rows = [row for row in read_summary(COMPARE_DIR / "summary.csv") if row.get("ok") == "yes"]
        fig, axes = plt.subplots(1, 2, figsize=(8.0, 4.0))
        markers = {("axion", 2.0): ("o", "#457b9d", "axion n=2"), ("axion", 3.0): ("s", "#1d3557", "axion n=3"), ("teds", 3.0): ("^", "#e76f51", "tEDS")}
        for key, (marker, color, label) in markers.items():
            model, n = key
            subset = [row for row in compare_rows if row["model"] == model and float(row["n"]) == n]
            axes[0].scatter([row["log10_z_peak"] for row in subset], [row["f_ede_peak"] for row in subset], marker=marker, color=color, s=52, label=label)
            axes[1].scatter([row["dm_mass_shift_max"] for row in subset], [row["f_ede_peak"] for row in subset], marker=marker, color=color, s=52, label=label)
        axes[0].axvline(math.log10(ZEQ), color="0.35", ls="--", lw=1)
        axes[0].axhspan(0.08, 0.14, color="#2a9d8f", alpha=0.12)
        axes[0].set_xlabel(r"$\log_{10} z_{\rm peak}$")
        axes[0].set_ylabel(r"$f_{\rm EDE,peak}$")
        axes[1].set_xlabel("max fractional DM mass shift")
        axes[1].set_ylabel(r"$f_{\rm EDE,peak}$")
        for ax in axes:
            ax.grid(True, alpha=0.25)
            ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        examples = [
            ("standard EDE, no coupling", "axion", 3.0, 0.3, 100.0, 0.0, "#6c757d", "-"),
            ("standard EDE + coupling", "axion", 3.0, 0.3, 100.0, 1.0, "#1d3557", "--"),
            ("standard EDE + coupling", "axion", 3.0, 1.0, 100.0, 0.5, "#457b9d", "-."),
            ("tEDS + coupling", "teds", 3.0, 0.05, 0.0, 0.68, "#e76f51", "-"),
        ]
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        for label, model, n, f, m, g, color, ls in examples:
            subset = [
                row for row in compare_rows
                if row["model"] == model and float(row["n"]) == n and float(row["f"]) == f
                and float(row["m"]) == m and float(row["g"]) == g and row.get("background")
            ]
            if not subset:
                continue
            bg = load_background(Path(str(subset[0]["background"])))
            ax.loglog(1.0 + bg["z"], bg["(.)Omega_scf"], color=color, ls=ls, lw=2.2, label=label)
        ax.axvline(1 + ZEQ, color="0.35", ls=":", lw=1.3, label="equality")
        ax.axhspan(0.08, 0.14, color="#2a9d8f", alpha=0.12, label="10% target band")
        ax.invert_xaxis()
        ax.set_xlabel(r"$1+z$")
        ax.set_ylabel(r"$f_{\rm EDE}=\Omega_\phi$")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def write_notebook() -> None:
    NOTEBOOK.parent.mkdir(parents=True, exist_ok=True)
    nb = nbf.v4.new_notebook()
    nb["cells"] = [
        nbf.v4.new_markdown_cell(
            "# tEDS coupling-trigger reproduction\n\n"
            "This notebook reproduces the background-level tEDS trigger behavior from arXiv:2212.08098 using this CLASS fork. "
            "It loads the generated CLASS background files, plots the plateau potential, the field release, the EDE fraction, "
            "the coupling-force dominance, and compares against axion-like potentials."
        ),
        nbf.v4.new_code_cell(
            "from pathlib import Path\n"
            "import csv, math\n"
            "import numpy as np\n"
            "import matplotlib.pyplot as plt\n\n"
            "REPO = Path.cwd().parent if Path.cwd().name == 'notebooks' else Path.cwd()\n"
            "TEDS_DIR = REPO / 'analysis' / 'teds_trigger_check'\n"
            "COMPARE_DIR = REPO / 'analysis' / 'potential_trigger_compare'\n"
            "ZEQ = 3400.0"
        ),
        nbf.v4.new_markdown_cell(
            "## Model definitions\n\n"
            "Standard axion-like EDE uses\n\n"
            "$$V_{\\rm axion}=m^2f^2\\,[1-\\cos(\\phi/f)]^n,$$\n\n"
            "with late-time averaged equation of state $w=(n-1)/(n+1)$.  Thus `n=3` has the desirable "
            "$w=1/2$ dilution, while `n=2` has $w=1/3$.\n\n"
            "The tEDS power-law plateau uses\n\n"
            "$$V_{\\rm tEDS}=V_0\\,\\frac{\\theta^6}{(1+\\theta^{6p})^{1/p}},\\qquad \\theta=\\phi/f.$$\n\n"
            "It is approximately flat at large $\\theta$ but behaves as $\\theta^6$ near the minimum, so it "
            "keeps the fast post-release dilution of an `n=3` axion-like field while separating stored energy "
            "from the bare slope.\n\n"
            "The DM coupling used here is\n\n"
            "$$m_{\\rm DM}(\\phi)=m_0(1+g\\phi^2),\\qquad \\partial_\\phi\\ln m_{\\rm DM}=\\frac{2g\\phi}{1+g\\phi^2}.$$\n\n"
            "In the background Klein-Gordon equation this adds the effective force "
            "$3\\rho_{\\rm DM}\\partial_\\phi\\ln m_{\\rm DM}$ on top of $V'$."
        ),
        nbf.v4.new_markdown_cell("## Optional: regenerate the CLASS backgrounds"),
        nbf.v4.new_code_cell(
            "# Run these if you want to regenerate the data from scratch.\n"
            "# %run ../scripts/teds_trigger_check.py\n"
            "# %run ../scripts/compare_trigger_potentials.py"
        ),
        nbf.v4.new_code_cell(
            "def read_summary(path):\n"
            "    with path.open() as handle:\n"
            "        rows = list(csv.DictReader(handle))\n"
            "    out = []\n"
            "    for row in rows:\n"
            "        clean = {}\n"
            "        for key, value in row.items():\n"
            "            try:\n"
            "                clean[key] = float(value)\n"
            "            except (TypeError, ValueError):\n"
            "                clean[key] = value\n"
            "        out.append(clean)\n"
            "    return out\n\n"
            "def load_background(path):\n"
            "    header = ''\n"
            "    with path.open() as handle:\n"
            "        for line in handle:\n"
            "            if line.startswith('#') and '1:z' in line:\n"
            "                header = line[1:]\n"
            "                break\n"
            "    names = [entry.split(':', 1)[1] for entry in header.split() if ':' in entry]\n"
            "    data = np.loadtxt(path)\n"
            "    return {name: data[:, i] for i, name in enumerate(names)}\n\n"
            "summary = read_summary(TEDS_DIR / 'summary.csv')\n"
            "summary"
        ),
        nbf.v4.new_markdown_cell("## Plateau potential"),
        nbf.v4.new_code_cell(
            "theta = np.logspace(-2, 1.2, 800)\n"
            "fig, ax = plt.subplots(figsize=(6.4, 4.2))\n"
            "for p in [2, 4, 8]:\n"
            "    v = theta**6 / (1 + theta**(6*p))**(1/p)\n"
            "    ax.loglog(theta, v, lw=2, label=f'p={p}')\n"
            "ax.axvspan(4, 8, color='gold', alpha=0.18, label='paper IC range')\n"
            "ax.set_xlabel(r'$\\theta=\\phi/f$')\n"
            "ax.set_ylabel(r'$V/V_0$')\n"
            "ax.grid(True, which='both', alpha=0.25)\n"
            "ax.legend(frameon=False);"
        ),
        nbf.v4.new_markdown_cell("## tEDS release histories"),
        nbf.v4.new_code_cell(
            "def bg_path(row):\n"
            "    theta = f\"{row['theta_i']:g}\".replace('.', 'p')\n"
            "    g = f\"{row['g']:g}\".replace('.', 'p')\n"
            "    return TEDS_DIR / f'theta{theta}_g{g}_coupled00_background.dat'\n\n"
            "data = [(row, load_background(bg_path(row))) for row in summary]\n"
            "colors = {8.0: '#264653', 6.0: '#2a9d8f', 4.0: '#e76f51'}\n"
            "fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.6), sharex=True)\n"
            "for row, bg in data:\n"
            "    c = colors[row['theta_i']]\n"
            "    x = 1 + bg['z']\n"
            "    label = fr'$\\theta_i={row[\"theta_i\"]:g}$, $g={row[\"g\"]:.2f}$'\n"
            "    axes[0,0].loglog(x, bg['(.)Omega_scf'], color=c, lw=2.2, label=label)\n"
            "    axes[0,1].semilogx(x, bg['phi_scf']/0.05, color=c, lw=2.2)\n"
            "    bare = np.abs(bg[\"V'_scf\"])\n"
            "    coupling = np.abs(bg['dVadd_contribution'])\n"
            "    axes[1,0].loglog(x, coupling/np.maximum(bare, 1e-300), color=c, lw=2.2)\n"
            "    mrel = 1 + row['g']*bg['phi_scf']**2\n"
            "    axes[1,1].semilogx(x, mrel/mrel[0]-1, color=c, lw=2.2)\n"
            "for ax in axes.flat:\n"
            "    ax.axvline(1+ZEQ, color='0.35', ls='--')\n"
            "    ax.invert_xaxis(); ax.grid(True, which='both', alpha=0.25)\n"
            "axes[0,0].set_ylabel(r'$\\Omega_\\phi$'); axes[0,1].set_ylabel(r'$\\phi/f$')\n"
            "axes[1,0].set_ylabel('coupling force / bare force'); axes[1,1].set_ylabel(r'$m_{DM}/m_{DM,ini}-1$')\n"
            "axes[1,0].set_xlabel(r'$1+z$'); axes[1,1].set_xlabel(r'$1+z$')\n"
            "axes[0,0].legend(frameon=False); fig.tight_layout()"
        ),
        nbf.v4.new_markdown_cell("## Potential comparison"),
        nbf.v4.new_code_cell(
            "rows = [row for row in read_summary(COMPARE_DIR / 'summary.csv') if row.get('ok') == 'yes']\n"
            "fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.9))\n"
            "styles = {('axion', 2.0): ('o', '#457b9d', 'axion n=2'), ('axion', 3.0): ('s', '#1d3557', 'axion n=3'), ('teds', 3.0): ('^', '#e76f51', 'tEDS')}\n"
            "for (model, n), (marker, color, label) in styles.items():\n"
            "    subset = [row for row in rows if row['model'] == model and float(row['n']) == n]\n"
            "    axes[0].scatter([r['log10_z_peak'] for r in subset], [r['f_ede_peak'] for r in subset], marker=marker, color=color, s=54, label=label)\n"
            "    axes[1].scatter([r['dm_mass_shift_max'] for r in subset], [r['f_ede_peak'] for r in subset], marker=marker, color=color, s=54, label=label)\n"
            "axes[0].axvline(math.log10(ZEQ), color='0.35', ls='--')\n"
            "axes[0].axhspan(0.08, 0.14, color='#2a9d8f', alpha=0.12)\n"
            "axes[0].set_xlabel(r'$\\log_{10} z_{peak}$'); axes[0].set_ylabel(r'$f_{EDE,peak}$')\n"
            "axes[1].set_xlabel('max fractional DM mass shift'); axes[1].set_ylabel(r'$f_{EDE,peak}$')\n"
            "for ax in axes: ax.grid(True, alpha=0.25); ax.legend(frameon=False)\n"
            "fig.tight_layout()"
        ),
        nbf.v4.new_markdown_cell("## Explicit standard-EDE comparison"),
        nbf.v4.new_code_cell(
            "examples = [\n"
            "    ('standard EDE, no coupling', 'axion', 3.0, 0.3, 100.0, 0.0, '#6c757d', '-'),\n"
            "    ('standard EDE + coupling', 'axion', 3.0, 0.3, 100.0, 1.0, '#1d3557', '--'),\n"
            "    ('standard EDE + coupling', 'axion', 3.0, 1.0, 100.0, 0.5, '#457b9d', '-.'),\n"
            "    ('tEDS + coupling', 'teds', 3.0, 0.05, 0.0, 0.68, '#e76f51', '-'),\n"
            "]\n"
            "fig, ax = plt.subplots(figsize=(7.2, 4.4))\n"
            "for label, model, n, f, m, g, color, ls in examples:\n"
            "    subset = [r for r in rows if r['model'] == model and float(r['n']) == n and float(r['f']) == f and float(r['m']) == m and float(r['g']) == g and r.get('background')]\n"
            "    if not subset:\n"
            "        continue\n"
            "    bg = load_background(Path(subset[0]['background']))\n"
            "    ax.loglog(1 + bg['z'], bg['(.)Omega_scf'], color=color, ls=ls, lw=2.2, label=label)\n"
            "ax.axvline(1+ZEQ, color='0.35', ls=':', lw=1.3, label='equality')\n"
            "ax.axhspan(0.08, 0.14, color='#2a9d8f', alpha=0.12, label='10% target band')\n"
            "ax.invert_xaxis(); ax.grid(True, which='both', alpha=0.25)\n"
            "ax.set_xlabel(r'$1+z$'); ax.set_ylabel(r'$f_{EDE}=\\Omega_\\phi$')\n"
            "ax.legend(frameon=False, fontsize=8); fig.tight_layout()"
        ),
    ]
    nbf.write(nb, NOTEBOOK)


def main() -> None:
    write_report()
    write_notebook()
    print(f"notebook = {NOTEBOOK}")
    print(f"report = {REPORT}")


if __name__ == "__main__":
    main()
