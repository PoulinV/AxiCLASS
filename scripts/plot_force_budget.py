#!/usr/bin/env python3
"""Plot bare-potential vs coupling contributions for representative coupled runs."""

from __future__ import annotations

import argparse
import configparser
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


@dataclass
class RunSpec:
    name: str
    background: Path
    ini: Path
    kind: str  # "axion" or "teds"


def read_ini_value(path: Path, key: str, cast=float):
    text = path.read_text().splitlines()
    for line in text:
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        lhs, rhs = line.split("=", 1)
        if lhs.strip() == key:
            return cast(rhs.strip())
    raise KeyError(f"{key} not found in {path}")


def load_background(path: Path):
    data = np.loadtxt(path, comments="#")
    header = path.read_text().splitlines()
    colline = next(line for line in header if line.startswith("#    1:z"))
    return data, colline


def compute_terms(spec: RunSpec):
    data, _ = load_background(spec.background)
    z = data[:, 0]
    H = data[:, 3]
    rho_idm = data[:, 17]
    phi = data[:, 25]
    dV = data[:, 29]

    g = float(read_ini_value(spec.ini, "g_idm_ede"))
    if spec.kind == "axion":
        f = float(read_ini_value(spec.ini, "f_axion"))
        dlnm = 2.0 * g * phi / (1.0 + g * phi**2)
        label = rf"$(1-\cos)^3$: $m/H_0={read_ini_value(spec.ini, 'm_axion'):.0f}$, $f={f:g}$, $g={g:g}$"
    else:
        dlnm = 2.0 * g * phi / (1.0 + g * phi**2)
        label = rf"tEDS: $\theta_i\approx 6$, $g={g:g}$"

    coupling = 3.0 * rho_idm * dlnm
    bare = dV
    total = bare + coupling
    scale = np.maximum(H**2, 1e-300)
    out = {
        "z": z,
        "bare": np.abs(bare) / scale,
        "coupling": np.abs(coupling) / scale,
        "total": np.abs(total) / scale,
        "fraction": np.abs(coupling) / np.maximum(np.abs(bare) + np.abs(coupling), 1e-300),
        "label": label,
    }
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        default="analysis/report_assets",
        help="Directory where the plot files will be written.",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    outdir = (root / args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    specs = [
        RunSpec(
            name="tEDS",
            background=root / "analysis/teds_trigger_check/theta6_g0p68_coupled00_background.dat",
            ini=root / "analysis/teds_trigger_check/theta6_g0p68_coupled.ini",
            kind="teds",
        ),
        RunSpec(
            name="Axion",
            background=root / "analysis/axion_natural_trigger_report/m30_f1_g100_background.dat",
            ini=root / "analysis/axion_natural_trigger_report/m30_f1_g1.ini",
            kind="axion",
        ),
    ]

    runs = [compute_terms(spec) for spec in specs]
    eq_z = 3400.0

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    fig.subplots_adjust(wspace=0.14, hspace=0.12)

    colors = {"bare": "#1f77b4", "coupling": "#d62728", "total": "#2ca02c"}
    for col, run in enumerate(runs):
        ax = axes[0, col]
        ax.plot(run["z"], run["bare"], color=colors["bare"], lw=2.0, label=r"$|V_{,\phi}|/H^2$")
        ax.plot(run["z"], run["coupling"], color=colors["coupling"], lw=2.0, label=r"$|3\rho_{\rm IDM}\,d\ln m/d\phi|/H^2$")
        ax.plot(run["z"], run["total"], color=colors["total"], lw=1.8, ls="--", label=r"$|V_{,\phi}+3\rho_{\rm IDM}d\ln m/d\phi|/H^2$")
        ax.axvline(eq_z, color="0.5", lw=1.0, ls=":")
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.grid(True, which="both", alpha=0.22)
        ax.set_title(run["label"], fontsize=11)
        if col == 0:
            ax.set_ylabel(r"Force terms / $H^2$")
        ax.legend(fontsize=8, loc="best", frameon=False)

        ax2 = axes[1, col]
        ax2.plot(run["z"], run["fraction"], color="#9467bd", lw=2.0)
        ax2.axhline(0.5, color="0.5", lw=1.0, ls=":")
        ax2.axvline(eq_z, color="0.5", lw=1.0, ls=":")
        ax2.set_ylim(0.0, 1.05)
        ax2.set_xscale("log")
        ax2.invert_xaxis()
        ax2.grid(True, which="both", alpha=0.22)
        if col == 0:
            ax2.set_ylabel(r"Coupling fraction")
        ax2.set_xlabel(r"$z$")

    fig.suptitle("Bare-slope and coupling contributions to the scalar-field force", y=0.98, fontsize=13)

    png = outdir / "force_budget_comparison.png"
    pdf = outdir / "force_budget_comparison.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")


if __name__ == "__main__":
    main()
