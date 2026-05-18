#!/usr/bin/env python3

import csv
import math
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
BASE = REPO / "analysis" / "type1_trigger_scan"
PLOTS = BASE / "plots"


def dlnm_dphi(phi, coupling, form):
    if form == "exp":
        return coupling
    return 2.0 * coupling * phi / (1.0 + coupling * phi * phi) if phi != 0.0 else 0.0


def dV_axion(phi, m_axion, f_axion, h0=72.81, n_axion=3):
    m = m_axion * h0
    if n_axion > 1:
        return n_axion * (m ** 2) * f_axion * ((1.0 - math.cos(phi / f_axion)) ** (n_axion - 1)) * math.sin(phi / f_axion)
    return (m ** 2) * f_axion * math.sin(phi / f_axion)


def phiN(row):
    a = 1.0 / (1.0 + row[0])
    H = row[3]
    return row[24] / (a * H) if a * H else float("nan")


def read_background(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            rows.append([float(x) for x in line.split()])
    return rows


def write_scatter(results):
    path = PLOTS / "scatter.dat"
    with open(path, "w") as handle:
        for row in results:
            if row["ok"] != "True":
                continue
            form_code = 0 if row["form"] == "poly" else 1
            handle.write(f"{row['z_peak']} {row['omega_peak']} {form_code}\n")
    return path


def write_history_files(results):
    ok = [row for row in results if row["ok"] == "True"]
    ok.sort(key=lambda row: float(row["score"]))
    top = ok[:4]
    paths = []
    for idx, row in enumerate(top, start=1):
        bg = Path(row["background_path"])
        rows = read_background(bg)
        out = PLOTS / f"history_{idx}.dat"
        form = row["form"]
        coupling = float(row["coupling"])
        m_axion = float(row["m_axion"])
        f_axion = float(row["f_axion"])
        with open(out, "w") as handle:
            for bg_row in rows:
                ratio_denom = dV_axion(bg_row[23], m_axion, f_axion)
                ratio = abs(3.0 * bg_row[15] * dlnm_dphi(bg_row[23], coupling, form) / ratio_denom) if ratio_denom != 0.0 else 0.0
                handle.write(f"{1.0 + bg_row[0]} {bg_row[18]} {phiN(bg_row)} {ratio}\n")
        paths.append((out, row))
    return paths


def write_gnuplot(history_files):
    script = PLOTS / "plots.gp"
    colors = ["'#1f77b4'", "'#d62728'", "'#2ca02c'", "'#9467bd'"]
    labels = []
    for _, row in history_files:
        labels.append(f"{row['form']}, th={float(row['theta_i']):.7f}, f={float(row['f_axion']):g}, m={float(row['m_axion']):g}, c={float(row['coupling']):g}")

    lines = [
        "set terminal pngcairo size 1200,800 enhanced font ',10'",
        f"set output '{(PLOTS / 'scan_scatter.png').as_posix()}'",
        "set logscale xy",
        "set xlabel 'z_peak(f_EDE)'",
        "set ylabel 'f_EDE,peak'",
        "set grid",
        "set key left top",
        "set arrow 1 from 3400, graph 0 to 3400, graph 1 nohead dt 2 lc rgb 'black'",
        "set arrow 2 from graph 0, first 0.1 to graph 1, first 0.1 nohead dt 2 lc rgb 'black'",
        "plot \\",
        f"  '{(PLOTS / 'scatter.dat').as_posix()}' using ($3==0?$1:1/0):2 with points pt 7 ps 1.5 lc rgb '#1f77b4' title 'poly', \\",
        f"  '{(PLOTS / 'scatter.dat').as_posix()}' using ($3==1?$1:1/0):2 with points pt 5 ps 1.5 lc rgb '#d62728' title 'exp'",
        "",
        "unset arrow 1",
        "unset arrow 2",
        "unset output",
        "set terminal pngcairo size 1200,1400 enhanced font ',10'",
        f"set output '{(PLOTS / 'top_histories.png').as_posix()}'",
        "set multiplot layout 3,1 title 'Top candidate histories'",
        "set logscale x",
        "set xrange [1e5:1]",
        "set key outside",
        "set grid",
        "set xlabel '1+z'",
        "set ylabel 'f_EDE'",
        "set arrow 11 from 3401, graph 0 to 3401, graph 1 nohead dt 2 lc rgb 'black'",
        "set arrow 12 from graph 0, first 0.1 to graph 1, first 0.1 nohead dt 2 lc rgb 'black'",
        "plot \\",
    ]
    for i, (path, _) in enumerate(history_files):
        sep = ", \\" if i < len(history_files) - 1 else ""
        lines.append(f"  '{path.as_posix()}' using 1:2 with lines lw 2 lc rgb {colors[i]} title '{labels[i]}'{sep}")
    lines += [
        "unset arrow 11",
        "unset arrow 12",
        "set ylabel \"phi'/(aH)\"",
        "set arrow 21 from 3401, graph 0 to 3401, graph 1 nohead dt 2 lc rgb 'black'",
        "plot \\",
    ]
    for i, (path, _) in enumerate(history_files):
        sep = ", \\" if i < len(history_files) - 1 else ""
        lines.append(f"  '{path.as_posix()}' using 1:3 with lines lw 2 lc rgb {colors[i]} title '{labels[i]}'{sep}")
    lines += [
        "unset arrow 21",
        "set ylabel '|3 rho dlnm/dphi / V_phi|'",
        "set arrow 31 from 3401, graph 0 to 3401, graph 1 nohead dt 2 lc rgb 'black'",
        "set arrow 32 from graph 0, first 1 to graph 1, first 1 nohead dt 2 lc rgb 'black'",
        "plot \\",
    ]
    for i, (path, _) in enumerate(history_files):
        sep = ", \\" if i < len(history_files) - 1 else ""
        lines.append(f"  '{path.as_posix()}' using 1:4 with lines lw 2 lc rgb {colors[i]} title '{labels[i]}'{sep}")
    lines += [
        "unset multiplot",
        "unset output",
    ]
    script.write_text("\n".join(lines) + "\n")
    return script


def main():
    PLOTS.mkdir(parents=True, exist_ok=True)
    results = list(csv.DictReader(open(BASE / "scan_results.csv")))
    write_scatter(results)
    history_files = write_history_files(results)
    gp = write_gnuplot(history_files)
    subprocess.run(["gnuplot", str(gp)], check=True, cwd=REPO)
    print("Wrote plots to", PLOTS)


if __name__ == "__main__":
    main()
