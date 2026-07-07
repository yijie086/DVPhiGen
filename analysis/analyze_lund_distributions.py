#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/dvphigen_mplconfig")

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_PROJECT = Path(__file__).resolve().parents[1]
PID_LABELS = {
    11: "e'",
    2212: "p'",
    321: "K+",
    -321: "K-",
}


def run(cmd):
    subprocess.run(cmd, check=True)


def compile_analyzer(project_dir, cpp_path, exe_path):
    src_dir = project_dir / "src"
    cmd = [
        "g++",
        "-std=c++17",
        "-O2",
        "-Wall",
        "-Wextra",
        "-I",
        str(src_dir),
        str(cpp_path),
        str(src_dir / "kinematics.cpp"),
        str(src_dir / "event.cpp"),
        "-o",
        str(exe_path),
    ]
    run(cmd)


def read_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return rows


def col(rows, key, dtype=float):
    return np.array([dtype(r[key]) for r in rows])


def padded_range(vals, frac=0.035):
    lo = float(np.nanmin(vals))
    hi = float(np.nanmax(vals))
    width = hi - lo
    if width <= 0.0:
        pad = max(abs(lo) * frac, frac)
    else:
        pad = width * frac
    return lo - pad, hi + pad


def plot_particle_level(rows, output_png):
    pids = [11, 2212, 321, -321]
    variables = [
        ("p", "p [GeV]"),
        ("theta_deg", r"$\theta$ [deg]"),
        ("phi_deg", r"$\phi$ [deg]"),
    ]

    fig, axes = plt.subplots(len(variables), len(pids), figsize=(15, 8.5))
    for j, pid in enumerate(pids):
        subset = [r for r in rows if int(r["pid"]) == pid]
        for i, (key, label) in enumerate(variables):
            ax = axes[i, j]
            vals = col(subset, key)
            ax.hist(vals, bins=80, histtype="stepfilled", alpha=0.85, color="#2f6f9f")
            if i == 0:
                ax.set_title(PID_LABELS[pid])
            ax.set_xlabel(label)
            ax.set_ylabel("counts")
            ax.grid(alpha=0.22)

    fig.suptitle("Particle-level LUND distributions")
    fig.tight_layout()
    fig.savefig(output_png, dpi=170)
    plt.close(fig)


def plot_event_level(rows, output_png):
    variables = [
        ("xB", r"$x_B$"),
        ("Q2", r"$Q^2$ [GeV$^2$]"),
        ("W", r"$W$ [GeV]"),
        ("t", r"$t$ [GeV$^2$]"),
        ("phi", r"$\phi$ [rad]"),
        ("cosThetaH", r"$\cos\theta_H$"),
        ("phiH", r"$\phi_H$ [rad]"),
        ("epsilon", r"$\epsilon$"),
        ("Mkk", r"$M(K^+K^-)$ [GeV]"),
    ]

    fig, axes = plt.subplots(3, 3, figsize=(13.5, 10))
    for ax, (key, label) in zip(axes.ravel(), variables):
        vals = col(rows, key)
        ax.hist(vals, bins=80, histtype="stepfilled", alpha=0.85, color="#8a5a2b")
        ax.set_xlabel(label)
        ax.set_ylabel("counts")
        ax.grid(alpha=0.22)

    fig.suptitle("Event-level reconstructed LUND distributions")
    fig.tight_layout()
    fig.savefig(output_png, dpi=170)
    plt.close(fig)


def plot_event_2d(rows, output_png):
    variables = [
        ("xB", r"$x_B$"),
        ("Q2", r"$Q^2$ [GeV$^2$]"),
        ("W", r"$W$ [GeV]"),
        ("minus_t", r"$-t$ [GeV$^2$]"),
    ]
    panels = [
        (variables[i], variables[j])
        for i in range(len(variables))
        for j in range(i + 1, len(variables))
    ]
    ranges = {key: padded_range(col(rows, key)) for key, _ in variables}

    fig, axes = plt.subplots(2, 3, figsize=(15, 8.6))
    for ax, ((xkey, xlabel), (ykey, ylabel)) in zip(axes.ravel(), panels):
        xvals = col(rows, xkey)
        yvals = col(rows, ykey)
        hist = ax.hist2d(
            xvals,
            yvals,
            bins=80,
            range=[ranges[xkey], ranges[ykey]],
            cmap="viridis",
        )
        ax.set_xlim(ranges[xkey])
        ax.set_ylim(ranges[ykey])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.18)
        cbar = fig.colorbar(hist[3], ax=ax)
        cbar.set_label("counts")

    fig.suptitle(r"Event-level pairwise 2D distributions: $x_B$, $Q^2$, $W$, $-t$")
    fig.tight_layout()
    fig.savefig(output_png, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze DVPhiGen LUND output and draw particle/event distributions."
    )
    parser.add_argument("lund", help="Input LUND file")
    parser.add_argument(
        "--project-dir",
        default=str(DEFAULT_PROJECT),
        help="DVPhiGen project directory containing src/kinematics.cpp",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory for CSV and PNG files",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Output file prefix. Default: input LUND stem",
    )
    parser.add_argument(
        "--eb",
        default="10.6",
        help="Fallback beam energy if the LUND header does not contain it",
    )
    args = parser.parse_args()

    lund_path = Path(args.lund).resolve()
    project_dir = Path(args.project_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or lund_path.stem

    script_dir = Path(__file__).resolve().parent
    cpp_path = script_dir / "analyze_lund.cpp"
    exe_path = outdir / f"{prefix}_lund_analyzer"
    particle_csv = outdir / f"{prefix}_particle_level.csv"
    event_csv = outdir / f"{prefix}_event_level.csv"
    particle_png = outdir / f"{prefix}_particle_level.png"
    event_png = outdir / f"{prefix}_event_level.png"
    event_2d_png = outdir / f"{prefix}_event_2d.png"

    compile_analyzer(project_dir, cpp_path, exe_path)
    run([str(exe_path), str(lund_path), str(particle_csv), str(event_csv), str(args.eb)])

    particle_rows = read_csv(particle_csv)
    event_rows = read_csv(event_csv)
    if not particle_rows:
        raise SystemExit("No particle rows were reconstructed.")
    if not event_rows:
        raise SystemExit("No event rows were reconstructed.")

    plot_particle_level(particle_rows, particle_png)
    plot_event_level(event_rows, event_png)
    plot_event_2d(event_rows, event_2d_png)

    print(f"particle CSV: {particle_csv}")
    print(f"event CSV:    {event_csv}")
    print(f"particle PNG: {particle_png}")
    print(f"event PNG:    {event_png}")
    print(f"event 2D PNG: {event_2d_png}")


if __name__ == "__main__":
    main()
