from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def plot_per_body_dashboard(metrics: dict[str, Any], output_path: Path) -> None:
    """
    One figure: four subplots — relative pos/vel and absolute ΔR (km) / ΔV (m/s) per body.
    Absolute-error subplots use log scale on Y when all values are positive.
    """
    rows = metrics.get("per_body") or []
    if not rows:
        return

    names = [r["name"] for r in rows]
    x = range(len(names))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n = len(names)
    fig_w = max(14, min(56, 6 + n * 0.28))
    tick_fs = max(5, min(9, 14 - n // 25))
    fig, axes = plt.subplots(2, 2, figsize=(fig_w, 10))

    ax = axes[0, 0]
    ax.bar(x, [r["rel_pos_error"] for r in rows], color="#1f77b4")
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=tick_fs)
    ax.set_ylabel("Relative position error")
    ax.set_title("rel_pos_error = |r_cpp−r_ref| / |r_ref|")
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[0, 1]
    ax.bar(x, [r["rel_vel_error"] for r in rows], color="#ff7f0e")
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=tick_fs)
    ax.set_ylabel("Relative velocity error")
    ax.set_title("rel_vel_error = |v_cpp−v_ref| / |v_ref|")
    ax.grid(True, axis="y", alpha=0.3)

    ax = axes[1, 0]
    abs_km = [max(r["abs_pos_error_km"], 1e-30) for r in rows]
    ax.bar(x, abs_km, color="#2ca02c")
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=tick_fs)
    ax.set_yscale("log")
    ax.set_ylabel("Absolute error (km)")
    ax.set_title("ΔR = |r_cpp − r_ref| (Euclidean)")
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:g}"))
    ax.grid(True, which="both", axis="y", alpha=0.3)

    ax = axes[1, 1]
    abs_v = [max(r["abs_vel_error"], 1e-30) for r in rows]
    ax.bar(x, abs_v, color="#d62728")
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=tick_fs)
    ax.set_yscale("log")
    ax.set_ylabel("Absolute ΔV (m/s)")
    ax.set_title("ΔV = |v_cpp − v_ref| (Euclidean)")
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:g}"))
    ax.grid(True, which="both", axis="y", alpha=0.3)

    fig.suptitle("C++ vs reference: per-body errors", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_abs_mean_vs_dt(
    dt_values: list[float],
    mean_abs_pos_err_km: list[float | None],
    output_path: Path,
) -> None:
    """Log–log: mean absolute position error (km) vs integration step dt (years)."""
    pairs = [(dt, e) for dt, e in zip(dt_values, mean_abs_pos_err_km) if e is not None]
    if not pairs:
        return
    valid_dt, valid_e = zip(*pairs)
    safe_e = [max(float(e), 1e-30) for e in valid_e]

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 6))
    plt.plot(valid_dt, safe_e, marker="o", markersize=8, linestyle="-", color="#2ca02c", linewidth=2)
    plt.xscale("log")
    plt.yscale("log")
    ax = plt.gca()
    ax.invert_xaxis()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f"{y:g}"))

    plt.xlabel("Integration step dt (years)", fontsize=12)
    plt.ylabel("Mean absolute position error (km)", fontsize=12)
    plt.title("Mean absolute ΔR (km) vs dt — C++ vs REBOUND", fontsize=13)
    plt.grid(True, which="major", ls="-", color="gray", alpha=0.5)
    plt.grid(True, which="minor", ls="--", color="gray", alpha=0.2)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
