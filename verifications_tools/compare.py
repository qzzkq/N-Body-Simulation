from __future__ import annotations
import json, math
from pathlib import Path
from typing import Any


# value: (pos_m, vel_ms, mass_kg)
_Snapshot = dict[str, tuple[tuple[float, float, float],
                            tuple[float, float, float],
                            float]]


def _parse(path: Path) -> _Snapshot:
    bodies: _Snapshot = {}
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < 12:
                continue
            bodies[p[0]] = (
                (float(p[3]), float(p[4]), float(p[5])),
                (float(p[6]), float(p[7]), float(p[8])),
                float(p[1]),
            )
    return bodies


def _com_pos(snap: _Snapshot) -> tuple[float, float, float]:
    mx = my = mz = mt = 0.0
    for pos, _vel, mass in snap.values():
        mx += mass * pos[0]
        my += mass * pos[1]
        mz += mass * pos[2]
        mt += mass
    return (mx / mt, my / mt, mz / mt) if mt else (0.0, 0.0, 0.0)


def _com_vel(snap: _Snapshot) -> tuple[float, float, float]:
    mvx = mvy = mvz = mt = 0.0
    for _pos, vel, mass in snap.values():
        mvx += mass * vel[0]
        mvy += mass * vel[1]
        mvz += mass * vel[2]
        mt += mass
    return (mvx / mt, mvy / mt, mvz / mt) if mt else (0.0, 0.0, 0.0)


def compare(cpp_path: Path, ref_path: Path, *, com_align: bool = True) -> dict[str, Any]:
    a = _parse(cpp_path)
    b = _parse(ref_path)
    names = sorted(set(a) & set(b))
    if not names:
        return {"error": "no common bodies", "n_common": 0}

    com_a = _com_pos(a) if com_align else (0.0, 0.0, 0.0)
    com_b = _com_pos(b) if com_align else (0.0, 0.0, 0.0)
    cvel_a = _com_vel(a)
    cvel_b = _com_vel(b)
    com_vel_drift_ms = math.sqrt(sum((x - y) ** 2 for x, y in zip(cvel_a, cvel_b)))
    com_pos_drift_km = math.sqrt(sum((x - y) ** 2 for x, y in zip(com_a, com_b))) / 1e3

    per_body: list[dict] = []
    pos_rel, vel_rel = [], []
    pos_abs_m, vel_abs_ms = [], []

    for name in names:
        pa, va, _ = a[name]
        pb, vb, _ = b[name]

        pa_s = (pa[0] - com_a[0], pa[1] - com_a[1], pa[2] - com_a[2])
        pb_s = (pb[0] - com_b[0], pb[1] - com_b[1], pb[2] - com_b[2])

        dpos = math.sqrt(sum((x - y) ** 2 for x, y in zip(pa_s, pb_s)))
        dvel = math.sqrt(sum((x - y) ** 2 for x, y in zip(va, vb)))
        rb   = math.sqrt(sum(x * x for x in pb_s))
        vbm  = math.sqrt(sum(x * x for x in vb))

        rp = dpos / max(rb, 1.0)
        rv = dvel / max(vbm, 1e-10)

        pos_rel.append(rp);  vel_rel.append(rv)
        pos_abs_m.append(dpos); vel_abs_ms.append(dvel)
        per_body.append({
            "name": name,
            "rel_pos_error":   rp,
            "rel_vel_error":   rv,
            "abs_pos_error_km": dpos / 1e3,
            "abs_vel_error_ms": dvel,
        })

    n = len(names)
    return {
        "n_common": n,
        "com_aligned": com_align,
        "com_pos_drift_km": com_pos_drift_km,
        "com_vel_drift_ms": com_vel_drift_ms,
        "mean_rel_pos_error": sum(pos_rel) / n,
        "max_rel_pos_error":  max(pos_rel),
        "mean_rel_vel_error": sum(vel_rel) / n,
        "max_rel_vel_error":  max(vel_rel),
        "mean_abs_pos_error_km": sum(pos_abs_m) / n / 1e3,
        "max_abs_pos_error_km":  max(pos_abs_m) / 1e3,
        "mean_abs_vel_error_ms": sum(vel_abs_ms) / n,
        "max_abs_vel_error_ms":  max(vel_abs_ms),
        "per_body": per_body,
    }


def write_log(metrics: dict[str, Any], path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    agg = {k: v for k, v in metrics.items() if k != "per_body"}
    lines = [
        "Verification metrics",
        "=" * 40,
        json.dumps(agg, indent=2),
        "",
        "Per body  (rel=dimensionless | pos in km | vel in m/s):",
        "-" * 60,
    ]
    for row in metrics.get("per_body", []):
        lines.append(
            f"  {row['name']:12s}  rel_pos={row['rel_pos_error']:.3e}"
            f"  rel_vel={row['rel_vel_error']:.3e}"
            f"  ΔR={row['abs_pos_error_km']:.3e} km"
            f"  ΔV={row['abs_vel_error_ms']:.3e} m/s"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_per_body(metrics: dict[str, Any], out: Path) -> None:
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("matplotlib not available, skipping plot")
        return

    rows = metrics.get("per_body") or []
    if not rows:
        return
    names = [r["name"] for r in rows]
    x = range(len(names))
    n = len(names)
    tick_fs = max(5, min(9, 14 - n // 25))
    fig_w = max(14, min(56, 6 + n * 0.28))

    fig, axes = plt.subplots(2, 2, figsize=(fig_w, 10))
    note = "  (COM-aligned)" if metrics.get("com_aligned") else ""

    for ax, key, label, color, title in [
        (axes[0, 0], "rel_pos_error",   "Relative pos error",   "#1f77b4", f"rel ΔR{note}"),
        (axes[0, 1], "rel_vel_error",   "Relative vel error",   "#ff7f0e", f"rel ΔV{note}"),
        (axes[1, 0], "abs_pos_error_km","Absolute ΔR (km)",     "#2ca02c", "ΔR (km)"),
        (axes[1, 1], "abs_vel_error_ms","Absolute ΔV (m/s)",    "#d62728", "ΔV (m/s)"),
    ]:
        vals = [max(r[key], 1e-30) for r in rows]
        ax.bar(x, vals, color=color)
        ax.set_xticks(list(x))
        ax.set_xticklabels(names, rotation=45, ha="right", fontsize=tick_fs)
        ax.set_ylabel(label)
        ax.set_title(title)
        if "abs" in key:
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:g}"))
        ax.grid(True, axis="y", alpha=0.3)

    drift = metrics.get("com_pos_drift_km", 0)
    fig.suptitle(f"C++ vs REBOUND — per-body errors  |  COM drift={drift:.2e} km", fontsize=13)
    fig.tight_layout()
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_dt_sweep(dt_values: list[float],
                  mean_abs_km: list[float | None],
                  out: Path) -> None:
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        return
    pairs = [(dt, e) for dt, e in zip(dt_values, mean_abs_km) if e is not None]
    if not pairs:
        return
    dts, errs = zip(*pairs)
    plt.figure(figsize=(10, 6))
    plt.plot(dts, [max(e, 1e-30) for e in errs], "o-", color="#2ca02c", linewidth=2)
    plt.xscale("log"); plt.yscale("log")
    plt.gca().invert_xaxis()
    plt.xlabel("dt (years)"); plt.ylabel("Mean |ΔR| (km)")
    plt.title("Accuracy vs integration step — C++ vs REBOUND")
    plt.grid(True, which="both", ls="--", alpha=0.4)
    plt.tight_layout()
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=300)
    plt.close()
