from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


def _parse(
    path: Path,
) -> dict[str, tuple[tuple[float, float, float], tuple[float, float, float]]]:
    bodies: dict[str, tuple[tuple[float, float, float], tuple[float, float, float]]] = {}
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < 12:
                continue
            name = p[0]
            pos = (float(p[3]), float(p[4]), float(p[5]))
            vel = (float(p[6]), float(p[7]), float(p[8]))
            bodies[name] = (pos, vel)
    return bodies


def compare(cpp_path: Path, ref_path: Path) -> dict[str, Any]:
    """
    Compare final snapshot TXT files (same format as SaveSystemToTextFile / REBOUND export).
    Positions are in meters, velocities in m/s.

    Returns relative errors (dimensionless) and absolute Euclidean errors:
    - ΔR = |r_cpp - r_ref|  (meters); reported in km for aggregates and per_body.
    - ΔV = |v_cpp - v_ref|  (m/s).
    """
    a = _parse(cpp_path)
    b = _parse(ref_path)
    names = sorted(set(a) & set(b))
    if not names:
        return {"error": "no common bodies", "n_common": 0}

    pos_rel: list[float] = []
    vel_rel: list[float] = []
    pos_abs_m: list[float] = []
    vel_abs_ms: list[float] = []
    per_body: list[dict[str, Any]] = []

    for n in names:
        pa, va = a[n]
        pb, vb = b[n]
        dpos = math.sqrt(sum((x - y) ** 2 for x, y in zip(pa, pb)))
        dvel = math.sqrt(sum((x - y) ** 2 for x, y in zip(va, vb)))
        rb = math.sqrt(sum(x * x for x in pb))
        vbm = math.sqrt(sum(x * x for x in vb))
        rp = dpos / max(rb, 1.0)
        rv = dvel / max(vbm, 1e-10)
        pos_rel.append(rp)
        vel_rel.append(rv)
        pos_abs_m.append(dpos)
        vel_abs_ms.append(dvel)
        per_body.append(
            {
                "name": n,
                "rel_pos_error": rp,
                "rel_vel_error": rv,
                "abs_pos_error_km": dpos / 1000.0,
                "abs_vel_error": dvel,
            }
        )

    n = len(names)
    return {
        "n_common": n,
        "mean_rel_pos_error": sum(pos_rel) / n,
        "max_rel_pos_error": max(pos_rel),
        "mean_rel_vel_error": sum(vel_rel) / n,
        "max_rel_vel_error": max(vel_rel),
        "mean_abs_pos_error": sum(pos_abs_m) / n / 1000.0,
        "max_abs_pos_error": max(pos_abs_m) / 1000.0,
        "mean_abs_vel_error": sum(vel_abs_ms) / n,
        "max_abs_vel_error": max(vel_abs_ms),
        "units": {
            "mean_abs_pos_error": "km",
            "max_abs_pos_error": "km",
            "mean_abs_vel_error": "m/s",
            "max_abs_vel_error": "m/s",
            "per_body.abs_pos_error_km": "km",
            "per_body.abs_vel_error": "m/s",
        },
        "per_body": per_body,
    }


def write_metrics_log(metrics: dict[str, Any], path: Path) -> None:
    """Human-readable log including per-body relative and absolute metrics."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "Verification metrics",
        "==================",
        json.dumps({k: v for k, v in metrics.items() if k != "per_body"}, indent=2),
        "",
        "Per body (rel = dimensionless vs ref | ΔR in km | ΔV in m/s):",
        "----------------------------------------------------------------",
    ]
    for row in metrics.get("per_body", []):
        lines.append(
            f"  {row['name']}: "
            f"rel_pos={row['rel_pos_error']:.6e}  rel_vel={row['rel_vel_error']:.6e}  |  "
            f"ΔR={row['abs_pos_error_km']:.6e} km  ΔV={row['abs_vel_error']:.6e} m/s"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
