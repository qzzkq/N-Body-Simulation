from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from verify_io import read_system_txt


def compare(cpp_txt: Path, ref_txt: Path) -> dict:
    cpp = {r.name: r for r in read_system_txt(cpp_txt)}
    ref = {r.name: r for r in read_system_txt(ref_txt)}
    common = sorted(set(cpp.keys()) & set(ref.keys()))
    if not common:
        raise ValueError("No common body names between snapshots.")

    pos_rel_errors = []
    vel_rel_errors = []
    per_body = {}
    for name in common:
        a = cpp[name]
        b = ref[name]
        dp = math.sqrt((a.px_m - b.px_m) ** 2 + (a.py_m - b.py_m) ** 2 + (a.pz_m - b.pz_m) ** 2)
        rp = max(1e-9, math.sqrt(b.px_m**2 + b.py_m**2 + b.pz_m**2))
        dv = math.sqrt((a.vx_ms - b.vx_ms) ** 2 + (a.vy_ms - b.vy_ms) ** 2 + (a.vz_ms - b.vz_ms) ** 2)
        rv = max(1e-9, math.sqrt(b.vx_ms**2 + b.vy_ms**2 + b.vz_ms**2))
        ep = (dp / rp) * 100.0
        ev = (dv / rv) * 100.0
        pos_rel_errors.append(ep)
        vel_rel_errors.append(ev)
        per_body[name] = {"pos_err_pct": ep, "vel_err_pct": ev}

    return {
        "matched_bodies": len(common),
        "mean_pos_err_pct": sum(pos_rel_errors) / len(pos_rel_errors),
        "max_pos_err_pct": max(pos_rel_errors),
        "mean_vel_err_pct": sum(vel_rel_errors) / len(vel_rel_errors),
        "max_vel_err_pct": max(vel_rel_errors),
        "per_body": per_body,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare C++ snapshot against reference.")
    parser.add_argument("--cpp", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--json-out", default=None)
    args = parser.parse_args()

    summary = compare(Path(args.cpp), Path(args.ref))
    print(json.dumps(summary, indent=2))
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(f"Saved report -> {args.json_out}")


if __name__ == "__main__":
    main()

