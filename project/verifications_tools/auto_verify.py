from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from compare_snapshots import compare, write_metrics_log
from verification_plots import plot_per_body_dashboard
from generate_globular_cluster import generate_m13_plummer
from generate_large_system import generate_disk_system
from run_cpp_simulation import run_cpp
from run_rebound_reference import run_rebound
from verify_io import sample_system, write_system_txt, read_system_txt


def run_solar_case(
    exe: Path,
    input_txt: Path,
    years: float,
    dt: float,
    algo: int,
    *,
    save_interval_steps: int,
) -> dict:
    out_base = f"cpp_solar_t{years}_dt{dt}"
    t0 = time.time()
    cpp_final = run_cpp(
        exe,
        input_txt,
        out_base,
        dt,
        years,
        render_mode=2,
        algo=algo,
        save_interval_steps=save_interval_steps,
    )
    cpp_elapsed = time.time() - t0

    ref_file = Path("data") / f"rebound_solar_t{years}.txt"
    t1 = time.time()
    run_rebound(input_txt, ref_file, years, "ias15")
    ref_elapsed = time.time() - t1

    metrics = compare(cpp_final, ref_file.resolve())
    return {
        "case": "solar",
        "save_interval_steps": save_interval_steps,
        "cpp_seconds": cpp_elapsed,
        "rebound_seconds": ref_elapsed,
        "cpp_final": str(cpp_final),
        "rebound_final": str(ref_file.resolve()),
        "metrics": metrics,
    }


def run_large_case(
    exe: Path,
    years: float,
    dt: float,
    count: int,
    sample_size: int,
    seed: int,
    algo: int,
    preset: str = "disk",
    *,
    save_interval_steps: int,
) -> dict:
    if preset == "m13":
        full_input = Path("data") / f"m13_large_{count}_seed{seed}.txt"
    else:
        full_input = Path("data") / f"large_{count}_seed{seed}.txt"
    if not full_input.exists():
        if preset == "m13":
            rows = generate_m13_plummer(count, seed)
            write_system_txt(
                full_input,
                rows,
                header=f"M13 Plummer model N={count} seed={seed}",
            )
        else:
            rows = generate_disk_system(count, seed)
            write_system_txt(full_input, rows, header=f"Large synthetic disk count={count}")
    full_base = f"cpp_large_full_n{count}_t{years}_dt{dt}"

    t0 = time.time()
    cpp_full = run_cpp(
        exe,
        full_input,
        full_base,
        dt,
        years,
        render_mode=2,
        algo=algo,
        save_interval_steps=save_interval_steps,
    )
    full_elapsed = time.time() - t0

    sample_input = Path("data") / f"large_sample_{sample_size}_from_{count}.txt"
    sample_rows = sample_system(read_system_txt(full_input), sample_size, seed)
    write_system_txt(sample_input, sample_rows, header=f"Sample {sample_size} from {count}")
    sample_base_cpp = f"cpp_large_sample_n{sample_size}_t{years}_dt{dt}"
    sample_cpp = run_cpp(
        exe,
        sample_input,
        sample_base_cpp,
        dt,
        years,
        render_mode=2,
        algo=algo,
        save_interval_steps=save_interval_steps,
    )

    rebound_out = Path("data") / f"rebound_large_sample_n{sample_size}_t{years}.txt"
    t1 = time.time()
    run_rebound(sample_input, rebound_out, years, "whfast")
    rebound_elapsed = time.time() - t1

    metrics = compare(sample_cpp, rebound_out.resolve())
    return {
        "case": "large",
        "preset": preset,
        "save_interval_steps": save_interval_steps,
        "full_n": count,
        "sample_n": sample_size,
        "cpp_full_seconds": full_elapsed,
        "rebound_sample_seconds": rebound_elapsed,
        "cpp_full_final": str(cpp_full),
        "cpp_sample_final": str(sample_cpp),
        "rebound_sample_final": str(rebound_out.resolve()),
        "metrics_sample": metrics,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Automatic C++ vs REBOUND verification. "
        "Default: solar-style run (--input). Use --large for big-N workflow."
    )
    parser.add_argument("--exe", required=True, help="Path to nBodySim / Simulate executable.")
    parser.add_argument(
        "--large",
        action="store_true",
        help="Large-N mode: generate or reuse disk/m13 TXT, full C++ run + sampled REBOUND compare.",
    )
    parser.add_argument("--input", default="data/solar_nasa.txt", help="Input TXT (solar path; default).")
    parser.add_argument("--years", type=float, default=10.0)
    parser.add_argument("--dt", type=float, default=1.0 / 3652.5)
    parser.add_argument("--algo", type=int, default=1, choices=[1, 2, 3], help="C++ algorithm menu id.")
    parser.add_argument(
        "--save-interval-steps",
        type=int,
        default=10,
        metavar="N",
        help="Passed to Simulate: save HDF5 frame every N integration steps (same as main.cpp prompt).",
    )
    parser.add_argument("--count", type=int, default=100000, help="--large only: bodies (without Sun for disk).")
    parser.add_argument(
        "--preset",
        choices=["disk", "m13"],
        default="disk",
        help="--large only: disk or M13 Plummer.",
    )
    parser.add_argument("--sample-size", type=int, default=2000, help="--large only.")
    parser.add_argument("--seed", type=int, default=42, help="--large only.")
    parser.add_argument("--report", default="data/verification_report.json")
    parser.add_argument(
        "--plot-dir",
        type=str,
        default="",
        help="If set, write per_body_comparison.png (and use with metrics log).",
    )
    parser.add_argument(
        "--metrics-log",
        type=str,
        default="",
        help="If set, write human-readable per-body metrics to this path.",
    )
    args = parser.parse_args()

    exe = Path(args.exe)
    save_n = max(1, int(args.save_interval_steps))

    if args.large:
        result = run_large_case(
            exe,
            args.years,
            args.dt,
            args.count,
            args.sample_size,
            args.seed,
            args.algo,
            preset=args.preset,
            save_interval_steps=save_n,
        )
    else:
        result = run_solar_case(
            exe,
            Path(args.input),
            args.years,
            args.dt,
            args.algo,
            save_interval_steps=save_n,
        )

    out_path = Path(args.report)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"Saved report -> {out_path}")

    metrics_key = "metrics_sample" if args.large else "metrics"
    m = result.get(metrics_key)
    if isinstance(m, dict) and m.get("per_body"):
        if args.metrics_log:
            write_metrics_log(m, Path(args.metrics_log))
            print(f"Saved metrics log -> {args.metrics_log}")
        if args.plot_dir:
            pdir = Path(args.plot_dir)
            pdir.mkdir(parents=True, exist_ok=True)
            plot_per_body_dashboard(m, pdir / "per_body_comparison.png")
            print(f"Saved plot -> {pdir / 'per_body_comparison.png'}")


if __name__ == "__main__":
    main()
