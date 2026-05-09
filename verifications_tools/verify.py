from __future__ import annotations
import argparse, json, time
from pathlib import Path

from snapshot import read, write, sample, Body
from sim import run_cpp, run_rebound
from compare import compare, write_log, plot_per_body, plot_dt_sweep
from generators import generate_disk, generate_m13


def _solar(exe, input_txt, years, dt, algo, save_n, com_align, plot_dir, log_path):
    out_base = f"cpp_solar_t{years}_dt{dt}"
    t0 = time.time()
    cpp_final = run_cpp(exe, input_txt, out_base, dt, years, algo=algo,
                        save_interval_steps=save_n)
    cpp_t = time.time() - t0

    ref = Path("data") / f"rebound_solar_t{years}.txt"
    t0 = time.time()
    run_rebound(input_txt, ref, years)
    reb_t = time.time() - t0

    m = compare(cpp_final, ref, com_align=com_align)

    if log_path:
        write_log(m, Path(log_path))
    if plot_dir and m.get("per_body"):
        plot_per_body(m, Path(plot_dir) / "per_body.png")

    return {"case": "solar", "cpp_seconds": cpp_t, "rebound_seconds": reb_t,
            "cpp_final": str(cpp_final), "rebound_final": str(ref), "metrics": m}


def _large(exe, years, dt, count, sample_size, seed, algo, preset,
           save_n, com_align, plot_dir, log_path):
    data_dir = Path("data")
    prefix = "m13" if preset == "m13" else "disk"
    full_input = data_dir / f"{prefix}_{count}_s{seed}.txt"
    if not full_input.exists():
        rows = generate_m13(count, seed) if preset == "m13" else generate_disk(count, seed)
        write(full_input, rows, header=f"{preset} N={count} seed={seed}")

    t0 = time.time()
    cpp_full = run_cpp(exe, full_input, f"cpp_{prefix}_n{count}_t{years}_dt{dt}",
                       dt, years, algo=algo, save_interval_steps=save_n)
    full_t = time.time() - t0

    samp_input = data_dir / f"{prefix}_sample_{sample_size}_from_{count}_s{seed}.txt"
    samp_rows = sample(read(full_input), sample_size, seed)
    write(samp_input, samp_rows, header=f"sample {sample_size} from {count}")

    samp_cpp = run_cpp(exe, samp_input, f"cpp_{prefix}_sample_n{sample_size}_t{years}_dt{dt}",
                       dt, years, algo=algo, save_interval_steps=save_n)

    ref = data_dir / f"rebound_{prefix}_sample_n{sample_size}_t{years}.txt"
    t0 = time.time()
    run_rebound(samp_input, ref, years, integrator="whfast")
    reb_t = time.time() - t0

    m = compare(samp_cpp, ref, com_align=com_align)

    if log_path:
        write_log(m, Path(log_path))
    if plot_dir and m.get("per_body"):
        plot_per_body(m, Path(plot_dir) / "per_body.png")

    return {"case": "large", "preset": preset, "full_n": count, "sample_n": sample_size,
            "cpp_full_seconds": full_t, "rebound_sample_seconds": reb_t,
            "cpp_full_final": str(cpp_full), "metrics_sample": m}


def _benchmark(exe, input_txt, years, algo, save_n, com_align, out_dir):
    dt_values = [1.0, 1/365.25, 1/3652.5, 1/36525.0, 1/365250.0]

    ref = Path("data") / f"rebound_bench_t{years}.txt"
    if not ref.exists():
        run_rebound(input_txt, ref, years)

    results, abs_km = [], []
    for dt in dt_values:
        out_base = f"bench_dt{dt:.2e}_t{years}"
        try:
            cpp_f = run_cpp(exe, input_txt, out_base, dt, years,
                            algo=algo, save_interval_steps=save_n)
            m = compare(cpp_f, ref, com_align=com_align)
            err = m.get("mean_abs_pos_error_km")
            abs_km.append(err)
            results.append({"dt": dt, "mean_abs_pos_error_km": err,
                            "mean_rel_pos_error": m.get("mean_rel_pos_error")})
            print(f"  dt={dt:.2e}  ΔR={err:.3g} km  rel={m.get('mean_rel_pos_error'):.3e}")
        except Exception as e:
            print(f"  dt={dt:.2e}  FAILED: {e}")
            abs_km.append(None)
            results.append({"dt": dt, "error": str(e)})

    out_dir = Path(out_dir or "data")
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dt_sweep(dt_values, abs_km, out_dir / "benchmark_dt_sweep.png")
    return {"case": "benchmark", "results": results}


def main() -> None:
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("mode", choices=["solar", "large", "benchmark"])
    ap.add_argument("--exe",     required=True, help="Path to Simulate binary")
    ap.add_argument("--input",   default="data/solar_nasa.txt")
    ap.add_argument("--years",   type=float, default=10.0)
    ap.add_argument("--dt",      type=float, default=1/3652.5)
    ap.add_argument("--algo",    type=int,   default=1, choices=[1, 2, 3, 4])
    ap.add_argument("--save-n",  type=int,   default=10, dest="save_n")
    ap.add_argument("--no-com-align", action="store_false", dest="com_align",
                    help="Disable COM-frame alignment before error computation")
    ap.add_argument("--preset",      choices=["disk", "m13"], default="disk")
    ap.add_argument("--count",       type=int, default=100000)
    ap.add_argument("--sample-size", type=int, default=3000, dest="sample_size")
    ap.add_argument("--seed",        type=int, default=42)
    ap.add_argument("--report",   default="data/report.json")
    ap.add_argument("--log",      default="",  help="Per-body metrics text log path")
    ap.add_argument("--plot-dir", default="",  dest="plot_dir")

    args = ap.parse_args()
    exe = Path(args.exe)

    if args.mode == "solar":
        result = _solar(exe, Path(args.input), args.years, args.dt, args.algo,
                        args.save_n, args.com_align,
                        args.plot_dir or None, args.log or None)
    elif args.mode == "large":
        result = _large(exe, args.years, args.dt, args.count, args.sample_size,
                        args.seed, args.algo, args.preset, args.save_n,
                        args.com_align, args.plot_dir or None, args.log or None)
    else:
        result = _benchmark(exe, Path(args.input), args.years, args.algo,
                            args.save_n, args.com_align, args.plot_dir or "data")

    out = Path(args.report)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"\nReport → {out}")


if __name__ == "__main__":
    main()
