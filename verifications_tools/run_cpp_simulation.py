from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def run_cpp(exe_path: Path, input_txt: Path, out_base: str, dt: float, years: float, render_mode: int = 1, algo: int = 1) -> Path:
    # menu answers for Simulate:
    # render mode, algorithm, dt, source=TXT(1), path, out name, realtime=0, target time
    answers = [
        str(render_mode),
        str(algo),
        f"{dt}",
        "1",
        str(input_txt),
        out_base,
        "0",
        f"{years}",
    ]
    proc = subprocess.run(
        [str(exe_path)],
        input="\n".join(answers) + "\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Simulate failed ({proc.returncode}).\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return Path("data") / f"{out_base}_final.txt"


def main() -> None:
    parser = argparse.ArgumentParser(description="Run C++ Simulate non-interactively.")
    parser.add_argument("--exe", required=True, help="Path to Simulate executable.")
    parser.add_argument("--input", required=True, help="Input TXT path.")
    parser.add_argument("--out-base", required=True, help="Output base name (without extension).")
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--years", type=float, required=True)
    parser.add_argument("--render-mode", type=int, default=1, choices=[0, 1, 2])
    parser.add_argument("--algo", type=int, default=1, choices=[1, 2, 3])
    args = parser.parse_args()

    out = run_cpp(Path(args.exe), Path(args.input), args.out_base, args.dt, args.years, args.render_mode, args.algo)
    print(f"Done. Final snapshot: {out}")


if __name__ == "__main__":
    main()

