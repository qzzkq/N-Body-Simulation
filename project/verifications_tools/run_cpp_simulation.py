from __future__ import annotations

import subprocess
from pathlib import Path


def run_cpp(
    exe: Path,
    input_txt: Path,
    out_base: str,
    dt: float,
    years: float,
    *,
    render_mode: int = 2,
    algo: int = 1,
    save_interval_steps: int = 10,
    timeout: float | None = None,
) -> Path:
    """
    Non-interactive run of Simulate: stdin lines must match project/src/main.cpp in order:
      1) render mode (0 cubes / 1 spheres / 2 points)
      2) algorithm (1 bruteforce / 2 Barnes-Hut / 3 CUDA if built)
      3) dt (years per step, string with full precision)
      4) initial conditions source: 1 = TXT
      5) path to input TXT
      6) output basename (no extension; used for data/<base>.h5 and data/<base>_final.txt)
      7) HDF5 frame save interval N (write every N-th integration step)
      8) realtime: 0 = bake to target time
      9) target simulation time (years)
    """
    exe = Path(exe).resolve()
    cwd = exe.parent
    input_abs = str(Path(input_txt).resolve())
    # Enough precision for dt like 1/3652.5
    dt_str = format(dt, ".17g")
    lines = [
        str(int(render_mode)),
        str(int(algo)),
        dt_str,
        "1",
        input_abs,
        out_base,
        str(max(1, int(save_interval_steps))),
        "0",
        str(float(years)),
    ]
    stdin_text = "\n".join(lines) + "\n"
    r = subprocess.run(
        [str(exe)],
        input=stdin_text,
        text=True,
        cwd=str(cwd),
        capture_output=True,
        timeout=timeout,
    )
    if r.returncode != 0:
        err = (r.stderr or "") + (r.stdout or "")
        raise RuntimeError(
            f"C++ simulator exited {r.returncode}.\n"
            f"--- stderr/stdout ---\n{err[-8000:]}"
        )
    final_path = cwd / "data" / f"{out_base}_final.txt"
    if not final_path.is_file():
        raise FileNotFoundError(f"Expected output missing: {final_path}")
    return final_path
