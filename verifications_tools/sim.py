from __future__ import annotations
import math, subprocess
from pathlib import Path
from snapshot import AU_METERS, SOLAR_MASS_KG, YEAR_SECONDS


def run_cpp(
    exe: Path,
    input_txt: Path,
    out_base: str,
    dt: float,
    years: float,
    *,
    algo: int = 1,
    save_interval_steps: int = 10,
    timeout: float | None = None,
) -> Path:
    exe = Path(exe).resolve()
    cwd = exe.parent
    cmd = [
        str(exe),
        "--realtime", "0",
        "--algo",     str(int(algo)),
        "--dt",       format(dt, ".17g"),
        "--source",   "1",
        "--txt",      str(Path(input_txt).resolve()),
        "--output",   out_base,
        "--save_every", str(max(1, int(save_interval_steps))),
        "--target",   str(float(years)),
    ]

    import sys
    proc = subprocess.Popen(cmd, cwd=str(cwd),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stderr_chunks: list[bytes] = []
    try:
        # C++ progress uses \r without \n — читаем мелкими чанками для live-вывода.
        while True:
            chunk = proc.stdout.read(64)
            if not chunk:
                break
            sys.stdout.buffer.write(chunk)
            sys.stdout.flush()
        proc.wait(timeout=timeout)
        stderr_chunks.append(proc.stderr.read())
        sys.stdout.write("\n")
        sys.stdout.flush()
    except subprocess.TimeoutExpired:
        proc.kill()
        raise
    if proc.returncode != 0:
        raise RuntimeError(
            f"Simulate exited {proc.returncode}.\n"
            + (b"".join(stderr_chunks).decode(errors="replace") or "")[-6000:]
        )
    final = cwd / "data" / f"{out_base}_final.txt"
    if not final.is_file():
        raise FileNotFoundError(f"Expected output missing: {final}")
    return final


def run_rebound(
    input_file: Path,
    output_file: Path,
    years: float,
    integrator: str = "ias15",
) -> None:
    import rebound

    input_file  = Path(input_file)
    output_file = Path(output_file)

    sim = rebound.Simulation()
    sim.G = 4.0 * math.pi ** 2
    sim.integrator = integrator

    names: list[str] = []
    extras: dict[str, tuple] = {}

    with input_file.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < 12:
                continue
            name = p[0]
            names.append(name)
            extras[name] = (float(p[2]), p[9], p[10], p[11])
            m  = float(p[1]) / SOLAR_MASS_KG
            x  = float(p[3]) / AU_METERS
            y  = float(p[4]) / AU_METERS
            z  = float(p[5]) / AU_METERS
            vx = float(p[6]) * (YEAR_SECONDS / AU_METERS)
            vy = float(p[7]) * (YEAR_SECONDS / AU_METERS)
            vz = float(p[8]) * (YEAR_SECONDS / AU_METERS)
            sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)

    com = sim.com()
    for p in sim.particles:
        p.vx -= com.vx
        p.vy -= com.vy
        p.vz -= com.vz

    sim.integrate(float(years))

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8") as f:
        f.write(f"# REBOUND {integrator} | T={years} yr\n")
        f.write("# Name Mass Density Px Py Pz Vx Vy Vz R G B\n")
        for i, part in enumerate(sim.particles):
            d, r, g, b = extras[names[i]]
            f.write(
                f"{names[i]} {part.m * SOLAR_MASS_KG:.15e} {d:.15e} "
                f"{part.x * AU_METERS:.15e} "
                f"{part.y * AU_METERS:.15e} "
                f"{part.z * AU_METERS:.15e} "
                f"{part.vx * (AU_METERS / YEAR_SECONDS):.15e} "
                f"{part.vy * (AU_METERS / YEAR_SECONDS):.15e} "
                f"{part.vz * (AU_METERS / YEAR_SECONDS):.15e} "
                f"{r} {g} {b}\n"
            )
