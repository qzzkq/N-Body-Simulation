from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

from verify_io import BodyRecord, write_system_txt


def generate_disk_system(count: int, seed: int) -> list[BodyRecord]:
    rng = random.Random(seed)
    rows: list[BodyRecord] = []
    rows.append(
        BodyRecord(
            name="Sun",
            mass_kg=1.98847e30,
            density=1410.0,
            px_m=0.0,
            py_m=0.0,
            pz_m=0.0,
            vx_ms=0.0,
            vy_ms=0.0,
            vz_ms=0.0,
            r=1.0,
            g=0.9,
            b=0.1,
        )
    )

    # Broad disk in meters: ~0.5 to 100 AU
    r_min = 0.5 * 1.495978707e11
    r_max = 100.0 * 1.495978707e11
    mu = 6.67430e-11 * 1.98847e30
    for i in range(count):
        rr = r_min + (r_max - r_min) * math.sqrt(rng.random())
        ang = rng.random() * 2.0 * math.pi
        height = (rng.random() - 0.5) * 0.02 * rr
        px, py, pz = rr * math.cos(ang), height, rr * math.sin(ang)
        vcirc = math.sqrt(mu / max(rr, 1.0))
        jitter = 0.95 + 0.1 * rng.random()
        vx, vy, vz = -vcirc * math.sin(ang) * jitter, 0.0, vcirc * math.cos(ang) * jitter
        rows.append(
            BodyRecord(
                name=f"Body_{i+1}",
                mass_kg=1e20 * (0.5 + rng.random()),
                density=1400.0,
                px_m=px,
                py_m=py,
                pz_m=pz,
                vx_ms=vx,
                vy_ms=vy,
                vz_ms=vz,
                r=0.7,
                g=0.8,
                b=1.0,
            )
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate large synthetic N-body disk.")
    parser.add_argument("--count", type=int, default=100000, help="Number of orbiting bodies (without Sun).")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="data/large_system.txt")
    args = parser.parse_args()

    rows = generate_disk_system(args.count, args.seed)
    write_system_txt(Path(args.out), rows, header=f"Synthetic large disk count={args.count}, seed={args.seed}")
    print(f"Saved {len(rows)} bodies -> {args.out}")


if __name__ == "__main__":
    main()

