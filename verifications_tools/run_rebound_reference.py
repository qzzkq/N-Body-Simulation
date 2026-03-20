from __future__ import annotations

import argparse
import math
from pathlib import Path

import rebound

from verify_io import AU_METERS, SOLAR_MASS_KG, YEAR_SECONDS, BodyRecord, read_system_txt, write_system_txt


def run_rebound(input_txt: Path, output_txt: Path, years: float, integrator: str = "ias15") -> None:
    sim = rebound.Simulation()
    sim.G = 4.0 * math.pi * math.pi
    sim.integrator = integrator

    rows = read_system_txt(input_txt)
    names: list[str] = []
    extra: dict[str, tuple[float, float, float, float]] = {}
    for r in rows:
        names.append(r.name)
        extra[r.name] = (r.density, r.r, r.g, r.b)
        sim.add(
            m=r.mass_kg / SOLAR_MASS_KG,
            x=r.px_m / AU_METERS,
            y=r.py_m / AU_METERS,
            z=r.pz_m / AU_METERS,
            vx=r.vx_ms * (YEAR_SECONDS / AU_METERS),
            vy=r.vy_ms * (YEAR_SECONDS / AU_METERS),
            vz=r.vz_ms * (YEAR_SECONDS / AU_METERS),
        )

    com = sim.com()
    for p in sim.particles:
        p.vx -= com.vx
        p.vy -= com.vy
        p.vz -= com.vz

    sim.integrate(float(years))

    out_rows: list[BodyRecord] = []
    for i, p in enumerate(sim.particles):
        n = names[i]
        dens, rr, gg, bb = extra[n]
        out_rows.append(
            BodyRecord(
                n,
                p.m * SOLAR_MASS_KG,
                dens,
                p.x * AU_METERS,
                p.y * AU_METERS,
                p.z * AU_METERS,
                p.vx * (AU_METERS / YEAR_SECONDS),
                p.vy * (AU_METERS / YEAR_SECONDS),
                p.vz * (AU_METERS / YEAR_SECONDS),
                rr,
                gg,
                bb,
            )
        )
    write_system_txt(output_txt, out_rows, header=f"REBOUND {integrator} T={years} years")


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute reference snapshot with REBOUND.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--years", type=float, required=True)
    parser.add_argument("--integrator", default="ias15")
    args = parser.parse_args()

    run_rebound(Path(args.input), Path(args.output), args.years, args.integrator)
    print(f"Saved reference -> {args.output}")


if __name__ == "__main__":
    main()

