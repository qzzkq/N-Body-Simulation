from __future__ import annotations

import math
from pathlib import Path

import rebound

SOLAR_MASS_KG = 1.98847e30
AU_METERS = 1.495978707e11
YEAR_SECONDS = 3.15576e7


def run_rebound(
    input_file: Path,
    output_file: Path,
    years_to_simulate: float,
    integrator: str = "ias15",
) -> None:
    """Same text format as rebound_sim.py; integrator e.g. ias15, whfast."""
    input_file = Path(input_file)
    output_file = Path(output_file)
    sim = rebound.Simulation()
    sim.G = 4.0 * math.pi**2
    sim.integrator = integrator

    extra_data: dict[str, tuple[float, str, str, str]] = {}
    particle_names: list[str] = []

    with input_file.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            name = parts[0]
            mass_kg = float(parts[1])
            density = float(parts[2])
            px, py, pz = float(parts[3]), float(parts[4]), float(parts[5])
            vx, vy, vz = float(parts[6]), float(parts[7]), float(parts[8])
            r, g, b = parts[9], parts[10], parts[11]

            extra_data[name] = (density, r, g, b)
            particle_names.append(name)

            m = mass_kg / SOLAR_MASS_KG
            x, y, z = px / AU_METERS, py / AU_METERS, pz / AU_METERS
            v_x = vx * (YEAR_SECONDS / AU_METERS)
            v_y = vy * (YEAR_SECONDS / AU_METERS)
            v_z = vz * (YEAR_SECONDS / AU_METERS)
            sim.add(m=m, x=x, y=y, z=z, vx=v_x, vy=v_y, vz=v_z)

    com = sim.com()
    for p in sim.particles:
        p.vx -= com.vx
        p.vy -= com.vy
        p.vz -= com.vz

    sim.integrate(float(years_to_simulate))

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8") as f:
        f.write(f"# REBOUND {integrator} | T = {years_to_simulate} years\n")
        f.write("# Name Mass Density Px Py Pz Vx Vy Vz R G B\n")
        for i, p in enumerate(sim.particles):
            name = particle_names[i]
            d, r, g, b = extra_data[name]
            f.write(
                f"{name} {p.m * SOLAR_MASS_KG:.15e} {d:.15e} "
                f"{p.x * AU_METERS:.15e} {p.y * AU_METERS:.15e} {p.z * AU_METERS:.15e} "
                f"{p.vx * (AU_METERS / YEAR_SECONDS):.15e} "
                f"{p.vy * (AU_METERS / YEAR_SECONDS):.15e} "
                f"{p.vz * (AU_METERS / YEAR_SECONDS):.15e} "
                f"{r} {g} {b}\n"
            )
