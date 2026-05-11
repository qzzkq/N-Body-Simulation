# TXT snapshot: Name Mass_kg Density_kg_m3 Px_m Py_m Pz_m Vx_ms Vy_ms Vz_ms R G B
from __future__ import annotations
import random
from dataclasses import dataclass
from pathlib import Path

SOLAR_MASS_KG = 1.98847e30
AU_METERS     = 1.495978707e11
YEAR_SECONDS  = 3.15576e7


@dataclass
class Body:
    name: str
    mass_kg: float
    density: float
    px_m: float; py_m: float; pz_m: float
    vx_ms: float; vy_ms: float; vz_ms: float
    r: float; g: float; b: float


def read(path: Path) -> list[Body]:
    rows: list[Body] = []
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) < 12:
                continue
            rows.append(Body(
                name=p[0],
                mass_kg=float(p[1]), density=float(p[2]),
                px_m=float(p[3]), py_m=float(p[4]), pz_m=float(p[5]),
                vx_ms=float(p[6]), vy_ms=float(p[7]), vz_ms=float(p[8]),
                r=float(p[9]), g=float(p[10]), b=float(p[11]),
            ))
    return rows


def write(path: Path, rows: list[Body], header: str = "") -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for line in (header or "").strip().split("\n"):
            if line.strip():
                f.write(line if line.startswith("#") else f"# {line}\n")
        f.write("# Name Mass Density Px Py Pz Vx Vy Vz R G B\n")
        for b in rows:
            f.write(
                f"{b.name} {b.mass_kg:.15e} {b.density:.15e} "
                f"{b.px_m:.15e} {b.py_m:.15e} {b.pz_m:.15e} "
                f"{b.vx_ms:.15e} {b.vy_ms:.15e} {b.vz_ms:.15e} "
                f"{b.r} {b.g} {b.b}\n"
            )


def sample(rows: list[Body], n: int, seed: int) -> list[Body]:
    if n >= len(rows):
        return list(rows)
    rng = random.Random(seed)
    idx = sorted(rng.sample(range(len(rows)), n))
    return [rows[i] for i in idx]
