from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path

SOLAR_MASS_KG = 1.98847e30
# Астрономическая единица (м), IAU-совместимое значение — как в rebound_sim / JPL.
AU_METERS = 1.495978707e11


@dataclass
class BodyRecord:
    name: str
    mass_kg: float
    density: float
    px_m: float
    py_m: float
    pz_m: float
    vx_ms: float
    vy_ms: float
    vz_ms: float
    r: float
    g: float
    b: float


def read_system_txt(path: Path) -> list[BodyRecord]:
    rows: list[BodyRecord] = []
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            rows.append(
                BodyRecord(
                    name=parts[0],
                    mass_kg=float(parts[1]),
                    density=float(parts[2]),
                    px_m=float(parts[3]),
                    py_m=float(parts[4]),
                    pz_m=float(parts[5]),
                    vx_ms=float(parts[6]),
                    vy_ms=float(parts[7]),
                    vz_ms=float(parts[8]),
                    r=float(parts[9]),
                    g=float(parts[10]),
                    b=float(parts[11]),
                )
            )
    return rows


def write_system_txt(path: Path, rows: list[BodyRecord], header: str = "") -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        if header.strip():
            for hline in header.strip().split("\n"):
                line = hline if hline.startswith("#") else f"# {hline}"
                f.write(f"{line}\n")
        f.write("# Name Mass Density Px Py Pz Vx Vy Vz R G B\n")
        for r in rows:
            f.write(
                f"{r.name} {r.mass_kg:.15e} {r.density:.15e} "
                f"{r.px_m:.15e} {r.py_m:.15e} {r.pz_m:.15e} "
                f"{r.vx_ms:.15e} {r.vy_ms:.15e} {r.vz_ms:.15e} "
                f"{r.r} {r.g} {r.b}\n"
            )


def sample_system(rows: list[BodyRecord], sample_size: int, seed: int) -> list[BodyRecord]:
    if sample_size >= len(rows):
        return list(rows)
    rng = random.Random(seed)
    idx = rng.sample(range(len(rows)), sample_size)
    idx.sort()
    return [rows[i] for i in idx]
