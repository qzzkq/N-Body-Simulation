from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import random
from typing import Iterable


SOLAR_MASS_KG = 1.98847e30
AU_METERS = 1.495978707e11
YEAR_SECONDS = 3.15576e7


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
    r: float = 1.0
    g: float = 1.0
    b: float = 1.0


def read_system_txt(path: str | Path) -> list[BodyRecord]:
    records: list[BodyRecord] = []
    with Path(path).open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 12:
                continue
            records.append(
                BodyRecord(
                    name=p[0],
                    mass_kg=float(p[1]),
                    density=float(p[2]),
                    px_m=float(p[3]),
                    py_m=float(p[4]),
                    pz_m=float(p[5]),
                    vx_ms=float(p[6]),
                    vy_ms=float(p[7]),
                    vz_ms=float(p[8]),
                    r=float(p[9]),
                    g=float(p[10]),
                    b=float(p[11]),
                )
            )
    return records


def write_system_txt(path: str | Path, records: Iterable[BodyRecord], header: str = "") -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        if header:
            f.write(f"# {header}\n")
        f.write("# Name Mass Density Px Py Pz Vx Vy Vz R G B\n")
        for rec in records:
            f.write(
                f"{rec.name} {rec.mass_kg:.15e} {rec.density:.15e} "
                f"{rec.px_m:.15e} {rec.py_m:.15e} {rec.pz_m:.15e} "
                f"{rec.vx_ms:.15e} {rec.vy_ms:.15e} {rec.vz_ms:.15e} "
                f"{rec.r:.6f} {rec.g:.6f} {rec.b:.6f}\n"
            )


def sample_system(records: list[BodyRecord], sample_size: int, seed: int = 42) -> list[BodyRecord]:
    if len(records) <= sample_size:
        return list(records)
    sun = [r for r in records if r.name.lower() == "sun"]
    others = [r for r in records if r.name.lower() != "sun"]
    rng = random.Random(seed)
    picked = rng.sample(others, k=max(0, sample_size - len(sun)))
    return sun + picked

