from __future__ import annotations

"""
Генерация большой системы по физически мотивированной модели реального шарового скопления.

Пресет `m13` — шаровое скопление Messier 13 (NGC 6205, «Великий шар Геркулеса»):
  ~300k–500k звёзд по наблюдениям; масса порядка ~6×10^5 M☉; расстояние ~7.4 kpc.
  Здесь: N тел с равными массами, изотропная модель Пламмера с масштабом a и массой M,
  согласованными с типичными оценками для M13 (см. обзоры по шаровым скоплениям / Baumgardt & Hilker).

Это не каталог Gaia по отдельным звёздам, а стандартная N-body постановка задачи
для плотного скопления — то, что обычно делают в статьях по прямому N-body.
"""

import argparse
import math
import random
from pathlib import Path

from verify_io import BodyRecord, SOLAR_MASS_KG, write_system_txt

G_SI = 6.67430e-11
PC_METERS = 3.085677581e16

# M13: типичные порядки величин (литература / обзоры GC)
M13_TOTAL_MASS_SOLAR = 6.0e5
# Plummer a: из r_h ≈ 3 pc для GC и связи r_h ≈ 0.766 a для модели Пламмера
M13_PLUMMER_A_PC = 4.0

# Разносим оттенки по кругу: соседние индексы не совпадают по цвету.
_GOLDEN = (math.sqrt(5.0) - 1.0) * 0.5


def _hsv_to_rgb(h: float, s: float, v: float) -> tuple[float, float, float]:
    h = h % 1.0
    i = int(h * 6.0)
    f = h * 6.0 - i
    p = v * (1.0 - s)
    q = v * (1.0 - f * s)
    t = v * (1.0 - (1.0 - f) * s)
    i = i % 6
    if i == 0:
        return v, t, p
    if i == 1:
        return q, v, p
    if i == 2:
        return p, v, t
    if i == 3:
        return p, q, v
    if i == 4:
        return t, p, v
    return v, p, q


def _star_rgb_from_index(i: int) -> tuple[float, float, float]:
    """Уникальный набор (r,g,b) на тело: hue по золотому сечению + лёгкий разброс S/V."""
    hue = (i * _GOLDEN) % 1.0
    sat = 0.42 + 0.55 * (0.5 + 0.5 * math.sin(i * 0.731))
    val = 0.72 + 0.26 * (0.5 + 0.5 * math.cos(i * 0.419))
    rr, gg, bb = _hsv_to_rgb(hue, min(1.0, sat), min(1.0, val))
    return max(0.05, rr), max(0.05, gg), max(0.05, bb)


def _gauss(rng: random.Random, sigma: float) -> float:
    u1 = max(rng.random(), 1e-300)
    u2 = rng.random()
    return sigma * math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def _sample_plummer_radius(rng: random.Random, a_m: float) -> float:
    u = rng.random()
    u = min(max(u, 1e-300), 1.0 - 1e-300)
    return a_m / math.sqrt(u ** (-2.0 / 3.0) - 1.0)


def _velocity_sigma_sq(r: float, a: float, m_total_kg: float) -> float:
    """Дисперсия одной декартовой компоненты скорости для изотропного Пламмера (равновесие)."""
    s = a * a + r * r
    return (G_SI * m_total_kg / (12.0 * math.sqrt(s))) * (a * a / s)


def generate_m13_plummer(count: int, seed: int) -> list[BodyRecord]:
    """
    N звёзд с равной массой, суммарная масса = M13_TOTAL_MASS_SOLAR M_sun,
    положения — Пламмер, скорости — изотропные с дисперсией из модели Пламмера.
    """
    if count < 1:
        raise ValueError("count must be >= 1")

    rng = random.Random(seed)
    m_total_kg = M13_TOTAL_MASS_SOLAR * SOLAR_MASS_KG
    m_star = m_total_kg / float(count)
    a_m = M13_PLUMMER_A_PC * PC_METERS

    px = [0.0] * count
    py = [0.0] * count
    pz = [0.0] * count
    vx = [0.0] * count
    vy = [0.0] * count
    vz = [0.0] * count

    for i in range(count):
        r = _sample_plummer_radius(rng, a_m)
        cos_theta = 2.0 * rng.random() - 1.0
        sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
        phi = rng.random() * 2.0 * math.pi
        x = r * sin_theta * math.cos(phi)
        y = r * sin_theta * math.sin(phi)
        z = r * cos_theta
        px[i], py[i], pz[i] = x, y, z

        sig2 = _velocity_sigma_sq(r, a_m, m_total_kg)
        sigma = math.sqrt(max(sig2, 0.0))
        vx[i] = _gauss(rng, sigma)
        vy[i] = _gauss(rng, sigma)
        vz[i] = _gauss(rng, sigma)

    # Центрируем по положению и импульсу (как в rebound_sim для солнечной системы)
    mx = sum(px) / count
    my = sum(py) / count
    mz = sum(pz) / count
    mvx = sum(vx) / count
    mvy = sum(vy) / count
    mvz = sum(vz) / count
    for i in range(count):
        px[i] -= mx
        py[i] -= my
        pz[i] -= mz
        vx[i] -= mvx
        vy[i] -= mvy
        vz[i] -= mvz

    rows: list[BodyRecord] = []
    for i in range(count):
        rr, gg, bb = _star_rgb_from_index(i)
        rows.append(
            BodyRecord(
                name=f"M13_{i+1:06d}",
                mass_kg=m_star,
                density=1410.0,
                px_m=px[i],
                py_m=py[i],
                pz_m=pz[i],
                vx_ms=vx[i],
                vy_ms=vy[i],
                vz_ms=vz[i],
                r=rr,
                g=gg,
                b=bb,
            )
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Сгенерировать N-body модель шарового скопления M13 (Пламмер, ~реальные M и a)."
    )
    parser.add_argument("--count", type=int, default=100_000, help="Число звёзд (равной массы).")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="data/m13_globular_100k.txt")
    args = parser.parse_args()

    rows = generate_m13_plummer(args.count, args.seed)
    hdr = (
        f"M13 (NGC 6205) Plummer model: N={args.count}, M_tot~{M13_TOTAL_MASS_SOLAR:.3e} M_sun, "
        f"a={M13_PLUMMER_A_PC} pc, seed={args.seed}"
    )
    write_system_txt(Path(args.out), rows, header=hdr)
    print(f"Saved {len(rows)} bodies -> {args.out}")


if __name__ == "__main__":
    main()
