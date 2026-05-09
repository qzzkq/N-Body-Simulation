from __future__ import annotations
import argparse, math, random
from pathlib import Path
from snapshot import Body, write, AU_METERS, SOLAR_MASS_KG, YEAR_SECONDS


_PLANETS = {
    "Sun":     {"id": "10",  "mass": 1.98847e30, "density": 1410, "color": "1.0 0.9 0.1"},
    "Mercury": {"id": "199", "mass": 3.3011e23,  "density": 5429, "color": "0.5 0.5 0.5"},
    "Venus":   {"id": "299", "mass": 4.8675e24,  "density": 5243, "color": "0.9 0.8 0.6"},
    "Earth":   {"id": "399", "mass": 5.9723e24,  "density": 5514, "color": "0.2 0.5 1.0"},
    "Mars":    {"id": "499", "mass": 6.4171e23,  "density": 3934, "color": "1.0 0.3 0.1"},
    "Jupiter": {"id": "599", "mass": 1.8982e27,  "density": 1326, "color": "0.8 0.7 0.6"},
    "Saturn":  {"id": "699", "mass": 5.6834e26,  "density": 687,  "color": "0.9 0.8 0.5"},
    "Uranus":  {"id": "799", "mass": 8.6810e25,  "density": 1270, "color": "0.6 0.8 0.9"},
    "Neptune": {"id": "899", "mass": 1.0241e26,  "density": 1638, "color": "0.2 0.3 0.8"},
}

def fetch_solar(date_str: str, out: Path) -> list[Body]:
    from astroquery.jplhorizons import Horizons
    from datetime import datetime, timedelta
    stop = (datetime.strptime(date_str, "%Y-%m-%d") + timedelta(days=1)).strftime("%Y-%m-%d")
    rows: list[Body] = []
    for name, meta in _PLANETS.items():
        r, g, b = (float(x) for x in meta["color"].split())
        if name == "Sun":
            rows.append(Body(name, meta["mass"], meta["density"],
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, r, g, b))
            continue
        obj = Horizons(id=meta["id"], location="@10",
                       epochs={"start": date_str, "stop": stop, "step": "1d"})
        v = obj.vectors()[0]
        AU_DAY_TO_MS = AU_METERS / 86400.0
        rows.append(Body(
            name=name, mass_kg=meta["mass"], density=meta["density"],
            px_m=float(v["x"]) * AU_METERS,
            py_m=float(v["y"]) * AU_METERS,
            pz_m=float(v["z"]) * AU_METERS,
            vx_ms=float(v["vx"]) * AU_DAY_TO_MS,
            vy_ms=float(v["vy"]) * AU_DAY_TO_MS,
            vz_ms=float(v["vz"]) * AU_DAY_TO_MS,
            r=r, g=g, b=b,
        ))
    write(out, rows, header=f"Solar system {date_str} (heliocentric, JPL Horizons)")
    return rows


def generate_disk(count: int, seed: int) -> list[Body]:
    rng = random.Random(seed)
    SUN_MASS = 1.98847e30
    EARTH_MASS = 5.9723e24
    G_SI = 6.674e-11
    rows: list[Body] = [
        Body("Sun", SUN_MASS, 1410.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.9, 0.1)
    ]
    for i in range(count):
        r_au  = rng.uniform(0.3, 40.0)
        r_m   = r_au * AU_METERS
        theta = rng.uniform(0, 2 * math.pi)
        phi   = rng.gauss(0, math.radians(3))
        x = r_m * math.cos(theta) * math.cos(phi)
        y = r_m * math.sin(theta) * math.cos(phi)
        z = r_m * math.sin(phi)
        v_circ = math.sqrt(G_SI * SUN_MASS / r_m)
        v_frac = rng.gauss(1.0, 0.02)
        vx = -v_circ * v_frac * math.sin(theta)
        vy =  v_circ * v_frac * math.cos(theta)
        vz =  rng.gauss(0, v_circ * 0.01)
        mass = EARTH_MASS * rng.uniform(0.1, 5.0)
        rows.append(Body(f"P{i}", mass, 3000.0, x, y, z, vx, vy, vz, 0.4, 0.6, 1.0))
    return rows


_M13_MASS_SOLAR  = 6e5
_M13_SCALE_AU    = 8000.0

def generate_m13(count: int, seed: int) -> list[Body]:
    rng = random.Random(seed)
    M_total = _M13_MASS_SOLAR * SOLAR_MASS_KG
    a = _M13_SCALE_AU * AU_METERS
    G_SI = 6.674e-11
    m_star = M_total / count
    density = 5000.0

    rows: list[Body] = []
    for i in range(count):
        # Plummer CDF: r = a / sqrt(u^{-2/3} - 1)
        u = rng.random()
        r = a / math.sqrt(u ** (-2 / 3) - 1.0)
        ct   = rng.uniform(-1, 1)
        st   = math.sqrt(max(0.0, 1 - ct * ct))
        phi  = rng.uniform(0, 2 * math.pi)
        x = r * st * math.cos(phi)
        y = r * st * math.sin(phi)
        z = r * ct

        # Aarseth, Henon & Wielen 1974: rejection sampling по q = v / v_esc.
        v_esc = math.sqrt(2.0 * G_SI * M_total / math.sqrt(r * r + a * a))
        while True:
            q = rng.random()
            g = q * q * (1 - q * q) ** 3.5
            if 0.1 * rng.random() < g:
                break
        v = q * v_esc
        ct2 = rng.uniform(-1, 1)
        st2 = math.sqrt(max(0.0, 1 - ct2 * ct2))
        psi = rng.uniform(0, 2 * math.pi)
        vx = v * st2 * math.cos(psi)
        vy = v * st2 * math.sin(psi)
        vz = v * ct2
        rows.append(Body(f"S{i}", m_star, density, x, y, z, vx, vy, vz, 1.0, 0.9, 0.7))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)

    s = sub.add_parser("solar")
    s.add_argument("--date", default="2025-03-01")
    s.add_argument("--out",  default="data/solar_nasa.txt")

    d = sub.add_parser("disk")
    d.add_argument("--count", type=int, default=10000)
    d.add_argument("--seed",  type=int, default=42)
    d.add_argument("--out",   default="data/disk.txt")

    m = sub.add_parser("m13")
    m.add_argument("--count", type=int, default=100000)
    m.add_argument("--seed",  type=int, default=42)
    m.add_argument("--out",   default="data/m13.txt")

    args = ap.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    if args.mode == "solar":
        fetch_solar(args.date, out)
    elif args.mode == "disk":
        rows = generate_disk(args.count, args.seed)
        write(out, rows, header=f"Synthetic disk count={args.count} seed={args.seed}")
    elif args.mode == "m13":
        rows = generate_m13(args.count, args.seed)
        write(out, rows, header=f"M13 Plummer N={args.count} seed={args.seed}")
    print(f"Written → {out}")

if __name__ == "__main__":
    main()
