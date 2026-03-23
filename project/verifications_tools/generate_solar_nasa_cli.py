from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
from pathlib import Path

from astroquery.jplhorizons import Horizons

from verify_io import BodyRecord, AU_METERS, write_system_txt


PLANETS = {
    "Sun": {"id": "10", "mass": 1.98847e30, "density": 1410, "color": (1.0, 0.9, 0.1)},
    "Mercury": {"id": "199", "mass": 3.3011e23, "density": 5429, "color": (0.5, 0.5, 0.5)},
    "Venus": {"id": "299", "mass": 4.8675e24, "density": 5243, "color": (0.9, 0.8, 0.6)},
    "Earth": {"id": "399", "mass": 5.9723e24, "density": 5514, "color": (0.2, 0.5, 1.0)},
    "Mars": {"id": "499", "mass": 6.4171e23, "density": 3934, "color": (1.0, 0.3, 0.1)},
    "Jupiter": {"id": "599", "mass": 1.8982e27, "density": 1326, "color": (0.8, 0.7, 0.6)},
    "Saturn": {"id": "699", "mass": 5.6834e26, "density": 687, "color": (0.9, 0.8, 0.5)},
    "Uranus": {"id": "799", "mass": 8.6810e25, "density": 1270, "color": (0.6, 0.8, 0.9)},
    "Neptune": {"id": "899", "mass": 1.0241e26, "density": 1638, "color": (0.2, 0.3, 0.8)},
}


def fetch_solar_snapshot(date_utc: str) -> list[BodyRecord]:
    start = date_utc
    stop = (datetime.strptime(date_utc, "%Y-%m-%d") + timedelta(days=1)).strftime("%Y-%m-%d")
    au_day_to_ms = AU_METERS / 86400.0
    out: list[BodyRecord] = []
    for name, data in PLANETS.items():
        vec = Horizons(id=data["id"], location="@10", epochs={"start": start, "stop": stop, "step": "1d"}).vectors()[0]
        px, py, pz = float(vec["x"]) * AU_METERS, float(vec["y"]) * AU_METERS, float(vec["z"]) * AU_METERS
        vx, vy, vz = float(vec["vx"]) * au_day_to_ms, float(vec["vy"]) * au_day_to_ms, float(vec["vz"]) * au_day_to_ms
        if name == "Sun":
            px = py = pz = vx = vy = vz = 0.0
        r, g, b = data["color"]
        out.append(BodyRecord(name, data["mass"], data["density"], px, py, pz, vx, vy, vz, r, g, b))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Download Solar System snapshot from JPL Horizons.")
    parser.add_argument("--date", default=None, help="UTC date YYYY-MM-DD. Defaults to now-365d.")
    parser.add_argument("--out", default="data/solar_nasa.txt", help="Output TXT path.")
    args = parser.parse_args()

    if args.date is None:
        args.date = (datetime.now(timezone.utc) - timedelta(days=365)).strftime("%Y-%m-%d")
    rows = fetch_solar_snapshot(args.date)
    write_system_txt(Path(args.out), rows, header=f"Solar snapshot from JPL on {args.date}")
    print(f"Saved {len(rows)} bodies -> {args.out}")


if __name__ == "__main__":
    main()

