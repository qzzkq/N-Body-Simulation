from astroquery.jplhorizons import Horizons
from datetime import datetime, timedelta, timezone

PLANETS = {
    "Sun":     {"id": "10",  "mass": 1.98847e30, "density": 1410, "color": "1.0 0.9 0.1"},
    "Mercury": {"id": "199", "mass": 3.3011e23,  "density": 5429, "color": "0.5 0.5 0.5"},
    "Venus":   {"id": "299", "mass": 4.8675e24,  "density": 5243, "color": "0.9 0.8 0.6"},
    "Earth":   {"id": "399", "mass": 5.9723e24,  "density": 5514, "color": "0.2 0.5 1.0"},
    "Mars":    {"id": "499", "mass": 6.4171e23,  "density": 3934, "color": "1.0 0.3 0.1"},
    "Jupiter": {"id": "599", "mass": 1.8982e27,  "density": 1326, "color": "0.8 0.7 0.6"},
    "Saturn":  {"id": "699", "mass": 5.6834e26,  "density": 687,  "color": "0.9 0.8 0.5"},
    "Uranus":  {"id": "799", "mass": 8.6810e25,  "density": 1270, "color": "0.6 0.8 0.9"},
    "Neptune": {"id": "899", "mass": 1.0241e26,  "density": 1638, "color": "0.2 0.3 0.8"}
}

now = datetime.now(timezone.utc)

past_start = (now - timedelta(days=365)).strftime('%Y-%m-%d')
past_stop  = (now - timedelta(days=364)).strftime('%Y-%m-%d')


def fetch_and_save(start_str, stop_str, filename):
    print(f"Скачивание данных на {start_str}...")
    
    with open(filename, 'w') as f:
        f.write(f"# Солнечная система на {start_str}\n")
        f.write("# Format: Name Mass Density Px Py Pz Vx Vy Vz R G B\n")
        
        for name, data in PLANETS.items():
            obj = Horizons(id=data["id"], location="@10", epochs={'start': start_str, 'stop': stop_str, 'step': '1d'})
            vec = obj.vectors()[0] 
            
            AU_TO_METERS = 1.495978707e11
            px = vec['x'] * AU_TO_METERS
            py = vec['y'] * AU_TO_METERS
            pz = vec['z'] * AU_TO_METERS
            
            AU_PER_DAY_TO_MS = (1.495978707e11) / (24 * 60 * 60)
            vx = vec['vx'] * AU_PER_DAY_TO_MS
            vy = vec['vy'] * AU_PER_DAY_TO_MS
            vz = vec['vz'] * AU_PER_DAY_TO_MS
            
            if name == "Sun":
                px, py, pz, vx, vy, vz = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

            line = f"{name} {data['mass']} {data['density']} {px:.3f} {py:.3f} {pz:.3f} {vx:.3f} {vy:.3f} {vz:.3f} {data['color']}\n"
            f.write(line)
            
    print(f"Сохранено в {filename}\n")

fetch_and_save(past_start, past_stop, "system_past.txt")
