import rebound
import math
import sys

SOLAR_MASS_KG = 1.98847e30
AU_METERS = 1.495978707e11
YEAR_SECONDS = 3.15576e7

def run_rebound_simulation(input_file, output_file, years_to_simulate):
    sim = rebound.Simulation()
    sim.G = 4.0 * math.pi**2 

    extra_data = {}
    particle_names = []

    try:
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 12: continue
                
                name = parts[0]
                mass_kg = float(parts[1])
                density = float(parts[2])
                px, py, pz = float(parts[3]), float(parts[4]), float(parts[5])
                vx, vy, vz = float(parts[6]), float(parts[7]), float(parts[8])
                r, g, b = parts[9], parts[10], parts[11]

                extra_data[name] = (density, r, g, b)
                particle_names.append(name)

                # Конвертация в IAU
                m = mass_kg / SOLAR_MASS_KG
                x, y, z = px / AU_METERS, py / AU_METERS, pz / AU_METERS
                v_x = vx * (YEAR_SECONDS / AU_METERS)
                v_y = vy * (YEAR_SECONDS / AU_METERS)
                v_z = vz * (YEAR_SECONDS / AU_METERS)

                sim.add(m=m, x=x, y=y, z=z, vx=v_x, vy=v_y, vz=v_z)
    except FileNotFoundError:
        print(f"Ошибка: Файл {input_file} не найден.")
        return

    print(f"Интегрируем {input_file} на {years_to_simulate} лет...")
    sim.integrate(float(years_to_simulate))

    with open(output_file, 'w') as f:
        f.write(f"# REBOUND IAS15 Output | T = {years_to_simulate} years\n")
        f.write("# Name Mass Density Px Py Pz Vx Vy Vz R G B\n")
        
        for i, p in enumerate(sim.particles):
            name = particle_names[i]
            d, r, g, b = extra_data[name]
            
            # Обратная конвертация
            f.write(f"{name} {p.m*SOLAR_MASS_KG:.15e} {d:.15e} "
                    f"{p.x*AU_METERS:.15e} {p.y*AU_METERS:.15e} {p.z*AU_METERS:.15e} "
                    f"{p.vx*(AU_METERS/YEAR_SECONDS):.15e} {p.vy*(AU_METERS/YEAR_SECONDS):.15e} {p.vz*(AU_METERS/YEAR_SECONDS):.15e} "
                    f"{r} {g} {b}\n")
    print(f"Результат сохранен в {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Использование: python rebound_sim.py <input.txt> <output.txt> <years>")
        print("Пример: python rebound_sim.py system_past.txt rebound_1y.txt 1.0")
    else:
        run_rebound_simulation(sys.argv[1], sys.argv[2], sys.argv[3])