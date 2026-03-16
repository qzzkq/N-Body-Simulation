import argparse
import subprocess
import sys
import os
import math
import time
import matplotlib.pyplot as plt

try:
    import nasa_api
    import rebound_sim
except ImportError as e:
    print(f"Ошибка импорта: {e}. Убедитесь, что nasa_api.py и rebound_sim.py лежат в этой же папке.")
    sys.exit(1)

def read_positions(filename):
    pos = {}
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 6:
                name = parts[0]
                pos[name] = (float(parts[3]), float(parts[4]), float(parts[5]))
    return pos

def calculate_mean_error(cpp_pos, rebound_pos):
    total_error = 0.0
    count = 0
    
    for name, r_coords in rebound_pos.items():
        if name in cpp_pos:
            c_coords = cpp_pos[name]
            diff_dist = math.sqrt(sum((r - c)**2 for r, c in zip(r_coords, c_coords)))
            ref_dist = math.sqrt(sum(r**2 for r in r_coords))
            
            if ref_dist > 0:
                total_error += (diff_dist / ref_dist) * 100.0
                count += 1
                
    return (total_error / count) if count > 0 else 0.0

def main():
    parser = argparse.ArgumentParser(description="Бенчмарк N-Body: С++ vs REBOUND")
    parser.add_argument("exe_dir", help="Путь к папке, где лежит исполняемый файл (Simulate.exe)")
    args = parser.parse_args()
    
    exe_path = os.path.join(args.exe_dir, "Simulate.exe")
    if not os.path.exists(exe_path):
        exe_path = os.path.join(args.exe_dir, "Simulate")
        if not os.path.exists(exe_path):
            print(f"Ошибка: Исполняемый файл не найден в {args.exe_dir}")
            sys.exit(1)

    os.makedirs("data", exist_ok=True)
    input_txt = "data/system_past.txt"
    rebound_txt = "data/rebound_truth.txt"
    target_years = 1.0

    print("=== 1. Получение данных от NASA ===")
    nasa_api.fetch_and_save(nasa_api.past_start, nasa_api.past_stop, input_txt)

    print(f"\n=== 2. Просчет эталона через REBOUND на {target_years} лет ===")
    rebound_sim.run_rebound_simulation(input_txt, rebound_txt, target_years)
    rebound_positions = read_positions(rebound_txt)
    print(f"Эталон просчитан. Загружено объектов: {len(rebound_positions)}")

    dt_values = [0.01, 0.005, 0.001, 0.0005, 0.0001]
    errors = []

    print("\n=== 3. Запуск C++ программы с разными dt ===")
    for dt in dt_values:
        out_base = f"test_dt_{dt}"
        cpp_out_txt = f"data/{out_base}_final.txt"

        answers = [
            "1",          
            "1",          
            str(dt),      
            "1",          
            input_txt,   
            out_base,    
            "0",          
            str(target_years) 
        ]
        
        input_str = "\n".join(answers) + "\n"
        print(f"Запуск симуляции (dt = {dt})... ", end="", flush=True)
        
        start_time = time.time()
        process = subprocess.Popen(
            [exe_path], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8'
        )
        process.communicate(input=input_str)
        
        if process.returncode == 0 and os.path.exists(cpp_out_txt):
            cpp_pos = read_positions(cpp_out_txt)
            err = calculate_mean_error(cpp_pos, rebound_positions)
            errors.append(err)
            print(f"Готово ({time.time() - start_time:.1f} сек). Ошибка: {err:.6f}%")
        else:
            print("ОШИБКА выполнения!")
            errors.append(None)

    valid_dt = [dt for dt, err in zip(dt_values, errors) if err is not None]
    valid_err = [err for err in errors if err is not None]

    print("\n=== 4. Построение графика ===")
    plt.figure(figsize=(9, 6))
    plt.plot(valid_dt, valid_err, marker='o', linestyle='-', color='b', linewidth=2)
    plt.xscale('log')
    plt.yscale('log')

    plt.gca().invert_xaxis()
    
    plt.xlabel("Шаг интеграции dt (годы)")
    plt.ylabel("Ср. относительная ошибка позиций (%)")
    plt.title("Зависимость точности C++ симуляции от шага dt\n(В сравнении с REBOUND IAS15)")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    
    chart_file = "accuracy_chart.png"
    plt.savefig(chart_file)
    print(f"График сохранен в файл: {chart_file}")
    plt.show()

if __name__ == "__main__":
    main()