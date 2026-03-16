import argparse
import subprocess
import sys
import os
import math
import time
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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
    target_years = 100.0

    print("=== 1. Получение данных от NASA ===")
    nasa_api.fetch_and_save(nasa_api.past_start, nasa_api.past_stop, input_txt)

    print(f"\n=== 2. Просчет эталона через REBOUND на {target_years} лет ===")
    rebound_sim.run_rebound_simulation(input_txt, rebound_txt, target_years)
    rebound_positions = read_positions(rebound_txt)
    print(f"Эталон просчитан. Загружено объектов: {len(rebound_positions)}")

    dt_values = [1.0, 1.0/365.25, 1.0/3652.5, 1.0/36525.0, 1.0/365250.0]
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

    print("\n=== 4. Построение графика ===")
    valid_dt = [dt for dt, err in zip(dt_values, errors) if err is not None]
    valid_err = [err for err in errors if err is not None]
    
    plt.figure(figsize=(10, 6))
    
    safe_err = [max(e, 1e-10) for e in valid_err]
    
    plt.plot(valid_dt, safe_err, marker='o', markersize=8, linestyle='-', color='#1f77b4', linewidth=2)
    
    plt.xscale('log')
    plt.yscale('log')
    ax = plt.gca()
    ax.invert_xaxis()
    
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f"{y:g}%"))
    
    dt_labels = []
    for dt in valid_dt:
        if abs(dt - 1.0) < 1e-6: 
            dt_labels.append("1.0")
        elif abs(dt - 1.0/365.0) < 1e-6: 
            dt_labels.append("1/365")
        elif abs(dt - 1.0/365.25) < 1e-6: 
            dt_labels.append("1/365.25")
        elif abs(dt - 1.0/3652.5) < 1e-6: 
            dt_labels.append("1/3652.5")
        elif abs(dt - 1.0/36525.0) < 1e-9: 
            dt_labels.append("1/36525")
        elif abs(dt - 1.0/365250.0) < 1e-10: 
            dt_labels.append("1/365250")
        else: 
            dt_labels.append(f"{dt:.1e}")

    plt.xticks(valid_dt, dt_labels, rotation=0, fontsize=10)

    for x, y in zip(valid_dt, safe_err):
        if y <= 1e-10:
            label_text = "0.0%"
        else:
            label_text = f"{y:.5f}%" if y > 0.00001 else f"{y:.1e}%"
            
        plt.annotate(label_text, (x, y), textcoords="offset points", xytext=(0, 10), 
                     ha='center', fontsize=9, weight='bold')

    plt.xlabel("Шаг интеграции dt (в годах)", fontsize=12)
    plt.ylabel("Ср. относительная ошибка позиций", fontsize=12)
    plt.title("Зависимость точности C++ симуляции от шага dt", fontsize=13)
    
    plt.grid(True, which="major", ls="-", color='gray', alpha=0.5)
    plt.grid(True, which="minor", ls="--", color='gray', alpha=0.2)
    

    plt.margins(y=0.2)
    
    plt.tight_layout()
    chart_file = "accuracy_chart.png"
    plt.savefig(chart_file, dpi=300)
    print(f"График сохранен в файл: {chart_file}")
    plt.show()  

if __name__ == "__main__":
    main()