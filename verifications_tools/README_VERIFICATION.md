# Verification Tools

Набор скриптов для автоматической проверки вашего симулятора против REBOUND, без изменений основного C++ кода.

## Скрипты

- `generate_solar_nasa_cli.py`  
  Скачивает снимок Солнечной системы из JPL Horizons (NASA API) и сохраняет в TXT.

- `generate_large_system.py`  
  Генерирует синтетическую большую систему (до сотен тысяч тел) в формате TXT, который понимает `Simulate`.

- `run_cpp_simulation.py`  
  Неинтерактивно запускает ваш `Simulate` с заданными параметрами (dt, years, algorithm, input).

- `run_rebound_reference.py`  
  Строит эталонный снимок через REBOUND из того же входного TXT.

- `compare_snapshots.py`  
  Сравнивает два финальных снимка: ошибки по позиции/скорости, mean/max, per-body.

- `auto_verify.py`  
  Общий оркестратор:
  - `solar`: C++ vs REBOUND на одной системе;
  - `large`: полный прогон C++ на большом N + сравнение C++ vs REBOUND на детерминированной подвыборке.

## Примеры

1) Скачать solar snapshot:

```bash
python generate_solar_nasa_cli.py --date 2025-03-01 --out data/solar_nasa.txt
```

2) Проверка solar:

```bash
python auto_verify.py \
  --exe ../nBodySim \
  --mode solar \
  --input data/solar_nasa.txt \
  --years 100 \
  --dt 0.0002737850787132101 \
  --algo 1 \
  --report data/solar_report.json
```

3) Проверка large N:

```bash
python auto_verify.py \
  --exe ../nBodySim \
  --mode large \
  --years 1.0 \
  --dt 0.001 \
  --algo 2 \
  --count 200000 \
  --sample-size 3000 \
  --seed 42 \
  --report data/large_report.json
```

## Зависимости Python

- `rebound`
- `astroquery` (только для `generate_solar_nasa_cli.py`)

Установка:

```bash
pip install rebound astroquery
```

