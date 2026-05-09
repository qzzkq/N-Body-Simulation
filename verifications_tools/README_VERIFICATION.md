# Verification Tools

5 модулей вместо 12. Запуск всегда из `verifications_tools/`.

## Структура

| Файл | Содержит |
|------|----------|
| `snapshot.py` | `Body` dataclass, `read/write/sample` — I/O примитивы |
| `generators.py` | `generate_disk`, `generate_m13`, `fetch_solar` + CLI |
| `sim.py` | `run_cpp`, `run_rebound` |
| `compare.py` | `compare` (с COM-коррекцией), `write_log`, `plot_per_body`, `plot_dt_sweep` |
| `verify.py` | Оркестратор: режимы `solar`, `large`, `benchmark` |

## COM-коррекция (включена по умолчанию)

Оба интегратора (C++ и REBOUND) переходят в систему центра масс через
вычитание COM-скорости, но из-за разницы в floating-point арифметике
они получают **чуть** разные COM-скорости. За 100 лет это даёт дрейф
COM-позиций на тысячи км — без коррекции он накладывается на реальную
ошибку и раздувает метрики.

`compare.py::compare(..., com_align=True)` перед вычислением ошибок
сдвигает каждый снимок в свой COM-кадр (`r_i → r_i - COM(snapshot)`).
Поле `com_pos_drift_km` в метриках показывает насколько далеко разошлись
COM-позиции двух снимков — это диагностика, а не ошибка.

Отключить: `--no-com-align`.

## Быстрый старт

```bash
cd verifications_tools

# 1. Сгенерировать входной файл (солнечная система)
python3 generators.py solar --date 2025-03-01 --out data/solar_nasa.txt

# 2. Прогнать сравнение (100 лет, алго 1=brute-force)
python3 verify.py solar \
  --exe ../project/nBodySim \
  --input data/solar_nasa.txt \
  --years 100 --dt 2.74e-4 --algo 1 \
  --log data/solar.log --plot-dir data/plots

# 3. Бенчмарк: точность vs dt
python verify.py benchmark \
  --exe ../project/build/Simulate \
  --input data/solar_nasa.txt \
  --years 100 --algo 1

# 4. Large-N (100k тел, M13 Plummer)
python verify.py large \
  --exe ../project/build/Simulate \
  --preset m13 --count 100000 --sample-size 3000 \
  --years 1 --algo 2
```

## Генераторы

```bash
python generators.py disk  --count 50000 --seed 42 --out data/disk.txt
python generators.py m13   --count 100000 --seed 42 --out data/m13.txt
python generators.py solar --date 2025-03-01 --out data/solar_nasa.txt
```

## Зависимости

```bash
pip install rebound astroquery matplotlib
```
