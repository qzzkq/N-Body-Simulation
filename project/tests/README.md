# Тесты N-Body Simulation

Два независимых набора тестов: C++ (Google Test) и Python (pytest).

---

## C++ тесты (Google Test)

### Структура

| Файл | ID | Описание |
|---|---|---|
| `test_phy.cpp` | PHY-01…08 | Физика: гравитационное ускорение, сохранение энергии/импульса, орбиты |
| `test_alg.cpp` | ALG-01…09 | Алгоритмы: BruteForce vs Barnes-Hut CPU vs CUDA, стабильность |
| `test_io.cpp`  | IO-01…10  | Ввод/вывод: HDF5 round-trip, TXT-парсинг, атрибуты файла симуляции |

### Сборка и запуск

```bash
cd project
cmake -B build_tests -DCMAKE_BUILD_TYPE=Release
cmake --build build_tests --target nbody_tests -j$(nproc)
cd build_tests && ctest --output-on-failure
```

Тесты собираются вместе с основным проектом (`add_subdirectory(tests)` уже добавлен в `CMakeLists.txt`).

### Сборка с поддержкой CUDA (ALG-02)

```bash
cmake -B build_tests -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=ON
cmake --build build_tests --target nbody_tests -j$(nproc)
```

Тест `ALG-02` (`BarnesHutCPU_vs_CUDA_50Bodies_100Steps`) компилируется только при `USE_CUDA=ON`. Без него все остальные тесты работают как обычно.

### Запуск одного теста

```bash
cd build_tests
./nbody_tests --gtest_filter="AlgTest.BruteForceVsBH_CPU_50Bodies_100Steps"
./nbody_tests --gtest_filter="PHY*"
./nbody_tests --gtest_filter="IOTest*"
```

### Зависимости

- CMake ≥ 3.21
- GCC/Clang с C++17
- HDF5 (libhdf5-dev)
- GLM (libglm-dev)
- Google Test — скачивается автоматически через FetchContent если не найден в системе

---

## Python тесты GUI (pytest)

### Структура

| Файл | ID | Описание |
|---|---|---|
| `gui/test_baker.py` | GUI-01…06, 10 | Логика CalculatePage: поиск бинарника, построение команды, мониторинг процесса |
| `gui/test_config.py`| GUI-07…09    | AppConfig: defaults, set/get, персистентность, битый JSON |
| `gui/conftest.py`   | —            | Фикстуры `app_config`, `tmp_gui_dir`; метка `@requires_display` |

Тесты, требующие Tkinter (`@requires_display`), пропускаются в headless-окружениях (нет `$DISPLAY` / `$WAYLAND_DISPLAY`).

### Запуск

```bash
# Использовать venv проекта (рекомендуется)
/home/qzxc/study/N-Body-Simulation/project/gui/venv/bin/python3 -m pytest \
    project/tests/gui/ -v

# Или установить зависимости в системный Python
pip install pytest customtkinter
pytest project/tests/gui/ -v
```

Из директории `project/tests/gui/` напрямую:

```bash
cd project/tests/gui
pytest -v
```

### Зависимости

```
pytest>=7.0
customtkinter>=5.0
```

Файл `gui/requirements.txt` содержит эти зависимости.

---

## Быстрый прогон всего

```bash
# C++ тесты
cmake -B build_tests && cmake --build build_tests --target nbody_tests -j$(nproc)
cd build_tests && ctest --output-on-failure && cd ..

# Python тесты
/home/qzxc/study/N-Body-Simulation/project/gui/venv/bin/python3 \
    -m pytest project/tests/gui/ -v
```
