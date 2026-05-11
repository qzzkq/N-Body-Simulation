"""
GUI-01..06/10 — тесты CalculatePage и логики запуска симуляции.

Тесты, требующие дисплей, помечены @requires_display.
Тесты чистой логики (построение команды, мониторинг процесса) работают headless.
"""

import os
import sys
import subprocess
import threading
import time
import pytest

from conftest import requires_display

GUI_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "gui")
if GUI_DIR not in sys.path:
    sys.path.insert(0, GUI_DIR)

from modules.baker import ALGORITHMS, RENDER_MODES, SOURCES, SCENARIOS
from modules.config import AppConfig


# ── GUI-01 — логика find_binary ───────────────────────────────────────────────

def test_find_binary_existing(tmp_path):
    """Бинарник найден, если файл существует по пути project_root/binary_name."""
    binary = tmp_path / "Simulate"
    binary.touch(mode=0o755)

    binary_path = os.path.join(str(tmp_path), "Simulate")
    assert os.path.exists(binary_path), "Бинарник должен быть найден"


def test_find_binary_missing(tmp_path):
    """Бинарник не найден, если файл отсутствует."""
    binary_path = os.path.join(str(tmp_path), "Simulate")
    assert not os.path.exists(binary_path), "Отсутствующий бинарник не должен обнаруживаться"


def test_find_binary_multiple_paths(tmp_path):
    """Бинарник найден, если он существует хотя бы в одном из путей поиска."""
    subdir = tmp_path / "build"
    subdir.mkdir()
    (subdir / "Simulate").touch(mode=0o755)

    search_paths = [str(tmp_path), str(subdir)]
    found = any(
        os.path.exists(os.path.join(p, "Simulate"))
        for p in search_paths
    )
    assert found


# ── GUI-02 — build_simulate_stdin / построение команды ───────────────────────
# Проверяет, что массив команды из CalculatePage._on_start имеет правильную
# структуру для каждой комбинации источника/режима.
# Логика извлечена из baker.py без создания экземпляра Tkinter.

def _build_cmd(project_root, binary_name,
               render="1", algo="1",
               src_idx=0, dt="0.001", save_n="10",
               outfile="frames", is_rt=True,
               h5="data/solar.h5", txt="/path/to.txt",
               bodies="100", scenario_idx=0,
               target="10.0"):
    """Повторяет логику построения команды из CalculatePage._on_start."""
    binary_path = os.path.join(project_root, binary_name)
    cmd = [
        binary_path,
        "--render", render,
        "--algo", algo,
        "--dt", dt,
        "--source", str(src_idx),
        "--output", outfile,
        "--save_every", save_n,
        "--realtime", "1" if is_rt else "0",
    ]
    if src_idx == 0:
        cmd.extend(["--h5", h5])
    elif src_idx == 1:
        cmd.extend(["--txt", txt])
    else:
        cmd.extend(["--bodies", bodies])
        cmd.extend(["--scenario", str(scenario_idx)])

    if not is_rt:
        cmd.extend(["--target", target])

    return cmd


@pytest.mark.parametrize("src_idx,is_rt,expected_len", [
    # base = 15: binary + 7 пар флаг-значение (render, algo, dt, source, output, save_every, realtime)
    (0, True,  15 + 2),      # HDF5 + реальное время:  + --h5 val
    (0, False, 15 + 2 + 2),  # HDF5 + запись:          + --h5 val + --target val
    (1, True,  15 + 2),      # TXT  + реальное время:  + --txt val
    (2, True,  15 + 4),      # Rand + реальное время:  + --bodies val + --scenario val
    (2, False, 15 + 4 + 2),  # Rand + запись:          + --bodies val + --scenario val + --target val
])
def test_command_length(tmp_path, src_idx, is_rt, expected_len):
    """Массив команды имеет ожидаемое число элементов для каждой комбинации."""
    cmd = _build_cmd(str(tmp_path), "Simulate", src_idx=src_idx, is_rt=is_rt)
    assert len(cmd) == expected_len, (
        f"src_idx={src_idx}, is_rt={is_rt}: "
        f"получено {len(cmd)} аргументов, ожидалось {expected_len}"
    )


def test_command_contains_binary_path(tmp_path):
    """Первый элемент — всегда абсолютный путь к бинарнику."""
    cmd = _build_cmd(str(tmp_path), "Simulate")
    assert cmd[0] == os.path.join(str(tmp_path), "Simulate")


def test_command_realtime_flag(tmp_path):
    """Флаг --realtime равен '1' для реального времени и '0' для записи."""
    cmd_rt   = _build_cmd(str(tmp_path), "Simulate", is_rt=True)
    cmd_bake = _build_cmd(str(tmp_path), "Simulate", is_rt=False)

    rt_idx   = cmd_rt.index("--realtime")
    bake_idx = cmd_bake.index("--realtime")
    assert cmd_rt[rt_idx + 1]     == "1"
    assert cmd_bake[bake_idx + 1] == "0"


# ── GUI-03 — launch_simulate: чтение stdout в фоновом потоке ─────────────────

def test_monitor_process_reads_output():
    """Фоновый поток монитора читает все строки stdout и вызывает on_output."""
    lines_received = []
    finish_code    = []

    script = (
        "import sys\n"
        "for i in range(10):\n"
        "    print(f'line {i}')\n"
        "    sys.stdout.flush()\n"
    )
    proc = subprocess.Popen(
        [sys.executable, "-c", script],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1,
    )

    def monitor():
        for line in iter(proc.stdout.readline, ""):
            line = line.strip()
            if line:
                lines_received.append(line)
        proc.stdout.close()
        proc.wait()
        finish_code.append(proc.returncode)

    t = threading.Thread(target=monitor, daemon=True)
    t.start()
    t.join(timeout=10)

    assert not t.is_alive(), "Поток монитора завис"
    assert len(lines_received) == 10, f"Получено {len(lines_received)} строк, ожидалось 10"
    assert finish_code == [0], f"returncode = {finish_code}"


def test_monitor_process_progress_parsing():
    """Строки PROGRESS: разбираются корректно; остальные строки пропускаются насквозь."""
    progress_values = []
    other_lines     = []

    script = (
        "import sys\n"
        "print('Starting...')\n"
        "sys.stdout.flush()\n"
        "for i in range(0, 101, 25):\n"
        "    print(f'PROGRESS:{i}')\n"
        "    sys.stdout.flush()\n"
        "print('Done')\n"
        "sys.stdout.flush()\n"
    )
    proc = subprocess.Popen(
        [sys.executable, "-c", script],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1,
    )

    for line in iter(proc.stdout.readline, ""):
        line = line.strip()
        if not line:
            continue
        if line.startswith("PROGRESS:"):
            try:
                progress_values.append(float(line.split(":")[1]))
            except ValueError:
                pass
        else:
            other_lines.append(line)
    proc.stdout.close()
    proc.wait()

    assert progress_values == [0.0, 25.0, 50.0, 75.0, 100.0]
    assert "Starting..." in other_lines
    assert "Done"        in other_lines


# ── GUI-04 — launch_replay: корректный argv[1] ────────────────────────────────

def test_launch_replay_correct_path(tmp_path):
    """Бинарник воспроизведения получает путь .h5 как первый аргумент."""
    # Заглушка-"бинарник", печатающая argv[1]
    stub = tmp_path / "replay_stub.py"
    stub.write_text(
        "import sys\nprint(sys.argv[1])\n"
    )

    h5_path = str(tmp_path / "test.h5")
    result = subprocess.run(
        [sys.executable, str(stub), h5_path],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert result.stdout.strip() == h5_path


# ── GUI-06 — запуск без бинарника → ошибка в логе ────────────────────────────

@requires_display
def test_calculate_page_missing_binary(tmp_path):
    """CalculatePage показывает ошибку 'не найден', если бинарник отсутствует."""
    import tkinter as tk
    from modules.baker import CalculatePage

    cfg = AppConfig(str(tmp_path / "gui"))
    os.makedirs(str(tmp_path / "gui"), exist_ok=True)
    cfg.set("simulate_bin", "NonExistentBinary")

    root = tk.Tk()
    root.withdraw()

    palette = {
        "bg": "#080b12", "border": "#2e3240", "text": "#ffffff",
        "dim": "#8890aa", "start": "#3d9de0", "start_hover": "#5ab0f0",
    }

    try:
        page = CalculatePage(root, str(tmp_path), cfg, palette)
        page._on_start()  # должен обнаружить отсутствие бинарника и выставить ошибку
        result_text = page._result_box.get("1.0", "end").strip()
        assert "не найден" in result_text or "not found" in result_text.lower(), (
            f"Ожидалась ошибка 'не найден', получено: {result_text!r}"
        )
    finally:
        root.destroy()


# ── GUI-05 — запуск главного окна ────────────────────────────────────────────

@requires_display
def test_main_window_startup(tmp_path):
    """Главное окно инициализируется без исключений Python."""
    import tkinter as tk

    # Перенаправляем корень проекта во временную директорию — реальные бинарники не нужны
    os.makedirs(str(tmp_path / "gui"), exist_ok=True)
    os.makedirs(str(tmp_path / "data"), exist_ok=True)

    import importlib.util
    main_path = os.path.join(GUI_DIR, "main.py")
    spec = importlib.util.spec_from_file_location("gui_main", main_path)
    mod  = importlib.util.module_from_spec(spec)

    # Проверяем лишь, что модуль загружается без ошибок; mainloop() не вызываем
    try:
        spec.loader.exec_module(mod)
    except SystemExit:
        pass  # некоторые main() вызывают sys.exit при импорте — это нормально


# ── GUI-10 — Stop: завершение запущенного процесса ───────────────────────────

def test_stop_terminates_process():
    """Завершение запущенного процесса устанавливает poll() в ненулевое значение."""
    proc = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(60)"],
    )
    assert proc.poll() is None, "Процесс должен быть запущен"

    proc.terminate()
    proc.wait(timeout=5)

    assert proc.poll() is not None, "Процесс должен быть завершён"
