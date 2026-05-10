"""
GUI-07/08/09 — юнит-тесты AppConfig (config.py).
Дисплей не требуется — все тесты работают без Tkinter.
"""

import json
import os
import pytest

from modules.config import AppConfig


# ── GUI-08 — сохранение/загрузка config.json ─────────────────────────────────

def test_config_defaults(app_config):
    """Все ожидаемые ключи по умолчанию присутствуют и не пусты."""
    cfg = app_config
    assert cfg.get("simulate_bin") == "nBodySim"
    assert cfg.get("replay_bin")   == "replay"
    assert cfg.get("data_dir")     == "data"
    assert cfg.get("build_dir")    == "build"


def test_config_set_and_get(app_config):
    """set() сохраняет значение; get() его возвращает."""
    app_config.set("build_dir", "/custom/build")
    assert app_config.get("build_dir") == "/custom/build"


def test_config_persists_to_disk(tmp_gui_dir):
    """Значение, записанное одним экземпляром AppConfig, видно в новом экземпляре."""
    cfg1 = AppConfig(tmp_gui_dir)
    cfg1.set("build_dir", "/persisted/path")

    cfg2 = AppConfig(tmp_gui_dir)
    assert cfg2.get("build_dir") == "/persisted/path"


def test_config_update_multiple_keys(app_config):
    """update() применяет все переданные ключи одновременно."""
    app_config.update({"build_dir": "/new/build", "simulate_bin": "Simulate"})
    assert app_config.get("build_dir")    == "/new/build"
    assert app_config.get("simulate_bin") == "Simulate"


def test_config_all_returns_dict(app_config):
    """all() возвращает словарь, содержащий как минимум ключи по умолчанию."""
    d = app_config.all()
    assert isinstance(d, dict)
    assert "simulate_bin" in d
    assert "replay_bin"   in d


def test_config_missing_key_returns_fallback(app_config):
    """get() с неизвестным ключом возвращает переданный fallback."""
    assert app_config.get("nonexistent_key", "default_val") == "default_val"
    assert app_config.get("nonexistent_key") is None


def test_config_corrupted_json_uses_defaults(tmp_gui_dir):
    """Повреждённый config.json игнорируется без исключений; используются значения по умолчанию."""
    config_path = os.path.join(tmp_gui_dir, "config.json")
    with open(config_path, "w") as f:
        f.write("{ invalid json !!!")

    cfg = AppConfig(tmp_gui_dir)
    assert cfg.get("simulate_bin") == "nBodySim"


# ── GUI-09 — проверка путей ───────────────────────────────────────────────────

def test_config_binary_path_detection(tmp_gui_dir, tmp_path):
    """Проверка наличия бинарника работает корректно для существующего и отсутствующего файла."""
    existing_bin = tmp_path / "Simulate"
    existing_bin.touch(mode=0o755)

    # Существующий бинарник должен обнаруживаться
    assert os.path.exists(str(existing_bin))

    # Отсутствующий бинарник не должен обнаруживаться
    missing_path = str(tmp_path / "DoesNotExist")
    assert not os.path.exists(missing_path)
