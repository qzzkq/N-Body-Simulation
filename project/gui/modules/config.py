"""
Простое хранилище настроек в gui/config.json
"""

import json
import os


_DEFAULTS = {
    "build_dir":        "build",
    "simulate_bin":     "nBodySim",
    "replay_bin":       "replay",
    "data_dir":         "data",
    "default_dt":       "0.000273785",
    "default_save_n":   "10",
    "default_out_file": "frames",
    "default_bodies":   "100",
    "default_target":   "10.0",
}


class AppConfig:
    def __init__(self, gui_dir: str):
        self._path = os.path.join(gui_dir, "config.json")
        self._data: dict = dict(_DEFAULTS)
        self._load()

    # ── публичный API ──────────────────────────────────────────
    def get(self, key: str, fallback=None):
        return self._data.get(key, fallback)

    def set(self, key: str, value):
        self._data[key] = value
        self._save()

    def update(self, mapping: dict):
        self._data.update(mapping)
        self._save()

    def all(self) -> dict:
        return dict(self._data)

    # ── внутренние ────────────────────────────────────────────
    def _load(self):
        if os.path.isfile(self._path):
            try:
                with open(self._path, "r", encoding="utf-8") as f:
                    saved = json.load(f)
                self._data.update(saved)
            except Exception:
                pass  # повреждённый config — используем defaults

    def _save(self):
        try:
            with open(self._path, "w", encoding="utf-8") as f:
                json.dump(self._data, f, indent=2, ensure_ascii=False)
        except Exception:
            pass
