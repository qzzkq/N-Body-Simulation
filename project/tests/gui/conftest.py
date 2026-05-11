"""
Pytest fixtures shared across GUI tests.
"""

import os
import sys
import json
import tempfile
import shutil
import pytest

# Make sure the gui modules are importable
GUI_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "gui")
if GUI_DIR not in sys.path:
    sys.path.insert(0, GUI_DIR)


@pytest.fixture
def tmp_gui_dir(tmp_path):
    """Temporary directory that acts as the gui/ folder (holds config.json)."""
    return str(tmp_path)


@pytest.fixture
def app_config(tmp_gui_dir):
    """Fresh AppConfig instance backed by a temp directory."""
    from modules.config import AppConfig
    return AppConfig(tmp_gui_dir)


# ── Tkinter availability ──────────────────────────────────────────────────────

def _has_display() -> bool:
    if sys.platform.startswith("linux"):
        return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
    return True  # macOS / Windows always have a display


requires_display = pytest.mark.skipif(
    not _has_display(),
    reason="No display available (headless environment)"
)
