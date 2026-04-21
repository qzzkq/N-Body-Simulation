#!/bin/bash
# ──────────────────────────────────────────────────────────────
#  N-Body GUI — venv creating
# ──────────────────────────────────────────────────────────────

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/venv"

echo "──────────────────────────────────────────"
echo "  N-Body Simulation GUI — setup venv"
echo "──────────────────────────────────────────"

if ! command -v python3 &>/dev/null; then
    echo "❌  python3 has not been found"
    exit 1
fi

PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
echo "✓  Python $PYTHON_VERSION найден"

if [ -d "$VENV_DIR" ]; then
    echo "✓  venvalready exists in: $VENV_DIR"
else
    echo "→ Create venv.."
    python3 -m venv "$VENV_DIR"
    echo "✓  venvhas been created"
fi

echo "→ Install dependencies..."
source "$VENV_DIR/bin/activate"
pip install --upgrade pip --quiet
pip install -r "$SCRIPT_DIR/requirements.txt"

echo ""
echo "──────────────────────────────────────────"
echo "  ✅  Ready!"
echo ""
echo "  Запуск GUI:"
echo "    ./launch.sh"
echo ""
echo "  Without .sh:"
echo "    source gui/venv/bin/activate"
echo "    python gui/main.py"
echo "──────────────────────────────────────────"
