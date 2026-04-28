#!/bin/bash
# ──────────────────────────────────────────────────────────────
#  N-Body GUI — запуск
# ──────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/venv"

if [ ! -d "$VENV_DIR" ]; then
    echo "❌  venv doesn't exists"
    exit 1
fi

source "$VENV_DIR/bin/activate"
python "$SCRIPT_DIR/main.py" "$@"
