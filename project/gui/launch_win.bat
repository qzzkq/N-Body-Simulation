@echo off
chcp 65001 > nul
set "VENV_DIR=%~dp0venv"

if not exist "%VENV_DIR%" (
    echo ❌ Виртуальное окружение не создано! Сначала запусти setup_venv.bat
    pause
    exit /b 1
)

echo → Запуск N-Body GUI...
call "%VENV_DIR%\Scripts\activate"
python "%~dp0main.py" %*