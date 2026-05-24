@echo off
chcp 65001 > nul
set "VENV_DIR=%~dp0venv"

echo ──────────────────────────────────────────
echo   N-Body Simulation GUI ─ setup venv (Win)
echo ──────────────────────────────────────────

where python >nul 2>nul
if %errorlevel% neq 0 (
    echo ❌ Python не найден в системе! Добавь его в PATH.
    pause
    exit /b 1
)

echo ✓ Python найден.
if exist "%VENV_DIR%" (
    echo ✓ Виртуальное окружение venv уже существует.
) else (
    echo → Создаю venv...
    python -m venv "%VENV_DIR%"
    echo ✓ venv успешно создан.
)

echo → Активация окружения и установка зависимостей...
call "%VENV_DIR%\Scripts\activate"

python -m pip install --upgrade pip --quiet
if exist "%~dp0requirements.txt" (
    pip install -r "%~dp0requirements.txt"
    echo ✓ Все зависимости из requirements.txt установлены.
) else (
    echo ❌ Файл requirements.txt не найден!
)

echo ──────────────────────────────────────────
echo   ✅ Всё готово!
echo   Для запуска используй launch.bat
echo ──────────────────────────────────────────
pause