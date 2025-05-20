@echo off
cd /d "%~dp0"

echo [INFO] Checking Python version . . .

python --version | find "3.12.7" >nul
if errorlevel 1 (
    echo.
    echo [ERROR] Python 3.12.7 is required.
    echo Please install the correct Python version before continuing.
    exit /b 1
)

echo.
if not exist ".venv" (
    echo [INFO] Creating virtual environment with Python 3.12.7 . . .
    python -m venv .venv
) else (
    echo [INFO] Virtual environment already exists.
)

if exist ".venv\Scripts\activate.bat" (
    call .venv\Scripts\activate.bat
) else (
    echo.
    echo [ERROR] Failed to activate virtual environment.
    exit /b 1
)

echo.
echo [INFO] Upgrading pip . . .
python -m pip install --upgrade pip

echo.
echo [INFO] Installing project in editable mode . . .
pip install -e .

echo.
echo [INFO] Creating the MFD . . .
python -c "from db.make import makeDB; makeDB()"

echo.
echo [SUCCESS] Setup complete!
pause