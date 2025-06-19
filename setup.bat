@echo off
cd /d "%~dp0"

echo [INFO] Checking Python version . . .

for /f "tokens=*" %%i in ('py -3.12 -V 2^>nul') do set PYVER=%%i

echo %PYVER% | find "3.12.7" >nul
if errorlevel 1 (
    echo.
    echo [ERROR] Python 3.12.7 is required.
    echo Please install Python 3.12.7 before continuing.
    pause
    exit /b 1
)

set "PYTHON=py -3.12"

echo.
if not exist ".venv" (
    echo [INFO] Creating virtual environment with Python 3.12.7 . . .
    %PYTHON% -m venv .venv
) else (
    echo [INFO] Virtual environment already exists.
)

if exist ".venv\Scripts\activate.bat" (
    call .venv\Scripts\activate.bat
) else (
    echo.
    echo [ERROR] Failed to activate virtual environment.
    pause
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
make-db

if errorlevel 1 (
    echo.
    echo [ERROR] Failed to create the MFD.
    pause
    exit /b 1
)

echo.
echo [SUCCESS] Setup complete!
pause