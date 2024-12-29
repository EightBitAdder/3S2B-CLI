@echo off
if not "%ERRORLEVEL%"=="0" exit /b %ERRORLEVEL%

cd /d "%~dp0"

python --version | find "3.12.7" >nul
if errorlevel 1 (
    echo Python 3.12.7 is required.
    exit /b 1
)

if not exist ".venv" (
    echo Creating virtual environment with Python 3.12.7...
    python -m venv .venv
) else (
    echo Virtual environment already exists.
)

call .venv\Script\activate

set PYTHONPATH=%cd%

echo Upgrading pip...
python -m pip install --upgrade pip

echo Installing dependencies...
pip intsall -r requirements.txt

echo Creating the MFD...
cd bin
python -c "from db import makeDB; makeDB()"
cd ..

echo Setup complete!