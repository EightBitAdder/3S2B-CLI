@echo off
cd /d "%~dp0"

if not exist ".venv\Scripts\activate.bat" (
	echo [ERROR] Virtual environment not found.
	echo Run setup.bat before continuing.
	exit /b 1
)

call .venv\Scripts\activate.bat
3s2b %*