@echo off
setlocal
cd /d "%~dp0"
rem Run from the Python environment containing UniDec's development dependencies.
python -m PyInstaller GUniDec.spec --noconfirm
if errorlevel 1 exit /b %errorlevel%
