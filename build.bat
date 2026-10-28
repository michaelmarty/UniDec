@echo off
setlocal
cd /d "%~dp0"
rem Run from the Python environment containing UniDec's development dependencies.
rem --clean prevents stale analysis data from omitting newly added package resources.
python -m PyInstaller GUniDec.spec --noconfirm --clean
if errorlevel 1 exit /b %errorlevel%
