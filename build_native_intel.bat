@echo off
setlocal

if not defined ONEAPI_ROOT set "ONEAPI_ROOT=C:\Program Files (x86)\Intel\oneAPI"
if not defined UNIDEC_BUILD_DIR set "UNIDEC_BUILD_DIR=%~dp0build\native-intel"

call "%ONEAPI_ROOT%\setvars-vcvarsall.bat"
if errorlevel 1 exit /b %errorlevel%
call "%ONEAPI_ROOT%\setvars.bat" intel64
if errorlevel 1 exit /b %errorlevel%

set "CC=icx-cl"
set "CXX=icx-cl"
cmake -S "%~dp0unidec\src" -B "%UNIDEC_BUILD_DIR%" -G Ninja -DCMAKE_BUILD_TYPE=Release -DUNIDEC_REQUIRE_INTELLLVM=ON
if errorlevel 1 exit /b %errorlevel%
cmake --build "%UNIDEC_BUILD_DIR%" --parallel
if errorlevel 1 exit /b %errorlevel%

if not exist "%~dp0unidec\bin\UniDec.exe" exit /b 1
if not exist "%~dp0unidec\bin\unideclib.dll" exit /b 1
echo Built %~dp0unidec\bin\UniDec.exe and unideclib.dll
