@ECHO OFF
pushd %~dp0
python generate_api_docs.py
if errorlevel 1 goto end
call make.bat html

:end
set EXIT_CODE=%ERRORLEVEL%
popd
exit /b %EXIT_CODE%
