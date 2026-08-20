@echo off
setlocal
pushd "%~dp0"

rem update version in public\UniDec\pyproject.toml and the public UniDec source

rem commit to dev git and release git

rem Run Tests on test_GUI and test_MUD and and ImporterTest

rem build pyinstaller upload

rem Docker image is built, tested, and published by .github/workflows/publish-to-docker-hub.yml

rem build wheel
    rem python -m build --outdir .\distpypi .\public\UniDec
rem test wheel
    rem create venv with: python -m venv .venv
    rem activate venv with: .venv\Scripts\activate
    rem install wheel with: pip install .\distpypi\UniDec-8.0.0-py3-none-any.whl
    rem test UniDec with: set PYTHONPATH=%CD%\public\UniDec
    rem test UniDec with: python -m unidec .\public\UniDec\unidec\bin\TestSpectra\test_1.txt
    rem deactivate venv with: deactivate
rem python -m twine upload --repository testpypi .\distpypi\* --config-file .pypirc (check that old wheels are deleted)
rem redo last line with pypi instead of testpypi
rem install with: pip install unidec and test in venv

rem Documentation is built and deployed by .github/workflows/deploy-documentation.yml

echo "Building"
"C:\Python313\Scripts\pyinstaller.exe" GUniDec.spec --clean --noconfirm
if errorlevel 1 goto end
call ".\dist\UniDec_Windows\GUI_UniDec.exe"

:end
set "EXIT_CODE=%ERRORLEVEL%"
popd
exit /b %EXIT_CODE%
