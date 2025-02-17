rem update version in pyproject.toml and enginebase

rem commit to svn and git

rem Run Tests on test_GUI and test_MUD and and ImporterTest

rem build pyinstaller upload

rem build docker and commit (see docker_commands.txt)

rem build wheel
    rem python -m build -o .\distpypi
rem test wheel
    rem create venv with: python -m venv .venv
    rem activate venv with: .venv\Scripts\activate
    rem install wheel with: pip install .\distpypi\UniDec-8.0.0-py3-none-any.whl
    rem test with: python -m unidec.IsoDec .\unidec\bin\TestSpectra\test_2.txt -precentroided
    rem test UniDec with: python -m unidec .\unidec\bin\TestSpectra\test_1.txt
rem python -m twine upload --repository testpypi .\distpypi\* --config-file .pypirc (check that old wheels are deleted)
rem redo last line with pypi instead of testpypi
rem install with: pip install unidec and test in venv

rem update docs with .\unidec_doc\make.bat html
rem paste docs into UniDecDocumentation and push to git

echo "Building"
rem C:\Python312\Scripts\pyinstaller.exe GUniDec.spec --noconfirm
C:\Users\MartyLabsOfficePC\Python\Scripts\pyinstaller.exe GUniDec.spec --noconfirm
rem call "C:\Python\UniDec3\dist\UniDec_Windows\GUI_UniDec.exe"