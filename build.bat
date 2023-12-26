rem update version in pyproject.toml and enginebase
rem commit to svn and git
rem Run Tests on test_GUI and test_MUD
rem build pyinstaller
rem build docker and commit
rem build wheel, test, and commit
rem python -m build -o .\distpypi
rem python -m twine upload --repository testpypi .\distpypi\* --config-file .pypirc (check that old wheels are deleted)
rem update docs with .\unidec_doc\make.bat html
rem paste docs into UniDecDocumentation and push to git

echo "Building"
C:\Python311\Scripts\pyinstaller.exe GUniDec.spec --noconfirm
rem call "C:\Python\UniDec3\dist\UniDec_Windows\GUI_UniDec.exe"