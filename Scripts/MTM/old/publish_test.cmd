@echo off
pushd "%~dp0"
python -m build --outdir .\distpypi .\public\UniDec
echo python -m twine upload --repository testpypi .\distpypi\*
popd
