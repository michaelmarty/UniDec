python -m build -o .\distpypi
echo python -m twine upload --repository testpypi .\distpypi\*

