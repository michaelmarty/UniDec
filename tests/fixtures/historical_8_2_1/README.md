# Historical UniDec 8.2.1 fixtures

Generated with upstream tag `v.8.2.1`, commit `fe848be8d0039aa838e72ce87313f8eb17d6c80e`. The historical writers identify their serialized format as `8.2.0`; the original writer outputs are retained byte-for-byte.

## Files

- `unidec-8.2.1-state.zip` - SHA-256: 1d2c75e148be3076513b7ade26666b78f51ea8830332d1f892c5b8a10615d099
- `metaunidec-8.2.1-project.hdf5` - SHA-256: 7a384ca3b597e6eecddc503c1c35b7b19202007a7c670fb12f4c0b506e0dbeb4

## Generation environment

- OS: Windows-11-10.0.26300-SP0
- Python: 3.13.5
- NumPy: 2.2.6
- h5py: 3.14.0
- UniDecImporter: 1.0.4
- Historical code loaded from the extracted pinned source: verified by `unidec.__file__` during generation.
- The generation script ran in a clean temporary directory using a venv created with `python -m venv --system-site-packages --without-pip`. [Generation recipe](../../../../../Plans/Completed/UNIDEC_HISTORICAL_FIXTURE_PLAN.md#generation-script)

### `python -m pip freeze` snapshot

The snapshot omits two workspace-specific editable and local installations.

```text
altgraph==0.17.5
annotated-types==0.8.0
anyio==4.16.0
attrdict==2.0.1
attrdict3==2.0.2
attrs==26.1.0
build==1.6.1
certifi==2026.10.23
cffi==2.1.1
charset-normalizer==3.5.2
click==8.5.0
clr_loader==0.2.7.post0
colorama==0.4.6
contourpy==1.3.3
coverage==7.17.0
cryptography==50.0.2
cycler==0.12.1
docutils==0.23
et_xmlfile==2.0.0
fonttools==4.59.0
ghp-import==2.1.0
h11==0.16.0
h5py==3.14.0
httpcore==1.0.9
httpcore2==2.16.0
httpx==0.28.1
httpx-sse==0.4.3
httpx2==2.16.0
id==1.6.1
idna==3.20
iniconfig==2.3.1
isodec==2.0.5
IsoSpecPy==2.6.0
jaraco.classes==3.4.0
jaraco.context==6.1.2
jaraco.functools==4.6.0
Jinja2==3.1.6
joblib==1.5.2
jsonschema==4.26.0
jsonschema-specifications==2025.9.1
keyring==25.7.0
kiwisolver==1.4.8
llvmlite==0.44.0
lxml==6.0.0
lxml_html_clean==0.4.2
Markdown==3.11.2
markdown-it-py==4.2.0
MarkupSafe==3.0.3
matplotlib==3.10.3
mcp==1.30.0
mcp-types==2.3.0
mdurl==0.1.2
mergedeep==1.3.4
mkdocs==1.6.1
mkdocs-get-deps==0.2.2
more-itertools==11.2.1
mpld3==0.5.12
narwhals==2.25.0
natsort==8.4.0
networkx==3.5
nh3==0.3.7
numba==0.61.2
numpy==2.2.6
openpyxl==3.2.0b1
opentelemetry-api==1.45.2
packaging==26.3
pandas==2.3.1
pandas-stubs==2.3.2.250926
pathspec==1.1.1
pefile==2024.8.26
pillow==12.0.0
pillow_heif==1.1.1
platformdirs==4.13.5
plotly==7.0.0
pluggy==1.7.0
psims==1.4.0
pubsub==0.1.2
pycparser==2.22
pydantic==2.14.2
pydantic-settings==2.16.2
pydantic_core==2.50.2
Pygments==2.21.0
pyimzML==1.5.5
pyinstaller==6.22.3
pyinstaller-hooks-contrib==2026.9
pyisogen==1.1.3
PyJWT==2.15.1
pymzml==2.5.11
pyparsing==3.2.3
pyproject_hooks==1.3.3
Pypubsub==4.0.3
pyteomics==5.0.2
pytest==9.1.1
python-dateutil==2.9.0.post0
python-dotenv==1.2.4
python-multipart==0.0.33
pythonnet==3.0.5
pytz==2025.2
pywin32==311
pywin32-ctypes==0.2.3
PyYAML==6.0.3
pyyaml_env_tag==1.1
readme_renderer==46.0
referencing==0.37.0
regex==2025.7.34
requests==2.34.2
requests-toolbelt==1.0.0
rfc3986==2.0.0
rich==15.0.0
rpds-py==2026.9.1
scikit-learn==1.7.2
scipy==1.16.1
scipy-stubs==1.16.3.0
setuptools==84.0.0
six==1.17.0
SQLAlchemy==2.1.6
sse-starlette==3.6.0
starlette==1.10.0
tenacity==9.1.4
threadpoolctl==3.6.0
truststore==0.10.5
twine==7.0.0
typing-inspection==0.4.4
typing_extensions==4.16.0
tzdata==2025.2
urllib3==2.8.0
uvicorn==0.54.0
watchdog==6.0.0
wheezy.template==3.2.3
wxPython==4.2.3
```

## Reopen results

Both files were copied by themselves into separate empty temporary directories and opened with the historical reader from the pinned source. Both checks passed. The ZIP restored `startz=4`, `endz=12`, and `procflag=1`; raw and processed arrays matched the acceptance values at `rtol=0, atol=1e-6`. The HDF5 project restored `startz=5`, `endz=15`, and `datanorm=0`; names were `low.txt`, `high.txt` in order and raw and processed arrays matched at `rtol=0, atol=1e-6`.

The current-reader regressions also passed against temporary copies of only the fixture under test. The ZIP has exactly the three expected root members and relative `input` and `output` config paths. The HDF5 `/config` has no `minratio`; the current reader supplies `minratio=0.0`.
