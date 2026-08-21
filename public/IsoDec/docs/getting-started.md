# Getting started

Install the wheel and load a two-column spectrum:

```shell
python -m pip install isodec
```

```python
from isodec import IsoDecRuntime

engine = IsoDecRuntime()
engine.process_file("spectrum.txt", assume_centroided=True)
print(len(engine.pks))
```

Built-in readers cover TXT, DAT, CSV, and NPZ files. Install
`isodec[unidec]` for raw vendor, mzML/mzXML, and I2MS formats. This optional
dependency is also required for workflows that consume UniDec reader objects.

`IsoDecRuntime.batch_process_spectrum` accepts an `N x 2` NumPy array directly
and is the preferred integration point for applications that already read
their own data.
