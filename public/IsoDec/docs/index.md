# IsoDec

IsoDec assigns charge states and deconvolves isotopically resolved mass
spectra. The standalone package combines a native charge-assignment core with
IsoGen isotope distributions.

```python
import numpy as np
from isodec import IsoDecRuntime

spectrum = np.loadtxt("spectrum.txt")
peaks = IsoDecRuntime().batch_process_spectrum(spectrum, centroided=True)
```

Continue with [Getting started](getting-started.md) or the [API](api.md).
