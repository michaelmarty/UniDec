"""Print the processed and deconvolved ion counts for one CD-MS NPZ file."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


# Allow this script to be run directly from PowerShell or PyCharm.
SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = next(
    parent
    for parent in SCRIPT_PATH.parents
    if (parent / "unidec" / "modules" / "CDEng.py").is_file()
)
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from unidec.modules.CDEng import UniDecCD


# Change this path to analyze a different CD-MS file.
DATA_FILE = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT\Protein Mix 21 Injections Set 1\merged"
    r"\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_21.npz"
)


def main() -> int:
    """Process and deconvolve DATA_FILE, then print both ion counts."""
    if not DATA_FILE.is_file():
        print(f"ERROR: File does not exist or is not accessible: {DATA_FILE}")
        return 1

    engine = UniDecCD()
    engine.open_file(str(DATA_FILE))

    if engine.process_data(transform=False) != 1:
        print("ERROR: CDEng returned an empty histogram during processing.")
        return 1
    before_deconvolution = len(engine.farray)

    # The data are already processed above, so do not process them a second time.
    if engine.run_deconvolution(process_data=False) == 0:
        print("ERROR: CDEng returned an empty histogram during deconvolution.")
        return 1

    # Deconvolution produces a weighted histogram rather than a per-ion array.
    # Its sum is the post-deconvolution ion count used elsewhere in UniDec.
    after_deconvolution = float(np.sum(engine.harray))

    print("\nCD-MS Ion Counts")
    print(f"File: {DATA_FILE}")
    print(f"Before deconvolution (processed): {before_deconvolution:,}")
    print(f"After deconvolution:             {after_deconvolution:,.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
