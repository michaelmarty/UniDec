"""Private UniDec development scripts."""

import sys
from pathlib import Path

PUBLIC_UNIDEC_ROOT = Path(__file__).resolve().parent.parent / "public" / "UniDec"
if not (PUBLIC_UNIDEC_ROOT / "unidec").is_dir():
    raise RuntimeError(f"UniDec subtree not found at {PUBLIC_UNIDEC_ROOT}")

public_root = str(PUBLIC_UNIDEC_ROOT)
if public_root not in sys.path:
    sys.path.insert(0, public_root)
