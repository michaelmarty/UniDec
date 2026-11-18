"""SI 8/9 using the speciated lipid name after ``|`` as the lipid identity."""

import os
import sys

import matplotlib

if __name__ == "__main__":
    matplotlib.use(os.environ.get("MPLBACKEND", "Agg" if "--no-show" in sys.argv else "WxAgg"))

import SI_8_9 as _original
from SI_8_9 import *


_original_recompute_adduct_networks = _original.recompute_adduct_networks
_original_prepare_strength = _original.prepare_strength


def _use_speciated_names(source):
    frame = source.copy()
    names = frame["Metabolite name"].fillna("").astype(str).str.split("|", n=1)
    speciated = names.str[1].fillna("").astype(str).str.strip()
    frame["Metabolite name"] = speciated.where(speciated.ne(""), names.str[0].astype(str).str.strip())
    return frame


def recompute_adduct_networks(source, graph_store=None, accepted_only=True):
    return _original_recompute_adduct_networks(
        _use_speciated_names(source), graph_store=graph_store, accepted_only=accepted_only)


def prepare_strength(source, all_ids=True, excluded_only=False, allow_empty=False):
    return _original_prepare_strength(
        _use_speciated_names(source), all_ids=all_ids, excluded_only=excluded_only, allow_empty=allow_empty)


def main(argv=None):
    original_doc = _original.__doc__
    _original.recompute_adduct_networks = recompute_adduct_networks
    _original.prepare_strength = prepare_strength
    _original.__doc__ = __doc__
    try:
        return _original.main(argv)
    finally:
        _original.recompute_adduct_networks = _original_recompute_adduct_networks
        _original.prepare_strength = _original_prepare_strength
        _original.__doc__ = original_doc


if __name__ == "__main__":
    main()
