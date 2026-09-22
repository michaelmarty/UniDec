"""Run from public/IsoDec: python tests/show_fragment_matches.py"""

from types import SimpleNamespace

from isodec import match_fragments
from isodec.fragment_view import show_fragment_matches
from test_fragment_matching import CA_MASSES, CA_SEQUENCE


peaks = SimpleNamespace(peaks=[SimpleNamespace(monoiso=mass, monoisos=[mass])
                               for mass in CA_MASSES])
match_fragments(peaks, CA_SEQUENCE, fragmentation_type="ETD", ppm_tolerance=20)
show_fragment_matches(CA_SEQUENCE, peaks)
