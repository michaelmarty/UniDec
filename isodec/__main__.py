"""Command-line interface for IsoDec."""

import argparse
import sys
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    __package__ = "isodec"

from .runtime import IsoDecRuntime


def build_parser():
    parser = argparse.ArgumentParser(
        prog="isodec",
        description="Assign charge states and deconvolve isotopically resolved mass spectra.",
    )
    parser.add_argument("spectrum", type=Path, help="Input spectrum (TXT, DAT, CSV, NPZ, or a UniDec-supported format)")
    parser.add_argument("-o", "--output", type=Path, help="Output TSV filename")
    parser.add_argument("--centroided", action="store_true", help="Treat the input as centroid data")
    parser.add_argument("--type", choices=("PEPTIDE", "RNA"), default="PEPTIDE", help="Analyte isotope model")
    parser.add_argument("--quiet", action="store_true", help="Suppress progress output")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if not args.spectrum.is_file():
        raise SystemExit(f"Spectrum does not exist: {args.spectrum}")

    engine = IsoDecRuntime(verbose=not args.quiet)
    engine.analyte_type = args.type
    engine.process_file(
        str(args.spectrum),
        assume_centroided=args.centroided,
        verbose=not args.quiet,
    )
    output = args.output or args.spectrum.with_name(f"{args.spectrum.stem}_isodec.tsv")
    engine.pks.export_tsv(output, report_multiple_monoisos=engine.config.report_multiple_monoisos)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
