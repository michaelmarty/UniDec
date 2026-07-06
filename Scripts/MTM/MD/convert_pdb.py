import argparse
import os
from pathlib import Path


filepath = r"C:\Users\mm96978\The University of Texas at Austin\MartyLab - General\Papers\LX_Review\Structures\defrnd.pdb"


def convert_pdb_to_mmcif(pdb_path, output_path=None, structure_id=None, preserve_atom_numbering=True):
    """Convert a PDB coordinate file to PDBx/mmCIF format."""
    try:
        from Bio.PDB import MMCIFIO, PDBParser
    except ImportError as exc:
        raise ImportError(
            "Biopython is required for PDB to PDBx/mmCIF conversion. "
            "Install it with: pip install biopython"
        ) from exc

    pdb_path = Path(os.path.expanduser(str(pdb_path))).resolve()
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    if output_path is None:
        output_path = pdb_path.with_suffix(".cif")
    else:
        output_path = Path(os.path.expanduser(str(output_path))).resolve()

    if structure_id is None:
        structure_id = pdb_path.stem

    output_path.parent.mkdir(parents=True, exist_ok=True)

    parser = PDBParser(QUIET=True, PERMISSIVE=True)
    structure = parser.get_structure(structure_id, str(pdb_path))

    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(output_path), preserve_atom_numbering=preserve_atom_numbering)

    return output_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a PDB coordinate file to PDBx/mmCIF format."
    )
    parser.add_argument(
        "pdb_path",
        nargs="?",
        default=filepath,
        help="Input PDB file. Defaults to the filepath variable in this script.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output PDBx/mmCIF path. Defaults to the input path with a .cif extension.",
    )
    parser.add_argument(
        "--structure-id",
        help="Structure ID to store in the mmCIF file. Defaults to the input file stem.",
    )
    parser.add_argument(
        "--renumber-atoms",
        action="store_true",
        help="Renumber atoms sequentially instead of preserving PDB atom serial numbers.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output_path = convert_pdb_to_mmcif(
        args.pdb_path,
        output_path=args.output,
        structure_id=args.structure_id,
        preserve_atom_numbering=not args.renumber_atoms,
    )
    print(f"Converted PDB to PDBx/mmCIF: {output_path}")
    return 0


if __name__ == "__main__":
    convert_pdb_to_mmcif(filepath)


    raise SystemExit(main())
