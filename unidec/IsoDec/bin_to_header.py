"""Convert a PyTorch model to a binary file and C header."""

from pathlib import Path

import torch


def pth_to_bin(pth_file, bin_file):
    """Flatten a PyTorch state dictionary into a raw float32 binary file."""
    pth_file = Path(pth_file)
    bin_file = Path(bin_file)

    state_dict = torch.load(pth_file, map_location="cpu", weights_only=True)
    weights = torch.cat([tensor.flatten() for tensor in state_dict.values()])
    weights.detach().cpu().numpy().tofile(bin_file)


def bin_to_header(bin_file, header_file, varname="model_weights"):
    """Write the bytes from ``bin_file`` to a C array in ``header_file``."""
    bin_file = Path(bin_file)
    header_file = Path(header_file)
    data = bin_file.read_bytes()
    guard = f"{varname.upper()}_H"

    with header_file.open("w", encoding="ascii", newline="\n") as outfile:
        outfile.write(f"#ifndef {guard}\n")
        outfile.write(f"#define {guard}\n\n")
        outfile.write(f"// Generated from {bin_file.name}\n")
        outfile.write(f"static const unsigned char {varname}[] = {{\n")

        for start in range(0, len(data), 12):
            chunk = data[start:start + 12]
            values = ", ".join(f"0x{byte:02x}" for byte in chunk)
            if start + 12 < len(data):
                values += ","
            outfile.write(f"  {values}\n")

        outfile.write("};\n")
        outfile.write(f"static const unsigned int {varname}_len = {len(data)};\n\n")
        outfile.write(f"#endif // {guard}\n")


if __name__ == "__main__":
    isodec_dir = Path(__file__).resolve().parent

    input_pth = isodec_dir / "phase_model_8.pth"
    output_bin = isodec_dir / "phase_model_8.bin"
    output_header = isodec_dir / "src_cmake" / "phase_model_8.h"
    variable_name = "phase_model_8_bin"

    pth_to_bin(input_pth, output_bin)
    print("Wrote binary model:", output_bin)

    bin_to_header(output_bin, output_header, variable_name)
    print("Wrote C header:", output_header)
