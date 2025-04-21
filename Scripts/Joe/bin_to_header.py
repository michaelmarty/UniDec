from unidec.IsoDec.IsoGen.isogen_base import *

def bin_to_header(bin_file, header_file, varname="model_weights"):
    with open(bin_file, "rb") as f:
        data = f.read()

    with open(header_file, "a") as f:
        f.write(f"#ifndef {varname.upper()}_H\n")
        f.write(f"#define {varname.upper()}_H\n\n")
        f.write(f"// Auto-generated weights file\n")
        f.write(f"static const unsigned char {varname}[] = {{\n")

        for i, byte in enumerate(data):
            if i % 12 == 0:
                f.write("  ")
            f.write(f"0x{byte:02x}")
            if i != len(data) - 1:
                f.write(", ")
            if (i + 1) % 12 == 0:
                f.write("\n")

        f.write("\n};\n")
        f.write(f"static const unsigned int {varname}_len = {len(data)};\n\n")
        f.write(f"#endif // {varname.upper()}_H\n")

if __name__ == "__main__":
    bin_to_header("C:\\Python\\UniDec3\\unidec\\IsoDec\\IsoGen\\isogenmass_model_8.bin", "C:\\Python\\UniDec3\\unidec\IsoDec\\src_cmake\\isogen_weights.h", varname="isogenmass_model_8_bin")