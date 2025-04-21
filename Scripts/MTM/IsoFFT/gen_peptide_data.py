import pyteomics
import numpy as np
from pyteomics import mass
from unidec.IsoDec.IsoGen.isogen_tools import rnas, rna_dict, dnas, dna_dict
import molmass as mm

def to_c_array_string(values, name="amino_acid_names", colcount=8):
    # apply formatting to each element
    values = ['"{}"'.format(v) for v in values]

    # split into rows with up to `colcount` elements per row
    rows = [values[i:i+colcount] for i in range(0, len(values), colcount)]

    # separate elements with commas, separate rows with newlines
    body = ',\n    '.join([', '.join(r) for r in rows])

    # assemble components into the complete string
    return 'char *{}[] = {{\n    {}}};'.format(name, body)

def format_array(values, formatter=str):
    return "{" + ",".join([formatter(int(v)) for v in values]) + "}"
def to_c_array(values, ctype="int", name="table", colcount=8):
    # apply formatting to each element
    values = [format_array(v) for v in values]

    # split into rows with up to `colcount` elements per row
    rows = [values[i:i+colcount] for i in range(0, len(values), colcount)]

    # separate elements with commas, separate rows with newlines
    body = ',\n    '.join([', '.join(r) for r in rows])

    # assemble components into the complete string
    return '{} {}[] = {{\n    {}}};'.format(ctype, name, body)

# amino_acids = list(mass.std_aa_mass.keys())
#
# print(amino_acids)
# print(to_c_array_string(amino_acids))
#
# vectors = []
#
# for a in amino_acids:
#     print(f'"{a}",')
#     atomic_composition = mass.std_aa_comp[a]
#
#     vector = np.zeros(5)
#     if 'H' in atomic_composition:
#         vector[0] = atomic_composition['H']
#     if 'C' in atomic_composition:
#         vector[1] = atomic_composition['C']
#     if 'N' in atomic_composition:
#         vector[2] = atomic_composition['N']
#     if 'O' in atomic_composition:
#         vector[3] = atomic_composition['O']
#     if 'S' in atomic_composition:
#         vector[4] = atomic_composition['S']
#     vectors.append(vector)
#
# print(to_c_array(vectors, ctype="int", name="amino_acid_vectors", colcount=5))
# print(len(vectors))

print(rnas)
print(to_c_array_string(rnas, name="rna_names"))

vectors = []
for r in rnas:
    print(f'"{r}",')
    vector = np.zeros(5)
    atomic_composition = rna_dict[r]
    print(atomic_composition)
    formula = mm.Formula(atomic_composition)
    comp = formula.composition()
    if 'H' in comp:
        vector[0] = comp['H'].count
    if 'C' in comp:
        vector[1] = comp['C'].count
    if 'N' in comp:
        vector[2] = comp['N'].count
    if 'O' in comp:
        vector[3] = comp['O'].count

    vectors.append(vector)

print(to_c_array(vectors, ctype="int", name="rna_vectors", colcount=5))
print(len(vectors))

print(dnas)
print(to_c_array_string(dnas, name="dna_names"))

vectors = []
for d in dnas:
    print(f'"{d}",')
    vector = np.zeros(5)
    atomic_composition = dna_dict[d]
    print(atomic_composition)
    formula = mm.Formula(atomic_composition)
    comp = formula.composition()
    if 'H' in comp:
        vector[0] = comp['H'].count
    if 'C' in comp:
        vector[1] = comp['C'].count
    if 'N' in comp:
        vector[2] = comp['N'].count
    if 'O' in comp:
        vector[3] = comp['O'].count

    vectors.append(vector)

print(to_c_array(vectors, ctype="int", name="dna_vectors", colcount=5))



