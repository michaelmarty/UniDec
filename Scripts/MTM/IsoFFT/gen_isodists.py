import molmass as mm
from molmass import ELEMENTS, Element
import numpy as np

def to_c_array(values, ctype="float", name="table", formatter=str, colcount=8):
    # apply formatting to each element
    values = [formatter(v) for v in values]

    # split into rows with up to `colcount` elements per row
    rows = [values[i:i+colcount] for i in range(0, len(values), colcount)]

    # separate elements with commas, separate rows with newlines
    body = ',\n    '.join([', '.join(r) for r in rows])

    # assemble components into the complete string
    return '{} {}[] = {{\n    {}}};'.format(ctype, name, body)

def to_c_array_string(values, name="element_names", colcount=8):
    # apply formatting to each element
    values = ['"{}"'.format(v) for v in values]

    # split into rows with up to `colcount` elements per row
    rows = [values[i:i+colcount] for i in range(0, len(values), colcount)]

    # separate elements with commas, separate rows with newlines
    body = ',\n    '.join([', '.join(r) for r in rows])

    # assemble components into the complete string
    return 'char *{}[] = {{\n    {}}};'.format(name, body)

def complex_formatter(value):
    return '{{{:.10f}, {:.10f}}}'.format(value.real, value.imag)

def complex_to_c_array(values, name="fftarray", colcount=8):
    # apply formatting to each element
    values = [complex_formatter(v) for v in values]

    # split into rows with up to `colcount` elements per row
    rows = [values[i:i+colcount] for i in range(0, len(values), colcount)]

    # separate elements with commas, separate rows with newlines
    body = ',\n    '.join([', '.join(r) for r in rows])

    # assemble components into the complete string
    return 'fftw_complex {}[] = {{\n    {}}};'.format(name, body)


def create_fft(abundance, length=64):
    isodist = np.zeros(length)
    for i, a in enumerate(abundance):
        isodist[i] = a
    fft_dist = np.fft.rfft(isodist)
    return fft_dist

isotope_data = []

ffts = []

for i, e in enumerate(ELEMENTS):
    isotopes = e.isotopes
    abundances = []
    name = e.symbol
    previousmassnumber = -1
    for iso in isotopes:
        isodict = isotopes[iso]
        mass = isodict.mass
        abundance = isodict.abundance
        massnumber = isodict.massnumber
        if previousmassnumber == -1:
            previousmassnumber = massnumber
        else:
            if massnumber != previousmassnumber + 1:
                for j in range(previousmassnumber + 1, massnumber):
                    isotope_data.append([e.symbol, i+1, j, 0])
                    abundances.append(0)

        abundances.append(abundance)
        isotope_data.append([e.symbol, i+1, mass, abundance])
        previousmassnumber = massnumber
    fft_dist = create_fft(abundances)
    # if name == "S":
    #     print(isotopes)
    #     print(isotope_data[-10:])
    #     lower = name.lower()
    #     print(complex_to_c_array(fft_dist, name=lower+"fftarray"))
    #     exit()
    ffts.append(fft_dist)

# elements = [e.symbol for e in ELEMENTS]
# print(len(elements))
# print(to_c_array_string(elements, name="element_names"))
# exit()

ffts = np.array(ffts)
print(ffts.shape)
# Flatten ffts
ffts = ffts.flatten()
print(ffts.shape)
print(complex_to_c_array(ffts, name="fftarray", colcount=4))
exit()
isotope_data = np.array(isotope_data)

numbers = isotope_data[:,1].astype(int)
masses = isotope_data[:,2].astype(float)
abundances = isotope_data[:,3].astype(float)

print(to_c_array(numbers, ctype="int", name="isotope_numbers"))
print(to_c_array(masses, ctype="double", name="isotope_masses"))
print(to_c_array(abundances, ctype="double", name="isotope_abundances", formatter=lambda x: "{:.10f}".format(x)))
print(len(isotope_data))



