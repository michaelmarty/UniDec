import molmass as mm
import numpy as np
import os
import time
from isogenatom import elements
from isogenatom_trainingdata import parse_formulas


def create_formulas():
    formulas = []
    for a in elements:
        formulas.append(a)

    for a in elements:
        for b in elements:
            formulas.append(a + b)

    for a in elements:
        for b in elements:
            for c in elements:
                formulas.append(a + b + c)

    # for a in elements:
    #     for b in elements:
    #         for c in elements:
    #             for d in elements:
    #                 formulas.append(a + b + c + d)

    for a in elements:
        formulas.append(a + a + a + a)
        formulas.append(a + a + a + a + a)
        formulas.append(a + a + a + a + a + a)
        formulas.append(a + a + a + a + a + a + a)
        formulas.append(a + a + a + a + a + a + a + a)
        formulas.append(a + a + a + a + a + a + a + a + a)
        formulas.append(a + a + a + a + a + a + a + a + a + a)
        formulas.append(a + a + a + a + a + a + a + a + a + a + a)
        formulas.append(a + a + a + a + a + a + a + a + a + a + a + a)
        formulas.append(a + a + a + a + a + a + a + a + a + a + a + a + a)
        formulas.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        formulas.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a + a)
        formulas.append(a + a + a + a + a + a + a + a + a + a + a + a + a + a + a + a)

    return formulas

if __name__ == "__main__":
    # Search pubchem and get all atomic formulas
    os.chdir("Z:\Group Share\JGP\PubChem")
    starttime= time.perf_counter()

    formulas = create_formulas()
    print("Created Formulas:", len(formulas))
    dists, goodforms = parse_formulas(formulas)

    print("Parsed Formulas:", len(goodforms))
    np.savez_compressed("isodists_synthetic_" + str(len(dists)) + ".npz", dists=dists, formulas=goodforms)

    print("Time:", time.perf_counter() - starttime)

