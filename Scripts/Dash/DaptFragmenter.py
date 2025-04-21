import numpy as np
import molmass

std_aa_mass = {
    'G': 57.02146372057,
    'A': 71.03711378471,
    'S': 87.03202840427001,
    'P': 97.05276384885,
    'V': 99.06841391299,
    'T': 101.04767846841,
    'C': 103.00918478471,
    'L': 113.08406397713001,
    'I': 113.08406397713001,
    'J': 113.08406397713001,
    'N': 114.04292744114001,
    'D': 115.02694302383001,
    'Q': 128.05857750527997,
    'K': 128.09496301399997,
    'E': 129.04259308796998,
    'M': 131.04048491299,
    'H': 137.05891185845002,
    'F': 147.06841391298997,
    'U': 150.95363508471,
    'R': 156.10111102359997,
    'Y': 163.06332853254997,
    'W': 186.07931294985997,
    'O': 237.14772686284996}
"""A dictionary with monoisotopic masses of the twenty standard
amino acid residues, selenocysteine and pyrrolysine.
"""
watermass = molmass.Formula("H2O").monoisotopic_mass
nh2mass = molmass.Formula("NH2").monoisotopic_mass
comass = molmass.Formula("CO").monoisotopic_mass
hmass = molmass.Formula("H").monoisotopic_mass
daptform = "C72H101N17O26"
daptmass = molmass.Formula(daptform).monoisotopic_mass
kynmass = molmass.Formula("C10H12N2O3").monoisotopic_mass - watermass
kynsidechain = molmass.Formula("C8H8NO").monoisotopic_mass
print(kynsidechain)
print(daptmass)

massdict = {
    "D": std_aa_mass["D"],
    "W": std_aa_mass["W"],
    "N": std_aa_mass["N"],
    "Thr": std_aa_mass["T"],
    "G": std_aa_mass["G"],
    "Orn": molmass.Formula("C5H12N2O2").monoisotopic_mass - watermass,
    "Ser": std_aa_mass["S"],
    "DSer": std_aa_mass["S"],
    "Ala": std_aa_mass["A"],
    "DAla": std_aa_mass["A"],
    "Dec": molmass.Formula("C10H19O").monoisotopic_mass,
    "DLT(1)": std_aa_mass["T"] - hmass,
    "GluMe": std_aa_mass["E"] + molmass.Formula("CH2").monoisotopic_mass,
    "Kyn": kynmass,
    "(1)": 0
}
print(massdict)

dapt = "Dec-W-N-D-DLT(1)-G-Orn-D-DAla-D-G-DSer-GluMe-Kyn-(1)"

isosmiles = "CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)NC3C(OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"

smiles = "CCCCCCCCCC(=O)NC(CC1=CNC2=CC=CC=C21)C(=O)NC(CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC3C(OC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"

res = dapt.split("-")

fragmods = {

    "": 0,
    # "--Dec": -massdict["Dec"],
    # "+K": molmass.Formula("K").monoisotopic_mass-hmass,
    # "+Ca2+": molmass.Formula("Ca").monoisotopic_mass - hmass - hmass,
    # "--DecH": -massdict["Dec"] - hmass,
    # "--Dec+H": -massdict["Dec"] + hmass,
    # "--DecTrp": -341.22290318179995,
    # "--DecTrpH": -341.22290318179995 - hmass,
    # "--DecTrp+H": -341.22290318179995 + hmass,
    # "--DecTrpAsn": -455.26583062294,
    # "--DecTrpAsnH": -455.26583062294 - hmass,
    # "--DecTrpAsn+H": -455.26583062294 + hmass,
    "--KynSC": -kynsidechain,
    "--KynSC+H": -kynsidechain + hmass,
    "--KynSC-H": -kynsidechain - hmass,
    # "--TrpSC": -130.066,
    # "--TrpSC+H": -130.066 + hmass,
    # "--TrpSC-H": -130.066 - hmass,
    # "-Unk": -433.26583,
    # "-Unk2": -57,
    # "-Unk3": -12,
    # "+CarboxyMethyl": 58.005479,
}


def get_mass(res):
    mass = 0
    for r in res:
        if r in massdict.keys():
            mass += massdict[r]
            pass
        elif r == "||...||":
            pass
        else:
            print("Not Found:", r)
    return mass


topmass = get_mass(res)
print(topmass)


def makename(residues):
    name = residues[0]
    if len(residues) > 1:
        for r in residues[1:]:
            name += "-" + r
    return name


def fragment(residues, lterm=0., rterm=0., filtercyclic=True):
    fragments = []
    names = []
    for i in range(1, len(residues)):
        left = residues[:i]
        right = residues[i:]
        lmarker = "||"
        rmarker = "||"

        if "DLxiThr(1)" in left and "(1)" in right and filtercyclic:
            left = left + right
            right = []
            lmarker = ""

        lmass = get_mass(left)
        rmass = get_mass(right)

        if lmass != topmass and lmass != 0:
            lname = makename(left) + lmarker
            if lname not in names:
                fragments.append([lname, lmass + lterm, "l"])
                names.append(lname)
        if rmass != topmass and rmass != 0:
            rname = rmarker + makename(right)
            if rname not in names:
                fragments.append([rname, rmass + rterm, "r"])
                names.append(rname)

    return fragments


def fragment2(residues, lterm=0., rterm=0., filtercyclic=True):
    fragments = []
    names = []
    for i in range(1, len(residues) - 1):
        for j in range(i + 1, len(residues)):
            left = residues[:i]
            center = residues[i:j]
            right = residues[j:]
            lmarker = "||"
            cmarker = "||"
            rmarker = "||"

            if "DLxiThr(1)" in left and "(1)" in right and filtercyclic:
                left = left + ["||...||"] + right
                right = []
                lmarker = ""

            if "DLxiThr(1)" in center and "(1)" in right and filtercyclic:
                center = center + right
                right = []
                cmarker = ""

            # if "DLxiThr(1)" in left and "(1)" in center and filtercyclic:
            #    left = left + center
            #    center = []

            lmass = get_mass(left)
            cmass = get_mass(center)
            rmass = get_mass(right)

            if lmass != topmass and lmass != 0:
                lname = makename(left) + lmarker
                if lname not in names:
                    fragments.append([lname, lmass + lterm, "l"])
                    names.append(lname)

            if cmass != topmass and cmass != 0:
                cname = rmarker + makename(center) + cmarker
                if cname not in names:
                    fragments.append([cname, cmass + rterm + lterm, "c"])
                    names.append(cname)

            if rmass != topmass and rmass != 0:
                rname = rmarker + makename(right)
                if rname not in names:
                    fragments.append([rname, rmass + rterm, "r"])
                    names.append(rname)

    return fragments


def get_all_fragments(lterm=watermass - hmass, rterm=hmass):
    fragments = fragment(res, lterm=lterm, rterm=rterm, filtercyclic=True)
    fragments2 = fragment2(res, lterm=lterm, rterm=rterm)

    fragmentsall = fragments + fragments2
    fragmentsall = np.array(fragmentsall)
    fragmentsall = np.unique(fragmentsall, axis=0)
    fragmentsall = fragmentsall[fragmentsall[:, 1].astype(float).argsort()]

    print(len(fragments))
    print(len(fragments2))
    print(len(fragmentsall))
    return fragmentsall


def gen_mz(fragments):
    fragmasses = fragments[:, 1].astype(float)
    fragmz = fragmasses + hmass
    fragmz2 = (fragmasses / 2) + hmass
    allfragmz = np.concatenate([fragmz, fragmz2])
    return allfragmz, np.round(fragmz, 3), np.round(fragmz2, 3)


def gen_all_mz_fragments(outfile="daptfagall.csv"):
    fpepsin = get_all_fragments()
    fraglist, mz1, mz2 = gen_mz(fpepsin)

    # Generate all fragments from by
    fby = get_all_fragments(lterm=-hmass, rterm=hmass)
    fraglistby, mz1by, mz2by = gen_mz(fby)

    # Generate all fragments from cz
    fcz = get_all_fragments(lterm=nh2mass, rterm=-nh2mass)
    fraglistcz, mz1cz, mz2cz = gen_mz(fcz)

    # Generate all fragments from ay
    fay = get_all_fragments(lterm=-comass - hmass, rterm=+comass - hmass)
    fraglistay, mz1ay, mz2ay = gen_mz(fay)

    # Stack all fragment lists
    fragmentsall = np.hstack((fpepsin, np.array([mz1, mz2, mz1by, mz2by, mz1cz, mz2cz, mz1ay, mz2ay]).transpose()))
    fragmentsall = np.vstack(
        (fragmentsall, ["DAPT", daptmass, "P", daptmass + hmass, (daptmass / 2.) + hmass, 0, 0, 0, 0, 0, 0]))
    fragmentsall = np.vstack(
        (
            ["Name", "Neutral Mass", "Cleavage Type", "Pepsin +1", "Pepsin +2", "by +1", "by +2", "cz +1", "cz +2",
             "ay +1",
             "ay +2"], fragmentsall))

    # Save to file
    np.savetxt(outfile, fragmentsall, fmt="%s", delimiter=",")

    return fragmentsall


if __name__ == "__main__":
    gen_all_mz_fragments()
    exit()

    # fragments = fragment(res, lterm=-hmass, rterm=hmass, filtercyclic=True)
    # print(fragments)

    fragmentsall = get_all_fragments()
    for f in fragmentsall:
        print(f)

    fmt = "%1.3f"
    fraglist, mz1, mz2 = gen_mz(fragmentsall)
    # Add Fraglist to fragmentsall array
    print(mz1.shape, fragmentsall.shape)
    fragmentsall = np.hstack((fragmentsall, np.array([mz1, mz2]).transpose()))
    fragmentsall = np.vstack((fragmentsall, ["DAPT", daptmass, "P", daptmass + hmass, (daptmass / 2.) + hmass]))
    print(fragmentsall)
    np.savetxt("daptfragpepsin.csv", fragmentsall, fmt="%s", delimiter=",")

    fragmentsall = get_all_fragments(lterm=-hmass, rterm=hmass)
    for f in fragmentsall:
        print(f)

    fmt = "%1.3f"
    fraglist, mz1, mz2 = gen_mz(fragmentsall)
    # Add Fraglist to fragmentsall array
    print(mz1.shape, fragmentsall.shape)
    fragmentsall = np.hstack((fragmentsall, np.array([mz1, mz2]).transpose()))
    fragmentsall = np.vstack((fragmentsall, ["DAPT", daptmass, "P", daptmass + hmass, (daptmass / 2.) + hmass]))
    print(fragmentsall)
    np.savetxt("daptfragby.csv", fragmentsall, fmt="%s", delimiter=",")
    exit()

    fraglist = np.concatenate([fraglist, [daptmass + hmass], [(daptmass / 2.) + hmass]])
    fraglist = np.round(fraglist, 3)
    # Save fraglist to text
    np.savetxt("daptfraglistpepsin.txt", fraglist, fmt=fmt)

    exit()
    fragmentsall = get_all_fragments()
    for f in fragmentsall:
        print(f)

    fmt = "%1.3f"

    fraglist = gen_mz(fragmentsall)
    fraglist = np.concatenate([fraglist, [daptmass + hmass], [(daptmass / 2.) + hmass]])
    fraglist = np.round(fraglist, 3)
    # Save fraglist to text
    np.savetxt("daptfraglistpepsin.txt", fraglist, fmt=fmt)

    fragmentsall = get_all_fragments(lterm=-hmass, rterm=hmass)
    for f in fragmentsall:
        print(f)

    fraglist = gen_mz(fragmentsall)
    fraglist = np.concatenate([fraglist, [daptmass + hmass], [(daptmass / 2.) + hmass]])
    fraglist = np.round(fraglist, 3)
    # Save fraglist to text
    np.savetxt("daptfraglistby.txt", fraglist, fmt=fmt)
