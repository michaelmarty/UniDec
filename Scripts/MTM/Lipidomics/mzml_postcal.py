from pyopenms import *

kendrick_ref_mass = 14.01565
mass_h = 1.007276467

def kendrick_analysis(mz, charge):
    mass = np.abs(mz * charge) - charge * mass_h
    kmass = mass * 14 / kendrick_ref_mass
    kmd = kmass - np.round(kmass)
    return mass, kmass, kmd

def kendrick_analysis_adduct(mz, adduct):
    aname, charge = translate_adduct(adduct)
    return kendrick_analysis(mz, charge)


adduct_translator = np.array([["[M+H]+", "[M+H]", "+1", 1, "[M+H]1+"],
                              ["[M+NH4]+", "[M+NH4]", "+1", 1, "[M+NH4]1+"],
                              ["[M+2H]2+", "[M+2H]", "+2", 2, "[M+2H]2+"],
                              ["[M-H]-", "[M-H]", "-1", -1, "[M-H]1-"],
                              ["[M-HCOO]-", "[M-HCOO]", "-1", -1, "[M-HCOO]1-"],
                              ["[M-CH3OO]-", "[M-CH3OO]", "-1", -1, "[M-CH3OO]1-"],
                              ["[M-2H]2-", "[M-2H]", "-2", -2, "[M-2H]2-"]])

def translate_adduct(adduct):
    if adduct in adduct_translator[:, 0]:
        index = np.where(adduct_translator[:, 0] == adduct)[0][0]
        adduct_row = adduct_translator[index]
        aname = adduct_row[-1]
        charge = adduct_row[-2]
    else:
        aname = ""
        charge = 0
    return aname, charge


def adjust_mz(s, delta):
    mz, intensity = s.get_peaks()
    mz2 = mz - delta
    s.set_peaks([mz2, intensity])
    return s


def correct_precursor(s, delta):
    p = s.getPrecursors()[0]
    pmz = p.getMZ()
    pmz2 = pmz - delta
    p.setMZ(pmz2)
    s.setPrecursors([p])
    return s


def adjust_rt(s, m, b):
    rt = s.getRT()
    newrt = rt * m + b
    s.setRT(newrt)
    return s



def postcal(infile, outfile, delta):
    exp = MSExperiment()
    MzMLFile().load(infile, exp)
    filtered = MSExperiment()
    for i, s in enumerate(exp):
        if s.getMSLevel() == 2:
            s = adjust_mz(s, delta)
            s = correct_precursor(s, delta)
            filtered.addSpectrum(s)
        else:
            s = adjust_mz(s, delta)
            filtered.addSpectrum(s)

    MzMLFile().store(outfile, filtered)


def postcal_RT(infile, outfile, m, b):
    exp = MSExperiment()
    MzMLFile().load(infile, exp)
    filtered = MSExperiment()
    for i, s in enumerate(exp):
        s = adjust_rt(s, m, b)
        filtered.addSpectrum(s)
    MzMLFile().store(outfile, filtered)
