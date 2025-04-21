from unidec.IsoDec.IsoGen.isogenatom import *
from unidec.IsoDec.IsoGen.isogen_tools import *
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *

def cleanup_formulas(formulas):
    formulas = np.unique(formulas)
    print(len(formulas))
    # Rmove everything after + or -
    for i in range(len(formulas)):
        f = formulas[i]
        if "+" in f:
            formulas[i] = f.split("+")[0]
        if "-" in f:
            formulas[i] = f.split("-")[0]

    formulas = np.unique(formulas)
    return formulas

def parse_formulas(formulas, isolen=128, cutoff=0.001):
    dists = []
    goodforms = []
    for f in formulas:
        try:
            dist = get_dist_from_formula(f, isolen, cutoff)
            dists.append(dist)
            goodforms.append(f)
        except:
            pass

    dists = np.array(dists)
    goodforms = np.array(goodforms)
    return dists, goodforms

# Set backend to Agg
mpl.use('WxAgg')

eng = IsoGenAtomEngine(isolen=32)

# Search pubchem and get all atomic formulas
os.chdir("Z:\Group Share\JGP\PubChem")
fname = "CID-Mass.txt"



current_formulas = []

saved_files = 0

with open(fname, 'r') as f:
    for i, line in enumerate(f):
        current_formulas.append(line.split("\t")[1])

        if i % 100000 == 0 and i > 0:
            current_formulas = cleanup_formulas(current_formulas)
            dists, goodforms = parse_formulas(current_formulas)
            np.savez_compressed("pubchem_subfiles\\isodists_"+str(saved_files)+".npz", dists=dists, formulas=goodforms)
            saved_files += 1
            current_formulas = []

exit()


#Randomly select 10 files on which to assess the model
files = np.random.choice(os.listdir("pubchem_subfiles"), 10)

formulas = []
dists = []

for file in files:
    data = np.load("pubchem_subfiles\\"+file)
    formulas.extend(data["formulas"])
    dists.extend(data["dists"])


#Get the isotopic distributions from each formula
pred_dists = [eng.predict(f) for f in formulas]
true_dists = dists

#Calculate the pearson chi squared error for each distribution
chi_squareds = [calculate_pearson_chisquared(pred, true) for pred, true in zip(pred_dists, true_dists)]

print("Mean Chi Squared:", np.mean(chi_squareds))

