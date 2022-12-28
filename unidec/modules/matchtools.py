import pandas as pd
from unidec.modules.biopolymertools import *


def sort_sitematches_by(indexes, masses, probs, type="mass"):
    """
    Simple function for sorting by either mass or probability
    :param indexes: The indexes array
    :param masses: The masses array
    :param probs: The probability array
    :param type: Either "mass" (default) or "prob" to sort by mass or probability respectively
    :return: indexes, masses, probabilities (all sorted)
    """
    # sort everything by either prob or mass
    if type == "prob":
        sortindex = np.argsort(probs)[::-1]  # Sort by decending probability
    else:
        sortindex = np.argsort(masses)
    smasses = masses[sortindex]
    sindex = indexes[sortindex]
    sprobs = probs[sortindex]
    return sindex, smasses, sprobs


# Function to generate the brute force lists needed for matching
def get_sitematch_list(gdf, sites=None, probs_cutoff=1,
                       masscolumn="Monoisotopic mass", namecolumn="Glycan", percent=True, sort="mass"):
    """
    Function to generate a list of potential mass matches from a DataFrame with different potential sites.
    The dataframe should include three+ columns. One column is the name of the species, specified by namecolumn.
    A second column with the mass of the species, specified by masscolumn.
    A third or more columns specifying the probabilities of finding a particular species at that site.
    The names of the columns specifying each site are specified by the sites parameter.
    :param gdf: A pandas DataFrame with the species as rows and the mass, name, and probabilities for each site as columns
    :param sites: A list of the column names for each of the sites.
            Probabilites are provided in each column for each species row at that site.
    :param probs_cutoff: Cutoff to remove probabilities below a certain value. Specified in percent.
    :param masscolumn: The name of the mass
    :param namecolumn: The name of the column in the DF specifying the species. Default is "Glycan".
    :param percent: True/False, specifies whether to assume probabilities are provided in percent and thus to divide by 100.
    :param sort: None to leave lists unsorted, "mass" to sort by mass, "prob" to sort by probabilities
    :return: indexes, masses, probs, names. Lists, all the same length. Each row a possible combination of site mods.
    Indexes is the index from the possible mods of each site for combination.
    Masses is the mass of that possible combination.
    Probs is the probability of that possible combination based on the product of each site.
    Names is the list of names in each site.
    """
    # Make data frames for each site with only the glycans with non-zero probability
    if sites is None:
        sites = ["S1", "S1", "S2", "S2", "S3", "S3"]
    dfs = [gdf[gdf[s] > probs_cutoff] for s in sites]
    lens = [len(df) for df in dfs]
    lsites = len(sites)

    # Extract the probs and masses to make it faster
    if percent:
        pvals = [dfs[i][s].to_numpy() / 100 for i, s in enumerate(sites)]
    else:
        pvals = [dfs[i][s].to_numpy() for i, s in enumerate(sites)]
    mvals = [dfs[i][masscolumn].to_numpy() for i, s in enumerate(sites)]
    names = [dfs[i][namecolumn].to_numpy() for i, s in enumerate(sites)]

    print("Starting to brute force combine...", lens)
    # Get all possible combinations of each
    indexes = np.array(list(np.ndindex(tuple(lens))))
    print("Total combinations: ", len(indexes))
    # Loop through the sites to calculate the probabilities and masses for all indexes
    probs = np.empty((lsites, len(indexes)))
    masses = np.empty((lsites, len(indexes)))
    for i, s in enumerate(sites):
        inds = indexes[:, i]
        probs[i] = pvals[i][inds]
        masses[i] = mvals[i][inds]
    probs = np.product(probs, axis=0)
    masses = np.sum(masses, axis=0)

    if sort is not None:
        print("Sorting")
        indexes, masses, probs = sort_sitematches_by(indexes, masses, probs, type=sort)
    # print("Cutting off low probs")
    try:
        b1 = probs / np.amax(probs) > 0  # probs_cutoff / 100.
        indexes, masses, probs = indexes[b1], masses[b1], probs[b1]
    except:
        pass
    print("Finished creating list. Starting to match.")
    return indexes, masses, probs, names


def sitematch_to_target(targetmass, indexes, masses, probs, tolerance=5):
    """
    Function to match a particular target mass to combinations provided
    :param targetmass: Mass to match to
    :param indexes: Indexes from get_sitematch_list
    :param masses: Masses from get_sitematch_liist
    :param probs: Probabilities from get_sitematch_list
    :param tolerance: Tolerance (in mass units) of the match.
    :return: Selected indexes, masses, and probs arrays where all masses in the returned arrays are +/- tolerance away
    from the targetmass
    """
    b1 = np.abs(masses - targetmass) < tolerance
    return indexes[b1], masses[b1], probs[b1]


def sitematch_to_excel(indexes, masses, probs, names, peakmasses, protmass, sites, outfile):
    """
    Write out an Excel file with the potential matches for each peak
    :param indexes: Indexes list from get_sitematch_list
    :param masses: Masses list from get_sitematch_list
    :param probs: Probs list from get_sitematch_list
    :param names: Names list from get_sitematch_list
    :param peakmasses: List of peak masses to match against. Each will have its own sheet of potential masses in the excel file
    :param protmass: Mass of the constant protein to add to each possible mod to
    :param sites: List of site names to output. Doesn't need to match the get_sitematch_list above.
    :param outfile: Excel file path to write to.
    :return: None
    """
    with pd.ExcelWriter(outfile) as writer:
        # Loop through all peaks
        for mass in peakmasses:
            # Calculate the difference in mass between the measured peak and the protein mass
            peakdelta = mass - protmass
            # Match the glycan masses to the delta mass
            matchindexes, matchmasses, matchprobs = sitematch_to_target(peakdelta, indexes, masses, probs)

            # Create DataFrame and load with key data
            matchdf = pd.DataFrame()
            matchdf["Measured Delta Mass"] = np.zeros_like(matchmasses) + peakdelta
            matchdf["Match Mass"] = matchmasses
            matchdf["Error"] = matchmasses - peakdelta
            matchdf["Probability"] = matchprobs
            # Normalize the probabilities
            try:
                matchdf["Sum Norm Prob"] = matchprobs / np.sum(matchprobs)
                matchdf["Max Norm Prob"] = matchprobs / np.amax(matchprobs)
            except:
                matchdf["Sum Norm Prob"] = matchprobs * 0
                matchdf["Max Norm Prob"] = matchprobs * 0
            # Find the name of the glycan for each site from the index
            for i, s in enumerate(sites):
                inds = matchindexes[:, i].astype(int)
                if len(inds) > 0:
                    matchdf[s] = names[i][inds]
                else:
                    matchdf[s] = ""

            # print(matchdf)
            # Write to excel into a sheet with the peak mass as the name
            matchdf.to_excel(writer, sheet_name=str(mass))


def calc_seqmasses(row):
    """
    Function to calculate a sequence mass from a row
    :param row: Row from a df with Sequence N as a column header for each sequence you want to calculate
    :return: Sequences, Masses, Labels
    """
    seqs = []
    masses = []
    labels = []
    for k in row.keys():
        if "Sequence" in k:
            seq = row[k]
            if type(seq) is str:
                seqs.append(seq)
                mass = calc_pep_mass(seq)
                masses.append(mass)
                labels.append(k)
    return np.array(seqs), np.array(masses), np.array(labels)


def calc_pairs(row):
    """
    For use with UniDec Pharma Pipeline
    Calculate the potential pairs from a row.
    :param row: Row from a df with Sequence N in the column heading designating the sequence. Seq N + Seq M will look for a pair.
    :return: Masses of pairs, labels of potential pairs
    """
    labels = []
    pairs = []
    for i, item in enumerate(row):
        if type(item) is str:
            if "+" in item or "Seq" in item:
                label = row.keys()[i]
                labels.append(label)
                pairing = item
                pairing = pairing.replace("Seq", "")
                pairing = np.array(pairing.split("+"))
                pairing = pairing.astype(int)
                pairs.append(pairing)
    pmasses = []
    for i, pair in enumerate(pairs):
        masses = []
        for j, p in enumerate(pair):
            seqname = "Sequence " + str(p)
            for k in row.keys():
                if k == seqname:
                    seq = row[k]
                    if type(seq) is str:
                        mass = calc_pep_mass(seq)
                        masses.append(mass)
        pmass = np.sum(masses)
        pmasses.append(pmass)

    for k in row.keys():
        if "Sequence" in k:
            seq = row[k]
            if type(seq) is str:
                mass = calc_pep_mass(seq)
                pmasses.append(mass)
                labels.append(k)
    return pmasses, labels


def UPP_check_peaks(row, pks, tol, moddf):
    peakmasses = pks.masses
    peakheights = [p.height for p in pks.peaks]
    # seqs, seqmasses, seqlabels = calc_seqmasses(row)
    pmasses, plabels = calc_pairs(row)
    print(peakmasses)
    print(np.transpose([pmasses, plabels]))
    modmasses = moddf["Mass"].to_numpy()
    modlabels = moddf["Name"].to_numpy()
    pmassgrid = np.array([modmasses + p for p in pmasses])

    correctint = 0
    incorrectint = 0
    unmatchedint = 0
    totalint = np.sum(peakheights)
    matches = []
    matchstring = ""

    for i, peakmass in enumerate(peakmasses):
        closeindex2D = np.unravel_index(np.argmin(np.abs(pmassgrid - peakmass)), pmassgrid.shape)
        closeindex = closeindex2D[0]
        error = np.argmin(np.abs(pmassgrid - peakmass))
        label = plabels[closeindex] + "+" + modlabels[closeindex2D[1]]
        if error < tol:
            matches.append(label)
            matchstring += " " + label
            if "Correct" in plabels[closeindex]:
                print("Correct", peakmass, label)
                correctint += peakheights[i]
            else:
                print("Incorrect", peakmass, label)
                incorrectint += peakheights[i]
        else:
            matches.append("")
            print("Unmatched", peakmass)
            unmatchedint += peakheights[i]

    percents = np.array([correctint, incorrectint, unmatchedint]) / totalint * 100
    print(percents)
    if correctint + incorrectint > 0:
        percents2 = np.array([correctint, incorrectint]) / (correctint + incorrectint) * 100
    else:
        percents2 = np.array([0, 0])
    print(percents2)
    return percents, percents2, matches, matchstring
