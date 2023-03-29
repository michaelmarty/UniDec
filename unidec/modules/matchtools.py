import pandas as pd
from unidec.modules.biopolymertools import *
import time
from unidec.tools import nearestunsorted
import os
import math


def file_to_df(path):
    extension = os.path.splitext(path)[1]
    if extension == ".csv":
        df = pd.read_csv(path)
    elif extension == ".xlsx" or extension == ".xls":
        df = pd.read_excel(path)
    else:
        print('Extension Not Recognized', extension)
        df = None
    return df


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
def get_sitematch_list(gdf, sites=None, probs_cutoff=0,
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


def index_to_sname(index, names, sitenames):
    print(index)
    fullname = ""
    for i in range(0, len(sitenames)):
        nindex = index[i]
        name = names[i][nindex]
        fullname += "[" + sitenames[i] + ":" + name + "]"
    return fullname


def site_match(pks, oligomasslist, oligoindexes, onames, sitenames, tolerance=None):
    print("Starting Match")
    starttime = time.perf_counter()
    matches = []
    errors = []
    peaks = []
    names = []
    # print(len(oligomasslist), oligomasslist)
    # print(len(oligoindexes), oligoindexes)
    # print(onames)

    for i in range(0, pks.plen):
        p = pks.peaks[i]
        target = p.mass
        nearpt = nearestunsorted(oligomasslist, target)
        match = oligomasslist[nearpt]
        error = target - match

        if tolerance is None or np.abs(error) < tolerance:
            name = index_to_sname(oligoindexes[nearpt], onames, sitenames)
        else:
            name = ""

        p.label = name
        p.match = match
        p.matcherror = error
        matches.append(match)
        errors.append(error)
        peaks.append(target)
        names.append(name)

    matchlist = [peaks, matches, errors, names]
    endtime = time.perf_counter()
    print("Matched in: ", endtime - starttime, "s")
    return np.array(matchlist)


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
            elif type(seq) is float or type(seq) is int:
                if math.isnan(seq):
                    seq = 0
                seqs.append("")
                masses.append(float(seq))
                labels.append(k)
    return np.array(seqs), np.array(masses), np.array(labels)


def parse_fmoddf(fmoddf):
    # print(fmoddf)
    total_mod_mass = 0
    total_mod_name = ""
    if fmoddf is not None:
        for i, row in fmoddf.iterrows():
            try:
                n = float(row["Number"])
            except:
                n = 1
            modmass = row["Mass"]
            modname = row["Name"]
            total_mod_mass += n * modmass
            total_mod_name += "+" + str(n) + "[" + modname + "]"
        print("Fixed Mods:", total_mod_mass, total_mod_name)
    return total_mod_mass, total_mod_name


def parse_vmoddf(vmoddf):
    # Set up the mod masses and labaels
    if vmoddf is not None:
        modmasses = vmoddf["Mass"].to_numpy()
        modlabels = vmoddf["Name"].to_numpy()
    else:
        modmasses = [0]
        modlabels = [""]

    # print("Variable Mods:", modmasses, modlabels)
    return modmasses, modlabels


def parse_global_fixed_mod(row):
    # Get the global fixed mod
    globalfixedmod = 0
    if "Global Fixed Mod" in row.keys():
        try:
            val = float(row["Global Fixed Mod"])
            if not math.isnan(val):
                globalfixedmod = val
                print("Global Fixed Mod:", globalfixedmod)
        except Exception:
            globalfixedmod = 0
    return globalfixedmod


def check_string_for_seq(redstring, pair):
    if redstring == "All":
        return True
    if redstring == "None":
        return False

    code = "Seq" + str(pair)
    if code in redstring:
        return True
    else:
        return False


def calc_bispecific_correct(row, a_code="BsAb (Correct)", b_code="LC1 Mispair (Incorrect)",
                            c_code="LC2 Mispair (Incorrect)"):
    """
    Calculate the correct bispecific pairing from three columns of heights.
    Ref: http://dx.doi.org/10.1080/19420862.2016.1232217
    :param row: df row with the three columns
    :param a_code: the code for the first column with the correct pairing
    :param b_code: the code for the second column with the incorrect pairing
    :param c_code: the code for the third column with the incorrect pairing
    :return: x, y, the fraction of the a peak that is correctly paired and light chain scrambled, respectively.
    """
    if a_code not in row.keys() or b_code not in row.keys() or c_code not in row.keys():
        return row
    # Get Heights
    a = row[a_code.lower() + " Height"]
    b = row[b_code.lower() + " Height"]
    c = row[c_code.lower() + " Height"]

    # Normalize
    s = np.sum([a, b, c])
    if s == 0:
        return row
    a /= s
    b /= s
    c /= s

    # Calculate
    d = (a / 2.) ** 2 - (b * c)
    if d > 0:
        x = (a / 2.) + np.sqrt(d)
    else:
        x = (a / 2.)
    y = a - x

    # Set
    row["BsAb Pairing Calculated (%)"] = x * 100
    row["Light Chain Scrambled (%)"] = y * 100

    return row


known_labels = ["Correct", "Incorrect", "Ignore"]


def calc_pairs(row, include_seqs=False, remove_zeros=True, fmoddf=None, keywords=None):
    """
    For use with UniDec Processing Pipeline
    Calculate the potential pairs from a row.
    :param row: Row from a df with Sequence N in the column heading designating the sequence. Seq N + Seq M will look for a pair.
    :param include_seqs: Boolean to include the isolated sequences in the output. Default is False. If you want to include these sequences, set to True or include Seq1 as a separate incorrect column.
    :param remove_zeros: Boolean to remove the pairs with masses of 0 from the output. Default is True.
    :param fmoddf: DataFrame of the Fixed Modifications to use for the sequences. If None, will not use any.
    :return: Masses of pairs, labels of potential pairs
    """
    labels = []
    pairs = []

    # If keywords is None, use the known labels
    if keywords is None:
        keywords = known_labels

    redstring = ""
    # Try to get reduced things
    if "Disulfides Oxidized" in row.keys():
        redstring = str(row["Disulfides Oxidized"])
        # print(redstring)

    modstring = ""
    # Try to get reduced things
    if "Apply Fixed Mods" in row.keys():
        modstring = str(row["Apply Fixed Mods"])
        # print(modstring)

    total_mod_mass, total_mod_name = parse_fmoddf(fmoddf)

    # Loop through the rows and look for a known_label in the key
    for k in row.keys():
        for label in keywords:
            if label in k:
                labels.append(k)
                pairing = row[k]
                # If it has a sequence, calculate parse the sequence name into an array
                if type(pairing) is str and "Seq" in pairing:
                    pairing = pairing.replace("Seq", "")
                    pairing = np.array(pairing.split("+"))
                    # pairing = pairing.astype(int)
                    pairs.append(pairing)
                # If it is a float or int, add it to the array
                elif type(pairing) is float or type(pairing) is int:
                    if math.isnan(pairing):
                        pairing = 0
                    pairs.append([pairing])
                else:
                    # If it's a string without Seq or something else, try to make it a float
                    try:
                        p = float(pairing)
                        pairs.append([p])
                    except:
                        print("Error parsing pairing:", pairing, label)

    # Loop through the pairs and calculate their masses
    pmasses = []
    for i, pair in enumerate(pairs):
        masses = []
        for j, p in enumerate(pair):
            mass = 0
            # If it is a number, just use it
            if type(p) is float or type(p) is int:
                if math.isnan(p):
                    p = 0
                mass = float(p)
            # If it is a string, parse it
            else:
                # Find the sequence
                seqname = "Sequence " + str(p)
                reduced = check_string_for_seq(redstring, p)

                for k in row.keys():
                    if k == seqname:
                        seq = row[k]
                        # Calculate mass of the sequence
                        if type(seq) is str:
                            mass = calc_pep_mass(seq, all_cyst_ox=reduced)
                        elif type(seq) is float or type(seq) is int:
                            # If it's already a number, just use it
                            if math.isnan(seq):
                                seq = 0
                            mass = float(seq)
                        else:
                            # Something is wrong
                            print("Likely error in mass value: ", seq)
                            mass = 0

            applyfixedmod = check_string_for_seq(modstring, p)
            # Add the fixed mods
            if mass != 0 and (modstring == "" or applyfixedmod):
                mass += total_mod_mass

            # Add the mass to the list
            masses.append(mass)

        pmass = np.sum(masses)
        pmasses.append(pmass)

    # Add the individual sequences if desired
    if include_seqs:
        for k in row.keys():
            if "Sequence" in k:
                seq = row[k]
                if type(seq) is str:
                    mass = calc_pep_mass(seq)
                elif type(seq) is float or type(seq) is int:
                    if math.isnan(seq):
                        seq = 0
                    mass = float(seq)
                else:
                    print("Likely error in mass value: ", seq)
                    mass = 0

                # Add the fixed mods
                if mass != 0:
                    mass += total_mod_mass
                pmasses.append(mass)
                labels.append(k)

    # Convert to numpy arrays
    pmasses = np.array(pmasses)
    labels = np.array(labels)

    # Remove zeros if desired
    if remove_zeros:
        b1 = pmasses > 0
        pmasses, labels = pmasses[b1], labels[b1]
    return pmasses, labels


def UPP_check_peaks(row, pks, tol, vmoddf=None, fmoddf=None, favor="Closest"):
    # Get Peak Masses and Heights
    peakmasses = pks.masses
    peakheights = [p.height for p in pks.peaks]
    # seqs, seqmasses, seqlabels = calc_seqmasses(row)

    # Get the favored match
    if "Favored Match" in row.keys():
        favor = row["Favored Match"]
        print("Favoring:", favor)

    # Get the global fixed mod
    globalfixedmod = parse_global_fixed_mod(row)

    # Calculate the potential pairs
    pmasses, plabels = calc_pairs(row, fmoddf=fmoddf)
    # Add the global fixed mod
    pmasses += globalfixedmod
    # print(peakmasses)
    print(np.transpose([pmasses, plabels]))
    if len(pmasses) == 0:
        print("No masses found. Exiting Pair Mode.")
        return row
    plabelints = np.zeros(len(plabels))

    # Get Variable Mods
    modmasses, modlabels = parse_vmoddf(vmoddf)

    # Create a grid of all possible combinations of peak masses and mod masses
    pmassgrid = np.array([modmasses + p for p in pmasses])
    # TODO: Allow a flag to turn off mods to ignore.

    # Set up the arrays to store the results
    correctint = 0
    incorrectint = 0
    unmatchedint = 0
    ignoreheight = 0
    totalint = np.sum(peakheights)
    matches = []
    matchstring = ""
    # Loop through the peaks and find the closest match
    for i, p in enumerate(pks.peaks):
        # Find the minimum error
        peakmass = p.mass
        errors = peakmass - pmassgrid
        abserr = np.abs(errors)
        minerror = np.amin(abserr)

        # Get the number of matches
        b1 = abserr < tol
        # Set some values on peak object
        p.numberalts = np.sum(b1)

        # If there are more than one matches, find all possible matches within the tolerance
        true_indexes = np.argwhere(b1)
        true_errors = abserr[b1]
        labels = plabels[true_indexes[:, 0]]

        p.altmatches = np.array([plabels[x[0]] + "+" + modlabels[x[1]] for x in true_indexes])
        p.altmatcherrors = errors[b1]

        # If something is within the tolerance, add it to the list
        if minerror < tol:
            # If the goal is to favor the closest, find the closest match and use that
            if favor == "Closest":
                closeindex2D = np.unravel_index(np.argmin(np.abs(pmassgrid - peakmass)), pmassgrid.shape)
            else:
                # Otherwise, if there is a single match, use that
                if p.numberalts == 1:
                    closeindex2D = np.unravel_index(np.argmin(np.abs(pmassgrid - peakmass)), pmassgrid.shape)
                else:
                    # If there are any matches with labels that match the favored, use those
                    if favor == "Incorrect":
                        favored_indexes = np.array(
                            ["Correct" not in label and "Ignore" not in label for label in labels])
                    else:
                        favored_indexes = np.array([favor in label for label in labels])

                    if np.sum(favored_indexes) > 0:
                        # Select the closet match from the favored matches
                        true_indexes = true_indexes[favored_indexes]
                        true_errors = true_errors[favored_indexes]
                        closeindex2D = true_indexes[np.argmin(true_errors)]
                    else:
                        # Otherwise, default back to the closest match
                        closeindex2D = np.unravel_index(np.argmin(np.abs(pmassgrid - peakmass)), pmassgrid.shape)

            closeindex = closeindex2D[0]
            # print(pmassgrid.shape, closeindex2D)
            p.match = pmassgrid[tuple(closeindex2D)]
            p.matcherror = peakmass - p.match

            label = plabels[closeindex] + "+" + modlabels[closeindex2D[1]]

            # Remove label from altmatches
            if len(p.altmatches) > 0:
                p.altmatches = np.delete(p.altmatches, np.argwhere(p.altmatches == label))

            matches.append(label)
            matchstring += " " + label

            rowkey = plabels[closeindex]

            if "Correct" in rowkey:
                # print("Correct", peakmass, label)
                correctint += peakheights[i]
                p.color = [0, 1, 0]  # Green

            elif "Ignore" in rowkey:
                # print("Ignore", peakmass, label)
                p.color = [0, 0, 1]  # Blue
                totalint -= peakheights[i]  # Remove this from the total intensity]
                ignoreheight += peakheights[i]

            elif "Incorrect" in rowkey:
                # print("Incorrect", peakmass, label)
                incorrectint += peakheights[i]
                p.color = [1, 0, 0]  # Red

            else:
                print("Error: Something went wrong with row name parsing", rowkey)
                # print("Incorrect", peakmass, label)
                # incorrectint += peakheights[i]
                # p.color = [1, 0, 0]  # Red

            plabelints[closeindex] += peakheights[i]

        else:
            label = "Unknown"
            matches.append("")
            # print("Unmatched", peakmass)
            unmatchedint += peakheights[i]
            p.color = [1, 1, 0]  # Yellow
            p.match = 0
            p.matcherror = 0
        p.label = label

    for i, label in enumerate(plabels):
        newrowkey = label.lower() + " Height"
        row[newrowkey] = plabelints[i]

        newrowkey2 = label.lower() + " %"
        if totalint != 0:
            value = plabelints[i] / totalint * 100
        else:
            value = 0
        row[newrowkey2] = value
        # print(row[newrowkey2], newrowkey2)

    row["Total Height"] = totalint
    row["Total correct Height"] = correctint
    row["Total incorrect Height"] = incorrectint
    row["Total unmatched Height"] = unmatchedint
    row["Total ignored Height"] = ignoreheight

    # Calculate the percentages
    percents = np.array([correctint, incorrectint, unmatchedint]) / totalint * 100
    # print(percents)
    if correctint + incorrectint > 0:
        percents2 = np.array([correctint, incorrectint]) / (correctint + incorrectint) * 100
    else:
        percents2 = np.array([0, 0])
    # print(percents2)

    row["correct %"] = percents[0]
    row["incorrect %"] = percents[1]
    row["unmatched %"] = percents[2]
    row["correct % Matched Only"] = percents2[0]
    row["incorrect % Matched Only"] = percents2[1]
    row["Matches"] = matchstring

    row = calc_bispecific_correct(row)

    return row


def dar_calc(pks, protein_mass, drug_mass, min_drugs, max_drugs, tolerance):
    """
    Calculate the DAR for a given ADC. This is the weighted average of the number of drugs conjugated to the antibody.
    :param pks: The peaklist object
    :param protein_mass: The mass of the protein
    :param drug_mass: The mass of the drug
    :param min_drugs: The minimum number of drugs conjugated
    :param max_drugs: The maximum number of drugs conjugated
    :param tolerance: The tolerance for matching the ADC mass to a peak
    :return: The DAR
    """
    # Get the peaks
    peak_masses = np.array([p.mass for p in pks.peaks])
    peak_heights = np.array([p.height for p in pks.peaks])
    # Set all peaks to yellow
    for p in pks.peaks:
        p.color = [1, 1, 0]
        p.label = ""

    # Loop over all the possible numbers of drugs conjugated
    adc_heights = []
    adc_nums = []
    for i in range(min_drugs, max_drugs):
        # Calc the mass of the protein + drug
        adc_mass = protein_mass + drug_mass * i
        # Find the closest peak to the mass
        cindex = np.argmin(np.abs(peak_masses - adc_mass))
        closest_peak = peak_masses[cindex]
        # If the closest peak is within the tolerance, use it
        if np.abs(closest_peak - adc_mass) < tolerance:
            # Get the peak height
            peak_height = peak_heights[cindex]
            # Add the peak height to the list
            adc_heights.append(peak_height)
            adc_nums.append(i)

            # Format peak
            pks.peaks[cindex].color = [0, 1, 0]
            pks.peaks[cindex].label = "ADC+" + str(i) + " Drugs"
        else:
            # Otherwise, add a zero
            adc_heights.append(0)
            adc_nums.append(i)

    # Make a numpy array of the heights and nums
    adc_heights = np.array(adc_heights)
    adc_nums = np.array(adc_nums)
    # Get the total ADC intensity
    total_adc = np.sum(adc_heights)
    if total_adc == 0:
        print("Error in DAR calculation, no intensity matched to ADC")
        return -1
    else:
        # Get the weighted average number of drugs conjugated
        dar_val = np.sum(adc_nums * adc_heights) / total_adc
        # Get the standard deviation of the number of drugs conjugated

        print("DAR:", dar_val)
        return dar_val
