import json
import numpy as np
import os
import unidec

# Function to change aspects of config object
def navia2config(config_obj, json_dict):
    config_obj.minmz = json_dict['dtp']['mz_low'][0]
    config_obj.maxmz = json_dict['dtp']['mz_upp'][0]
    config_obj.subbuff = json_dict['dtp']['sub_value'][0]
    config_obj.smooth = json_dict['dtp']['gau_sigma'][0]
    config_obj.intthresh = json_dict['dtp']['intensity_threshold'][0]
    if (json_dict['dtp']['intensity_threshold'] == 'Substract Minimum'):
        config_obj.subtype = 0
    elif (json_dict['dtp']['intensity_threshold'] == 'Substract Line'):
        config_obj.subtype = 1
    elif (json_dict['dtp']['intensity_threshold'] == 'Substract Curved'):
        config_obj.subtype = 2
    config_obj.binsize = 0  # not sure here
    pass


def navia2manualfile(json_dict, mfilepath):
    man_list = []
    for i_ser in json_dict['ser_data'].keys():
        for i_peak in range(len(json_dict['ser_data'][i_ser]['x_low'])):
            #         print(str(i_ser) + ' ' + str(i_peak))
            x_max = 0.5 * (
                    json_dict['ser_data'][i_ser]['x_low'][i_peak] + json_dict['ser_data'][i_ser]['x_upp'][i_peak])
            x_diff = x_max - json_dict['ser_data'][i_ser]['x_low'][i_peak]
            x_charge = json_dict['ser_data'][i_ser]['charge'][i_peak]
            man_list.append([x_max, x_diff, x_charge])
    man_list = np.array(man_list)
    np.savetxt(mfilepath, man_list, fmt='%1.2f')
    return man_list


# Write oligomer file
def navia2ofile(json_dict, ofilepath):
    subunits = json_dict['subunits']
    masses = np.array(subunits["mass"])
    mins = np.array(subunits["min"])
    maxes = np.array(subunits["max"])
    names = np.array(subunits["name"])

    olist = []
    for i in range(0, len(masses)):
        m1 = masses[i]
        try:
            m2 = mins[i]
        except:
            m2 = 0
        try:
            m3 = maxes[i]
        except:
            m3 = 1
        try:
            name = names[i]
        except:
            name = ""
        olist.append([0, m1, m2, m3, name])
    olist = np.array(olist)
    np.savetxt(ofilepath, olist, fmt='%s')
    # I rewrote this part to avoid using pandas and have users required to install a new library
    # oligomer_frame = pd.DataFrame.from_dict(json_dict['subunits'])
    # oligomer_frame = oligomer_frame.reindex(columns=['mass', 'min', 'max', 'name'])
    # oligomer_frame.insert(0, 'offset', [0.0 for i in range(oligomer_frame.ndim)], True)
    # man_target = os.path.join(os.getcwd(), target_folder, file_name + '_ofile.dat')
    # oligomer_frame.to_csv('Test.csv', index=False, header=False, sep='\t')
    return olist


def navia_import(path):
    # Create a UniDec engine to help with things
    eng = unidec.UniDec()
    config_obj = eng.config
    try:
        # Where to find NaViA session and where to save UniDEC stuff
        with open(path, "r") as myfile:
            data = myfile.readlines()
        json_dict = json.loads(data[0])
        mz_raw = np.transpose([json_dict['raw_data']['mz'], json_dict['raw_data']['Intensity']])
        rawfile = os.path.splitext(path)[0] + ".txt"
        np.savetxt(rawfile, mz_raw)
        eng.open_file(rawfile)
        navia2config(config_obj, json_dict)
        config_obj.manuallist=navia2manualfile(json_dict, config_obj.manualfile)
        config_obj.oligomerlist=navia2ofile(json_dict, config_obj.ofile)
        if len(config_obj.manuallist) > 0:
            eng.config.manualfileflag = True
        eng.export_config()
        #eng.load_config(eng.config.confname)
        return rawfile
    except IOError:
        print("Cannot open file '%s'." % path)
        return None


if __name__ == '__main__':
    path = "C:\\Data\\Navia\\Archive\\GroEL.navia"
    path = "C:\\Data\\Navia\\Archive\\MthK-octamer.navia"
    navia_import(path)
