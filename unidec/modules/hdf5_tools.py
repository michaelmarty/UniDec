import h5py
import numpy as np


def replace_dataset(group, name, data):
    if name in list(group.keys()):
        del group[name]
    group.create_dataset(name, data=data, compression="gzip")
    pass


def replace_dataset_strings(group, name, data):
    if name in list(group.keys()):
        del group[name]
    asciilist = [n.encode("ascii", "ignore") for n in data]
    group.create_dataset(name, data=asciilist)
    pass


def get_dataset(group, name):
    try:
        if name in list(group.keys()):
            data = group.get(name)[:]
        else:
            data = np.array([])
    except:
        data = np.array([])
    return np.array(data)


def get_data(f, group, dataset):
    hdf = h5py.File(f)
    g = hdf.require_group(group)
    data = get_dataset(g, dataset)
    hdf.close()
    return data


def get_metadata(f, key):
    hdf = h5py.File(f)
    g = hdf.require_group("ms_dataset")
    num = g.attrs["num"]
    out = []
    for i in range(0, num):
        g2 = hdf.require_group("ms_dataset/" + str(i))
        out.append(g2.attrs[key])
    hdf.close()
    return np.array(out)


def get_num(f):
    hdf = h5py.File(f)
    g = hdf.require_group("ms_dataset")
    num = g.attrs["num"]
    hdf.close()
    return num


def read_attr(thing, string, config):
    if string in config.attrs.keys():
        val = config.attrs.get(string)
        if isinstance(val, np.ndarray):
            return val[0]
        else:
            return val
    else:
        return thing


def get_param(f, parameter):
    hdf = h5py.File(f)
    config_group = hdf.get("config")
    val = read_attr(None, parameter, config_group)
    return val
