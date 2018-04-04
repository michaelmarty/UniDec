import h5py
import numpy as np


def replace_dataset(group, name, data):
    if name in group.keys():
        del group[name]
    group.create_dataset(name, data=data)
    pass


def replace_dataset_strings(group, name, data):
    if name in group.keys():
        del group[name]
    ascii = [n.encode("ascii", "ignore") for n in data]
    group.create_dataset(name, data=ascii)
    pass


def get_dataset(group, name):
    if name in group.keys():
        data = group.get(name)[:]
    else:
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
    num=g.attrs["num"]
    out=[]
    for i in xrange(0,num):
        g2=hdf.require_group("ms_dataset/"+str(i))
        out.append(g2.attrs[key])
    hdf.close()
    return np.array(out)

def get_num(f):
    hdf = h5py.File(f)
    g = hdf.require_group("ms_dataset")
    num = g.attrs["num"]
    hdf.close()
    return num