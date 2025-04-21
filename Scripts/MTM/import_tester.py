import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from Scripts.old.data_reader import DataImporter as DI
from unidec.modules.thermo_reader.ThermoImporter import ThermoDataImporter as TI
matplotlib.use("TkAgg")

path = "C:\\Data\\20210126_ADHclean.RAW"
#path = "C:\\Data\\223321.RAW"

ti = TI(path)


scanrange=[8, 9]


array = ti.msrun.GetCentroidArray()
print(array)
#data = ti.get_data(time_range=scanrange)
exit()
importer = DI(path)
filt=np.array(importer.msrun.filters())
print(len(filt))


data = DI(path).get_data(time_range=scanrange)

data2 = TI(path)

#print(np.mean(np.diff(data[:,1]-data2[:,1])))
plt.plot(data2[:,0], data2[:,1], color="r")
plt.plot(data[:,0], data[:,1], color = "k")

plt.show()

'''
from unidec.modules.thermo_reader.ThermoImporter import ThermoDataImporter as TI
path = "C:\\Data\\file.raw"
ti = TI(path)
data = ti.get_data()


exit()
importer = DI(path)
filt=np.array(importer.msrun.filters())
print(len(filt))


data = DI(path).get_data(time_range=scanrange)

'''