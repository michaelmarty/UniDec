from annika_functions import calc_percent


path = "Z:\\Group Share\\Annika\\LXMS, Fall 2024\\Discs August\\2024_08_06 Post Exchange Discs\\2024_08_06_ALS_02_(+)1_PostEx_EqD7_APCIpos.raw"

cmass = 369.36
d7cmass = 376.3
popc_mass = 760.60
d7pc_mass = 753.5

startime = 6.5
endtime= 8.5

starttime2 = 6.5
endtime2 = 8.5

starttime3 = 6.5
endtime3 = 8.5

data = calc_percent(path, cmass, d7cmass, popc_mass, d7pc_mass, startime, endtime, starttime2, endtime2, starttime3, endtime3, 10)