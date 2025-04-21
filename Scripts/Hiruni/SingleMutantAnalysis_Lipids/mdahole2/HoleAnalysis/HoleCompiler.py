import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from HoleFunctions import *

#dir = "Z:\\Group Share/Group/VMD/AqpZMutants/WT_AA_TOCL2_2/namd\\hole"
#dir = "Z:\\Group Share/Group/VMD/AqpZMutants/WT_AA_TOCL2_2/namd\\hole"
#dir = "Z:\\Group Share/Group/VMD/AqpZMutants/WT_AA_POCL2_1/namd\\hole"
#dir = "Z:\\Group Share/Group/VMD/AqpZMutants/WT_AA_POPG_1/namd\\hole"
#dir = "Z:\\Group Share/Group/VMD/AqpZMutants/WT_AA_PE_PG_CL_1/namd\\hole"
#dir = "Z:\\Group Share/Group/VMD/AqpZMutants/WT_AA_PEPGCL_310K_1/namd\\hole"
dir = "Z:\\Group Share/Group/VMD/AqpZMutants/WT_AA_PEPG_1/namd\\hole"


os.chdir(dir)

avgs, all, datas = process_multiholes()

df = pd.DataFrame(avgs, columns=['rxn', 'avg', 'std'])
df.to_excel("profile.xlsx", index=False)