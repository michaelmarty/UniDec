import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
from unidec.modules.isotopetools import *
from unidec.IsoDec.plots import *
from unidec.IsoDec.match import *
from unidec.IsoDec.encoding import *
from matplotlib import rcParams
from unidec.IsoDec.models import *

# Set backend to Agg
mpl.use('WxAgg')

rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['font.family'] = 'Arial'

