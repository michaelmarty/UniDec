import pandas as pd
import os

topdir = r"Z:\Group Share\Shivam\Reducing_FP\Final_results"
os.chdir(topdir)

transition_list_path = r"Combined_TransitionList.csv"
isolation_list_path = r"HEK_IsoList.csv"
output_path = r"Combined_TransitionList_updated_RT.csv"


trans = pd.read_csv(transition_list_path)
iso = pd.read_csv(isolation_list_path)

# Loop through transition list and update RT and RT Window based on isolation list

