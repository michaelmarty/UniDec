import pybliometrics.scopus as sc
import pandas as pd
import time
from string import digits
import re

def get_digits(string):
    try:
        year = re.findall(r'\d+', string)[0]
    except:
        print(string)
        year = 2021
    return year

filepath = "C:\\Users\\marty\\OneDrive - University of Arizona\\Service\\APR\\PubData\\GradStudentDatabase2.xlsx"
df = pd.read_excel(filepath)

for index, row in df.iterrows():
    ssem = row["Start_Sem"]
    esem = row["End_Sem"]

    syear = get_digits(ssem)
    eyear = get_digits(esem)

    name = row["Name"]
    splits = name.split(",")
    last_name = splits[0]
    first_name = splits[1]
    row["First Name"] = first_name
    row["Last Name"] = last_name

    df.loc[index, "First Name"] = first_name
    df.loc[index, "Last Name"] = last_name
    df.loc[index, "Start"] = syear
    df.loc[index, "End"] = eyear


df.to_excel(filepath, index=False)