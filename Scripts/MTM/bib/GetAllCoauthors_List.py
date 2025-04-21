from GetAllCoAuthors import *
import pandas as pd
import os

file = "C:\\Users\\marty\\OneDrive - University of Arizona\\Grants\\Collaborations\Scopus\\Mentors_ Scopus.xlsx"

# cd to the directory of the file
os.chdir(os.path.dirname(file))

df = pd.read_excel(file)

for i, row in df.iterrows():
    name = row['Name']
    id = row['Scopus']
    modid = "9-s2.0-" + str(id)
    print(name, id)
    try:
        coauthors = get_coauthors(modid)
        if len(coauthors) > 0:
            coauthors.to_excel(name + ".xlsx", index=False)
    except Exception as e:
        print("ERROR:", e)
        print("Could not get coauthors for", name, id)
        continue