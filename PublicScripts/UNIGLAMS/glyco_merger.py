import pandas as pd

file2 = "C:\Data\Luis Genentech\human_glycans_formatted2.csv"
file = "C:\Data\Luis Genentech\Complete Glycan List.csv"
outfile = "C:\Data\Luis Genentech\Merged Glycan List.csv"

df = pd.read_csv(file)
df2 = pd.read_csv(file2)

df3 = pd.concat([df, df2])

df3 = df3.drop_duplicates(subset="Glycan")
df3.reset_index(inplace=True)
print(df3)
df3.to_csv(outfile)