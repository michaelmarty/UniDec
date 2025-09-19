from unidec.IsoDec.runtime import IsoDecRuntime
import os

eng = IsoDecRuntime()

# Find all files in the directory that are csv but don't have an associated centroided file
path = "C:\\Data\\Yuri\\"
os.chdir(path)

files = [f for f in os.listdir(path) if f.endswith(".csv") and not os.path.exists(os.path.splitext(f)[0] + "_1_centroided.txt")]

print("Files to centroid:", files)

eng.batch_centroid(files)