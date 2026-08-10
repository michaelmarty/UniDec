import os

topdir = "Z:\Group Share\JGP"

# count the number of .raw files in this directory recursively

count = 0
for root, dirs, files in os.walk(topdir):
    for file in files:
        if file.endswith(".raw"):
            count += 1

print(count)