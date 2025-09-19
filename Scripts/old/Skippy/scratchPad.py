import numpy as np
x = [1.1, 2.2, 3.1, 4.3, 5.2, 6.5, 7.8, 8.2]
noiseRanges = [(0, 6), (12.5, 15), (12.5, 36)]
print(noiseRanges[0][1])
#print(x[x[:] > 3])

# print index of first value greater than 3
lbound = np.argmax(np.array(x) > 3)
hbound = np.argmax(np.array(x) > 6)

print(x[lbound:hbound])

# print all values of x greater than 3 without using a loop
#print(np.array(x[np.array(x) > 3]))


