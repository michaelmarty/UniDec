import os
import numpy as np
import matplotlib.pyplot as plt

def create_hadamard_matrix(bit_depth):
    if bit_depth == 1:
        return np.array([[1]])
    else:
        H_prev = create_hadamard_matrix(bit_depth - 1)
        top = np.hstack((H_prev, H_prev))
        bottom = np.hstack((H_prev, -H_prev))
        H = np.vstack((top, bottom))
        return H

def create_my_matrix(bit_depth, flipped=False):
    length = 2 ** bit_depth
    output = []
    ones = [1] * length
    output.append(ones)
    for i in range(bit_depth):
        row = []
        block_size = 2 ** i
        toggle = True
        count = 0
        for j in range(length):
            row.append(1 if toggle else 0)
            count += 1
            if count == block_size:
                toggle = not toggle
                count = 0
        output.append(row)

        if flipped:
            # Flip the row
            flipped = np.logical_not(row).astype(int).tolist()
            output.append(flipped)
        for j in range(block_size-1):
            # Rotate the output matrix
            row = np.roll(row, 1)
            output.append(row)

            if flipped:
                flipped = np.roll(flipped, 1)
                output.append(flipped)
    output = np.array(output)
    print(output)
    return output

m = create_my_matrix(10)

# m4 = [[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
#       [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],
#       [0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0],
#       [1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1],
#       [1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0],
#       [0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1],
#       [0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0],
#       [1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1],
#       [1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0],
#       [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1],
#       [0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0],
#       [1,0,0,1,1,0,0,1,1,0,0,1,0,1,0,1],
#       [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],
#       [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]]
#
# m2 = [[1,1,0,0], [0,0,1,1], [1,0,0,1], [0,1,1,0],[1,0,1,0], [0,1,0,1]]
#
# m3 = [[1,1,1,1,0,0,0,0],
#       [0,0,0,0,1,1,1,1],
#       [0,1,1,1,1,0,0,0],
#       [1,0,0,0,0,1,1,1],
#       [0,0,1,1,1,1,0,0],
#       [1,1,0,0,0,0,1,1],
#       [0,0,0,1,1,1,1,0],
#       [1,1,1,0,0,0,0,1],
#       [1,1,0,0,1,1,0,0],
#       [0,0,1,1,0,0,1,1],
#       [0,1,1,0,0,1,1,0],
#       [1,0,0,1,1,0,0,1],
#       [1,0,1,0,1,0,1,0],
#       [0,1,0,1,0,1,0,1]]
#
# m = np.array(m3)
# m = create_hadamard_matrix(8)

# Find inverse of m
minv = np.linalg.pinv(m)

# print("M:")
# print(m)
# print("M Inverse:")
# print(minv)

input = np.random.rand(m.shape[1]) * 10

# Remove all but the most intense input value
# input[np.argsort(input)[:-1]] = 0

output = np.dot(input, m.T)
decon = np.dot(minv, output)
print("Input:", input.shape)
print("Output:", output.shape)
print("Deconvolved:", decon.shape)

if not np.allclose(input, decon, atol=1e-5):
    print("Decon failed!")
else:
    print("Deconvolution successful!")

exit()

plt.bar(np.arange(len(input))-0.2, input, width=0.3, label='Input')
# plt.bar(np.arange(len(output)), output, width=0.2, label='Output')
plt.bar(np.arange(len(decon))+0.2, decon, width=0.3, label='Deconvolved')
plt.legend()
plt.show()

