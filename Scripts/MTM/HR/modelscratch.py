import torch
import os
import numpy as np
from unidec.IsoDec.encoding import encode_phase

example = np.array([[5.66785531e+02, 1.47770838e+06],
                    [5.67057354e+02, 1.54980838e+06],
                    [5.67507468e+02, 5.21600520e+07],
                    [5.67708173e+02, 8.35557760e+07],
                    [5.67908401e+02, 7.28264240e+07],
                    [5.68060254e+02, 1.87337225e+06],
                    [5.68108674e+02, 4.35435520e+07],
                    [5.68239256e+02, 3.88155375e+06],
                    [5.68309390e+02, 2.05468060e+07],
                    [5.68509951e+02, 7.18109250e+06],
                    [5.68707871e+02, 2.30373500e+06],
                    [5.69150563e+02, 1.57598062e+06],
                    [5.69243121e+02, 1.96390440e+07],
                    [5.69334393e+02, 6.82677120e+07],
                    [5.69425337e+02, 1.22867432e+08],
                    [5.69516492e+02, 1.45702336e+08],
                    [5.69607541e+02, 1.20801936e+08],
                    [5.69698595e+02, 1.06786072e+08],
                    [5.69789906e+02, 6.56232960e+07],
                    [5.69881208e+02, 3.41013880e+07],
                    [5.69972168e+02, 1.70930360e+07],
                    [5.70063432e+02, 9.17621100e+06],
                    [5.70699369e+02, 1.96462650e+06]])
example.astype(np.float32)
phases = encode_phase(example)



model = "C:\\Python\\UniDec3\\unidec\\IsoDec\\phase_model_8.pth"

outfile = "C:\\Data\\IsoNN\\phase_model_8.bin"

model = torch.load(model, weights_only=True)

output = []

#binfile = np.fromfile(outfile, dtype=np.float32)
#print(binfile)
#print(len(binfile))
models = [m for m in model]

w1 = model[models[0]]
#print(w1)
b1 = model[models[1]]
#print(b1)
w2 = model[models[2]]
print(w2)
b2 = model[models[3]]
print(b2)

o = torch.nn.Linear(w1.shape[1], w1.shape[0])
o.weight = torch.nn.Parameter(w1)
o.bias = torch.nn.Parameter(b1)

o2 = torch.nn.Linear(w2.shape[1], w2.shape[0])
o2.weight = torch.nn.Parameter(w2)
o2.bias = torch.nn.Parameter(b2)

inp = phases.flatten()
input = torch.tensor(inp, dtype=torch.float32).to(o.weight.device)
#print(inp[:4])

h1 = o(input)
h11 = torch.nn.functional.relu(h1)
h2 = o2(h11)

print(w2.cpu().numpy()[1])
print(w2.shape)
#exit()
#print(h1[:10])
#print(h11[:10])
#print(h2[:20])
#print(torch.argmax(h2))

#print(h11.sum())
#print(h2.sum())
exit()
for m in model:
    flat = model[m].flatten()
    print(flat)
    output.append(flat)

output = torch.cat(output).cpu().detach().numpy()

output.tofile(outfile)
