import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from unidec.IsoDec.encoding import charge_phase_calculator, encode_phase
from unidec.IsoDec.match import create_isodist

mpl.use('WxAgg')

mass_diff_c = 1.0033
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
#example = example[:10]

# example = create_isodist(569, 16, example)
x = 1000.250 / mass_diff_c
x = np.arange(0,1,0.1)
#x = 0.068182
i = 11
rescale = x * 2 * np.pi * i
y = np.sin(rescale)
z = np.cos(rescale)
phase = (np.arctan2(z, y) / (2 * np.pi)) % 1

phase2 = (-x * i) % 1
print(x, phase, phase2)


print(charge_phase_calculator(example))

x = example[:, 0]

phases = encode_phase(example)
print(example)
print(phases)
exit()
plt.subplot(122)
plt.imshow(phases, aspect='auto')
plt.subplot(121)
plt.plot(range(1, len(phases) + 1), np.amax(phases, axis=1))
plt.show()
