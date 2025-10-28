import os
import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud

# import file 20231202 JDS Bgal groEL bit5 zp7 inj4s cyc1m_2023-12-07-03-46-56_TIC.txt from path Z:\Group Share\Skippy\Projects\HT\Example data for MTM\20231202 JDS Bgal groEL bit5 zp7 inj4s cyc1m_2023-12-07-03-46-56_unidecfiles
# import text file
path = 'Z:\\Group Share\\Skippy\\Projects\\HT\\Example data for MTM\\20231202 JDS Bgal groEL bit5 zp7 inj4s cyc1m_2023-12-07-03-46-56_unidecfiles\\20231202 JDS Bgal groEL bit5 zp7 inj4s cyc1m_2023-12-07-03-46-56_TIC.txt'
data = np.loadtxt(path, delimiter = ' ')
x = data[:,0]
y = data[:,1]

# Calculate SNR of data for x < 8
# calculate signal to noise ratio
# signal is the mean of the data
# noise is the standard deviation of the data
signal = np.mean(y[x < 8])
noise = np.std(y[x < 8])

#print(noise)

#Get the maximum value of the data for x > 9 and x < 11
# calculate maximum value of data
max_value = np.max(y[(x > 9) & (x < 11)])
#print(max_value)

# apply a savitzky-golay filter to the data
# smooth data
ySmth = ud.savgol(y, 25, 3)
#plt.plot(x, ySmth)

y_base = y - ySmth
# plot y_base
plt.plot(x[(x > 9) & (x < 11)], y_base[(x > 9) & (x < 11)])

# calculate standard deviation of y_base after x = 8
# calculate standard deviation of y_base
y_base_std = np.std(y_base[(x > 9) & (x < 11)])

#y_base_std = np.std(y_base)
print(y_base_std)

# plot data

#plt.plot(x, y)
plt.show()


