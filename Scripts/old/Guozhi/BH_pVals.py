
import os
import matplotlib.pyplot as plt
path = "C:\\Users\\fishaska\\Desktop\\In progress\\Lipidomix\\Test"
os.chdir(path)


X, Y = [], []
for line in open('p value.txt', 'r'):
  values = [float(s) for s in line.split()]
  X.append(values[0])
  Y.append(values[1])

#x=[0,30]
#y=[0,0.05]
#plt.plot(x, y)

plt.hlines([0.05], 0,30, linestyle="--", color="black")

plt.plot(X, Y, "ro")
plt.xlabel('Rank (k)', size=20)
plt.ylabel('p-values', size=20)
plt.title("BH method", size=20)
plt.show()