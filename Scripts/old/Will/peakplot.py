import sys

sys.path.append("C:\\Python\\UniDec")

from matplotlib import rcParams

#rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=0.5
rcParams['axes.linewidth']=0.5
#rcParams['font.size']=18
# matplotlib.rc('font', family='sans-serif')
# matplotlib.rc('font', serif='Helvetica')

fint=[]
fname = 'Z:\wresager\\Bleaching\\Figs\\RAW fig\\Extracts.txt'
a = np.loadtxt(fname, skiprows=1)
time= a[:,0]
APO= a[:,1]
POPE= a[:,3]
POPC= a[:,2]
print(time)


fit1, afitdat=exp_fit(time,APO,bguess=0)
plt.plot(time,afitdat,color='green',label="APO")
plt.scatter(time,APO,color='green')

fit2, efitdat=exp_fit(time,POPE,bguess=0)
plt.plot(time,efitdat, color='purple',label="POPE")
plt.scatter(time,POPE, color='purple')

fit3, cfitdat=exp_fit(time,POPC,bguess=0)
plt.plot(time,cfitdat,color="orange", label="POPC")
plt.scatter(time,POPC,color="orange")


#plt.title('Peak Aggregation Over Time')
plt.xlabel("Time(min)")
plt.ylabel("Normalized Peak Intesity")
axes=plt.gca()
axes.set_ylim(0,1)
axes.set_xlim(2,15)
plt.legend()

plt.figure()
ax=plt.bar([1,2,3],[fit1[0],fit3[0],fit2[0]])
ax=plt.gca()
ax.set_xticklabels(["APO","POPC","POPE"])

plt.show()
