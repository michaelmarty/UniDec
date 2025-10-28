import sys
sys.path.append("C:\\Python\\UniDec")

directory = "C:\\Users\\michael.marty\\Dropbox\\Mix ND manuscript\\Data\\fft result PG"
os.chdir(directory)
comps = [0,25,50,75,100]
reps = [1,2,3]


for c in comps:
    ws=[]
    for r in reps:
        fname = str(c)+"-"+str(r)+".txt"
        data = np.loadtxt(fname)
        data[:,1] = data[:,1]/np.amax(data[:,1])
        wstd=weighted_std(data[:,0],data[:,1])
        wfit=isolated_peak_fit(data[:,0],data[:,1],0)[0][0,0]
        ws.append(wfit)
        #plt.plot(c,wfit, marker="o")
        #plt.plot(data[:,0],data[:,1])
    plt.errorbar(c,np.mean(ws), marker = "s", yerr = np.std(ws))
plt.show()