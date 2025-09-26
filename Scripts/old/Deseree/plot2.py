import sys
sys.path.append("C:\\Python\\UniDec")

dir="Z:\Group Share\Deseree\MSPs\\201705112_DR_MSP1D1T0m_0_unidecfiles"
file="201705112_DR_MSP1D1T0m_0_mass.txt"

path=os.path.join(dir,file)

dir="Z:\Group Share\Deseree\MSPs\\201705112_DR_MSP1D1G1_0_unidecfiles"
file="201705112_DR_MSP1D1G1_0_mass.txt"

path2=os.path.join(dir,file)

dir="Z:\Group Share\Deseree\MSPs\\20170529_DJR_T1_0_unidecfiles"
file="20170529_DJR_T1_0_mass.txt"

path3=os.path.join(dir,file)

dir="Z:\Group Share\Deseree\MSPs\\201705112_DR_MSP1D1R1_0_unidecfiles"
file="201705112_DR_MSP1D1R1_0_mass.txt"

path4=os.path.join(dir,file)

dir="Z:\Group Share\Deseree\MSPs\\201705112_DR_MSP1D1T2_0_unidecfiles"
file="201705112_DR_MSP1D1T2_0_mass.txt"

path5=os.path.join(dir,file)

paths=[path, path2, path3,path4, path5]
colors=["black","Blue","purple","orange","green",]
for i,path in enumerate(paths):
    print(path)

    data=np.loadtxt(path)
    data=datachop(data,22000,22350)
    plt.plot(data[:,0] /1000.0, data[:,1]/np.amax(data[:,1])-i*0.00,color=colors[i])
plt.show()

plt.savefig('Z:\\Group Share\\Deseree\\BeltsShift.pdf')
