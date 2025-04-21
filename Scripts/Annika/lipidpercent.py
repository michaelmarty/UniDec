#calculate the percent lipid in a nanodisc by the peak area

POPCratio = 0.876

cholratio = 2


peakareaPOPC = 133178.52

peakareaD7PC = 6214.98
peakareachol = 1079.78

peakareaD7chol = 2391.89


uMPOPC = (peakareaPOPC * (POPCratio/peakareaD7PC))

uMchol = (peakareachol * (cholratio/peakareaD7chol))

print(uMPOPC)
print(uMchol)

cholpercent = (uMchol/(uMchol+uMPOPC))*100

print(cholpercent)