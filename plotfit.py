from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize

def fit_func(ldata, param0,param1,param2,param3,param4):
    psir=param0+param1/(2*ldata-1)/(2*ldata+3)+param2/(2*ldata-3)/(2*ldata-1)/(2*ldata+3)/(2*ldata+5)+param3/(2*ldata-5)/(2*ldata-3)/(2*ldata-1)/(2*ldata+1)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)+param4/(2*ldata-5)/(2*ldata-3)/(2*ldata-1)/(2*ldata+1)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)/(2*ldata-7)/(2*ldata+9)
    return psir
    


columnoffset=5
timecolumn=0
nummodes=6
datatable =np.loadtxt("/mnt/data/sdorsher/Fortranp9.9e0.1n32/psir_l.asc", skiprows=1)
t0=200.
interporder=4
interpkind='cubic'
tstoredlist=[]
lstoredlist=[]

for modenum in range(1,31):
    tnearest=0.0
    indexnearest=0
    lbest=0

    tstored = list(np.zeros(interporder))
    lstored=list(np.zeros(interporder))
    lbest=list(np.zeros(interporder))
    for ii in range(len(datatable[:,timecolumn])):
        if datatable[ii,timecolumn]<t0:
            #print(datatable[:,timecolumn])
            tnearest=datatable[ii,timecolumn]
            indexnearest=ii
    for ii in range(interporder):
        print(modenum,t0, tnearest, indexnearest)
        tstored[ii]=datatable[indexnearest-(interporder-1)/2+ii,timecolumn]
        print(tstored[ii])
        lstored[ii]=datatable[indexnearest-(interporder-1)/2+ii, columnoffset+modenum]
    func=interp1d(tstored,lstored,kind=interpkind)
    lbest=func(t0)
    lstoredlist.append(lbest)
    
            
llist=np.array(range(1,31))
print(len(llist),len(lstoredlist))
psir=np.array(lstoredlist)

paramopt, paramcov = optimize.curve_fit(fit_func, llist,psir)

plt.plot(llist,psir,'x')
plt.plot(llist,fit_func(llist,paramopt[0],paramopt[1], paramopt[2],paramopt[3],paramopt[4]),'-')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.ylabel('Re(dpsi/dr)')
plt.xlabel('l mode')
plt.title('DG order = 32, t=200')
plt.show()
