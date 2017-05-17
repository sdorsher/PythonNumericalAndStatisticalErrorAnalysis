from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

columnoffset=5
timecolumn=0
nummodes=6
datatable =np.loadtxt("/mnt/data/sdorsher/Fortranp9.9e0.1n40/psir_l.asc", skiprows=1)
t0=100.
interporder=4
interpkind='cubic'

for modenum in range(1,31):
    tnearest=0.0
    indexnearest=0
    lbest=0
    tstoredlist=[]
    lstoredlist=[]

    tstored = list(np.zeros(interporder))
    lstored=list(np.zeros(interporder))
    lbest=list(np.zeros(interporder))
    for ii in range(len(datatable[:,timecolumn])):
        if datatable[ii,timecolumn]<t0:
            #print(datatable[:,timecolumn])
            tnearest=datatable[ii,timecolumn]
            indexnearest=ii
    for ii in range(interporder):
        print(modenum,tnearest, indexnearest)
        tstored[ii]=datatable[indexnearest-(interporder-1)/2+ii,timecolumn]
        print(tstored[ii])
        lstored[ii]=datatable[indexnearest-(interporder-1)/2+ii, columnoffset+modenum]
        func=interp1d(tstored,lstored,kind=interpkind)
        lbest=func(t0)
        lstoredlist.append[lbest]
        
            
ls=np.array(range(1,31))
print(len(ls),len(lstoredlist))
        
plt.plot(np.array(range(1,31)),np.array(lstoredlist),'x')
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.ylabel('Re(dpsi/dr)')
plt.xlabel('l mode')
plt.title('DG order = 40')
plt.show()
