from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize

def fit_func(ldata, param0,param1,param2,param3,param4):
    psir=param0+param1/(2*ldata-1)/(2*ldata+3)+param2/(2*ldata-3)/(2*ldata-1)/(2*ldata+3)/(2*ldata+5)+param3/(2*ldata-5)/(2*ldata-3)/(2*ldata-1)/(2*ldata+1)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)+param4/(2*ldata-5)/(2*ldata-3)/(2*ldata-1)/(2*ldata+1)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)/(2*ldata-7)/(2*ldata+9)
    return psir
    
def test_func(ldata,C,A):
    return C*ldata**-A

finfcolumn=4
lcolumn=1
nummodes=31
t0=570
startmode=2
datatable =np.loadtxt("coeffsbyl570.csv", skiprows=startmode)
            
llist=datatable[:,lcolumn]
psir=datatable[:,finfcolumn]

#errscale=zeros(len(llist))
#for ii in range(len(llist)):
    

paramopt, paramcov = optimize.curve_fit(fit_func, llist,psir)
linearparam, linearcov = optimize.curve_fit(test_func, llist, psir)



plt.plot(llist,psir,'x',label='Extrapolated data')
plt.plot(llist,fit_func(llist,paramopt[0],paramopt[1], paramopt[2],paramopt[3],paramopt[4]),'-', label='Fit using Hesthaven et al')
plt.plot(llist,test_func(llist,linearparam[0], linearparam[1]), '--', label='Power law fit, exp=-'+str(linearparam[1]))
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.legend(loc='upper right')
plt.ylabel('Re(dpsi/dr)')
plt.xlabel('l mode')
plt.title('Radial self force at infinite DG order, t=' + str(t0)+ ', starting from l=' + str(startmode))
plt.show()
