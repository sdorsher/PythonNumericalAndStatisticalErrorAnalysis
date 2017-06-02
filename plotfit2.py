from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize

def fit_func(ldata, param1,param2,param3,param4):
    psir=param1/(2*ldata-1)/(2*ldata+3)+param2/(2*ldata-3)/(2*ldata-1)/(2*ldata+3)/(2*ldata+5)+param3/(2*ldata-5)/(2*ldata-3)/(2*ldata-1)/(2*ldata+1)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)+param4/(2*ldata-5)/(2*ldata-3)/(2*ldata-1)/(2*ldata+1)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)/(2*ldata-7)/(2*ldata+9)
    return psir

def fit_func2(ldata,param1,param2,param3):
    psir=param1/(2*ldata-1)/(2*ldata+3)+param2/(2*ldata-3)/(2*ldata-1)/(2*ldata+3)/(2*ldata+5)+param3/(2*ldata-5)/(2*ldata-3)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)
    return psir
    

finfcolumn=4
lcolumn=1
nummodes=31

startmode=14
datatable =np.loadtxt("coeffsbyl570.csv", skiprows=startmode)

t0 = datatable[0,0]
llist=datatable[:,lcolumn]
psir=datatable[:,finfcolumn]

errscale=np.zeros(len(llist))
for ii in range(len(llist)):
    errscale[ii]=llist[ii]**-1.
errscale2=np.zeros(len(llist))
for ii in range(len(llist)):
    errscale2[ii]=llist[ii]**-2.
    
paramopt, paramcov = optimize.curve_fit(fit_func2, llist,psir,sigma=errscale)
paramopt2, paramcov2=optimize.curve_fit(fit_func2,llist,psir)
paramopt3, paramcov3 = optimize.curve_fit(fit_func2, llist,psir,sigma=errscale2)


residual1=np.zeros(len(llist))
residual2=np.zeros(len(llist))
residual3=np.zeros(len(llist))
for ii in range(len(llist)):
    residual1[ii]=psir[ii]-fit_func2(llist[ii],paramopt[0],paramopt[1],paramopt[2])
    residual2[ii]=psir[ii]-fit_func2(llist[ii],paramopt2[0],paramopt2[1],paramopt2[2])
    residual3[ii]=psir[ii]-fit_func2(llist[ii],paramopt3[0],paramopt3[1],paramopt3[2])

#plt.plot(llist,psir,'x',label='Data, extrapolated to infinite DG order')
#plt.plot(llist,fit_func2(llist,paramopt[0],paramopt[1],paramopt[2]),'-', label='3 term fit, without scaled weignts')
#plt.plot(llist,fit_func2(llist,paramopt2[0],paramopt2[1],paramopt2[2]),'-', label='3 term fit, with sigma~l^-1')
#plt.plot(llist,fit_func2(llist,paramopt3[0],paramopt3[1],paramopt3[2]),'-', label='3 term fit, with sigma~l^-2')
plt.plot(llist,0*llist,'-')
plt.plot(llist,residual1,'x-', label="3 term fit residual, sigma~1")
plt.plot(llist,residual2,'o-', label="3 term fit residual, sigma~l^-1")
plt.plot(llist,residual3,'+-', label="3 term fit residual, sigma~l^-2" )

ax=plt.gca()
#plt.axis([0.1,100,1e-7,1e-4])
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.legend(loc='upper right')
plt.ylabel('Re(dpsi/dr)')
plt.xlabel('l mode')
plt.title('Radial self force fit residuals, t=' + str(t0))
          #', starting from l=' + str(startmode))
plt.show()
