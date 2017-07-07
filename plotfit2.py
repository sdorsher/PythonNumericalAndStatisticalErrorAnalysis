from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize



def fit_func1(ldata,param1):
    psir=param1/(2.*ldata-1.)/(2.*ldata+3.)
    return psir

def fit_func2(ldata,param1,param2):
    psir=param1/(2.*ldata-1.)/(2.*ldata+3.)+param2/(2.*ldata-1.)/(2.*ldata+3.)/(2.*ldata-3.)/(2.*ldata+5.)
    return psir

def fit_func3(ldata,param1,param2,param3):
    psir=param1/(2.*ldata-1.)/(2.*ldata+3.)+param2/(2.*ldata-3.)/(2.*ldata-1.)/(2.*ldata+3.)/(2.*ldata+5.)+param3/(2.*ldata-5.)/(2.*ldata-3.)/(2.*ldata-1.)/(2.*ldata+3.)/(2.*ldata+5.)/(2.*ldata+7.)
    return psir

def fit_func4(ldata, param1,param2,param3,param4):
    psir=param1/(2*ldata-1)/(2*ldata+3)+param2/(2*ldata-3)/(2*ldata-1)/(2*ldata+3)/(2*ldata+5)+param3/(2*ldata-5)/(2*ldata-3)/(2*ldata-1)/(2*ldata+1)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)+param4/(2*ldata-5)/(2*ldata-3)/(2*ldata-1)/(2*ldata+1)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)/(2*ldata-7)/(2*ldata+9)
    return psir

    
terms = 2
fit_func=fit_func2
if terms==1:
    fit_func=fit_func1
if terms==2:
    fit_func=fit_func2
if terms==3:
    fit_func=fit_func3
if terms==4:
    fit_func=fit_func4
#finfcolumn=4
finfcolumn=2
#lcolumn=1
lcolumn=0
nummodes=29

startmode=18
#datatable =np.loadtxt("coeffsbyl610_12_16_20.csv", skiprows=startmode)
datatable =np.loadtxt("bestinfoverinitorder.csv")

#t0 = datatable[0,0]
t0=610
llist=datatable[startmode:nummodes,lcolumn]
psir=datatable[startmode:nummodes,finfcolumn]

errscale=np.zeros(len(llist))
for ii in range(len(llist)):
    errscale[ii]=llist[ii]**-1.
errscale2=np.zeros(len(llist))
for ii in range(len(llist)):
    errscale2[ii]=llist[ii]**-2.
    
paramopt, paramcov = optimize.curve_fit(fit_func, llist,psir,sigma=errscale)
paramopt2, paramcov2=optimize.curve_fit(fit_func,llist,psir)
paramopt3, paramcov3 = optimize.curve_fit(fit_func, llist,psir,sigma=errscale2)
print psir[0], psir[len(psir)-1]

residual1=np.zeros(len(llist))
residual2=np.zeros(len(llist))
residual3=np.zeros(len(llist))
for ii in range(len(llist)):
    if terms==1:
        residual1[ii]=psir[ii]-fit_func(llist[ii],paramopt[0])
        residual2[ii]=psir[ii]-fit_func(llist[ii],paramopt2[0])
        residual3[ii]=psir[ii]-fit_func(llist[ii],paramopt3[0])
for ii in range(len(llist)):
    if terms==2:
        residual1[ii]=psir[ii]-fit_func(llist[ii],paramopt[0],paramopt[1])
        residual2[ii]=psir[ii]-fit_func(llist[ii],paramopt2[0],paramopt2[1])
        residual3[ii]=psir[ii]-fit_func(llist[ii],paramopt3[0],paramopt3[1])
for ii in range(len(llist)):
    if terms==3:
        residual1[ii]=psir[ii]-fit_func(llist[ii],paramopt[0],paramopt[1],paramopt[2])
        residual2[ii]=psir[ii]-fit_func(llist[ii],paramopt2[0],paramopt2[1],paramopt2[2])
        residual3[ii]=psir[ii]-fit_func(llist[ii],paramopt3[0],paramopt3[1],paramopt3[2])
    if terms ==4:
        residual1[ii]=psir[ii]-fit_func(llist[ii],paramopt[0],paramopt[1],paramopt[2], paramopt[3])
        residual2[ii]=psir[ii]-fit_func(llist[ii],paramopt2[0],paramopt2[1],paramopt2[2], paramopt2[3])
        residual3[ii]=psir[ii]-fit_func(llist[ii],paramopt3[0],paramopt3[1],paramopt3[2], paramopt3[3])
#plt.plot(llist,psir,'x',label='Data, extrapolated to infinite DG order')
#plt.plot(llist,fit_func2(llist,paramopt[0],paramopt[1],paramopt[2]),'-', label='3 term fit, without scaled weignts')
#plt.plot(llist,fit_func2(llist,paramopt2[0],paramopt2[1],paramopt2[2]),'-', label='3 term fit, with sigma~l^-1')
#plt.plot(llist,fit_func2(llist,paramopt3[0],paramopt3[1],paramopt3[2]),'-', label='3 term fit, with sigma~l^-2')
plt.plot(llist,0*llist,'-')
plt.plot(llist,residual1,'x-', label=str(terms)+" term fit residual, sigma~1")
plt.plot(llist,residual2,'o-', label=str(terms)+" term fit residual, sigma~l^-1")
plt.plot(llist,residual3,'+-', label=str(terms)+" term fit residual, sigma~l^-2" )

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
