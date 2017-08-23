from math import *
import numpy as np
import scipy.optimize as optimization
import matplotlib.pyplot as plt
import csv
import pylab

def zerooffsetexpdecay(x,alpha,coeff):
    return coeff*exp(-alpha*x)

datatablefiniteorder=np.loadtxt('genrawsummodes2.csv')
datatablefinf=np.loadtxt('selfForceOverTime_fit3.csv')
summedfinfcolumn = 4 #column for finf with unextrapolated sum
tcolumn=0
t0=390
psirfinfsummed=0
tindex=0
#datatablefiniteorder=np.loadtxt('genrawdatamodes390.csv')
#datatablefinf=np.loadtxt('bestinfoverinitorder390.csv')
#summedfinfcolumn = 7 #column for finf with unextrapolated sum
#tcolumn=0
#t0=390
#psirfinfsummed=0
#tindex=0
#assume times match in both files
for index in range(len(datatablefinf[:,summedfinfcolumn])):
    #print index, abs(datatablefinf[index,tcolumn]-t0)
    if abs(datatablefinf[index,tcolumn]-t0)<0.01:
        #print index, tcolumn, t0
        tindex=index
        psirfinfsummed=datatablefinf[tindex,summedfinfcolumn]
#print t0, tindex, psirfinfsummed


absoluteerr=np.zeros([8])
#for order in range(len(datatablefiniteorder[tindex,1:])):
for order in range(len(datatablefiniteorder[tindex,1:])):
    absoluteerr[order]=abs(psirfinfsummed-datatablefiniteorder[index,order+1])
    print absoluteerr[order]
abserr=absoluteerr.reshape(1,8)
#logabserr=np.log10(abserr)
orders=np.array([range(16,48,4)])
#print absoluteerr.shape
#print orders.shape


orders2=orders.flatten()
abserr2=abserr.flatten()
#popt,pcov=optimization.curve_fit(zerooffsetexpdecay, orders2, abserr2)

plt.xlabel("DG order")
plt.ylabel("Absolute error")
plt.title("|Psir(DG order)-Psir(Finf)|, no l-mode fit")
#plt.semilogy('log')
plt.plot(orders2,abserr2,'o')
#plt.plot(orders,popt[0]*np.exp(-popt[1]*orders),'x')
#t = np.arange(0.01, 20.0, 0.01)
#t = np.arange(16, 44, 4)
#plt.semilogy(t, np.exp(-t/500.))
#print orders.flatten()
#plt.semilogy(orders2, np.exp(-orders2/500.))
plt.grid(True)
plt.show()

