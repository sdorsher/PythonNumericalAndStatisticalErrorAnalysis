from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.optimize as optimization
import csv
import pylab

#def truncfunc(Flclalpha, xdata, ydata):
#    Fl=Flclalpha[0]
#    cl=Flclalpha[1]
#    alpha=Fclapha[2]
#    return(ydata-cl*exp(-alpha*xdata))

def ratiofunc(alpha, n1, n2, n3, yratio):
    ratio=(1-exp(-alpha*(n2-n1)))/(1-exp(-alpha*(n3-n1)))-yratio
    print ratio
    return ratio
    
def ratiofuncprime(alpha, n1, n2, n3, yratio):
    numerator=exp(alpha*(n3-n2))*(n3-n2+exp(alpha*(n3-n1))*(n2-n1)
                                  +exp(alpha*(n2-n1))*(n1-n3))
    denominator=(1-exp(alpha*(n3-n1)))**2
    ratioderiv=numerator/denominator
    print ratioderiv
    return ratioderiv
                                                          
fio=open("coeffsbyl590_24_28_32.csv","a")
columnoffset=5
timecolumn=0
#modenum=6
#modenum=30
nummodes = 31
datatable =np.loadtxt("/mnt/data/sdorsher/Fortranp9.9e0.1n40/psir_l.asc", skiprows=1)
#orders=[8,16,24,32, 40,56,36,44,28,52,20]#48 bad
orders=[20,24,28,32,33, 36,40,44]#48 bad
#[8,16]
#orderspred=[25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56]
orderspred=[36,44,52, 28,20, 48] #20 defective? #52 still running, changedir
t0=200.
#t0=786.7
interporder=4
interpkind='cubic'

i1=2
i2=i1+1
i3=i1+2


for modenum in range(4,6):
     datatablelist=list(np.zeros(0))
     tstoredlist=list(np.zeros(0))
     lbestarr=np.zeros([len(orders)])
     lstoredlist=list(np.zeros(0))
     for count in range(len(orders)):
         tnearest=0.0
         indexnearest=0
         lbest=0
         tstored = list(np.zeros(interporder))
         lstored=list(np.zeros(interporder))
         loadstring=""
         if((orders[count]==48) or (orders[count]==28) or (orders[count]==20)):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_restart/psir_l.asc"
         elif((orders[count]==36) or (orders[count]==44)):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_2_restart/psir_l.asc"
         elif(orders[count]==52):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_3_rerestart/psir_l.asc"
         else:
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"/psir_l.asc"
         datatable=np.loadtxt(loadstring,skiprows=1)
         for ii in range(len(datatable[:,timecolumn])):
             if datatable[ii,timecolumn]<t0:
                 tnearest=datatable[ii,timecolumn]
                 indexnearest=ii
         print orders[count], indexnearest, tnearest, interporder, len(datatable[:,timecolumn])
         for ii in range(interporder):
             tstored[ii]=datatable[indexnearest-(interporder-1)/2+ii,timecolumn]
             lstored[ii]=datatable[indexnearest-(interporder-1)/2+ii, columnoffset+modenum]
         func=interp1d(tstored,lstored,kind=interpkind)
         lbest=func(t0)
         lbestarr[count]=lbest
             
         tstoredlist.append(tstored)
         lstoredlist.append(lstored)

     finf=0

     lbestnew=np.zeros(len(lbestarr))
     for ii in range(len(lbestarr)):
         lbestnew[ii]=abs(lbestarr[ii]-finf)

     
     plt.plot(1./np.array(orders),lbestnew, 'o')
     ax=plt.gca()

     plt.ylabel('Re(dpsi/dr)')
     plt.xlabel('DG order')
     plt.title('Mode l='+str(modenum))
     plt.show()
     
fio.close()
    #fig=plt.plot(orders,abs(lbestarr-bestvals[0]), 'o-')
#ax=plt.gca()
#ax.set_yscale('log')
#plt.plot(orders, abs(-bestvals[1]*exp(-orders),'--')
#plt.show()

