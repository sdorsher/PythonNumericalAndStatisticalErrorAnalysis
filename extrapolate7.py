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

def func(n, alpha, ccoeff, finf):
    return finf-ccoeff*np.exp(-alpha*n)

start=7
step = 1
fio=open("coeffsbyl.csv","a")
columnoffset=5
timecolumn=0
#modenum=6
#modenum=30
nummodes = 6
datatable =np.loadtxt("/mnt/data/sdorsher/Fortranp9.9e0.1n40/psir_l.asc", skiprows=1)
#orders=[8,16,24,32, 40,56,36,44,28,48,52,20]
orders=[20, 24,28, 32,36,40,44] #not 48,33
#[8,16]
#orderspred=[25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56
#orderspred=[36,44, 28,20, 48] #20 defective? #52 still running, changedir
orderspred=[32,33,36,40,44] #20 defective? #52 still running, changedir
#t0=472.721
#t0=654
#t0=500
t0=480
#t0=786.7
interporder=4
interpkind='cubic'

i1=2
i2=i1+1
i3=i1+2


for modenum in range(start,31,step):
     print "modenum = ", modenum
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
         elif((orders[count]==36) or (orders[count]==44) or (orders[count]==48)):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_2_restart/psir_l.asc"
         elif(orders[count]==52):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_3_rerestart/psir_l.asc"
         else:
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"/psir_l.asc"
         datatable=np.loadtxt(loadstring,skiprows=1)
         print loadstring
         for ii in range(len(datatable[:,timecolumn])):
             if datatable[ii,timecolumn]<t0:
                 #print(datatable[:,timecolumn])
                 tnearest=datatable[ii,timecolumn]
                 indexnearest=ii
         print orders[count], indexnearest, interporder, len(datatable[:,timecolumn])
         for ii in range(interporder):
             tstored[ii]=datatable[indexnearest-(interporder-1)/2+ii,timecolumn]
             lstored[ii]=datatable[indexnearest-(interporder-1)/2+ii, columnoffset+modenum]
         func=interp1d(tstored,lstored,kind=interpkind)
         lbest=func(t0)
         lbestarr[count]=lbest
             
         tstoredlist.append(tstored)
         lstoredlist.append(lstored)
     #print(t0)
     #print(lbestarr)
     #print(tstoredlist)
     #print(lstoredlist)
     
     #print(orders)
     #print(lbestarr)
     
     yratio=(lbestarr[2]-lbestarr[3])/(lbestarr[2]-lbestarr[4])
     print lbestarr[2], lbestarr[3], lbestarr[4], yratio
     alpha0=0.2
     
     #pylab.plot(alpha,ratiofunc(alpha,n1,n2,n3))
     #pylab.show()
     alpha =optimization.newton(ratiofunc,alpha0,fprime=ratiofuncprime, tol=1e-14,args=(orders[i1],orders[i2],orders[i3],yratio),fprime2=None)
     #alpha =optimization.bisect(ratiofunc,alphaupper,alphalower,args=(orders[i1],orders[i2],orders[i3],yratio))
     
     print "alpha= ", alpha
     
     
     ccoeff = (lbestarr[i1]-lbestarr[i2])/(exp(-alpha*orders[i1])-exp(-alpha*orders[i2]))
     finf = lbestarr[i3]-ccoeff*exp(-alpha*orders[i3])
     
     print alpha, ccoeff, finf

     sigm=np.zeros(len(orders))
     for ii in range(len(orders)):
         if(alpha>0):
             sigm[ii]=exp(-alpha*orders[ii])
             #best guess at orders scaling given initial offset and scaling
             #from extrapolation
         else:
            sigm[ii]=1.

     #fit_params,fit_cov=optimization.curve_fit(func,orders,lbestarr)


                                               #p0=np.array([alpha, ccoeff, finf]),sigma=sigm)

     #alpha=fit_params[0]
     #ccoeff=fit_params[1]
     #finf=fit_params[2]

     
     csvwriter=csv.writer(fio,delimiter=' ')
     csvwriter.writerow([t0, modenum,alpha,ccoeff,finf])
     lpred= np.zeros(len(orderspred))
     lbestnew=np.zeros(len(lbestarr))
     for ii in range(len(lbestarr)):
         lbestnew[ii]=abs(lbestarr[ii]-finf)
     for ii in range(len(orderspred)):
         lpred[ii]=abs(ccoeff*exp(-alpha*orderspred[ii]))
        
    #ordertot=np.concatenate(orders,orderspred)
    #ltot=np.concatenate(lbestarr,lpred)
     
     
     #plt.plot(orders,lbestnew, 'o')
     #plt.plot(orderspred, lpred, 'x')
     #ax=plt.gca()
     #ax.set_yscale('log')
     #plt.ylabel('Re(dpsi/dr)')
     #plt.xlabel('DG order')
     #plt.title('Mode l='+str(modenum)+' time = ' + str(t0))
     #plt.show()
     
fio.close()


