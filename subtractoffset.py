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
    
def ratiofuncprime(alpha, n1, n2, n3, yrato):
    numerator=exp(alpha*(n3-n2))*(n3-n2+exp(alpha*(n3-n1))*(n2-n1)
                                  +exp(alpha*(n2-n1))*(n1-n3))
    denominator=(1-exp(alpha*(n3-n1)))**2
    ratioderiv=numerator/denominator
    print ratioderiv
    return ratioderiv

def errorfunc(order,finf,cl,alpha):
    return finf-cl*np.exp(-alpha*order)

t0=470
orders=[12,16,20, 24,28, 32,36,40,44] #not 48
minmode=3
maxmode=minmode+1
i1=0
i2=i1+1
i3=i1+2

fio=open("coeffsbyl"+str(t0)+"_"+str(orders[i1])+"_"+str(orders[i2])+"_"+str(orders[i3])+".csv","a")
columnoffset=5
timecolumn=0
#modenum=6
#modenum=30
nummodes = 6
datatable =np.loadtxt("/mnt/data/sdorsher/Fortranp9.9e0.1n40/psir_l.asc", skiprows=1)
#orders=[8,16,24,32, 40,56,36,44,28,48,52,20]
#[8,16]
#orderspred=[25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56
#orderspred=[36,44, 28,20, 48] #20 defective? #52 still running, changedir
#orderspred=[32,33,36,40,44] #20 defective? #52 still running, changedir
#t0=472.721
#t0=786.7
interporder=4
interpkind='cubic'





for modenum in range(minmode,maxmode):
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
         if(orders[count]==28):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_restart/psir_l.asc"
         elif(orders[count]==24):
            loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_3_restart/psir_l.asc"
         elif(orders[count]==16 or orders[count]==36 or orders[count]==44):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_2_restart/psir_l.asc"
         elif(orders[count]==20):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_restart/psir_l.asc"
         elif(orders[count]==12 or orders[count]==32 or  orders[count]==40):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"/psir_l.asc"
         datatable=np.loadtxt(loadstring,skiprows=1)
         print loadstring
         for ii in range(len(datatable[:,timecolumn])):
             if datatable[ii,timecolumn]<t0:
                 #print(datatable[:,timecolumn])
                 tnearest=datatable[ii,timecolumn]
                 indexnearest=ii
         #print orders[count], indexnearest, interporder, len(datatable[:,timecolumn])
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
     
     yratio=(lbestarr[i1]-lbestarr[i2])/(lbestarr[i1]-lbestarr[i3])
     print lbestarr[i1], lbestarr[i2], lbestarr[i3], yratio
             
     alphaupper=1.0-yratio
     alphalower=0.0-yratio
     alpha0=0.2

     #print alphaupper, alphalower, yratio
     
     #pylab.plot(alpha,ratiofunc(alpha,n1,n2,n3))
     #pylab.show()
     alpha =optimization.newton(ratiofunc,alpha0,fprime=ratiofuncprime, tol=1e-14,args=(orders[i1],orders[i2],orders[i3],yratio),fprime2=None)
     

     #alpha =optimization.bisect(ratiofunc,alphaupper,alphalower,args=(orders[i1],orders[i2],orders[i3],yratio))
     
     print "alpha= ", alpha
     
     
     ccoeff = (lbestarr[i1]-lbestarr[i2])/(exp(-alpha*orders[i1])-exp(-alpha*orders[i2]))
     finf = lbestarr[i3]-ccoeff*exp(-alpha*orders[i3])
     
     print alpha, ccoeff, finf


     #opt_params, param_cov=optimization.curve_fit(errorfunc, orders, lbestarr, p0=np.array([alpha, ccoeff, finf]))

     #alpha=opt_params[0]
     #ccoeff=opt_params[1]
     #finf=opt_params[2]

     #print alpha, ccoeff, finf

     fio3=open("modesbyorder"+str(t0)+"_l"+str(modenum)+"_i"+str(i1)+".csv","a")
     for ii in range(len(orders)):
         csvwriter=csv.writer(fio3,delimiter=' ')
         csvwriter.writerow([orders[ii],lbestarr[ii]])
     fio3.close()
     


     
     csvwriter=csv.writer(fio,delimiter=' ')
     csvwriter.writerow([modenum,alpha,ccoeff,finf])
     #lpred= np.zeros(len(orderspred))
     lbestnew=np.zeros(len(lbestarr))
     for ii in range(len(lbestarr)):
         lbestnew[ii]=abs(lbestarr[ii]-finf)
     #for ii in range(len(orderspred)):
     #    lpred[ii]=abs(ccoeff*exp(-alpha*orderspred[ii]))
     
     #ordertot=np.concatenate(orders,orderspred)
     #ltot=np.concatenate(lbestarr,lpred)
     
     
     plt.plot(orders,lbestnew, 'o-')
     #plt.plot(orderspred, lpred, 'x')
     ax=plt.gca()
     ax.set_yscale('log')
     #ax.set_xscale('log')
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

