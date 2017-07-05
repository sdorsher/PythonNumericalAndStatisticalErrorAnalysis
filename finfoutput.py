# outputs a table of finf in order of increasing mode along the rows and
# and in order of increasing order along the columns. order number is first
# column

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
    return (exp(-alpha*n2)-exp(-n1*alpha))/(exp(-n3*alpha)-exp(-alpha*n1))-yratio
    
def ratiofuncprime(alpha, n1, n2, n3, yratio):
    print alpha, n1, n2, n3
    return exp(-alpha*(-n1+n2+n3))*(exp(alpha*n3)*(-n1+n2)+exp(alpha*n2)*(n1-n3)+exp(alpha*n1)*(-n2+n3))/(-1+exp(alpha*(n1-n3)))**2
    
   
def func(n, alpha, ccoeff, finf):
    return finf-ccoeff*np.exp(-alpha*n)
t0=590
orders=[12,16, 20, 24,28, 32,36,40,44] #not 48,33
start=0
step = 1
fio2=open("finfoutput590.csv","a")

columnoffset=5
timecolumn=0
#modenum=6
#modenum=30
nummodes = 6
datatable =np.loadtxt("/mnt/data/sdorsher/Fortranp9.9e0.1n40/psir_l.asc", skiprows=1)


#orders=[ 24,28, 32,36,40,44] #not 48,33
#[8,16]
#orderspred=[36,40,44] #20 defective? #52 still running, changedir

#t0=472.721
#t0=654
#t0=500
#t0=786.7
interporder=4
interpkind='cubic'

finfi=np.zeros([len(orders)-2,31-start])


for modenum in range(start,31,step):
     print "modenum = ", modenum
     datatablelist=list(np.zeros(0))
     tstoredlist=list(np.zeros(0))
     lbestarr=np.zeros([len(orders)])
     psirarr=np.zeros([len(orders),31-start])
     lstoredlist=list(np.zeros(0))
     for count in range(len(orders)):
         tnearest=0.0
         indexnearest=0
         lbest=0
         tstored = list(np.zeros(interporder))
         lstored=list(np.zeros(interporder))
         loadstring=""
         if((orders[count]==48) or (orders[count]==28)):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_restart/psir_l.asc"
         elif((orders[count]==36) or (orders[count]==44) or (orders[count]==48)):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_2_restart/psir_l.asc"
         elif(orders[count]==52):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_3_rerestart/psir_l.asc"
         elif(orders[count]==16):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_2/psir_l.asc"
         elif(orders[count]==20):
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"_restart/psir_l.asc"
         else:
             loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"/psir_l.asc"

         print loadstring


         
         datatable=np.loadtxt(loadstring,skiprows=1)
         for ii in range(len(datatable[:,timecolumn])):
             if datatable[ii,timecolumn]<t0:
                 tnearest=datatable[ii,timecolumn]
                 indexnearest=ii
         #print orders[count], indexnearest, interporder, len(datatable[:,timecolumn])
         for ii in range(interporder):
             tstored[ii]=datatable[indexnearest-(interporder-1)/2+ii,timecolumn]
             lstored[ii]=datatable[indexnearest-(interporder-1)/2+ii, columnoffset+modenum]
         func=interp1d(tstored,lstored,kind=interpkind)
         lbest=func(t0)
         lbestarr[count]=lbest
         psirarr[count,modenum-start]=lbest    
         tstoredlist.append(tstored)
         lstoredlist.append(lstored)
         
     for i1 in range(0,len(orders)-3):
         i2=i1+1
         i3=i1+2

         yratio=(lbestarr[i1]-lbestarr[i2])/(lbestarr[i1]-lbestarr[i3])
         print lbestarr[i1], lbestarr[i2], lbestarr[i3], yratio
         alpha0=0.5
     
         alphamax=alpha0
         alphamin=1e-12
         ratiofnreturn=-1.
         while (ratiofnreturn<0.):
             alphamax*=1.5
             ratiofnreturn=ratiofunc(alphamax,orders[i1],orders[i2],orders[i3],yratio)
             print ratiofnreturn, alphamax, yratio, orders[i1], orders[i2], orders[i3]

         rationfnreturn=1.
         while(ratiofnreturn>0.):
             alphamin/=2
             ratiofnreturn=ratiofunc(alphamin,orders[i1],orders[i2],orders[i3],yratio)
             print ratiofnreturn, alphamin, yratio
         
         alpha =optimization.bisect(ratiofunc,alphamin,alphamax,args=(orders[i1],orders[i2],orders[i3],yratio))
     
     
         print "alpha= ", alpha
     
     
         ccoeff = (lbestarr[i1]-lbestarr[i2])/(exp(-alpha*orders[i1])-exp(-alpha*orders[i2]))
         finf = lbestarr[i3]+ccoeff*exp(-alpha*orders[i3])
     
         print alpha, ccoeff, finf
         finfi[i1,modenum]=finf

for i1 in range(len(orders)-2):         
     csvwriter=csv.writer(fio2,delimiter=' ')
     csvwriter.writerow(np.append([orders[i1]],finfi[i1,:]))
fio2.close()

