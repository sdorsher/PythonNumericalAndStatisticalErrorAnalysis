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
    return exp(4*alpha)/(1+exp(4*alpha))-yratio
    #return exp(-alpha*n2)*(exp(alpha*n1)-exp(n2*alpha))/(exp(alpha*(n1-n3))-1)-yratio
    
   
#def func(n, alpha, ccoeff, finf):
#    return finf-ccoeff*np.exp(-alpha*n)
t0=610
orders=[12,16, 20, 24,28, 32,36,40,44] #not 48,33
start=0
stop=31
i10=0

step = 1
#t0=472.721
#t0=654
#t0=500
#t0=786.7
interporder=4
interpkind='cubic'


columnoffset=5
timecolumn=0
#modenum=6
#modenum=30
nummodes = 6
datatable=[]
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
    elif(orders[count]==12 or  orders[count]==32 or orders[count]==40):
        loadstring="/mnt/data/sdorsher/Fortranp9.9e0.1n"+str(orders[count])+"/psir_l.asc"
        
    print loadstring
        
    dat=np.loadtxt(loadstring,skiprows=1)
    datatable.append(dat)

    
#orders=[ 24,28, 32,36,40,44] #not 48,33
#[8,16]
#orderspred=[36,40,44] #20 defective? #52 still running, changedir



for modenum in range(start,stop,step):
    finfarr=np.zeros(len(orders)-2)
    for i1 in range(0,len(orders)-2):
         i2=i1+1
         i3=i1+2
         fio=open("coeffsbestfinf"+str(t0)+"_"+str(orders[i1])+"_"+str(orders[i2])+"_"+str(orders[i3])+".csv","a")
         print "modenum = ", modenum
         orderspred=orders[i3+1:] #20 defective? #52 still running, changedir
         datatablelist=list(np.zeros(0))
         tstoredlist=list(np.zeros(0))
         lbestarr=np.zeros([len(orders)])
         psirarr=np.zeros([len(orders),stop-start])
         lstoredlist=list(np.zeros(0))
         for count in range(len(orders)):
             
             for ii in range(len(datatable[count][:,timecolumn])):
                 if datatable[count][ii,timecolumn]<t0:
                     tnearest=datatable[count][ii,timecolumn]
                     indexnearest=ii
                     #print orders[count], indexnearest, interporder, len(datatable[:,timecolumn])
             for ii in range(interporder):
                tstored[ii]=datatable[count][indexnearest-(interporder-1)/2+ii,timecolumn]
                lstored[ii]=datatable[count][indexnearest-(interporder-1)/2+ii, columnoffset+modenum]
             func=interp1d(tstored,lstored,kind=interpkind)
             lbest=func(t0)
             lbestarr[count]=lbest
             psirarr[count,modenum-start]=lbest    
             tstoredlist.append(tstored)
             lstoredlist.append(lstored)
        
         yratio=(lbestarr[i1]-lbestarr[i2])/(lbestarr[i1]-lbestarr[i3])
         print "yratio=", lbestarr[i1], lbestarr[i2], lbestarr[i3], yratio
         alpha0=0.5
          
         alphamax=alpha0
         alphamin=1.e-12
         ratiofnreturn=-1.
         if(yratio>0.5 and yratio<1.0):
             while (ratiofnreturn<0.):
                 alphamax*=1.5
                 ratiofnreturn=ratiofunc(alphamax,orders[i1],orders[i2],orders[i3],yratio)
             print ratiofnreturn, alphamax, yratio, orders[i1], orders[i2], orders[i3]
     
             while (ratiofnreturn>0.):
                 alphamin/=2.
                 ratiofnreturn=ratiofunc(alphamin,orders[i1],orders[i2],orders[i3],yratio)
             print ratiofnreturn, alphamin, yratio
                 
             alpha =optimization.bisect(ratiofunc,alphamin,alphamax,args=(orders[i1],orders[i2],orders[i3],yratio))
             #fprime=ratiofuncprime, tol=1e-14,args=(orders[i1],orders[i2],orders[i3],yratio),fprime2=None
            
             print "alpha= ", alpha
            
          
             ccoeff = (lbestarr[i1]-lbestarr[i2])/(exp(-alpha*orders[i1])-exp(-alpha*orders[i2]))
             finf = lbestarr[i3]-ccoeff*exp(-alpha*orders[i3])
             print alpha, ccoeff, finf
         else:
             alpha=np.nan
             ccoeff=np.nan
             finf=np.nan
             print "Mode failed!"
         finfarr[i1]=finf
         csvwriter=csv.writer(fio,delimiter=' ')
         csvwriter.writerow([t0, modenum,alpha,ccoeff,finf])
         fio.close()
    plt.plot(orders[:len(orders)-2],finfarr,'o-')
    ax=plt.gca()
    ax.set_xlabel("Starting order of extrapolation")
    ax.set_ylabel("Finf")
    plt.title("Infinite order self force for l="+str(modenum)+", t="+str(t0))
    plt.show()
 

