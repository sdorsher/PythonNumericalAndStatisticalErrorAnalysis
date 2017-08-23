from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.optimize as optimization
import csv
import pylab
import math

#def truncfunc(Flclalpha, xdata, ydata):
#    Fl=Flclalpha[0]
#    cl=Flclalpha[1]
#    alpha=Fclapha[2]
#    return(ydata-cl*exp(-alpha*xdata))

def ratiofunc(alpha, n1, n2, n3, hratio):

    return (exp(alpha*(n2-n1))-1)/(1-exp(alpha*(n2-n3)))-hratio
    #return exp(4*alpha)/(1+exp(4*alpha))-hratio
    
    
   
#def func(n, alpha, ccoeff, finf):
#    return finf-ccoeff*np.exp(-alpha*n)
t0=450
#steps of four
#orders=[12,16, 20, 24,28, 32,36,40,44] #not 48,33
orders=[16, 20, 24,28, 32,36,40,44] #not 48,33
start=3
stop=start+1
i10=0
i20=i10+1
i30=i10+2
step = 1
fio=open("coeffsbyl"+str(t0)+"_"+str(orders[i10])+"_"+str(orders[i20])+"_"+str(orders[i30])+".csv","a")
#fio2=open("selfforceallmodes.csv","a")

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

i1=i10
i2=i20
i3=i30


for modenum in range(start,stop,step):
     print "modenum = ", modenum
     if (modenum<0):
         i1=i10+1
         i2=i20+1
         i3=i30+1
     orderspred=orders[i3+1:] #20 defective? #52 still running, changedir
     datatablelist=list(np.zeros(0))
     tstoredlist=list(np.zeros(0))
     lbestarr=np.zeros([len(orders)])
     psirarr=np.zeros([len(orders),stop-start])
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
         print orders[count]
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
     
     hratio=(lbestarr[i1]-lbestarr[i2])/(lbestarr[i2]-lbestarr[i3])
     print "hratio=", lbestarr[i1], lbestarr[i2], lbestarr[i3], hratio

     #alpha=0.25*log(hratio)
     alpha0=0.25*log(hratio)
     
     alphamax=alpha0
     alphamin=alpha0
     ratiofnreturn=-1.
     if(hratio>-1.0):
     #if(hratio>0.5 and hratio<1.0): 
         while (ratiofnreturn<0.):
             alphamax*=1.5
             #print alphamax, ratiofnreturn
             ratiofnreturn=ratiofunc(alphamax,orders[i1],orders[i2],orders[i3],hratio)
         print ratiofnreturn, alphamax, hratio, orders[i1], orders[i2], orders[i3]

         while (ratiofnreturn>0.):
             alphamin/=2.
             ratiofnreturn=ratiofunc(alphamin,orders[i1],orders[i2],orders[i3],hratio)
         print ratiofnreturn, alphamin, hratio
     
         alpha =optimization.bisect(ratiofunc,alphamin,alphamax,args=(orders[i1],orders[i2],orders[i3],hratio))
         
     
         print "alpha= ", alpha
    
         
         ccoeff = (lbestarr[i1]-lbestarr[i2])/(exp(-alpha*orders[i1])-exp(-alpha*orders[i2]))
         finf = lbestarr[i3]-ccoeff*exp(-alpha*orders[i3])
         
         print alpha, ccoeff, finf
     else:
         finf=math.nan
         #lbestarr[len(orders)-1]-8*10**-14
         print "Mode failed!"
     csvwriter=csv.writer(fio,delimiter=' ')
     if(hratio>0.5 and hratio<1.0):
         csvwriter.writerow([t0, modenum,alpha,ccoeff,finf])
     lpred= np.zeros(len(orderspred))
     lbestnew=np.zeros(len(lbestarr))
     for ii in range(len(lbestarr)):
         lbestnew[ii]=abs(lbestarr[ii]-finf)
     if(hratio>0.5 and hratio<1.0):
         for ii in range(len(orderspred)):
             lpred[ii]=abs(ccoeff*exp(-alpha*orderspred[ii]))

     #plt.semilogy(orders,lbestnew,marker='o',label='Data')
     #if(hratio>0.5 and hratio<1.0):
         #plt.plot(orderspred,lpred,marker='^',label='Predicted')
     #plt.semilogy(orders[i1:i3+1],lbestnew[i1:i3+1],'ro',label="Points used in extrapolation")
     #ax=plt.gca()
     #ax.set_yscale('log')
     #ax.legend(loc='lower left')
     #plt.xlabel('DG order')
     #plt.ylabel('Radial self force minus Finf')
     #plt.title('l='+str(modenum)+", extrapolated from orders "+str(orders[i1])+", "+str(orders[i2])+", and "+str(orders[i3]))
     #plt.show()    

fio.close()
#for count in range(len(orders)):
#    for modenum in range(start,stop,step):
        #csvwriter2=csv.writer(fio2,delimiter=' ')
        #csvwriter2.writerow(np.append([orders[count]],psirarr[count,:]))
        
#fio2.close()

