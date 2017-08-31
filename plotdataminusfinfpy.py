import matplotlib.pyplot as plt
import numpy as np
from math import *
import pylab

def plotdataminusfinf(orders,lbestarr,finf,i1,i2,i3,modenum):
    lbestnew=np.zeros(len(lbestarr))
    for ii in range(len(lbestarr)):
        lbestnew[ii]=abs(lbestarr[ii]-finf)
        #if(hratio>1.):
        #for ii in range(len(orderspred)):
        #lpred[ii]=abs(ccoeff*exp(-alpha*orderspred[ii]))
        
    plt.semilogy(orders,lbestnew,marker='o',label='Data')
    #if(hratio>0.5 and hratio<1.0):
    #plt.plot(orderspred,lpred,marker='^',label='Predicted')
    plt.semilogy(orders[i1:i3+1],lbestnew[i1:i3+1],'ro',label="Points used in extrapolation")
    ax=plt.gca()
    #ax.set_yscale('log')
    ax.legend(loc='lower left')
    plt.xlabel('DG order')
    plt.ylabel('Radial self force minus Finf')
    plt.title('l='+str(modenum)+", extrapolated from orders "+str(orders[i1])+", "+str(orders[i2])+", and "+str(orders[i3]))
    plt.show()    
