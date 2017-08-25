from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.optimize as optimization
import csv
import pylab
import math
import sys

#def truncfunc(Flclalpha, xdata, ydata):
#    Fl=Flclalpha[0]
#    cl=Flclalpha[1]
#    alpha=Fclapha[2]
#    return(ydata-cl*exp(-alpha*xdata))

def ratiofunc(alpha, n1, n2, n3, hratio):

    return (exp(alpha*(n2-n1))-1)/(1-exp(alpha*(n2-n3)))-hratio
    #return exp(4*alpha)/(1+exp(4*alpha))-hratio
    
def main(t0in,startorder,showplot,writedata,start,stop):
 
    finf =0
    ccoeff =0
    alpha =0

    #either show plot or write file    
    
    t0=t0in
    orders=[16, 20, 24,28, 32,36,40,44] #not 48,33, 12, 56
    i10=orders.index(startorder)
    i20=i10+1
    i30=i10+2
    step = 1
    fio=open("placeholder.csv","a")
    if(writedata==1):
        print "opening output file"
        fio=open("coeffsbyl"+str(t0)+"_"+str(orders[i10])+"_"+str(orders[i20])+"_"+str(orders[i30])+".csv","a")
    csvwriter=csv.writer(fio,delimiter=' ')
    print "defining csvwriter"    
    columnoffset=5
    timecolumn=0
    nummodes = 6
    datatable =np.loadtxt("/mnt/data/sdorsher/Fortranp9.9e0.1n40/psir_l.asc", skiprows=1)


    i1=i10
    i2=i20
    i3=i30

    loadstring=""
    loadstring="genrawdata"+str(t0)+".csv"
    print loadstring
    
    datatable=np.loadtxt(loadstring)
    print "genrawdata loaded"
    psirarr=datatable
    
    outputarr = np.zeros([datatable.shape[0],5])
    for modenum in range(start,stop,step):
        lbestarr=datatable[modenum,1:]
        print "modenum = ", modenum
        orderspred=orders[i3+1:]
        datatablelist=list(np.zeros(0))
        tstoredlist=list(np.zeros(0))
        #lbestarr=np.zeros([len(orders)])
        #psirarr=np.zeros([len(orders),stop-start])
        lstoredlist=list(np.zeros(0))

        #psirarr=datatable    
        
        
        
        hratio=(lbestarr[i1]-lbestarr[i2])/(lbestarr[i2]-lbestarr[i3])
        print "hratio=", lbestarr[i1], lbestarr[i2], lbestarr[i3], hratio

        alpha0=0
        if(hratio>0.):
            alpha0=0.25*log(hratio)
        else:
            alpha0=1.
        alphamax=alpha0
        alphamin=alpha0
        ratiofnreturn=-1.
        if(hratio>1.0):
            while (ratiofnreturn<0.):
                alphamax*=1.5
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
            finf=np.nan
            ccoeff=np.nan
            alpha=np.nan
            print "Mode and order failed!"
        outputarr[modenum,0]=t0
        outputarr[modenum,1]=modenum
        outputarr[modenum,2]=alpha
        outputarr[modenum,3]=ccoeff
        outputarr[modenum,4]=finf
            
        if (writedata==1):
            if(hratio>1.):
                csvwriter.writerow([t0, modenum,alpha,ccoeff,finf])
                print "writing to file"
        lpred= np.zeros(len(orderspred))
        lbestnew=np.zeros(len(lbestarr))
        for ii in range(len(lbestarr)):
            lbestnew[ii]=abs(lbestarr[ii]-finf)
        if(hratio>0.5 and hratio<1.0):
            for ii in range(len(orderspred)):
                lpred[ii]=abs(ccoeff*exp(-alpha*orderspred[ii]))
        if(showplot==1):
            plt.semilogy(orders,lbestnew,marker='o',label='Data')
            if(hratio>0.5 and hratio<1.0):
                plt.plot(orderspred,lpred,marker='^',label='Predicted')
            plt.semilogy(orders[i1:i3+1],lbestnew[i1:i3+1],'ro',label="Points used in extrapolation")
            ax=plt.gca()
            #ax.set_yscale('log')
            ax.legend(loc='lower left')
            plt.xlabel('DG order')
            plt.ylabel('Radial self force minus Finf')
            plt.title('l='+str(modenum)+", extrapolated from orders "+str(orders[i1])+", "+str(orders[i2])+", and "+str(orders[i3]))
            plt.show()    
   
    fio.close()
    #for count in range(len(orders)):
    #    for modenum in range(start,stop,step):
           #csvwriter2=csv.writer(fio2,delimiter=' ')
           #csvwriter2.writerow(np.append([orders[count]],psirarr[count,:]))
           
           #fio2.close()
    if(stop-start==1):
        return alpha, ccoeff, finf
    else:
        return outputarr
   
if __name__=="__main__":
    if (len(sys.argv)<5):
        print "Usage python extrapolate7.py t0 startorder showplot writedata lchosen(opt)"
        print "t0 is time"
        print "startorder is the starting order of the extrapolation"
        print "showplot is 1 to produce an exponentially decaying plot"
        print "of DG order versus psir-finf and is 0 to write data for all modes"
        print "writedata is 1 to write data to a csv file coeffsbyl*_*_*_*.csv"
        print "lchosen is the specificmode chosen (optional)." 
        print "by default runs from 0 to lmax"
        print "requires interpolated raw data generated by rawdatagen.py"
        exit()
    t0in=int(sys.argv[1])
    startorder=int(sys.argv[2])
    showplot = int(sys.argv[3])
    writedata = int(sys.argv[4])
    lchosen=0
    if(len(sys.argv)==6):
        lchosen=int(sys.argv[5])
    start = 0
    stop = 31
    if (len(sys.argv)==6): #select chosen l if 
        start=lchosen
        stop=lchosen+1

    main(t0in,startorder,showplot,writedata,start,stop)
