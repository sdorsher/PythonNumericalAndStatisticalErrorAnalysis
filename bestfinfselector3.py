from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.optimize as optimization
import csv
import pylab
import sys, getopt, cmath, os
from extrapolate7 import main as extrap7main
import numpy.ma as ma
from plotdataminusfinfpy import plotdataminusfinf

#masked array

#def truncfunc(Flclalpha, xdata, ydata):
#    Fl=Flclalpha[0]
#    cl=Flclalpha[1]
#    alpha=Fclapha[2]
#    return(ydata-cl*exp(-alpha*xdata))


def case5(best_start_sub, best_end_sub, fdata,forders):
    best_i1=0
    deriv=1
    secondderiv=-1
    jj=best_start_sub+1
    while ((deriv*secondderiv<0) and jj<best_end_sub):
        #if(~isnan(finfarr[ii]) and ~isnan(finfarr[ii-1]) and ~isnan(finfarr[ii+1])):
        deriv=(fdata[jj+1]-fdata[jj-1])/(forders[jj+1]-forders[jj-1])*2
        secondderiv=(fdata[jj+1]-2*fdata[jj]+fdata[jj-1])*4/(forders[jj+1]-forders[jj-1])**2
        best_i1=jj
        jj+=1
        print deriv, secondderiv
        return best_i1
def case4(best_start_sub, best_end_sub, fdata):
    
    
    print "case 4"
    datatoaverage=fdata[best_start_sub:best_end_sub+1]
    avg1 = np.average(datatoaverage)
    std1=np.std(datatoaverage)
    print avg1, std1
    mymask = (np.abs(datatoaverage-avg1)>1*std1)
    datavetoedtoavg=ma.masked_array(data=datatoaverage, mask=mymask)
    avg2=datavetoedtoavg.mean()
    #find nearest value see bookmarked page
    #count on the fact that finf will not be order 1 and fill with 1's
    vetoeddata=datavetoedtoavg.filled([1.])
    print vetoeddata
    #select best index based on smallest deviation from average based on vetoed data
    best_i1=(np.abs(vetoeddata-avg2)).argmin()+best_start_sub
    print avg2, best_i1
    return best_i1

def ratiofunc(alpha, n1, n2, n3, yratio):

    return (exp(alpha*(n2-n1))-1)/(1-exp(alpha*(n2-n3)))-yratio
    #return exp(4*alpha)/(1+exp(4*alpha))-yratio
    #return exp(-alpha*n2)*(exp(alpha*n1)-exp(n2*alpha))/(exp(alpha*(n1-n3))-1)-yratio
    
def main(argv):
    
    #def func(n, alpha, ccoeff, finf):
    #    return finf-ccoeff*np.exp(-alpha*n)

    if(len(sys.argv)!=4):
        print "Usage bestfinfselector.py t0 showplot1 showplot2"
        print "\tt0 is the start time"
        print "Produces a table of mode vs order chosen and Finf selected"
        print "showplot1: plots finf vs dg order"
        print "showplot2: plots |Psir-finf| vs dg order"
        exit()

    print sys.argv[0], sys.argv[1]
    t0=int(sys.argv[1])
    showplot1=int(sys.argv[2])
    showplot2=int(sys.argv[3])
    print t0
    orders=np.array([16, 20, 24,28, 32,36,40,44]) #not 48,33
    start=0
    stop=4
    i10=0
    
    step = 1
    fio3=open("bestinfoverinitorder"+str(t0)+".csv","a")
    loadstring = "genrawdata"+str(t0)+".csv"

    print loadstring
    datatable=np.loadtxt(loadstring)
    if(len(orders)!=datatable.shape[1]-1):
        print "wrong orders included in genrawdata"
        print datatable.shape
        print len(orders)
        exit()
    
    
    columnoffset=5
    timecolumn=0
    nummodes = 6
    for modenum in range(datatable.shape[0]):
        lbestarr=datatable[modenum,1:]
        finfarr=np.zeros(len(orders)-2)
        overall_best_chisq_dof=0
        overall_best_i1=0
        finfs=np.zeros(len(orders)-2)
        for i1 in range(0,len(orders)-2):
            i2=i1+1
            i3=i1+2
            #loadtxt="coeffsbyl"+str(t0)+"_"+str(orders[i1])+"_"+str(orders[i2])+"_"+str(orders[i3])+".csv","a"
            alpha, ccoeff, finf = extrap7main(t0,orders[i1],0,0,modenum,modenum+1)
            print "modenum = ", modenum, " start order = ", i1
            #orderspred=orders[(i3+1):] #20 defective? #52 still running, changedir
 
            finfarr[i1]=finf

        #find the starting index of the longest contiguous sub array
        start_i_sub=[]
        end_i_sub=[]
        truish = ~np.isnan(finfarr)
        print truish
        print finfarr
        for ii in range(len(finfarr)-1):
            if(ii==0 and truish[0]):
                start_i_sub.append(0)
            if ((ii==len(finfarr)-2) and truish[ii+1]):
                end_i_sub.append(ii+1)
            if( (not truish[ii]) and truish[ii+1]):
                start_i_sub.append(ii+1)
            if (truish[ii] and  (not truish[ii+1])):
                end_i_sub.append(ii)

        #print start_i_sub
        #print end_i_sub
        print start_i_sub
        print end_i_sub

        
        if (len(start_i_sub)!=len(end_i_sub)):
            print len(start_i_sub), len(end_i_sub), "start and end don't match"
            exit()

        #now that I've found the pairs of sub array indices, find the longest such pair

        len_sub=0
        best_start_sub=0
        best_end_sub=0
        for ii in range(len(start_i_sub)):
            len_this_sub=end_i_sub[ii]-start_i_sub[ii]+1
            if(len_this_sub>=len_sub):
                len_sub=len_this_sub
                best_start_sub=start_i_sub[ii]
                best_end_sub=end_i_sub[ii]
        


        

        fdata=finfarr.flatten()
        orders2=np.array(orders)
        orders3=orders2.flatten()        
        forders=orders2[1:len(orders)-1]

        print best_start_sub, best_end_sub, len_sub
        nonepassed=False
        if (len_sub >= 5):
            print "case 5 or more"
            best_i1=case5(best_start_sub, best_end_sub, fdata,forders)
            print best_i1
            print len(finfarr)
            if (best_i1-1<len(finfarr)/2):
                best_i1=case4(best_start_sub, best_end_sub, fdata)
        elif (len_sub==1):
            print "case 1"
            best_i1=best_end_sub
        elif(len_sub==2):
            print "case 2"
            best_i1=best_end_sub
        elif(len_sub==3):
            print "case 3"
            best_i1=case4(best_start_sub, best_end_sub,fdata)
            #best_i1=(best_end_sub+best_start_sub)/2
        elif(len_sub==4):
            print "case 4"
            #find average value and take value closest to average (this is case like plateau usually)
            best_i1=case4(best_start_sub, best_end_sub, fdata)
        else:
            print "no order passed for this mode, l=" + str(modenum)+ ", t=" +str(t0)
            best_i1=-1
            nonepassed=True
        print best_i1
        if(showplot1==1):
            plt.plot(orders[:len(orders)-2],finfarr,'o-')
            ax=plt.gca()
            ax.set_xlabel("Starting order of extrapolation")
            ax.set_ylabel("Finf")
            plt.title("Infinite order self force for l="+str(modenum)+", t="+str(t0))
            plt.show()


        if((showplot2==1) and (not nonepassed)):
           
            bestfinf=finfarr[best_i1]
            plotdataminusfinf(orders,lbestarr,bestfinf,best_i1,best_i1+1,best_i1+2,modenum)
        
        #ydata=np.array(np.abs(datatable[modenum,1:]-finfarr[best_i1-1]))
        #ydata2=ydata.flatten()
        #plt.semilogy(orders3,ydata2,'o-',label="Best choice Psir-Finf")
        #plt.semilogy(orders3[best_i1-1:best_i1+2],ydata2[best_i1-1:best_i1+2],'ro',label="Points used in extrapolation")
        #ax=plt.gca()
        #ax.set_xlabel("DG order")
        #ax.set_ylabel("|Psir-Finf|")
        #ax.legend(loc='lower left')
        #plt.title("Radial self force with Finf subtracted for l="+str(modenum)+", t="+str(t0))
        #plt.show()


        finfout='NaN'
        if(not nonepassed):
            finfout=finfarr[best_i1]
        print modenum, len_sub, best_i1, finfout
        csvwriter3=csv.writer(fio3,delimiter=' ')
        csvwriter3.writerow([modenum, len_sub, best_i1, finfout])
    fio3.close()
    
if __name__=="__main__":
    main(sys.argv[:1])
