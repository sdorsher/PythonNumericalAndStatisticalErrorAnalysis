from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.optimize as optimization
import csv
import pylab
import sys, getopt, cmath, os
from extrapolate7 import main as extrap7main


#def truncfunc(Flclalpha, xdata, ydata):
#    Fl=Flclalpha[0]
#    cl=Flclalpha[1]
#    alpha=Fclapha[2]
#    return(ydata-cl*exp(-alpha*xdata))

def ratiofunc(alpha, n1, n2, n3, yratio):

    return (exp(alpha*(n2-n1))-1)/(1-exp(alpha*(n2-n3)))-yratio
    #return exp(4*alpha)/(1+exp(4*alpha))-yratio
    #return exp(-alpha*n2)*(exp(alpha*n1)-exp(n2*alpha))/(exp(alpha*(n1-n3))-1)-yratio
    
def main(argv):
    
    #def func(n, alpha, ccoeff, finf):
    #    return finf-ccoeff*np.exp(-alpha*n)

    if(len(sys.argv)!=2):
        print "Usage bestfinfselector.py t0"
        print "\tt0 is the start time"
        print "Produces a table of mode vs order chosen and Finf selected"
        print "Chooses a mode using the highest order with valid Finf values"
        exit()

    print sys.argv[0], sys.argv[1]
    t0=int(sys.argv[1])
    print t0
    orders=[16, 20, 24,28, 32,36,40,44] #not 48,33
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
        for i1 in range(0,len(orders)-2):
            i2=i1+1
            i3=i1+2
            #loadtxt="coeffsbyl"+str(t0)+"_"+str(orders[i1])+"_"+str(orders[i2])+"_"+str(orders[i3])+".csv","a"
            alpha, ccoeff, finf = extrap7main(t0,orders[i1],0,0,modenum,modenum+1)
            print "modenum = ", modenum, " start order = ", i1
            #orderspred=orders[(i3+1):] #20 defective? #52 still running, changedir
 

            best_chisq_dof=0
            best_start=0
            best_end=0
            if(~np.isnan(finf)):
                for lbeststartindex in range(len(lbestarr)-2):
                    for lbestendindex in range(lbeststartindex+3,len(lbestarr)+1):
                        orderstrunc=orders[lbeststartindex:lbestendindex]
                        lbestarrtrunc=lbestarr[lbeststartindex:lbestendindex]
                        loglbestarrtrunc=np.log(abs(lbestarrtrunc-finf))
                         
                        forders,residuals,rank,sing,rcond=np.polyfit(orderstrunc,loglbestarrtrunc,1,full=True)
                        #plt.xlabel("DG order")
                        #plt.ylabel("Log(Radial self force)")
                        #plt.title("Determining the best starting index by fitting line segments")
                        #plt.plot(orderstrunc,loglbestarrtrunc, marker='o',label="Data")
                        #plt.plot(orderstrunc,np.polyval(forders,orderstrunc), label="Linear least squares fit")
                        #plt.legend(loc='lower left')
                        #plt.show()
                        chisq_dof=residuals/(len(orderstrunc)-2)
                        #print chisq_dof
                        if (abs(best_chisq_dof-1))>(abs(chisq_dof-1)):
                            best_chisq_dof=chisq_dof
                            best_start=lbeststartindex
                            best_end=lbestendindex
                            #print best_chisq_dof, best_start,best_end
            if (abs(overall_best_chisq_dof-1))>(abs(best_chisq_dof-1)):
                overall_best_chisq_dof=best_chisq_dof
                best_i1=i1
                
            finfarr[i1]=finf
            #csvwriter=csv.writer(fio,delimiter=' ')
            #csvwriter.writerow([t0, modenum,alpha,ccoeff,finf])
            #fio.close()
        #plt.plot(orders[:len(orders)-2],finfarr,'o-')
        #ax=plt.gca()
        #ax.set_xlabel("Starting order of extrapolation")
        #ax.set_ylabel("Finf")
        #plt.title("Infinite order self force for l="+str(modenum)+", t="+str(t0))
        #plt.show()

        ydata=np.array(np.abs(datatable[modenum,1:]-finfarr[best_i1]))
        print ydata.shape
        ydata2=ydata.flatten()
        print ydata2.shape
        orders2=np.array(orders)
        print orders2.shape
        orders3=orders2.flatten()
        print orders3.shape
        plt.semilogy(orders3,ydata2,'o-',label="Best choice Psir-Finf")
        plt.semilogy(orders3[best_i1:best_i1+3],ydata2[best_i1:best_i1+3],'ro',label="Points used in extrapolation")
        ax=plt.gca()
        ax.set_xlabel("DG order")
        ax.set_ylabel("|Psir-Finf|")
        ax.legend(loc='lower left')
        plt.title("Radial self force with Finf subtracted for l="+str(modenum)+", t="+str(t0))
        plt.show()


        
        ii=0
        finfsum=0.
        sumcount=0
        #take median of array:
        
        finfs=finfarr[~np.isnan(finfarr)]
        finfssort=np.sort(finfs)
        lenfinf=len(finfssort)
        mid=int(int(lenfinf)/int(2))
        oddfinf=lenfinf % 2
        bestfinf=0.
        if(oddfinf==1):
            bestfinf=finfssort[mid]
        else:
            bestfinf=(finfssort[mid]+finfssort[mid-1])/2.
            
        print modenum, min(finfssort), max(finfssort), bestfinf, best_i1, finfarr[best_i1]
        csvwriter3=csv.writer(fio3,delimiter=' ')
        csvwriter3.writerow([modenum, min(finfssort), max(finfssort), bestfinf, best_i1, finfarr[best_i1]])
    fio3.close()
    
if __name__=="__main__":
    main(sys.argv[:1])
