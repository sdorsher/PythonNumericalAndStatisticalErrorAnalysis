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
        exit()

    print sys.argv[0], sys.argv[1]
    t0=int(sys.argv[1])
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

        print finfarr
        #find the starting index of the longest contiguous sub array
        start_i_sub=[]
        end_i_sub=[]
        truish = ~np.isnan(finfarr)
        print truish
        for ii in range(len(finfarr)-1):
            if(ii==0 and truish[1]):
                start_i_sub.append(0)
            elif ((ii==len(finfarr)-2) and truish[ii+1]):
                end_i_sub.append(ii+1)
            elif( (not truish[ii]) and truish[ii+1]):
                start_i_sub.append(ii+1)
            elif (truish[ii] and  (not truish[ii+1])):
                end_i_sub.append(ii)

        #print start_i_sub
        #print end_i_sub
                
        if (len(start_i_sub)!=len(end_i_sub)):
            print len(start_i_sub), len(end_i_sub), "start and end don't match"
            exit()

        #now that I've found the pairs of sub array indices, find the longest such pair
        len_sub=0
        best_start_sub=0
        best_end_sub=0
        for ii in range(len(start_i_sub)):
            len_this_sub=end_i_sub[ii]-start_i_sub[ii]+1
            if(len_this_sub>len_sub):
                len_sub=len_this_sub
                best_start_sub=start_i_sub
                best_end_sub=end_i_sub
        

        deriv=1
        secondderiv=-1
        best_i1=0
        jj=best_start_sub
        nanbit=False

        ydata=np.array(np.abs(datatable[modenum,1:]-finfarr[best_i1]))
        ydata2=ydata.flatten()
        orders2=np.array(orders)
        orders3=orders2.flatten()
        
        print best_start_sub, best_end_sub
        while ((deriv*secondderiv<0) and jj<=best_end_sub):
            #if(~isnan(finfarr[ii]) and ~isnan(finfarr[ii-1]) and ~isnan(finfarr[ii+1])):
            print ydata2.shape, orders3.shape, ydata2[jj]
            deriv=(ydata2[jj]-ydata2[jj-1])/(orders3[jj]-orders3[jj-1])
            secondderiv=(ydata2[jj+1]-2*ydata2[jj]+ydata2[jj-1])/(orders3[jj]-orders3[jj-1])**2
            best_i1=jj
            jj+=1
                    
        #plt.plot(orders[:len(orders)-2],finfarr,'o-')
        #ax=plt.gca()
        #ax.set_xlabel("Starting order of extrapolation")
        #ax.set_ylabel("Finf")
        #plt.title("Infinite order self force for l="+str(modenum)+", t="+str(t0))
        #plt.show()


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
