from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.optimize as optimization
import csv
import pylab
import sys

def main(argv):
     if(len(sys.argv)!=2):
          print "Usage python interpolateandprint.py t0"
          print "t0 is the time to which the data should be interpolated"
          exit()
     t0=int(sys.argv[1])
     #def func(n, alpha, ccoeff, finf):
     #    return finf-ccoeff*np.exp(-alpha*n)
     orders=[12,16, 20, 24,28, 32,36,40,44] #not 48,33
     start=0
     stop=31
     #i10=1
     #i20=i10+1
     #i30=i10+2
     step = 1
     #fio=open("coeffsbyl"+str(t0)+"_"+str(orders[i10])+"_"+str(orders[i20])+"_"+str(orders[i30])+".csv","a")
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
     
     #i1=i10
     #i2=i20
     #i3=i30
     
     
     
     fio=open("rawinterpdatat"+str(t0)+".csv", "a")
     
     for modenum in range(start,stop,step):
          print "modenum = ", modenum
          #if (modenum<0):
              #i1=i10+1
              #i2=i20+1
              #i3=i30+1
          #orderspred=orders[i3+1:] #20 defective? #52 still running, changedir
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
              print datatable.shape, t0
              for ii in range(len(datatable[:,timecolumn])):
                   if (datatable[ii,timecolumn]<t0):
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
     
         
          csvwriter=csv.writer(fio,delimiter=' ')
          csvwriter.writerow(np.append(np.array([modenum]),lbestarr))    
     
if __name__=="__main__":
    main(sys.argv[:1])
