from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.optimize as optimization
import csv
import pylab




t0=380
orders=[12,16, 20, 24,28, 32,36,40,44] #not 48,33
stop=30
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

lbestarr=np.zeros([len(orders),stop+1])
psirarr=np.zeros([len(orders),stop+1])

datatablelist=list(np.zeros(0))
for count in range(0,len(orders)):
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

    for modenum in range(0,stop,step):
        print "modenum = ", modenum
        tnearest=0.0
        indexnearest=0
        lbest=0
        tstored = list(np.zeros(interporder))
        lstored=list(np.zeros(interporder))
        loadstring=""

        
        for ii in range(len(datatable[:,timecolumn])):
            if datatable[ii,timecolumn]<t0:
                tnearest=datatable[ii,timecolumn]
                indexnearest=ii
               
        for ii in range(interporder):
            tstored[ii]=datatable[indexnearest-(interporder-1)/2+ii,timecolumn]
            lstored[ii]=datatable[indexnearest-(interporder-1)/2+ii, columnoffset+modenum]
        func=interp1d(tstored,lstored,kind=interpkind)
        lbest=func(t0)
        lbestarr[count,modenum]=lbest
        psirarr[count,modenum]=lbest    

sumOverModes=np.sum(psirarr[count,:])

#datatowrite=np.append([modenum],psirarr[:,modenum])
datatowrite=psirarr
datatowrite.tofile("genrawdata"+str(t0)+".csv",sep=' ', format='%10.5f')


print sumOverModes
