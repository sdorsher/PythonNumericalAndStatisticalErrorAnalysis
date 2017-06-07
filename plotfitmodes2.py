
from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D
import csv

#def fit_func(ldata, param1,param2,param3,param4):
#    psir=param1/(2*ldata-1)/(2*ldata+3)+param2/(2*ldata-3)/(2*ldata-1)/(2*ldata+3)/(2*ldata+5)+param3/(2*ldata-5)/(2*ldata-3)/(2*ldata-1)/(2*ldata+1)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)+param4/(2*ldata-5)/(2*ldata-3)/(2*ldata-1)/(2*ldata+1)/(2*ldata+3)/(2*ldata+5)/(2*ldata+7)/(2*ldata-7)/(2*ldata+9)
#    return psir

def fit_func1(ldata,param1):
    psir=param1/(2.*ldata-1.)/(2.*ldata+3.)
    return psir

def fit_func2(ldata,param1,param2):
    psir=param1/(2.*ldata-1.)/(2.*ldata+3.)+param2/(2.*ldata-1.)/(2.*ldata+3.)/(2.*ldata-3.)/(2.*ldata+5.)
    return psir

def fit_func3(ldata,param1,param2,param3):
    psir=param1/(2.*ldata-1.)/(2.*ldata+3.)+param2/(2.*ldata-3.)/(2.*ldata-1.)/(2.*ldata+3.)/(2.*ldata+5.)+param3/(2.*ldata-5.)/(2.*ldata-3.)/(2.*ldata+3.)/(2.*ldata+5.)/(2.*ldata+7.)
    return psir


def sum_func1(lmin):
    return lmin/(4.*lmin**2.-1.)

def sum_func2(lmin):
    return lmin/3./(9.-40.*lmin**2.+16.*lmin**4.)

def sum_func3(lmin):
    return lmin/5./(2.*lmin-5.)/(2.*lmin-3.)/(3.*lmin-1.)/(2.*lmin+1.)/(2.*lmin+3.)/(2.*lmin+5.)


def fit_lmode_and_sum(fitindex,llist,psir,maxmodefit,sigmaweights,startindex):

    
    fit_func=fit_func1
    if fitindex==1:
        fit_func=fit_func1
    elif fitindex==2:
        fit_func=fit_func2
    elif fitindex==3:
        fit_func=fit_func3
    paramopt, paramcov = optimize.curve_fit(fit_func, llist,psir,sigma=sigmaweights)
        
    partialsum=[]
    partialsum2=[]
    psirtosum=np.sum(psir[0:maxmodefit+1])
    unextrapolatedsum=np.sum(psirtosum)
    extrapolatedsum1=sum_func1(maxmodefit+1)*paramopt[0]
    partialsum.append(extrapolatedsum1)
    extrapolatedsum2 = 0
    extrapolatedsum3 = 0
    if fitindex>1:
        extrapolatedsum2 = sum_func2(maxmodefit+1)*paramopt[1]
        partialsum.append(extrapolatedsum2)
    if fitindex==3:
        extrapolatedsum3 = sum_func3(maxmodefit+1)*paramopt[2]
        partialsum.append(extrapolatedsum3)
            
    sumtot=unextrapolatedsum+extrapolatedsum1+extrapolatedsum2+extrapolatedsum3

       
    print sumtot,paramopt[0],paramopt[1],unextrapolatedsum,extrapolatedsum1,extrapolatedsum2,extrapolatedsum3           

    return sumtot

            

finfcolumn=4
lcolumn=1
nummodes=31

timecolumn=0
interporder=4
columnoffset=5
interpkind='cubic'
#startmode=14
maxmodefit=30
#initially, sum up to maxmodefit then use fit parameters from there (because there's no reason you'd run with more modes than you intend to fit)
#also should try summing up to start mode and using fit from there



t0 = 570
allindices=np.array(range(0,31,1))
startindeces=np.array(range(14,15,1))
finalindices=np.array(range(26,27,1))
fitindices=np.array(range(2,3,1))
sumtotalarr=np.zeros([len(fitindices)*len(startindeces)*len(finalindices)])
startx=np.zeros(len(sumtotalarr))
finaly=np.zeros(len(sumtotalarr))
orders=[12,16,20,24,28, 32,36,40,44,0] #not 48,33
i1=3
i2=i1+1
i3=i2+1
#order 400 is actually infinite
usesigma=False

lbestarr=np.zeros([len(allindices),len(orders)])

bestselfforcearr=np.zeros(len(orders))



count =0
for order in range(len(orders)):
    skprows=1
    if(orders[count]==0):
        loadstring = "coeffsbyl"+str(t0)+"_"+str(orders[i1])+"_"+str(orders[i2])+"_"+str(orders[i3])+".csv" #extrapolate from 24, 28, 32
        skprows=0
    elif((orders[count]==48) or (orders[count]==28) or (orders[count]==20)):
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

        
    datatable=np.loadtxt(loadstring,skiprows=skprows)
    print loadstring
    if(orders[count]!=0):
        for ii in range(len(datatable[:,timecolumn])):
            if datatable[ii,timecolumn]<t0:
                tnearest=datatable[ii,timecolumn]
                indexnearest=ii
        for modenum in allindices:
            tstored=np.zeros(interporder)
            lstored=np.zeros(interporder)
            for ii in range(interporder):
                tstored[ii]=datatable[indexnearest-(interporder-1)/2+ii,timecolumn]
                lstored[ii]=datatable[indexnearest-(interporder-1)/2+ii, columnoffset+modenum]
        
            func=interp1d(tstored,lstored,kind=interpkind)
            lbest=func(t0)
            lbestarr[modenum,count]=lbest
                

    fiti=0
    for fitindex in fitindices:
        starti=0
        for startindex in startindeces:
            modei=0
            for maxmodefit in finalindices:
                llist=np.zeros(maxmodefit+1-startindex)
                psir=np.zeros(len(llist))
                if (orders[count]==0):
                    llist=datatable[startindex:maxmodefit+1,lcolumn]
                    psir=datatable[startindex:maxmodefit+1,finfcolumn]
                else:
                    llist=np.array(range(startindex,maxmodefit+1,1))
                    psir=lbestarr[startindex:maxmodefit+1,count]
                    sigmaweights=np.ones(len(llist))
                    
                if(usesigma):
                    for ii in range(len(llist)):
                        sigmaweights[ii]=llist[ii]**-2.
                    
                    
                sumtotal=fit_lmode_and_sum(fitindex,llist,psir,maxmodefit,sigmaweights,startindex)
                print sumtotal
                sumtotalarr[fiti*len(finalindices)*len(startindeces)+starti*(len(finalindices))+modei]=sumtotal
                startx[fiti*len(finalindices)*len(startindeces)+starti*(len(finalindices))+modei]=startindex
                finaly[fiti*len(finalindices)*len(startindeces)+starti*(len(finalindices))+modei]=maxmodefit
                
                modei+=1
            
            starti+=1
        fiti+=1
        #bestselfforce=np.max(sumtotalarr)
        #bestselfforce2=np.max(sumtotalarr2)
        #bestindexselfforce=np.argmax(sumtotalarr)
        #bestindexselfforce2=np.argmax(sumtotalarr2)
        #beststartx=startx[bestindexselfforce]
        #beststartx2=startx[bestindexselfforce2]
        #beststarty=finaly[bestindexselfforce]
        #beststarty2=finaly[bestindexselfforce2]
        bestselfforcearr[count]=np.max(sumtotalarr)

        
        print orders[count], np.max(sumtotalarr)
#        print orders[count], bestselfforce2, beststartx2, beststarty2, bestindexselfforce2
        
    count+=1
        #plt.plot(llist,0*llist,'-')
            #plt.plot(llist,residual, label="start l="+str(startindex)+" end l="+str(maxmodefit))
            #ax=plt.gca()
            #plt.legend(loc='upper right')
            #plt.ylabel('Re(dpsi/dr)')
            #plt.xlabel('l mode')
            #plt.title('Radial self force fit residuals, t=' + str(t0))
            #', starting from l=' + str(startmode))

psiradjusted=np.zeros(len(orders)-1)
for ii in range(len(orders)-1):
    psiradjusted[ii]=abs(bestselfforcearr[len(orders)-1]-bestselfforcearr[ii])
print orders
print bestselfforcearr
print psiradjusted
    #plt.plot(orders[:len(orders)-1],bestselfforcearr[:len(orders)-1],'x-',label='Finite DG order')
#plt.plot(orders[:len(orders)-1],bestselfforcearr[len(orders)-1]*np.ones(len(orders)-1),'-',label='Infinite DG order')
plt.plot(orders[:len(orders)-1],psiradjusted,'o',label='Finite minus infinte order self force')
ax=plt.gca()
ax.set_yscale('log')
plt.legend(loc='upper right')
plt.ylabel('Summed radial self force')
plt.xlabel('DG order')
plt.title('lmin=14, lmax=26, Finf extrapolated from order '+str(orders[i1])+", "+str(orders[i2])+", "+str(orders[i3]))
plt.show()

            
#fig=plt.figure()
#ax=fig.gca(projection='3d')
#ax.scatter(startx, finaly,sumtotalarr, c='b', marker='o')
#ax.scatter(startx, finaly, sumtotalarr2, c='r',marker='^')
#ax=plt.gca()


#print startx
#ax.set_zlabel("")
#ax.set_xlim(min(startx),max(startx))
#ax.set_ylim(min(finaly),max(finaly))
#ax.set_zlim(min(sumtotalarr),max(sumtotalarr))
#plt.xlabel("Start l mode")
#plt.ylabel("Final l mode")
#plt.zlabel("Summed radial self force")
#plt.title("Variation of total radial self force with start and end ponits of fit")
#plt.show()

