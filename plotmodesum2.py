
from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D


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


finfcolumn=4
lcolumn=1
nummodes=31

#startmode=14
maxmodefit=30
#initially, sum up to maxmodefit then use fit parameters from there (because there's no reason you'd run with more modes than you intend to fit)
#also should try summing up to start mode and using fit from there

datatable =np.loadtxt("coeffsbyl570.csv")

t0 = datatable[0,0]

startindeces=np.array(range(14,22,1))
finalindices=np.array(range(26,31,1))
fitindices=np.array(range(2,4,1))
sumtotalarr=np.zeros([len(fitindices)*len(startindeces)*len(finalindices)])
sumtotalarr2=np.zeros([len(fitindices)*len(startindeces)*len(finalindices)])
startx=np.zeros(len(sumtotalarr))
finaly=np.zeros(len(sumtotalarr))



fiti=0
for fitindex in fitindices:
    starti=0
    for startindex in startindeces:
        modei=0
        for maxmodefit in finalindices:
            llist=datatable[startindex:maxmodefit,lcolumn]
            psir=datatable[startindex:maxmodefit,finfcolumn]

            sigmaweights=np.zeros(len(llist))
            for ii in range(len(llist)):
                sigmaweights[ii]=llist[ii]**-2.
    
            fit_func=fit_func1
            if fitindex==1:
                fit_func=fit_func1
            elif fitindex==2:
                fit_func=fit_func2
            elif fitindex==3:
                fit_func=fit_func3
            paramopt, paramcov = optimize.curve_fit(fit_func, llist,psir, sigma=sigmaweights)
            paramopt2, paramcov2 = optimize.curve_fit(fit_func, llist,psir)
    
            partialsum=[]
            partialsum2=[]
            psirtosum=np.sum(datatable[0:maxmodefit,finfcolumn])
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

            extrapolatedsum1_2=sum_func1(maxmodefit+1)*paramopt2[0]
            partialsum.append(extrapolatedsum1)
            extrapolatedsum2_2 = 0
            extrapolatedsum3_2 = 0
            if fitindex>1:
                extrapolatedsum2_2 = sum_func2(maxmodefit+1)*paramopt2[1]
                partialsum.append(extrapolatedsum2)
            if fitindex==3:
                extrapolatedsum3_2 = sum_func3(maxmodefit+1)*paramopt2[2]
                partialsum.append(extrapolatedsum3)
            
            sumtotal=unextrapolatedsum+extrapolatedsum1+extrapolatedsum2+extrapolatedsum3
            sumtotal2=unextrapolatedsum+extrapolatedsum1+extrapolatedsum2+extrapolatedsum3
            sumtotalarr[fiti*len(finalindices)*len(startindeces)+starti*(len(finalindices))+modei]=sumtotal
            sumtotalarr2[fiti*len(finalindices)*len(startindeces)+starti*(len(finalindices))+modei]=sumtotal
            startx[fiti*len(finalindices)*len(startindeces)+starti*(len(finalindices))+modei]=startindex
            finaly[fiti*len(finalindices)*len(startindeces)+starti*(len(finalindices))+modei]=maxmodefit
            
            print startindex, maxmodefit, sumtotal
            #print "params=", paramopt
            #print "partialsum=",partialsum


            residual=np.zeros(len(llist))
            for ii in range(len(llist)):
                if fitindex==1:
                    residual[ii]=psir[ii]-fit_func1(llist[ii],paramopt[0])
                if fitindex==2:
                    residual[ii]=psir[ii]-fit_func2(llist[ii],paramopt[0],paramopt[1])
                if fitindex==3:
                    residual[ii]=psir[ii]-fit_func3(llist[ii],paramopt[0],paramopt[1],paramopt[2])

            #plt.plot(llist,0*llist,'-')
            #plt.plot(llist,residual, label="start l="+str(startindex)+" end l="+str(maxmodefit))
            #ax=plt.gca()
            #plt.legend(loc='upper right')
            #plt.ylabel('Re(dpsi/dr)')
            #plt.xlabel('l mode')
            #plt.title('Radial self force fit residuals, t=' + str(t0))
            #', starting from l=' + str(startmode))
            modei+=1
        
        starti+=1
    fiti+=1
fig=plt.figure()
ax=fig.gca(projection='3d')
ax.scatter(startx, finaly,sumtotalarr, c='b')
ax.scatter(startx, finaly, sumtotalarr2, c='r')
ax=plt.gca()

print startx
ax.set_zlabel("")
ax.set_xlim(min(startx),max(startx))
ax.set_ylim(min(finaly),max(finaly))
ax.set_zlim(min(sumtotalarr),max(sumtotalarr))
plt.xlabel("Start l mode")
plt.ylabel("Final l mode")
#plt.zlabel("Summed radial self force")
plt.title("Variation of total radial self force with start and end ponits of fit")
plt.show()
