#Works fine as long as you don't loop over fitindices. Produces a wire grid plot (a surface plot) of the summed self force as a function of lmin and lmax

from mpl_toolkits.mplot3d import Axes3D
from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import csv
from matplotlib import ticker
import sys, getopt, cmath, os


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
    psir=param1/(2.*ldata-1.)/(2.*ldata+3.)+param2/(2.*ldata-3.)/(2.*ldata-1.)/(2.*ldata+3.)/(2.*ldata+5.)+param3/(2.*ldata-5.)/(2.*ldata-3.)/(2.*ldata-1.)/(2.*ldata+3.)/(2.*ldata+5.)/(2.*ldata+7.)
    return psir

def fit_func4(ldata,param1,param2,param3,param4):
    psir=param1/(2.*ldata-1.)/(2.*ldata+3.)+param2/(2.*ldata-3.)/(2.*ldata-1.)/(2.*ldata+3.)/(2.*ldata+5.)+param3/(2.*ldata-5.)/(2.*ldata-3.)/(2.*ldata-1.)/(2.*ldata+3.)/(2.*ldata+5.)/(2.*ldata+7.)+param4/(2.*ldata-7.)/(2.*ldata-5.)/(2.*ldata-3.)/(2.*ldata-1.)/(2.*ldata+3.)/(2.*ldata+5.)/(2.*ldata+7.)/(2.*ldata+9.)
    return psir
    
def sum_func1(lmin):
    return lmin/(4.*lmin**2.-1.)

def sum_func2(lmin):
    return lmin/3./(9.-40.*lmin**2.+16.*lmin**4.)

def sum_func3(lmin):
    return lmin/5./(2.*lmin-5.)/(2.*lmin-3.)/(2.*lmin-1.)/(2.*lmin+1.)/(2.*lmin+3.)/(2.*lmin+5.)

def sum_func4(lmin):
    return lmin/7./(2.*lmin-7.)/(2.*lmin-5.)/(2.*lmin-3.)/(2.*lmin-1.)/(2.*lmin+1.)/(2.*lmin+3.)/(2.*lmin+5.)/(2.*lmin+7.)

def main(argv):

    if (len(sys.argv)<7):
        print "Usage fitcoefftable.py t0 useBestFinf useReducedRange showSurfacePlot useAvg termsChosen startOrder(opt)"
        print "\tt0 is the initial time"
        print "\tuseBestFinf is 0 to use Finf based on a certain set of three \n\t\t extrapolation points starting with startOrder.\n\t\tFor useBestFinf=0 run extrapolate7.py first"
        print "\tuseBestFinf is 1 to use Finf selected based on the highest \n\t\t order for which Finf was defined.\n\t\tFor useBestFinf=1 run bestfinfselector.py first."
        print "\tuseReducedRange=0 chooses lmin=14-19 lmax=24-30."
        print "\tuseReducedrange=1 chooses lmin=14-17, lmax=22-25, which avoids roundoff error near periastron"
        print "\tshowSurfacePlot=1 shows the surface plot, showSurfacePlot=0 hides it"
        print "\tuseAvg=1 uses an average over the termsChosen surface to\n\t\testimate the total self force for that time."
        print "\tuseAvg=0 uses the central point of the terms chosen surface\n\\t\t to estimate the total self force for that time."
        exit()
    t0=int(sys.argv[1])
    useBestFinf=(int(sys.argv[2])==1)
    useReducedRange=(int(sys.argv[3])==1)
    showPlot=(int(sys.argv[4])==1)
    useAvg=(int(sys.argv[5])==1)
    termsChosen=int(sys.argv[6])
    if(len(sys.argv)==8):
        startOrder=int(sys.argv[7])
    minplotnum=2
    finfcolumn=4
    lcolumn=1
    if(useBestFinf):
        finfcolumn=2
        lcolumn=0
    else:
        finfcolumn=4
        lcolumn=1
    nummodes=31

    #startmode=14
    maxmodefit=30
    #initially, sum up to maxmodefit then use fit parameters from there (because there's no reason you'd run with more modes than you intend to fit)
    #also should try summing up to start mode and using fit from there

    #datatable =np.loadtxt("coeffsbyl570_28_32_36.csv")
    #datatable =np.loadtxt("coeffsbyl590_24_28_32.csv")
    #datatable =np.loadtxt("coeffsbyl610_24_28_32.csv")
    datatable=[]
    if(useBestFinf):
        datatable =np.loadtxt("bestinfoverinitorder"+str(t0)+".csv")
    else:
        datatable =np.loadtxt("coeffsbyl"+str(t0)+"_"+str(startOrder)+"_"+str(startOrder+4)+"_"+str(startOrder+8)+".csv")
    

    #t0 = datatable[0,0]
    #t0 = 630

    plotnosigma=True
    plotsigma=False
    startindeces=[]
    finalindices=[]
    if(useReducedRange):
        startindeces=np.array(range(14,18,1))
        finalindices=np.array(range(22,26,1))
    else:
        startindeces=np.array(range(14,20,1))
        finalindices=np.array(range(24,31,1))
    fitindices=[2,3,4]
    sumtotalarr=np.zeros([len(fitindices),len(startindeces)*len(finalindices)])
    sumtotalarr2=np.zeros([len(fitindices),len(startindeces)*len(finalindices)])
    unextrapolatedarr=np.zeros([len(fitindices),len(startindeces)*len(finalindices)])
    startx=np.zeros(len(sumtotalarr[0,:]))
    finaly=np.zeros(len(sumtotalarr[0,:]))



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
                elif fitindex==4:
                    fit_func=fit_func4
                #print fitindex, startindex, maxmodefit,"nosigma"
                paramopt, paramcov = optimize.curve_fit(fit_func, llist,psir)
                #print fitindex, startindex, maxmodefit, "with sigma"
                paramopt2, paramcov2 = optimize.curve_fit(fit_func, llist,psir,sigma=sigmaweights)
                partialsum=[]
                partialsum2=[]
                psirtosum=np.sum(datatable[0:maxmodefit+1,finfcolumn])
                unextrapolatedsum=np.sum(psirtosum)
                extrapolatedsum1=sum_func1(maxmodefit+1)*paramopt[0]
                partialsum.append(extrapolatedsum1)
                extrapolatedsum2 = 0
                extrapolatedsum3 = 0
                extrapolatedsum4 = 0
                if fitindex>1:
                    extrapolatedsum2 = sum_func2(maxmodefit+1)*paramopt[1]
                    partialsum.append(extrapolatedsum2)
                if fitindex>2:
                    extrapolatedsum3 = sum_func3(maxmodefit+1)*paramopt[2]
                    partialsum.append(extrapolatedsum3)
                if fitindex>3:
                    extrapolatedsum4 = sum_func4(maxmodefit+1)*paramopt[3]
                    partialsum.append(extrapolatedsum4)
                extrapolatedsum1_2=sum_func1(maxmodefit+1)*paramopt2[0]
                partialsum2.append(extrapolatedsum1)
                extrapolatedsum2_2 = 0
                extrapolatedsum3_2 = 0
                extrapolatedsum4_2 = 0
                if fitindex>1:
                    extrapolatedsum2_2 = sum_func2(maxmodefit+1)*paramopt2[1]
                    partialsum2.append(extrapolatedsum2)
                if fitindex>2:
                    extrapolatedsum3_2 = sum_func3(maxmodefit+1)*paramopt2[2]
                    partialsum2.append(extrapolatedsum3)
                if fitindex>3:
                    extrapolatedsum4_2 = sum_func4(maxmodefit+1)*paramopt2[3]
                    partialsum2.append(extrapolatedsum3)
                sumtotal=unextrapolatedsum+extrapolatedsum1+extrapolatedsum2+extrapolatedsum3+extrapolatedsum4
                sumtotal2=unextrapolatedsum+extrapolatedsum1_2+extrapolatedsum2_2+extrapolatedsum3_2+extrapolatedsum4_2
                unextrapolatedarr[fiti,starti*(len(finalindices))+modei]=unextrapolatedsum
                sumtotalarr[fiti,starti*(len(finalindices))+modei]=sumtotal
                sumtotalarr2[fiti,starti*(len(finalindices))+modei]=sumtotal2
                if(fiti==0):
                    startx[starti*(len(finalindices))+modei]=startindex
                    finaly[starti*(len(finalindices))+modei]=maxmodefit
            
                temp1 = 0
                temp2=0
                temp3=0
                if fitindex>1:
                    temp1 = paramopt[1]
                if fitindex>2:
                    temp2 = paramopt[2]
                if fitindex>3:
                    temp3 = paramopt[3]
                #with open('parametertable.dat', 'a') as file:
                #   file.write(str(fitindex)+","+str(startindex)+","+str(maxmodefit)+","+str(sumtotal)+","+str(paramopt[0])+","+str(temp1)+","+str(temp2)+","+str(temp3)+","+str(unextrapolatedsum)+","+str(extrapolatedsum1)+","+str(extrapolatedsum2)+","+str(extrapolatedsum3))
                #    file.write("\n")
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
                    if fitindex==4:     
                        residual[ii]=psir[ii]-fit_func4(llist[ii],paramopt[0],paramopt[1],paramopt[2],paramopt[3])

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
    ax=plt.axes(projection='3d')
    #ax.scatter(startx, finaly,sumtotalarr, c='b', marker='o')
    #ax.scatter(startx, finaly, sumtotalarr2, c='r',marker='^')
    startx2 = np.reshape(startx,(len(startindeces),len(finalindices)))
    finaly2 = np.reshape(finaly,(len(startindeces),len(finalindices)))
    termChosenIndex=fitindices.index(termsChosen)
    if 1 in fitindices:
        oneindex=fitindices.index(1)
        zvals1 = np.reshape(sumtotalarr[oneindex,:],(len(startindeces),len(finalindices)))
        zvals1_2 = np.reshape(sumtotalarr2[oneindex,:],(len(startindeces),len(finalindices)))
    if 2 in fitindices:
        twoindex=fitindices.index(2)
        #print twoindex
        zvals2 =np.reshape(sumtotalarr[twoindex,:],(len(startindeces),len(finalindices)))
        zvals2_2 =np.reshape(sumtotalarr2[twoindex,:],(len(startindeces),len(finalindices)))
    if 3 in fitindices:
        threeindex=fitindices.index(3)
        zvals3 = np.reshape(sumtotalarr[threeindex,:],(len(startindeces),len(finalindices)))
        zvals3_2 = np.reshape(sumtotalarr2[threeindex,:],(len(startindeces),len(finalindices)))

    if 4 in fitindices:
        fourindex=fitindices.index(4)
        zvals4 = np.reshape(sumtotalarr[fourindex,:],(len(startindeces),len(finalindices)))
        zvals4_2 = np.reshape(sumtotalarr2[fourindex,:],(len(startindeces),len(finalindices)))

    selfForce=0
    selfForce=0
    unextrapSelfForce=0
    if(useAvg):
        selfForce=np.average(sumtotalarr[termChosenIndex,:])
        selfForce2=np.average(sumtotalarr2[termChosenIndex,:])
        unextrapSelfForce=np.average(unextrapolatedarr[termChosenIndex,:])
    else:
        zvalschosen=np.reshape(sumtotalarr[termChosenIndex,:],(len(startindeces),len(finalindices)))
        zvalschosen2=np.reshape(sumtotalarr2[termChosenIndex,:],(len(startindeces),len(finalindices)))
        zvalschosenunextrap=np.reshape(unextrapolatedarr[termChosenIndex,:],(len(startindeces),len(finalindices)))
        midx=len(startindeces)/2-1
        midy=len(finalindices)/2-1
        selfForce=zvalschosen[midx,midy]
        selfForce2=zvalschosen2[midx,midy]
        unextrapSelfForce=zvalschosenunextrap[midx,midy]

    if 1 in fitindices:
        if plotnosigma:
            ax.plot_wireframe(startx2,finaly2,zvals1,rstride=1,cstride=1,color="red",label='1 term')
        if plotsigma:
            ax.plot_wireframe(startx2,finaly2,zvals1_2,rstride=1,cstride=1,color="black",label='weights')
    if 2 in fitindices:
        if plotnosigma:
            ax.plot_wireframe(startx2,finaly2,zvals2,rstride=1,cstride=1,color="blue",label='2 term')
        if plotsigma:
            ax.plot_wireframe(startx2,finaly2,zvals2_2,rstride=1,cstride=1,color="orange",label='weights')
    if 3 in fitindices:
        if plotnosigma:
            ax.plot_wireframe(startx2,finaly2,zvals3,rstride=1,cstride=1,color="green",label='3 term')
        if plotsigma:
            ax.plot_wireframe(startx2,finaly2,zvals3_2,rstride=1,cstride=1,color="purple",label='weights')
            ax.legend(loc='lower left')
    if 4 in fitindices:
        if plotnosigma:
            ax.plot_wireframe(startx2,finaly2,zvals4,rstride=1,cstride=1,color="orange",label='4 term')
        if plotsigma:
            ax.plot_wireframe(startx2,finaly2,zvals4_2,rstride=1,cstride=1,color="green",label='weights')
            ax.legend(loc='lower left')
    #ax=plt.gca()

    #print sumtotalarr
    ax.set_zlabel("")
    ax.set_xlim(min(startindeces),max(startindeces))
    ax.set_ylim(min(finalindices),max(finalindices))
    #ax.set_zlim(min(min(zvals1),min(zvals2),min(zvals3)),max(max(zvals1),max(zvals2),max(zvals3)))

    formatter=ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    ax.w_zaxis.set_major_formatter(formatter)
    ax.ticklabel_format(axis="z",style="sci",scilimits=(0,0))
    plt.xlabel("Start l mode")
    plt.ylabel("Final l mode")
    #plt.zlabel("Total radial self force")
    plt.title("Total radial self force, using DG error extrapolation per l-mode, t="+str(t0))
    if(showPlot):
        plt.show()
    print t0, selfForce, selfForce2, unextrapSelfForce
    return t0, selfForce, selfForce2, unextrapSelfForce
    
if __name__=="__main__":
    main(sys.argv[1:])
