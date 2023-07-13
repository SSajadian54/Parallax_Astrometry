import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
################################################################################

yr=float(364.0);
f0=open("./files/MONTB/IBH_MONT1b.dat","r")
nr=sum(1 for line in f0)
par=np.zeros(( nr , 53 ))
par=np.loadtxt("./files/MONTB/IBH_MONT1b.dat") 




for i in range(nr):   
    icon, lat, lon= int(par[i,0]), par[i,1], par[i,2]
    strucl, Ml, Dl, vl= par[i,3], par[i,4], par[i,5], par[i,6]
    strucs, cl, mass, Ds,Tstar, Rstar, logl= par[i,7], par[i,8], par[i,9], par[i,10], par[i,11], par[i,12], par[i,13]
    types, col, vs, MW149, Mk, mW149, mK=    par[i,14], par[i,15],par[i,16], par[i,17], par[i,18], par[i,19], par[i,20]
    mbs1, mbs2, fb1, fb2,Nb1, Nb2, Ex1, Ex2= par[i,21], par[i,22],par[i,23], par[i,24], par[i,25], par[i,26], par[i,27], par[i,28]
    tE, RE, t0, mul, Vt, u0, opd,ros,tetE=par[i,29],par[i,30],par[i,31], par[i,32], par[i,33], par[i,34], par[i,35],par[i,36], par[i,37]
    flagf,flagD, dchi, ndw, li, mus1, mus2=par[i,38], par[i,39],par[i,40], par[i,41], par[i,42], par[i,43], par[i,44]
    xi, mul1, mul2, piE=  par[i,45], par[i,46],par[i,47], par[i,48]
    ampM, errM, ampA, errA, dchip, dchia= par[i,49], par[i,50],par[i,51], par[i,52], par[i,53], par[i,54]
    
    pirel= float(1.0/Dl-1.0/Ds)## marcs
    print("Parameters: i,  icon, tE, tetE, piE :   ", i,  icon,  tE,    tetE,    piE)
    
    
    if(flagf>0 and icon>117): 
        f1=open("./files/MONTB/datB{0:d}.dat".format(icon),"r")
        nd=sum(1 for line in f1)  
        dat=np.zeros((nd,10)); 
        dat= np.loadtxt("./files/MONTB/datB{0:d}.dat".format(icon))
        
        
        f2=open("./files/MONTB/magB{0:d}.dat".format(icon),"r")
        nm=sum(1 for line in f2)  
        mag=np.zeros((nm,15));   
        mag= np.loadtxt("./files/MONTB/magB{0:d}.dat".format(icon))
        
        dat1=np.zeros((nd,10));   d1=0         
        dat2=np.zeros((nd,10));   d2=0
        for j in range(nd): 
            if(dat[j,8]<0.5):
                dat1[d1,:]= dat[j,:]#ELT extra observation
                d1+=1 
            if(dat[j,8]>0.5):
                dat2[d2,:]= dat[j,:]#Roman common observation
                d2+=1     
        delx= abs( np.max(mag[:,9]) -  np.min(mag[:,9]) )*0.1
        dely= abs( np.max(mag[:,10]) - np.min(mag[:,10]))*0.1         
        ######################################################################## 
        '''
        ymin=0.0;  ymax=0.0
        ymin= np.min(np.concatenate((dat1[:d1,3],dat2[:d2,3],mag[:,3]),axis=0))- np.max(np.concatenate((dat1[:d1,4],dat2[:d2,4]), axis=0))
        ymax= np.max(np.concatenate((dat1[:d1,3],dat2[:d2,3],mag[:,3]),axis=0))+ np.max(np.concatenate((dat1[:d1,4],dat2[:d2,4]), axis=0))
        dy=float(ymax-ymin)/50.0
        plt.clf()
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        plt.errorbar(dat1[:d1,0]/yr,dat1[:d1,3],yerr=dat1[:d1,4],fmt=".",markersize=7.5,color='darkred',ecolor='red',
         elinewidth=0.1,capsize=0, alpha=0.7)
        plt.errorbar(dat2[:d2,0]/yr,dat2[:d2,3],yerr=dat2[:d2,4],fmt=".",markersize=0.9,color='g',ecolor='#C1FFC1',
        elinewidth=0.01,capsize=0, alpha=0.3)
        plt.plot(mag[:,0]/yr, mag[:,3], "k-", lw=1.5, label=r"$\rm{Magnification}$")
        plt.plot(mag[:,0]/yr, mag[:,4], "b--",lw=1.5, label=r"$\rm{Magnification}+\rm{parallax}$")
        plt.xlabel(r"$\rm{time(yrs)}$", fontsize=18)
        plt.ylabel(r"$\rm{Magnification}$", fontsize=18)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        plt.xlim([0.0,5.0])
        plt.ylim([ymin-dy,ymax+dy])
        plt.title(r"$t_{\rm{E}}\rm{(days)}=$"+str(round(tE,1))+ r"$,~\theta_{\rm E}\rm{(mas)}=$"+str(round(tetE,2))+
          r"$,~\pi_{\rm E}=$"+str(round(piE,3))+r"$,~m_{\rm{bg},~W149}(\rm{mag})=$"+str(round(mbs1,1)),fontsize=16, color="k")
        plt.legend()
        ax1.legend(prop={"size":15.5})
        fig=plt.gcf()
        fig.savefig("./lightc/2/lcB{0:d}.jpg".format(icon),dpi=200)
        
        ########################################################################  
         
        
        plt.clf()
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        plt.errorbar(dat1[:d1,5],dat1[:d1,6],yerr=dat1[:d1,7],xerr=dat1[:d1,7],fmt=".",markersize=6.5,color='darkred',ecolor='red',     
        elinewidth=0.1,capsize=0,alpha=0.6)
        plt.errorbar(dat2[:d2,5],dat2[:d2,6],yerr=dat2[:d2,7],xerr=dat2[:d2,7],fmt=".",markersize=0.9,color='g',ecolor='#C1FFC1',
        elinewidth=0.01,capsize=0,alpha=0.3)
        plt.plot(mag[:,9],  mag[:,10], "-", color='k',lw=1.5, label=r"$\rm{Deflection}$", alpha=1.0)
        plt.plot(mag[:,13], mag[:,14],"--", color='b',lw=1.5, label=r"$\rm{Deflection}+\rm{Parallax}$", alpha=1.0)
    
        plt.xlim([np.min(mag[:,9])-delx,  np.max(mag[:,9]) +delx])
        plt.ylim([np.min(mag[:,10])-dely, np.max(mag[:,10])+dely])
        plt.xlabel(r"$\Delta x\rm{(mas)}$", fontsize=17)
        plt.ylabel(r"$\Delta y\rm{(mas)}$", fontsize=17)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        
        plt.title(r"$M_{\rm{l}}(M_{\odot})=$"+str(round(Ml,1))+ r"$,~D_{\rm{l}}\rm{(kpc)}=$"+str(round(Dl,1))+r"$,~\pi_{\rm{rel}}=$"+str(round(pirel,3))+r"$,~m_{\rm{bg},~K}(\rm{mag})=$"+str(round(mbs2,1)),fontsize=16, color="k")
        ax1.legend(prop={"size":15.5})
        fig=plt.gcf()
        fig.savefig("./lightc/2/dfB{0:d}.jpg".format(icon),dpi=200)
        '''
        ########################################################################
                
        plt.clf()
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        plt.plot(mag[:,9],  mag[:,10],"-",  color='k', lw=1.5, label=r"$\rm{Astrometric}~\rm{Deflection}$", alpha=1.0)
        plt.plot(mag[:,13], mag[:,14],"--", color='b',lw=1.5, label=r"$\rm{Deflection}+\rm{Parallax}$",  alpha=1.0)
        plt.plot(mag[:,11], mag[:,12],":",  color='m', lw=1.5, alpha=1.0)
        plt.xlim([np.min(mag[:,9])-delx,  np.max(mag[:,9]) +delx])
        plt.ylim([np.min(mag[:,10])-dely, np.max(mag[:,10])+dely])
        plt.xlabel(r"$\Delta x\rm{(mas)}$", fontsize=20)
        plt.ylabel(r"$\Delta y\rm{(mas)}$", fontsize=20)
        plt.xticks(fontsize=21, rotation=0)
        plt.yticks(fontsize=21, rotation=0)
        plt.title(r"$M_{\rm{l}}(M_{\odot})=$"+str(round(Ml,1))+ r"$,~D_{\rm{l}}\rm{(kpc)}=$"+str(round(Dl,1))+r"$,~\pi_{\rm{rel}}=$"+str(round(pirel,3))+r"$,~\theta_{\rm E}\rm{(mas)}=$"+str(round(tetE,2)),fontsize=19.5, color="k")
        ax1.legend(prop={"size":17.0})
        ax1.grid("True")
        ax1.grid(linestyle='dashed')
        fig=plt.gcf()
        plt.subplots_adjust(hspace=.0)
        fig.tight_layout(pad=0.1)
        fig.savefig("./lightc/2/dftB{0:d}.jpg".format(icon),dpi=200)
        print("********* Lightcurve & deflection was plotted ****, No:  ", icon)
        
        ########################################################################
    
      
