import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.figure import Figure
from matplotlib import colors
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
################################################################################

cm=colors.ListedColormap(['purple', 'blue', 'darkgreen','yellowgreen', 'orange', 'red'])
yr=float(364.0);
f0=open("./files/MONTC/IBH_MONT1c.dat","r")
nr=sum(1 for line in f0)
par=np.zeros(( nr , 49 ))
par=np.loadtxt("./files/MONTC/IBH_MONT1c.dat") 

################################################################################

for i in range(nr):
    icon, lat, lon, right, decli=     int(par[i,0]), par[i,1], par[i,2], par[i,3], par[i,4]
    strucl, Ml, Dl, vl   =               par[i,5], par[i,6], par[i,7], par[i,8]
    strucs, cl, Ds,logT, logL, vs=       par[i,9], par[i,10], par[i,11], par[i,12], par[i,13], par[i,14]
    mab1, mab2, map1, map2,mbs1=        par[i,15],par[i,16], par[i,17], par[i,18], par[i,19]
    mbs2,fb1, fb1, Nb1, Nb2, Ex1, Ex2=  par[i,20], par[i,21], par[i,22],par[i,23], par[i,24], par[i,25], par[i,26]
    tE, RE, t0, mul, Vt, u0, opd,tetE=   par[i,27], par[i,28], par[i,29],par[i,30],par[i,31], par[i,32], par[i,33], par[i,34]
    flagf,flagdet, dchi, ndw, mus1, mus2=par[i,35],par[i,36], par[i,37], int(par[i,38]), par[i,39],par[i,40]
    xi, mul1, mul2, piE=                 par[i,41], par[i,42], par[i,43],  par[i,44]
    ampM, errM, ampA, errA= par[i,45], par[i,46],par[i,47], par[i,48]
    
    pirel= float(1.0/Dl-1.0/Ds)
    print("Parameters: i,  icon, tE, tetE, piE :   ", i,  icon,  tE,    tetE,    piE)
    
    
    if(flagf>0 ): 
        f1=open("./files/MONTC/datC{0:d}.dat".format(icon),"r")
        nd=sum(1 for line in f1)  
        dat=np.zeros((nd,11)); 
        dat= np.loadtxt("./files/MONTC/datC{0:d}.dat".format(icon))
        
        
        f2=open("./files/MONTC/magC{0:d}.dat".format(icon),"r")
        nm=sum(1 for line in f2)  
        mag=np.zeros((nm,27))  
        mag= np.loadtxt("./files/MONTC/magC{0:d}.dat".format(icon))
        
   
        dat1=np.zeros((nd,11));   d1=0         
        dat2=np.zeros((nd,11));   d2=0
        for j in range(nd): 
            if(dat[j,8]<0.5):
                dat1[d1,:]= dat[j,:]## ELT extra observation
                d1+=1 
            if(dat[j,8]>0.5):
                dat2[d2,:]=dat[j,:]## LSST common observation
                d2+=1    
                if(dat[j,10]<0 or dat[j,10]>5): 
                    print("Error   , ", j, icon,   dat[j, 10])
                    input("Enter a number ")
                    
        delx= abs( np.max(mag[:,9]) -  np.min(mag[:,9]) )*0.1
        dely= abs( np.max(mag[:,10]) - np.min(mag[:,10]))*0.1             
        ######################################################################## 
        
        ymin=0.0;  ymax=0.0; 
        ymin= np.min(np.concatenate((dat1[:d1,3],dat2[:d2,3],mag[:,3]),axis=0))- np.max(np.concatenate((dat1[:d1,4],dat2[:d2,4]), axis=0))
        ymax= np.max(np.concatenate((dat1[:d1,3],dat2[:d2,3],mag[:,3]),axis=0))+ np.max(np.concatenate((dat1[:d1,4],dat2[:d2,4]), axis=0))
        dy=float(ymax-ymin)/100.0
        plt.clf()
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        for j in range(d2): 
            nc=abs(int(dat2[j,10]))
            plt.errorbar(dat2[j,0]/yr,dat2[j,3],yerr=dat2[j,4], fmt=".", markersize=6.2, color=cm(nc),ecolor=cm(nc),
            elinewidth=0.6, capsize=0,alpha=0.3)
        plt.errorbar(dat1[:d1,0]/yr,dat1[:d1,3],yerr=dat1[:d1,4],fmt=".",markersize=6.2,color='darkred',ecolor='r',elinewidth=0.6,
            capsize=0,alpha=0.7)
        plt.plot(mag[:,0]/yr, mag[:,3], "k-", lw=1.5, label=r"$\rm{Magnification}$")
        plt.plot(mag[:,0]/yr, mag[:,4], "b--",lw=1.5, label=r"$\rm{Magnification}+\rm{parallax}$")
        plt.xlabel(r"$\rm{time(yrs)}$", fontsize=18)
        plt.ylabel(r"$\rm{Magnification}$", fontsize=18)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        plt.xlim([0.0,10.0])
        plt.ylim([ymin-dy,ymax+dy])
        plt.title(r"$t_{\rm{E}}\rm{(days)}=$"+str(round(tE,1))+ r"$,~\theta_{\rm E}\rm{(mas)}=$"+str(round(tetE,2))+
          r"$,~\pi_{\rm E}=$"+str(round(piE,3))+r"$,~m_{\rm{bg},~r}(\rm{mag})=$"+str(round(mbs1,1)),fontsize=16, color="k")
        plt.legend()
        ax1.legend(prop={"size":15.5})
        fig=plt.gcf()
        fig.savefig("./lightc/3/lcC{0:d}.jpg".format(icon),dpi=200)
        
        ######################################################################## 
        
       
        plt.clf()
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        
        for j in range(d2): 
            nc=abs(int(dat2[j,10]))
            plt.errorbar(dat2[j,5],dat2[j,6],yerr=dat2[j,7],xerr=dat2[j,7], fmt=".", markersize=6.2, color=cm(nc),ecolor=cm(nc),
            elinewidth=0.6, capsize=0,alpha=0.3)
        plt.errorbar(dat1[:d1,5],  dat1[:d1,6],yerr=dat1[:d1,7],xerr=dat1[:d1,7], fmt=".", markersize=6.2, color='darkred', ecolor='r',
        elinewidth=0.6, capsize=0,alpha=0.7)
       
        #plt.errorbar(dat1[:d1,5],dat1[:d1,6],yerr=dat1[:d1,7],xerr=dat1[:d1,7],fmt=".",markersize=6.5,color='b',ecolor='#97FFFF',     
        #elinewidth=0.1,capsize=0,alpha=0.6)
        #plt.errorbar(dat2[:d2,5],dat2[:d2,6],yerr=dat2[:d2,7],xerr=dat2[:d2,7],fmt=".",markersize=1.5,color='g',ecolor='#C1FFC1',
        #elinewidth=0.1,capsize=0,alpha=0.4)
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
        fig.savefig("./lightc/3/dfC{0:d}.jpg".format(icon),dpi=200)
        
        ########################################################################
             
        plt.clf()
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        plt.plot(mag[:,9],  mag[:,10],"-", color='k', lw=1.5, label=r"$\rm{Astrometric}~\rm{Deflection}$", alpha=1.0)
        plt.plot(mag[:,13], mag[:,14],"--",color='b', lw=1.5, label=r"$\rm{Deflection}+\rm{Parallax}$",  alpha=1.0)
        plt.plot(mag[:,11], mag[:,12],":", color='m', lw=1.5, alpha=1.0)
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
        fig.savefig("./lightc/3/dftC{0:d}.jpg".format(icon),dpi=200)
        print("********* Lightcurve & deflection was plotted ****, No:  ", icon)
        
        ########################################################################
        '''
        plt.clf()
        fig=plt.figure(figsize=(8,6))
        ax1=fig.add_subplot(111)
        plt.errorbar(dat1[:d1,0]/yr,dat1[:d1,1],yerr=dat1[:d1,2], fmt=".", markersize=6.2,color='darkred',ecolor='red',elinewidth=1.1, capsize=0, alpha=0.7)
        for j in range(d2): 
            nc=abs(int(dat2[j,10]) )
            plt.errorbar(dat2[j,0]/yr,dat2[j,1],yerr=dat2[j,2], fmt=".", markersize=6.2, color=cm(nc),ecolor=cm(nc),elinewidth=1.1, capsize=0,alpha=0.7)
        plt.plot(mag[:,0]/yr, mag[:,1], "m:", lw=0.9, label=r"$u-\rm{LSST}$")
        plt.plot(mag[:,0]/yr, mag[:,2], "m--",lw=0.9)
        plt.plot(mag[:,0]/yr, mag[:,15],':',  color="blue",    lw=0.9)
        plt.plot(mag[:,0]/yr, mag[:,16],'--',  color="blue",   lw=0.9, label=r"$g-\rm{LSST}$")
        plt.plot(mag[:,0]/yr, mag[:,17],':',  color="green",   lw=0.9)
        plt.plot(mag[:,0]/yr, mag[:,18],'--', color="green",   lw=0.9, label=r"$r-\rm{LSST}$")
        plt.plot(mag[:,0]/yr, mag[:,19],':',  color="lime",    lw=0.9)
        plt.plot(mag[:,0]/yr, mag[:,20],'--', color="lime",    lw=0.9, label=r"$i-\rm{LSST}$")
        plt.plot(mag[:,0]/yr, mag[:,21],':',  color="orange",  lw=0.9)
        plt.plot(mag[:,0]/yr, mag[:,22],'--', color="orange",  lw=0.9, label=r"$z-\rm{LSST}$")
        plt.plot(mag[:,0]/yr, mag[:,23],':',  color="red",     lw=0.9)
        plt.plot(mag[:,0]/yr, mag[:,24],'--', color="red",     lw=0.9, label=r"$y-\rm{LSST}$")
        plt.plot(mag[:,0]/yr, mag[:,25],':',  color="darkred", lw=0.9)
        plt.plot(mag[:,0]/yr, mag[:,26],'--', color="darkred", lw=0.9, label=r"$K-\rm{ELT}$")
        plt.xlabel(r"$\rm{time(yrs)}$", fontsize=18)
        plt.ylabel(r"$\rm{Magnitude}(\rm{mag})$", fontsize=18)
        plt.xticks(fontsize=17, rotation=0)
        plt.yticks(fontsize=17, rotation=0)
        plt.xlim([0.0,10.0])
        plt.title(r"$t_{\rm E}\rm{(days)}=$"+str(round(tE,1))+ r"$,~~\theta_{\rm E}\rm{(mas)}=$"+str(round(tetE,2))+r"$,~~\pi_{\rm E}=$"+str(round(piE,3))+r"$,~m_{\rm{base},~r}=$"+str(round(mbs1,1))+r"$,~m_{\rm{base},~K}=$"+str(round(mbs2,1)),fontsize=14, color="k")
        plt.gca().invert_yaxis()
        plt.legend()
        ax1.legend(prop={"size":15.5})
        fig=plt.gcf()
        fig.savefig("./lightc/3/lctC{0:d}.jpg".format(icon),dpi=200)
        '''
        ########################################################################
        
        
        
        
