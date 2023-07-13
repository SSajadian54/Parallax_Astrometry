import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.figure import Figure
import warnings
warnings.filterwarnings("ignore")
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
###############################################################
yr=float(364.0)
pp=float(100.0)
###############################################################
f1=open("./files/MONTC/IBH_MONT1c.dat","r")
f2=open("./files/MONTC/IBH_MONT2c.dat","r")
f3=open("./files/MONTC/IBH_MONT3c.dat","r")
nm= sum(1 for line in f1) 
par=np.zeros((nm,51)) 
dat=np.zeros((nm,38)) 
mat=np.zeros((nm,14)) 
par= np.loadtxt("./files/MONTC/IBH_MONT1c.dat")
dat= np.loadtxt("./files/MONTC/IBH_MONT2c.dat")
mat= np.loadtxt("./files/MONTC/IBH_MONT3c.dat")

fij=open("./files/HistoMCC/resultMCC1.dat","w")
fij.close(); 
fij=open("./files/HistoMCC/resultMCC2.dat","w")
fij.close(); 

##############################################################

nd=7;  na=4;  nq=17; nv=11
F = np.zeros([nd , nd]);  
G = np.zeros([na , na])
Era=np.zeros((nd)); 
Erb=np.zeros((na));
resu= np.zeros((nq))
numc= np.zeros(nq+1)
numa= np.zeros(nq+1)
numb= np.zeros(nq+1)
hisd0= np.zeros((nm,10))
hisd2= np.zeros((nm,10))
hisd3= np.zeros((nm,10))
hisd4= np.zeros((nm,10))
epsi= np.zeros(3)
arry=np.zeros((nm, nq+4))
parp=np.zeros((nm,9))
chis=np.zeros((nm,5))
l=0;  nc=0; 
xar=np.zeros((nv,4))
xar[:,0]=np.array([2.0, 4.8 , 8.0, 14.0, 20.0, 28.0, 33.0, 38.0, 42.5, 46.5, 50.0])
xar[:,1]=np.array([ 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
xar[:,2]=np.array([0.1,2.0,4.5,6.0,7.0, 7.7, 8.5, 9.2, 10.0, 11.0, 12.5])
xar[:,3]=np.array([15.5,16.4,17.3,18.2,19.1,20.0,20.9,21.8, 22.7, 23.6, 24.5])
##############################################################
def func(para,cc): 
    i1=-1
    for j in range(nv-1):  
        if( float((para-xar[j,cc])*(para-xar[j+1,cc]))<0.0  or  para==xar[j,cc]):
            i1=j;  
            break
    return(i1)        
##############################################################
nam0=["counter","b", "l","Right", "Declination", #5
   "l.struc", "l.mass", "l.Dl", "l.vl", #9
    "s.struc", "s.cl", "s.Ds", "s.logT", "s.logL","vs", r"$Mr_{LSST}$", "Mk_{ELT}", r"$mr_{LSST}$", r"$mk_{ELT}$", #19
    "s.magb_{r-LSST}", "s.magb_{k-ELT}",r"$fb_{r-LSST}$", r"$fb_{k-ELT}$","Nb_{r-LSST}",r"$Nb_{k-ELT}$", "Ext[r-LSST]", "Ext[k-ELT]", #27
    "l.tE", "l.RE(AU)", "l.t0", "l.mul", "l.Vt", "l.u0", "s.opt(e6)", "tetE" , #35
    "flagf","flagdet", "dchi", "ndw", "mus1", "mus2", "xi","mul1", "mul2", "piE", "ampM", "errM", "ampA", "errA", "dchip_ph_N",
    "dchi_as_N"] #51
   
for i in range(51):
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    ax= plt.gca()              
    plt.hist(par[:,i],30,histtype='bar',ec='darkgreen',facecolor='green',alpha=0.7, rwidth=1.5)
    y_vals = ax.get_yticks()
    ax.set_yticks(y_vals)
    ax.set_yticklabels(['{:.2f}'.format(float(1.0*x*(1.0/nm))) for x in y_vals]) 
    y_vals = ax.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
    ax.set_xlabel(str(nam0[i]),fontsize=19,labelpad=0.1)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig=plt.gcf()
    fig.savefig("./files/HistoMCC/hsimC{0:d}.jpg".format(i),dpi=200)
print ("****  All histos are plotted *****************************" )   

#############################################################
nph=int(0.0);  
nas=int(0.0);
ngo=int(0.0)
jf=int(0)
eps= float(0.00000000005463263454624313654)
mapp=np.zeros((nv,nv,2)); 
halo=0.0; 
for i in range(nm): 
    coun1, lat, lon, right, decli=   int(par[i,0]), par[i,1], par[i,2], par[i,3], par[i,4]
    strucl, Ml, Dl, vl   =               par[i,5], par[i,6], par[i,7], par[i,8]
    strucs, cl, Ds,logT, logL, vs=       par[i,9], par[i,10], par[i,11], par[i,12], par[i,13], par[i,14]
    mab1, mab2, map1, map2,magb1=        par[i,15],par[i,16], par[i,17], par[i,18], par[i,19]
    magb2,fb1, fb1, Nb1, Nb2, Ex1, Ex2=  par[i,20], par[i,21], par[i,22],par[i,23], par[i,24], par[i,25], par[i,26]
    tE, RE, t0, mul, Vt, u0, opd,tetE=   par[i,27], par[i,28], par[i,29],par[i,30],par[i,31], par[i,32], par[i,33], par[i,34]
    flagf,flagdet, dchi, ndw, mus1, mus2=par[i,35],par[i,36], par[i,37], int(par[i,38]), par[i,39],par[i,40]
    xi, mul1, mul2, piE=                 par[i,41], par[i,42], par[i,43],  par[i,44]
    ampM, errM, ampA, errA, dchip, dchia=par[i,45], par[i,46],par[i,47], par[i,48], par[i,49], par[i,50]
    if(strucl<4): halo+=1.0   
       
    lpir=np.log10(1.0/Dl-1.0/Ds)
    mus=np.sqrt(mus1*mus1+mus2*mus2)
    mul=np.sqrt(mul1*mul1+mul2*mul2)
    xls=float(Dl/Ds)
    dd= float(ampA)##/float(ampM) )
    
    
    #if(not(dd>0.0 and dchia<1.0)):
    parp[jf,:]= dd, Ml, Dl, lpir, np.log10(tetE) , np.log10(piE), strucl, np.log10(tE), u0
    jf+=1 
            
    coun2,t0b,u0b,tEb,xib,fb1,mbs1,fb2,mbs2,piEb=dat[i,0],dat[i,1],dat[i,2],dat[i,3],dat[i,4],dat[i,5],dat[i,6],dat[i,7],dat[i,8],dat[i,9]
    F[0,0]=dat[i,10]
    F[1,0], F[1,1]=dat[i,11], dat[i,12];
    F[2,0], F[2,1], F[2,2]= dat[i,13], dat[i,14], dat[i,15];
    F[3,0], F[3,1], F[3,2], F[3,3]=dat[i,16], dat[i,17], dat[i,18], dat[i,19];
    F[4,0], F[4,1], F[4,2], F[4,3], F[4,4]= dat[i,20], dat[i,21], dat[i,22], dat[i,23], dat[i,24];
    F[5,0], F[5,1], F[5,2], F[5,3], F[5,4], F[5,5]=dat[i,25], dat[i,26], dat[i,27], dat[i,28], dat[i,29], dat[i,30];
    F[6,0], F[6,1], F[6,2], F[6,3], F[6,4], F[6,5], F[6,6]= dat[i,31], dat[i,32], dat[i,33], dat[i,34], dat[i,35], dat[i,36], dat[i,37]   
    
    
    coun3, tetEb, mus1b, mus2b = mat[i,0], mat[i,1], mat[i,2], mat[i,3]
    G[0][0]= mat[i,4]
    G[1][0], G[1][1]= mat[i,5], mat[i,6]
    G[2][0], G[2][1], G[2][2]= mat[i,7], mat[i,8], mat[i,9]
    G[3][0], G[3][1], G[3][2], G[3,3]= mat[i,10], mat[i,11], mat[i,12], mat[i,13]
    
    
    if(coun1!=coun2 or coun2!=coun3 or coun3!=coun1 or abs(t0b-t0)>0.001 or abs(u0b-u0)>0.001 or abs(tEb-tE)>0.001 or abs(xib-xi)>0.001 
    or abs(fb1-fb1)>0.4 or abs(magb1-mbs1)>1.0 or t0<=0.0 or  abs(piEb-piE)>0.001 or abs(tetEb-tetE)>0.001 or abs(mus1-mus1b)>0.001 or 
    abs(mus2-mus2b)>0.001 or Dl>Ds or u0>1.0or u0<0.0 or t0>(10.0*yr) or tE<0.0 or Ml<2.0  or Ml>600.0 or flagdet<1 
    or dchi<800.0  or tetE<0.0  or vl<0.0  or Vt<0.0 or ndw<1.0 or  ampM<0.0 or ampA<0.0  or errM<0.0  or errA<0.0 or dchia<0 or dchip<0): 
        print("Error:  ",  coun1, coun2, coun3, t0, t0b,  u0, u0b,  tE, tEb,  xi, xib, fb1, fb2,  magb1, mbs1, mbs2)
        print("Error ",  piE, piEb, tetE, tetEb, mus1, mus1b, mus2, mus2b,  Dl, Ds, Ml )
        print("Error ",  flagdet, dchi, vl, vs, Vt, Ml, ndw )
        print("Amplitudes, errors:  ", ampM,  ampA,  errM,  errA)
        print("Delta chi2s:  ",  dchip,   dchia, u0)
        input("Enter a number !!!")
    
    mbs=float(mbs1+mbs2)*0.5
    fb= float( fb1+ fb2)*0.5
    
    for k in range(nd): 
        F[k,(k+1):]= F[(k+1):,k]
    for k in range(na): 
        G[k,(k+1):]= G[(k+1):,k]
    
       
    
    corr1=0.0;
    inverF= np.mat(F).I
    inverG= np.mat(G).I
    for k in range(nd):  
        Era[k]=np.sqrt(abs(inverF[k,k]))  
    corr1= inverF[2,3]/(Era[2]*Era[3])## correlation between tE, xi         
    for k in range(na):  
        Erb[k]=np.sqrt(abs(inverG[k,k]))     
    if(all(Era)>0.0 and all(Erb)>0.0):    
        resu[0]=abs(Era[0]*1.0/(t0+eps))  
        resu[1]=abs(Era[1]*1.0/(u0+eps))
        resu[2]=abs(Era[2]*1.0/(tE+eps))  
        resu[3]=abs(Era[3]*1.0/(xi+eps))
        resu[4]=abs(Era[4]*1.0/(fb+eps))   
        resu[5]=abs(Era[5]*1.0/(mbs+eps))  
        resu[6]=abs(Era[6]*1.0/(piE+eps))
        resu[7]=abs(Erb[0]*1.0/(tetE+eps))
        resu[8]=abs(Erb[1]*1.0/(mus1+eps))
        resu[9]=abs(Erb[2]*1.0/(mus2+eps))
        resu[10]=np.sqrt((resu[7])**2.0 + (resu[6])**2.0 )### relative error in mass 
        resu[11]=abs(resu[10]*(Ds-Dl)/Ds)## relative error in Dl 
        f1= abs(resu[7]**2.0 + resu[2]**2.0 + (Era[3]*np.tan(xi))**2.0 -2.0*resu[2]*Era[3]*np.tan(xi)*corr1 )
        f2= abs(resu[7]**2.0 + resu[2]**2.0 + (Era[3]/np.tan(xi))**2.0 -2.0*resu[2]*Era[3]/np.tan(xi)*corr1 )
        sig_mul1=np.sqrt( Erb[1]**2.0 + mul**2.0 *(np.cos(xi))**2.0 *f1)
        sig_mul2=np.sqrt( Erb[2]**2.0 + mul**2.0 *(np.sin(xi))**2.0 *f2)  
        resu[12]= abs(sig_mul1)/abs(mul1+eps)
        resu[13]= abs(sig_mul2)/abs(mul2+eps)
        resu[14]=np.sqrt(sig_mul1*sig_mul1 + sig_mul2*sig_mul2)/(mul+eps)#relative error in Lens proper motion
        resu[15]=np.sqrt(Erb[1]**2.0 + Erb[2]**2.0)/(mus+eps)#            relative error in source proper motion
        resu[16]=abs(Erb[3]*1.0/(piE+eps))## relative error in piE through astrometry 
        resu[:]=resu[:]*pp### in persent 
        
        chis[i, :]= piE, np.log10(dchip), np.log10(resu[6]) , np.log10(dchia), np.log10(resu[16])
        
        '''        
        print ("******************** counter:  ", i )
        print ("Err(t0)*100.0/t0:  ",     resu[0] )
        print ("Err(u0)*100.0/u0:  ",     resu[1] ) 
        print ("Err(tE)*100.0/tE:  ",     resu[2]  )
        print ("Err(ksi)*100.0/ksi:  ",   resu[3] )
        print ("Err(fb)*100.0/fb:  ",     resu[4]   )          
        print ("Err(mbs)*100.0/mbs:  ",   resu[5] )
        print ("Err(piE)*100.0/piE:  ",   resu[6]   )
        print ("Err(tetE)*100.0/tetE:  ", resu[7] )
        print ("Err(mus1)*100.0/mus1:  ", resu[8])
        print ("Err(mus2)*100.0/mus2:  ", resu[9])
        print ("Err(Mass)*100.0/Mass: " , resu[10]) 
        print ("Err(Dl)*100.0/Dl: " ,     resu[11] )
        print ("Err(mul1)*100.0/mul1: " , resu[12] )
        print ("Err(mul2)*100.0/mul2: " , resu[13] ) 
        print ("Err(mul)*100.0/mul : " ,  resu[14] ) 
        print ("Err(mus)*100.0/mus : " ,  resu[15] ) 
        ''' 
        fij=open("./files/HistoMCC/resultMCC1.dat","a")
        np.savetxt(fij,resu.reshape((1,nq)),fmt="%.7e %.7e  %.7e  %.7e  %.7e  %.7e  %.7e  %.7e %.7e  %.7e  %.7e   %.7e   %.7e  %.7e  %.7e  %.7e  %.7e") 
        fij.close() 
        
        
        if(np.isfinite(resu).all and all(resu)>0.0):     
            for j in range(nq): 
                if(resu[j]< 7.0):  numc[j]+= 1.0
                if(resu[j]< 4.0):  numb[j]+= 1.0
                if(resu[j]< 1.0):  numa[j]+= 1.0
                if(resu[j]==0.0):  resu[j]= +1.0e-10
 
            if(resu[6]<7.0 or resu[16]<7.0):     numc[nq]+=1.0      
            if(resu[6]<4.0 or resu[16]<4.0):     numb[nq]+=1.0      
            if(resu[6]<1.0 or resu[16]<1.0):     numa[nq]+=1.0      
            
            if(resu[7]<1.0  and (resu[6]<1.0  or resu[16]<1.0)  and resu[2]<1.0):   epsi[0]+=1.0
            if(resu[7]<4.0  and (resu[6]<4.0  or resu[16]<4.0)  and resu[2]<4.0):   epsi[1]+=1.0
            if(resu[7]<7.0 and  (resu[6]<7.0  or resu[16]<7.0)  and resu[2]<7.0):   epsi[2]+=1.0
                
            arry[l,:nq]=resu
            arry[l, nq]=   func(Ml,0)
            arry[l, nq+1]= func(float(Dl/Ds),1)
            arry[l, nq+2]= func(Ds,2)
            arry[l, nq+3]= func(mbs,3)
            
            hisd0[l,:]= np.array([np.log10(tE), mbs, np.log10(piE), np.log10(tetE), t0/yr,  np.log10(mus), np.log10(mul), u0, Dl,lpir])
            l+=1
            if(resu[6]<4.0 or resu[16]<4.0):## parallax measuring in either photometry or astrometry
                hisd2[nph,:]= np.array([np.log10(tE),mbs, np.log10(piE), np.log10(tetE), t0/yr,  np.log10(mus), np.log10(mul), u0, Dl,lpir])
                nph+=1    
            if(resu[16]>4.0 and resu[6]<4.0):## parallax measuring in only photometry     
                hisd3[nas,:]= np.array([np.log10(tE),mbs, np.log10(piE), np.log10(tetE), t0/yr,  np.log10(mus), np.log10(mul), u0, Dl,lpir])
                nas+=1   
            if(resu[6]>4.0 and resu[16]<4.0): ##parallax measuring in only astrometry 
                hisd4[ngo,:]= np.array([np.log10(tE),mbs, np.log10(piE), np.log10(tetE), t0/yr,  np.log10(mus), np.log10(mul), u0, Dl,lpir])
                ngo+=1 

print ("fraction of halo-lensing events:  ",  float(halo*100.0/l)  )
print ("*****************************************************" )
print ("tot_number:  ", nm, " accepted number:  ",  l, "detected number:  ",  nc) 
print ("Efficiency for detectign parallax in photometry(only): ",   float(nas*100.0/l) )
print ("Efficiency for detectign parallax in astrometry(only): ",   float(ngo*100.0/l) )
print ("Efficiency for detectign parallax in either astro or photo: ",   float(nph*100.0/l) )
print ("*****************************************************" )
######################################################################

xdet=[r"$\log_{10}[t_{\rm E}(days)]$",  r"$m_{\rm{base}}(\rm{mag})$", r"$\log_{10}[\pi_{\rm{E}}]$",  r"$\log_{10}[\theta_{\rm{E}}(\rm{mas})]$",  r"$t_{0}(years)$", r"$\log_{10}[\mu_{s}(\rm{mas/days})]$", r"$\log_{10}[\mu_{\rm{l}}(\rm{mas/days})]$",  r"$u_{0}$",  r"$D_{\rm{l}}(\rm{kpc})$", r"$\log_{10}[\pi_{\rm{rel}}(\rm{mas})]$" ]# r"$D_{\rm{s}}(\rm{kpc})$"]

xd1=[2.0, 15.5, -4.0,-1.0, 0.0, -2.5,-2.1, 0.0, 0.001, -5.0]
xd2=[3.45,23.5, -0.5, 3.0, 10.0,-2.2,-0.0, 1.0, 50.0, 2.0]

lab1=[np.mean(np.power(10,hisd0[:l,0])), np.mean(hisd0[:l,1]), np.mean(np.power(10,hisd0[:l,2])), 
      np.mean(np.power(10,hisd0[:l,3])), np.mean(hisd0[:l,4]), np.mean(np.power(10,hisd0[:l,5])), 
      np.mean(np.power(10,hisd0[:l,6])), np.mean(hisd0[:l,7]), np.mean(hisd0[:l,8]), np.mean(hisd0[:l,9])]

lab3=[np.mean(np.power(10,hisd2[:nph,0])), np.mean(hisd2[:nph,1]), np.mean(np.power(10,hisd2[:nph,2])), 
      np.mean(np.power(10,hisd2[:nph,3])), np.mean(hisd2[:nph,4]), np.mean(np.power(10,hisd2[:nph,5])), 
      np.mean(np.power(10,hisd2[:nph,6])), np.mean(hisd2[:nph,7]), np.mean(hisd2[:nph,8]),         np.mean(hisd2[:nph,9])]

lab4=[np.mean(np.power(10,hisd3[:nas,0])), np.mean(hisd3[:nas,1]), np.mean(np.power(10,hisd3[:nas,2])), 
      np.mean(np.power(10,hisd3[:nas,3])), np.mean(hisd3[:nas,4]), np.mean(np.power(10,hisd3[:nas,5])), 
      np.mean(np.power(10,hisd3[:nas,6])), np.mean(hisd3[:nas,7]), np.mean(hisd3[:nas,8]), np.mean(hisd3[:nas,9])]

lab5=[np.mean(np.power(10,hisd4[:ngo,0])), np.mean(hisd4[:ngo,1]), np.mean(np.power(10,hisd4[:ngo,2])), 
      np.mean(np.power(10,hisd4[:ngo,3])), np.mean(hisd4[:ngo,4]), np.mean(np.power(10,hisd4[:ngo,5])), 
      np.mean(np.power(10,hisd4[:ngo,6])), np.mean(hisd4[:ngo,7]), np.mean(hisd4[:ngo,8]), np.mean(hisd4[:ngo,9])]


for i in range(10):
    plt.clf()
    plt.clf()
    fig=plt.figure(figsize=(8,6))
    ax1=plt.gca()
    plt.hist(hisd0[:l,i],  25,histtype='bar', ec='darkgreen',facecolor='g',alpha=0.4,rwidth=1.5,label=r"$\rm{Mean=}$"+str(round(lab1[i],2)))
    plt.hist(hisd2[:nph,i],25,histtype='step',color='r',alpha=0.9, lw=2.1,label= r"$\rm{Mean=}$"+str(round(lab3[i],2))   )
    plt.hist(hisd3[:nas,i],25,histtype='step',color='b',alpha=0.9, lw=2.1,label= r"$\rm{Mean=}$"+str(round(lab4[i],2))   )
    plt.hist(hisd4[:ngo,i],9, histtype='step',color='k',alpha=0.9, lw=2.1,label= r"$\rm{Mean=}$"+str(round(lab5[i],2))   )
    y_vals=ax1.get_yticks()
    ax1.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/len(hisd0[:l,i]))) for x in y_vals]) 
    y_vals = ax1.get_yticks()
    ax1.set_ylim([ 0.5, np.max(y_vals) ])
    plt.yscale('log')
    plt.yticks(fontsize=20, rotation=0)
    plt.xticks(fontsize=20, rotation=0)
    ax1.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=20,labelpad=0.1)
    ax1.set_xlabel(str(xdet[i]) , fontsize=20,  labelpad=0.1)
    plt.xlim([ xd1[i],  xd2[i] ])
    plt.grid("True")
    plt.grid(linestyle='dashed')
    plt.legend()
    plt.legend(loc='best',fancybox=True, shadow=True)
    plt.legend(prop={"size":19})
    fig= plt.gcf()
    plt.subplots_adjust(hspace=.0)
    fig.tight_layout(pad=0.5)
    fig.savefig("./files/HistoMCC/hdetC{0:d}.jpg".format(i),dpi=200)
print ("****  All histo_detected_parameters are plotted **********************")   

################################################################################        
 
nam=[r"$\log_{10}[\sigma_{t_{0}}\big/t_{0}(\%) ]$", 
     r"$\log_{10}[\sigma_{u_{0}}\big/u_{0}(\%)]$", 
     r"$\log_{10}[\sigma_{t_{\rm E}}\big/t_{\rm E}(\%)]$", 
     r"$\log_{10}[\sigma_{\xi}\big/ \xi(\%)]$",     
     r"$\log_{10}[\sigma_{f_{b}}\big/f_{b}(\%)]$", 
     r"$\log_{10}[\sigma_{m_{\rm base}}\big/m_{\rm base}(\%)]$", 
     r"$\log_{10}[\sigma_{\pi_{\rm E}}\big/ \pi_{\rm E}(\%)]$",  
     r"$\log_{10}[\sigma_{\theta_{\rm E}}\big/ \theta_{\rm E}(\%)]$",  
     r"$\log_{10}[\sigma_{\mu_{\rm s, n1}}\big/\mu_{\rm s, n1}(\%)]$", 
     r"$\log_{10}[\sigma_{\mu_{\rm s, n2}}\big/ \mu_{\rm s, n2}(\%)]$", 
     r"$\log_{10}[\sigma_{M_{\rm l}}\big/ M_{\rm l}(\%)]$", 
     r"$\log_{10}[\sigma_{D_{\rm l}}\big/ D_{\rm l}(\%)]$", 
     r"$\log_{10}[\sigma_{\mu_{\rm l, n1}}\big/ \mu_{\rm l, n1}(\%)]$", 
     r"$\log_{10}[\sigma_{\mu_{\rm l, n2}}\big/ \mu_{\rm l, n2}(\%)]$",
     r"$\log_{10}[\sigma_{\mu_{\rm l}}\big/ \mu_{\rm l}(\%)]$", 
     r"$\log_{10}[\sigma_{\mu_{\rm s}}\big/ \mu_{\rm s}(\%)]$",
     r"$\log_{10}[\sigma_{\pi_{\rm E}}\big/ \pi_{\rm E}(\%)]$"]     
     
################################################################################ 

nb1=0.0; nb2=0.0
for i in range(l): 
    if(arry[i,6]< 4.0):  
        nb1+=1;
for i in range(l): 
    if(arry[i,16]<4.0):  
        nb2+=1;        
effm=round(float(nb1*100.0/l),2)
effa=round(float(nb2*100.0/l),2)
print ("Effis:  ",  effm,   effa)
plt.clf()
fig= plt.figure(figsize=(8,6))
ax= plt.gca()  
plt.hist(np.log10(arry[:l,6]),40,histtype='bar',ec='darkgreen',facecolor='g',alpha=0.85, rwidth=1.5,label=r"$\rm{Photometry}$")
plt.hist(np.log10(arry[:l,16]),40,histtype='bar',ec='b',  facecolor='b',alpha=0.3, rwidth=1.5,label=r"$\rm{Astrometry}$")
y_vals = ax.get_yticks()
ax.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/l)) for x in y_vals]) 
y_vals = ax.get_yticks()
plt.ylim([np.min(y_vals), np.max(y_vals)])
plt.xlim([ -2, 6.0])
plt.title(r"$\rm{Simulation}~(\rm{C})$", fontsize=19)
plt.axvline(x=np.log10(4.0), color='k', linestyle='--', lw=2.)
plt.text(-1.95, np.max(y_vals)*0.85, r"$\epsilon_{\rm{m}}[\%]$="+str(effm), color="g", fontsize=18)
plt.text(-1.95, np.max(y_vals)*0.75, r"$\epsilon_{\rm{a}}[\%]$="+str(effa), color="k", fontsize=18)
plt.xlabel(str(nam[6]), fontsize=19,labelpad=0.15)
ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
plt.xticks(fontsize=19, rotation=0)
plt.yticks(fontsize=19, rotation=0)
plt.legend()
plt.legend(loc='upper right',fancybox=True, shadow=True)
plt.legend(prop={"size":18})
fig=plt.gcf()
fig.savefig("./files/HistoMCC/herrpiC.jpg",dpi=200)
print("Effi parallax detection in photometry,  in astrometry:  ",   effm,   effa)
################################################################################     
        
x1=[-3.0,-1.5,-1.0, -2.0,-1.0, -3.0, -2.0,-3.5, -2.0, -2.5, -1.5, -1.5, -1.5 ,  -1.0, -1.0,-2.0,-2.0]
x2=[ 3.0, 3.0, 3.0,  6.0, 3.0,  1.0,  6.0, 2.1,  1.5,  1.0,  5.0 , 3.0 , 5.0,    5.0,  5.0, 1.1, 6.0]    
for i in range(nq): 
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    ax= plt.gca()  
    plt.hist(np.log10(arry[:l,i]),40, histtype='bar',ec='darkgreen',facecolor='green',alpha=0.55, rwidth=1.5)
    y_vals = ax.get_yticks()
    ax.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/l)) for x in y_vals]) 
    y_vals = ax.get_yticks()
    plt.ylim([np.min(y_vals)+0.01, np.max(y_vals)*0.95])
    plt.xlim([x1[i], x2[i]])
    plt.axvline(x=0.6989 , color='b', linestyle='-', lw=1.6)
    plt.xlabel(str(nam[i]), fontsize=19,labelpad=0.15)
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    fig=plt.gcf()
    fig.savefig("./files/HistoMCC/hErrorC{0:d}.jpg".format(i),dpi=200)
    numa[i]=float(numa[i]*100.0/l)
    numb[i]=float(numb[i]*100.0/l)
    numc[i]=float(numc[i]*100.0/l)
numa[nq]=float(numa[nq]*100.0/l)    
numb[nq]=float(numb[nq]*100.0/l)    
numc[nq]=float(numc[nq]*100.0/l)            
epsi[0]=float(epsi[0]*100.0/l)
epsi[1]=float(epsi[1]*100.0/l)
epsi[2]=float(epsi[2]*100.0/l)      
print ("****  All histo_Error_VS_parameters are plotted **********************")       
################################################################################      
        
fii=open("./files/HistoMCC/resultMCC2.dat","a")
param=np.array([numa[7], numa[nq], numa[6], numa[16], numa[10], numa[11], numa[2], epsi[0] ]) 
np.savetxt(fii,param.reshape((1,8)),fmt="$%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$\n") 

param=np.array([numb[7], numb[nq], numb[6], numb[16], numb[10], numb[11], numb[2], epsi[1] ]) 
np.savetxt(fii,param.reshape((1,8)),fmt="$%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$\n") 

param=np.array([numc[7], numc[nq], numc[6], numc[16], numc[10], numc[11], numc[2], epsi[2] ]) 
np.savetxt(fii,param.reshape((1,8)),fmt="$%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$\n") 
fii.close()         
      
      
        
################################################################################
nla=[r"$M_{\rm l}(M_{\odot})$", r"$D_{\rm{l}}(\rm{kpc})$", r"$\log_{10}[\pi_{\rm{rel}}(\rm{mas})]$",
r"$\log_{10}[\theta_{\rm E}(\rm{mas})]$", r"$\log_{10}[\pi_{\rm E}]$" , r"$struct$" , r"$\log_{10}[t_{\rm E}(days)]$" , r"$u_{0}$"]
parp1=np.zeros((nm, 9));  n1=0
parp2=np.zeros((nm, 9));  n2=0
for i in range(nm): 
    if(parp[i,6]>3):#LMC 
        parp1[n1, :]=parp[i,:]
        n1+=1
    else:#GH
        parp2[n2,:]=parp[i,:]
        n2+=1  
x1=[0.2,0.0,-4.0,-1.5,-4.5,-0.1, 1.8,0.0]
x2=[2.8,50.0, 1.0, 2.8, 0.0,6.1,  3.4,1.0]        
for i in range(8):  
    plt.clf()
    plt.cla()
    fig=plt.figure(figsize=(8,6))
    plt.scatter(parp1[:n1,i+1], parp1[:n1,0]*1.0, marker='^',c='b', s=9.0, alpha=1.0)
    plt.scatter(parp2[:n2,i+1], parp2[:n2,0]*1.0, marker='o',c='b', s=9.0, alpha=1.0)
    #plt.xlim([x1[i], x2[i]])    
    plt.yscale('log')
    plt.xlabel(str(nla[i]), fontsize=21)
    plt.ylabel(r"$\delta\theta_{\rm{RMS}}(\rm{mas}), \rm{Parallax}-\rm{induced}$", fontsize=21)
    plt.xticks(fontsize=21,  rotation=0)
    plt.yticks(fontsize=21, rotation=0)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    plt.subplots_adjust(hspace=.0)
    fig=plt.gcf()
    fig.tight_layout(pad=0.1)
    fig.savefig("./files/HistoMCC/ErrsC{0:d}.jpg".format(i),dpi=200)
print ("****  pin/piM scatter  are plotted ************************************")               

################################################################################
'''
#piE, dchip, resu[6] , dchia, resu[16]
plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
plt.scatter(chis[:,0], chis[:,1], marker='o',c='b', s=5.0, label=r"$\rm{Photometry}$")
plt.scatter(chis[:,0], chis[:,3], marker='o',c='r', s=5.0, label=r"$\rm{Astrometry}$")
plt.xlabel(r"$\pi_{\rm{E}}$", fontsize=18)
plt.ylabel(r"$\log_{10}[\Delta \chi^{2}]$", fontsize=18)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xscale('log')
plt.grid("True")
plt.grid(linestyle='dashed')
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":18})
fig=plt.gcf()
fig.savefig("./files/HistoMCC/Scatter1.jpg",dpi=200)
################################################################################
plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
plt.scatter(chis[:,0], chis[:,2], marker='o',c='b', s=5.0, label=r"$\rm{Photometry}$")
plt.scatter(chis[:,0], chis[:,4], marker='o',c='r', s=5.0, label=r"$\rm{Astrometry}$")
plt.xlabel(r"$\pi_{\rm{E}}$", fontsize=18)
plt.ylabel(r"$\log_{10}[\sigma\pi_{\rm{E}}/\pi_{\rm{E}}]$", fontsize=18)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xscale('log')
plt.grid("True")
plt.grid(linestyle='dashed')
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":18})
fig=plt.gcf()
fig.savefig("./files/HistoMCC/Scatter2.jpg",dpi=200)


################################################################################
plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
plt.scatter(chis[:,1], chis[:,2], marker='o',c='b', s=5.0, label=r"$\rm{Photometry}$")
plt.xlabel(r"$\log_{10}[\Delta \chi^{2}_{n}]$", fontsize=18)
plt.ylabel(r"$\log_{10}[\sigma\pi_{\rm{E}}/\pi_{\rm{E},~P}]$", fontsize=18)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":18})
fig=plt.gcf()
fig.savefig("./files/HistoMCC/Scatter3.jpg",dpi=200)
################################################################################
plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
plt.scatter(chis[:,3], chis[:,4], marker='o',c='b', s=5.0, label=r"$\rm{Astrometry}$")
plt.xlabel(r"$\log_{10}[\Delta \chi^{2}_{n}]$", fontsize=18)
plt.ylabel(r"$\log_{10}[\sigma\pi_{\rm{E}}/\pi_{\rm{E},~A}]$", fontsize=18)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":18})
fig=plt.gcf()
fig.savefig("./files/HistoMCC/Scatter4.jpg",dpi=200)
################################################################################
'''



