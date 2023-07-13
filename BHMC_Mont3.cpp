#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>

using namespace std;

#define Nx 7
#define Ny 4

const int Num=10000;
const double MaxD=65.0;///kpc
const double step=MaxD/(double)Num/1.0;///step in kpc

const double pi= M_PI;
const double RA=180.0/M_PI;
const double Hp=6.62607004*pow(10.0,-34);
const double KP=3.08568025*pow(10.,19); // in meter.
const double G=6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity= 3.0*pow(10.0,8.0);//velosity of light
const double Msun=1.98892*pow(10., 30); //in [kg].
const double Rsun= 6.957*pow(10.0,8.0); ///solar radius [meter]
const double vro_sun=226.0;
const double AU=1.4960*pow(10.0,11.0);
const double year=float(364.0);
const double binary_fraction=double(2.0/3.0);
const double VSunR =11.1;
const double VSunT =vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ =7.25;
const double Dsun=8.0; 

///============================ Besancon constant ==============///

const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[7]={3.1,2.5,3.1,3.1,3.1,3.1,2.5};
const int    N1=27004, N2=36558, N3=2612, N4=491;///CMD_BESANCON, thinD, bulge, thickD, halo

///===============LSST constant===============================///
const int M=6;///No. of filter  ugrizy  of LSST
const double Tobs=10.0*year;///LSST observational time 10 years
const double texp=30.0;///in second be sure about Dcm[????
const double delta2=0.005;///systematic errors
const double gama[M]=  {0.037,0.038,0.039,0.039,0.040,0.040};///a function of sky brightness and airmass in different wavelengths.
const double seeing[M]={0.77,0.73,0.70,0.67,0.65,0.63};///seeing
const double msky[M]=  {22.9,22.3,21.2,20.5,19.6,18.6};
const double Cm[M]=    {22.92,24.29,24.33,24.20,24.07,23.69};
const double Dci[M]=   {0.67,0.21,0.11,0.08,0.05,0.04};
const double km[M]=    {0.451,0.163,0.087,0.065,0.043,0.138};///sky extinction
const double sigma[M+1]={0.022,0.02,0.017,0.017,0.027,0.027,0.04};
const double wav[M+1]= {0.3671,0.4827,0.6223,0.7546,0.8691,0.9712,2.1739};///in[micrometer] u g r i z y K
const double AlAv[M+1]={1.48,  1.05,  0.8, 0.6,  0.45, 0.4, 0.118}; 
const double thre[M+1]={23.4,24.6,24.3,23.6,22.9,21.7,23.0};///depth of single visit in ugrizy+  K-bad
const double satu[M+1]={15.2,16.3,16.0,15.3,14.6,13.4,12.0};///saturation limit of single visit in ugrizy+ K-band
const double FWHM[M+1]={1.22087,1.10136,0.993103,0.967076,0.951766,0.936578, 5.0*0.004};//LSST [arcsec] ugrizy + K-band


const int Nm=82;//rows in sigmaM_LSST.txt   https://arxiv.org/abs/1201.5140
const int Na=96;//rows in sigmaA_LSST.txt
const int Nl=20000;//filter LSST rows
const int YZ=3578;//No.yzma.txt rows
const int met= 70;//No rows metal.txt  

const double cade1=3.0;//in days 
const double cade2=10.0;//ELT cadence 
const double dt=cade1;
const int coun=int(Tobs/cade1 +10) + int(Tobs/cade2+10);

////========================== LMC =============================
const double bLMC= -32.9; 
const double lLMC= 280.5; 
const double RaLMC= 80.89375;
const double DecLMC=-68.2438888888889;  
const double sLMC= 5.0;///degree (half of size of LMC)
const double DLMC= 49.97;///KPC 

const int nn=int(455);
const double dec1= -72.5; 
const double dec2= -64.8; 
const double rig1=float((4.0+27.0/60.0)*15.0 );
const double rig2=float((5.0+50.0/60.0)*15.0 );
const double ddec=5.0*fabs(dec2-dec1)/nn/1.0;
const double drig=5.0*fabs(rig2-rig1)/nn/1.0;
const int nrd=10000; 
const int Nel=int(7);//https://academic.oup.com/mnras/article/494/3/4413/5813442

///================================================================
struct source{
    int nums,struc, cl;
    double Ds,TET,FI, lat, lon;
    double right, decli;
    double od_disk,od_ThD,od_bulge,od_halo,opt;///optical  depth
    double od_dlmc, od_blmc, od_hlmc; 
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double rho_hlmc[Num], rho_dlmc[Num],rho_blmc[Num]; 
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs, Romins,nstart, nstarti;
    double Nblend[M+1], blend[M+1], Fluxb[M+1], magb[M+1], Ai[M+1], Mab[M+1], Map[M+1];
    double logT, logL, vs;
    double Romax, deltao;
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
    double mus1, mus2, mul1, mul2;
    double xv, yv, zv; 
    double pos1, pos2,pos10, pos20, Astar, xi, ux, uy;
    double  def1a,  def2a, def1b,  def2b,  def1c,  def2c; 
    double fb[2], mbs[2];
    double ut, ut0; 
    double ampM, ampA, errM, errA;  
};
struct lens{
    int numl,struc;
    double Ml, Dl, vl , Vt, xls;
    double rhomaxl,tE,RE;
    double mul,u0;
    double t0;
    double tetE;
    double piE, pirel; 
};
struct astromet{
   double tetp;
   double omegae;
   double vearth;
   double fact;
   double Ve_n1, Ve_n2;
   double ue_n1, ue_n2;
};
struct CMD{
    double logT_d[N1], logL_d[N1], Mab_d[M][N1], Mk_d[N1]; int cl_d[N1];  ///thin disk
    double logT_b[N2], logL_b[N2], Mab_b[M][N2], Mk_b[N2]; int cl_b[N2];  /// bulge
    double logT_t[N3], logL_t[N3], Mab_t[M][N3], Mk_t[N3]; int cl_t[N3];  ///thick disk
    double logT_h[N4], logL_h[N4], Mab_h[M][N4], Mk_h[N4]; int cl_h[N4];  /// halo
};
struct covarian{
   double FishM[Nx][Nx];
   double FishA[Ny][Ny];
};
struct lsst{
     double magl[Nm], errl[Nm];
     double maga[Na], erra[Na];
     double mage[Nel],ere1[Nel], ere2[Nel];  
     double airm[Nl];
     int  filter[Nl];         
};
struct galactic{
   double l[nrd], b[nrd], RA[nrd], Dec[nrd]; 
};
struct yfilter{
    double Age[YZ]; 
    double   B[YZ];  
    double   M[YZ];   
    double  mm[YZ]; 
    int   number[met];
    int    count[met];   
    double Metal[met];
}; 
///===================== FUNCTION ==============================
void read_cmd(CMD & cm,  yfilter & yf);
void func_source(source & s, CMD & cm);
void func_lens(lens & l, source & s);
void vrel(source & s , lens & l);
void Disk_model(source & s, int );
void optical_depth(source & s);
double RandN(double sigma,double Nn);
double RandR(double down, double up);
void lightcurve(source & s, lens & l, astromet & as, double tim); 
double ylsst(yfilter & yf, double metals, double ages);  
double  errlsst(lsst & ls, double ghadr, int flag);
double   errELT(lsst & ls, double ghadr, int flag);
double errlsst2(double maga, int fi, double airm);
//////////////////////////////////////////////////////
    time_t _timeNow;
      unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;
///////////////////////////////////////////////////////
///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//
int main()
{
///****************************************************************************
//gettimeofday(&newTime,NULL);
//ftime(&tp);
	time(&_timeNow);
	_randStream = fopen("/dev/urandom", "r");
	_dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
    _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
///****************************************************************************
     source s;
     lens l;
     CMD cm;
     lsst ls; 
     astromet as;
     galactic ga;
     covarian co; 
     yfilter yf;      
     
///===========================================================
    FILE * lsstf;
    FILE* convert; 
    FILE *meta; 
    FILE *yze; 
    FILE* magg;  
    FILE* data3; 
    FILE * elt;
    elt=fopen("./files/sigma_ELT.txt","r");
    if(!elt){cout<<"cannot read sigma_ELT.txt: "<<endl; int yye;  cin>>yye;}
    for(int i=0; i<Nel;++i){
    fscanf(elt,"%lf   %lf   %lf\n",&ls.mage[i], &ls.ere1[i], &ls.ere2[i] );
    ls.ere2[i]=ls.ere2[i]*0.001;//mili arcsec
    cout<<"mage:  "<<ls.mage[i]<<"\t err_magni(mag): "<<ls.ere1[i]<<"\t err_astro(mili arcsec):  "<<ls.ere2[i]<<endl;}
    cout<<"**** File sigma_ELT.txt was read ****"<<endl;
    fclose(elt);
    
    
    /*  we calculate the LSST error using the formula (directly) 
    lsstf=fopen("./files/sigmaM_LSST.txt","r");
    if(!lsstf){cout<<"cannot read  sigmaM_LSST.txt: "<<endl; int yye;  cin>>yye;}
    for(int i=0; i<Nm;++i){
    fscanf(lsstf,"%lf   %lf \n",&ls.magl[i], &ls.errl[i] );
    ls.errl[i]=ls.errl[i]*0.001; }//LSST photometry accuracy [mag]
    cout<<"**** File sigmaM_LSST.txt was read ****"<<endl;
    fclose(lsstf); */
    
    

    
    lsstf=fopen("./files/sigmaA_LSST.txt","r");
    if(!lsstf){cout<<"cannot read sigmaA_LSST.txt: "<<endl; int yye;  cin>>yye;}
    for(int i=0; i<Na;++i){
    fscanf(lsstf,"%lf   %lf\n",&ls.maga[i], &ls.erra[i]); }//astrometry LSST accuracy in [marcse]
    cout<<"**** File sigmaA_LSST.txt was read ****"<<endl;
    fclose(lsstf);
    
   
   
 
    lsstf=fopen("./files/filterLSST.txt","r");
    if(!lsstf){cout<<"cannot read filterLSST.txt: "<<endl; int yye;  cin>>yye;}
    for(int i=0; i<Nl; ++i){
    fscanf(lsstf,"%lf   %d\n",&ls.airm[i],&ls.filter[i]);  
    if(ls.airm[i]<0.0 or ls.filter[i]>=6){ 
    cout<<"ERROR airmass: "<<ls.airm[i]<<"\t filter: "<<ls.filter[i]<<endl;  int yye; cin>>yye; }}
    cout<<"**** File filterLSST.txt was read ****"<<endl;    
    fclose(lsstf);
   
   
  
    convert=fopen("./files/convert_coordinate_2.dat","r");
    if(!convert){cout<<"cannot read convert_coordinate.txt: "<<endl; int yye;  cin>>yye;}
    int coo=0; 
    for(int i=0; i<100; ++i){
    for(int j=0; j<100; ++j){
    fscanf(convert,"%lf  %lf   %lf  %lf \n",&ga.RA[coo], &ga.Dec[coo], &ga.l[coo], &ga.b[coo]);  coo++; }}
    cout<<"**** File convert_coordinate_2.dat was read ****"<<endl;   
    fclose(convert); 
   
   
  

    meta=fopen("./files/metal.txt","r"); 
    if(!meta){cout<<"cannot read metal.txt: "<<endl; int yye;  cin>>yye;}
    for(int i=0; i<met; ++i){
    fscanf(meta,"%lf   %d    %d\n",&yf.Metal[i],&yf.count[i],&yf.number[i]);}
    //cout<<"metal:  "<<yf.Metal[i]<<"\t count:  "<<yf.count[i]<<"\t number "<<yf.number[i]<<endl;
    cout<<"**** File metal.txt was read ****"<<endl;
    fclose(meta); 
    
    

    yze=fopen("./files/yzma.txt", "r"); 
    if(!yze){cout<<"cannot read yzma.txt: "<<endl; int yye;  cin>>yye;}
    for(int i=0; i<YZ; ++i){
    fscanf(yze,"%lf   %lf   %lf  %lf\n",&yf.Age[i],&yf.mm[i],&yf.B[i],&yf.M[i]);}
    cout<<"**** File yzma.txt was read ****"<<endl;
    fclose(yze); 
    
    read_cmd(cm, yf);                            
    
   
///===========================================================
   FILE *fil0;
   fil0=fopen("./files/density/right.txt", "w"); 
   fclose(fil0);

    FILE* fil1;
    fil1=fopen("./files/MONT1/IBH_MONT1c.dat","w");
    fclose(fil1);

    FILE* fil2;
    fil2=fopen("./files/MONT1/IBH_MONT2c.dat","w");
    fclose(fil2);
    
    FILE* fil3;
    fil3=fopen("./files/MONT1/IBH_MONT3c.dat","w");
    fclose(fil3);

///===========================================================

    double obs[10][2]={0.0};
    obs[0][0]=0.0*year+0.0;    obs[0][1]=0.0*year+float(7.0*30.5);
    obs[1][0]=1.0*year+0.0;    obs[1][1]=1.0*year+float(7.0*30.5);
    obs[2][0]=2.0*year+0.0;    obs[2][1]=2.0*year+float(7.0*30.5);
    obs[3][0]=3.0*year+0.0;    obs[3][1]=3.0*year+float(7.0*30.5);
    obs[4][0]=4.0*year+0.0;    obs[4][1]=4.0*year+float(7.0*30.5);
    obs[5][0]=5.0*year+0.0;    obs[5][1]=5.0*year+float(7.0*30.5);
    obs[6][0]=6.0*year+0.0;    obs[6][1]=6.0*year+float(7.0*30.5);
    obs[7][0]=7.0*year+0.0;    obs[7][1]=7.0*year+float(7.0*30.5);
    obs[8][0]=8.0*year+0.0;    obs[8][1]=8.0*year+float(7.0*30.5);
    obs[9][0]=9.0*year+0.0;    obs[9][1]=9.0*year+float(7.0*30.5);




    int    counter=0, sea, filt;
    int    flag_det, hh, nml;
    int    flagf, ndw, fi, datf1, datf2, tt;
    char   filnam1[40], filnam2[40];
    double errs, errg, mindd, fel, amp;
    double magnio, test, dist, deltaA;
    double magni[M], magni0[M], Astar0;
    double timp1,timp2, lonn, chi, chi1, chi2, dchi2, dchi;
    double chi1a, chi2a, dchia, def, def0, sil; ///*** 
    double flag0, flag1, flag2;
    double lonn0, lat0, flago, prob1, prob2;
    double magw, derm1[2], derm2[2], dera1[2], derb1[2], dera2[2], derb2[2]; 
    double derm1f, derm2f, dera1f, derb1f, dera2f, derb2f, numc; 
    double Delta1[Nx], Delta2[Ny];
    float  sig[2]={+1.0, -1.0};  
   
    
    double *timn=new double[coun];
    double *magn=new double[coun];
    double *soux=new double[coun];
    double *souy=new double[coun];
    double *errm=new double[coun];
    double *erra=new double[coun];
    int    *tele=new int [coun];  
   

    as.tetp=double(M_PI/3.0);
    as.omegae=double(2.0*M_PI/year); ///radian per day
    as.vearth= as.omegae;//*0.001/(24.0*3600.0); ///km/s  on the Earth surface 
    as.fact=double(1000.0*3600.0*24.0/AU);




for (int icon=0;  icon<100; ++icon){
    
    
    nml=0;
    s.right= RaLMC + RandR(-2.0,2.0);  
    s.decli= DecLMC+ RandR(-2.0,2.0);   
    hh=-1; mindd=10000.0; 
    for(int i=0; i<10000; ++i){
    dist= sqrt( (s.right- ga.RA[i])*(s.right- ga.RA[i]) + (s.decli- ga.Dec[i])*(s.decli- ga.Dec[i]) ); 
    if(dist<=mindd){ mindd=dist; hh= i; }}
    if(hh<0){ cout<<"ERROR:right:    "<<s.right<<"\t s.dec:  "<<s.decli<<endl;   int uue;  cin>>uue; }
    s.lon=double(ga.l[hh]);   
    s.lat=double(ga.b[hh]);
    s.TET=(360.0-s.lon)/RA;
    s.FI=  s.lat/RA; 
    
    
    Disk_model(s, icon); 
   
 
    do{
    do{
    func_source(s, cm);
    func_lens(l, s);
    }while(l.tE<=0.1 or l.tE>4000.0); // which is around 10 years
    optical_depth(s);
    
    s.mus1= double(s.SV_n1- s.VSun_n1)*1000.0*3600.0*24.0/(s.Ds*AU);///mas/days
    s.mus2= double(s.SV_n2- s.VSun_n2)*1000.0*3600.0*24.0/(s.Ds*AU);///mas/days
    s.mul1= double(s.LV_n1- s.VSun_n1)*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days
    s.mul2= double(s.LV_n2- s.VSun_n2)*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days 
    l.piE=  double(l.pirel)/l.tetE; 
    
    
    test=(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    if(test<=s.blend[2] and s.magb[2]<=thre[2]){
        
   
    counter+=1;
    flagf=0;
    if(counter<200){
    flagf=1;
    sprintf(filnam1,"./files/MONT1/%c%c%c%c%d.dat",'m','a','g','C', counter);
    sprintf(filnam2,"./files/MONT1/%c%c%c%c%d.dat",'d','a','t','C', counter);
    magg=fopen(filnam1,"w");    
    data3=fopen(filnam2,"w");}
     
          
          
     for(int i=0; i<coun; ++i){
     timn[i]=0.0; magn[i]=0.0; soux[i]=0.0;  souy[i]=0.0;  errm[i]=0.0;  erra[i]=0.0; tele[i]=-1;}
     filt=int(RandR(1.0, Nl-2.0)); 
     flag0=flag1=flag2=0.0;
     ndw=0;  flag_det=0;
     chi=0.0;    chi1=0.0;  chi2=0.0;  chi1a=0.0;  chi2a=0.0;  
     timp1=0.0;  timp2=0.0; 
     fel=0.0;
     s.ampM=0.0; s.ampA=0.0; s.errM=1364540.0;  s.errA=146238290.0;  
     numc=0.0;  
     
    
     for(double tim=float(-3.0*year);  tim<float(12.0*year); tim=tim+dt){
     
     
     
     lightcurve(s,l,as,tim);  
     Astar0= double(s.ut0*s.ut0+2.0)/sqrt(s.ut0*s.ut0*(s.ut0*s.ut0+4.0));
     s.Astar=double(s.ut*s.ut + 2.0)/sqrt(s.ut*s.ut*(s.ut*s.ut+4.0));
     if(Astar0<1.0) {cout<<"Error Astar0:  "<<Astar0<<endl; int uue; cin>>uue;    Astar0=1.0; }
     if(s.Astar<1.0){cout<<"Error Astar:  "<<s.Astar<<endl; int uue; cin>>uue;    s.Astar=1.0; }
     
     for(int i=0;i<int(M+1); ++i){
     magni0[i] =s.magb[i]-2.5*log10( Astar0*s.blend[i] + 1.0 - s.blend[i]);
     magni[i]  =s.magb[i]-2.5*log10(s.Astar*s.blend[i] + 1.0 - s.blend[i]);}
     def= sqrt( s.pos1*s.pos1 + s.pos2*s.pos2 );  
     def0=sqrt(s.pos10*s.pos10 + s.pos20*s.pos20 ); 
     
     if(flagf>0) 
     fprintf(magg,"%.4lf  %.4lf   %.4lf  %.4lf    %.4lf   %.4lf   %.4lf  %.4lf   %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf   %.4lf  %.4lf    %.4lf %.4lf   %.4lf %.4lf   %.4lf %.4lf   %.4lf %.4lf   %.4lf %.4lf   %.4lf\n",
     tim,magni0[0], magni[0], Astar0, s.Astar, s.pos10, s.pos20, s.pos1, s.pos2, s.def1a, s.def2a, s.def1b, s.def2b, s.def1c, s.def2c,
     magni0[1], magni[1],magni0[2], magni[2],magni0[3], magni[3],magni0[4], magni[4],magni0[5], magni[5],magni0[6], magni[6]);//27
     
     
     
     
     if( fabs(tim-l.t0)<1.5*l.tE ){///j,جه شود که به دلیل تغییر زمان کمترین فاصله این قسمت قابل اعتماد نیست
     s.ampA+= float((s.def1a-s.def1c)*(s.def1a-s.def1c) + (s.def2a-s.def2c)*(s.def2a-s.def2c));///(s.def1a*s.def1a+s.def2a*s.def2a); 
     s.ampM+= float((s.Astar- Astar0)*(s.Astar-Astar0));///Astar0/Astar0;///}///Ampl in magnification 
     numc+=1.0;}
    
    
    
         
     
     flago=0.0;   
     sea=-1; 
     for(int kl=0; kl<10; kl=kl+1){
     if(float((obs[kl][0]-tim)*(obs[kl][1]-tim))<=0.0 and tim>0.0 and tim<Tobs){flago=1.0;  sea=kl;   break;}}
     if(flago>0.5){
     timp1 += dt; 
     timp2 += dt; 
     
     
     datf1=0;//LSST
     prob1=RandR(0.0,100.0);//Bad weather
     prob2=RandR(0.0,100.0);//Uniform observation
     if(timp1>=cade1){
     timp1=timp1- cade1;
     if(prob1>10.0 and prob2>=10.0) datf1=1;}
     
     
     datf2=0;//ELT
     prob1=RandR(0.0,100.0);//Bad weather
     prob2=RandR(0.0,100.0);//Uniform observation
     if(timp2>=cade2){
     timp2=timp2-cade2;
     if(prob1>10.0 and prob2>=10.0) datf2=1;}
     fel=0.0;
     if((tim<l.t0 and s.Astar>1.3) or tim>l.t0) fel=1.0; 
      
      
     ///LSST DATA **************************************
     if(datf1>0){
     filt+=1;  
     if(filt>(Nl-2)) filt=0;
     fi=int(ls.filter[filt]);
     if(magni[fi]>=satu[fi] and magni[fi]<=thre[fi]){
     errg= errlsst2(magni[fi], fi, ls.airm[filt]);
     errs= errlsst(ls, magni[2], 1); 
     magnio=magni[fi] + RandN(errg,3.0);
     chi  +=fabs( (magnio- s.magb[fi])*(magnio- s.magb[fi])/(errg*errg));
     chi1 +=fabs( (magnio-  magni[fi])*(magnio-  magni[fi])/(errg*errg));   
     chi2 +=fabs( (magnio- magni0[fi])*(magnio- magni0[fi])/(errg*errg)); 
     sil=  RandN(errs*sqrt(2.0), 3.0) ;
     chi1a+= fabs( (def-def0)*(def-def0)/(errs*errs*2.0) );
     chi2a+= fabs( (def+sil-def)*(def+sil-def)/(errs*errs*2.0) );
     timn[ndw]=tim;
     magn[ndw]=magni[2];///scale to r-LSST band, but the errors are different 
     errm[ndw]=errg;
     soux[ndw]=s.pos1;
     souy[ndw]=s.pos2;
     erra[ndw]=errs;
     tele[ndw]=0; 
     flag2=0.0;
     if(fabs(magnio-s.magb[fi])>fabs(4.0*errg))    flag2=1.0;
     if(ndw>2 and float(flag0+flag1+flag2)>2.0)   flag_det=1;
     flag0=flag1;
     flag1=flag2;
     deltaA=fabs(pow(10.0,-0.4*errg)-1.0)*s.Astar;
     if(deltaA<s.errM)  s.errM=deltaA; 
     if(errs< s.errA)   s.errA=errs; 
     //s.errM += deltaA; 
     //s.errA += errs; 
     if(flagf>0) 
     fprintf(data3,"%.4lf  %.4lf  %.6lf  %.6lf  %.4lf  %.4lf   %.4lf  %.6lf  %.1f  %d   %d\n",
     tim,magnio,errg,s.Astar+RandN(deltaA,3.0),deltaA, s.def1c+RandN(errs,3.0), s.def2c+RandN(errs,3.0), errs, flago, sea, int(fi) );
     if(sea<0 or sea>9 or flago<0.5 or filt<0 or filt>int(Nl-1) or fi<0 or fi>5 or ndw>=coun){
     cout<<"Error(LSST):  filter:  "<<fi<<"\t counter_filter:  "<<filt<<endl;
     cout<<" season:  "<<sea<<"\t flago:  "<<flago<<"\t ndw:  "<<ndw<<endl;  int yys; cin>>yys;}
     ndw+=1; } }
     
     
     ///ELT data *********************************************************
     if(datf2>0 and fel>0.0 and magni[6]>=satu[6] and magni[6]<=thre[6]){//ELT data
     errg= errELT(ls, magni[6], 0); 
     errs= errELT(ls, magni[6], 1);      
     magnio=magni[6] + RandN(errg ,3.0);
     
     chi  +=fabs( (magnio- s.magb[6])*(magnio- s.magb[6])/(errg*errg));
     chi1 +=fabs( (magnio-  magni[6])*(magnio-  magni[6])/(errg*errg));
     chi2 +=fabs( (magnio- magni0[6])*(magnio- magni0[6])/(errg*errg)); 
     sil=  RandN(errs*sqrt(2.0), 3.0) ;
     chi1a+= fabs( (def-def0)*(def-def0)/(errs*errs*2.0) );
     chi2a+= fabs( (def+sil-def)*(def+sil-def)/(errs*errs*2.0) );
     timn[ndw]=tim;
     magn[ndw]=magni[6];
     errm[ndw]=errg;
     soux[ndw]=s.pos1;
     souy[ndw]=s.pos2;
     erra[ndw]=errs;
     tele[ndw]=1;
     deltaA=fabs(pow(10.0,-0.4*errg)-1.0)*s.Astar; 
     if(deltaA<s.errM)  s.errM=deltaA; 
     if(errs< s.errA)   s.errA=errs; 
     //s.errM += deltaA; 
     //s.errA += errs; 
     if(flagf>0) fprintf(data3,"%.4lf  %.4lf  %.6lf  %.6lf  %.4lf  %.4lf   %.4lf  %.6lf %.1f  %d  %d\n",
     tim,magnio,errg,s.Astar+RandN(deltaA,3.0),deltaA, s.def1c+RandN(errs,3.0), s.def2c+RandN(errs,3.0), errs, -flago, sea, -1);
     ndw+=1;
     if(sea<0 or sea>9 or flago<0.5 or ndw>=coun){
     cout<<"Error(ELT) sea:   "<<sea<<"\t flago:  "<<flago<<endl;
     cout<<"Error ndw:  "<<ndw<<"\t coun:  "<<coun<<endl;  int uue; cin>>uue;}} 
     }
     }//end of loop time
     if(flagf>0) {fclose(data3); fclose(magg);}
   ///************ WFIRST *********************************************************************
    dchi=fabs(chi-chi1);
    dchi2=fabs(chi1-chi2)/ndw;
    dchia=fabs(chi1a-chi2a)/ndw;
    s.ampM=sqrt(s.ampM/numc);//RMS
    s.ampA=sqrt(s.ampA/numc);//RMS
    
    //s.errM =double(s.errM/ndw); //<Error_Magnification>
    //s.errA =double(s.errA/ndw); //<Error_Astrometry>
    if(dchi>800.0 and flag_det>0 and ndw>5){
    nml+=1;
    fil1=fopen("./files/MONT1/IBH_MONT1c.dat","a+");
    fprintf(fil1,
   "%d  %.5lf   %.5lf   %.5lf  %.5lf  "///5
   "%d  %.5lf   %.5lf   %.5lf  "///9
   "%d   %d  %.5lf   %.5lf  %.4lf  %.4lf  %.4lf %.4lf  %.4lf   %.4lf  "   //19
   "%.5lf  %.5lf  %.5lf  %.5lf   %.1lf  %.1lf   %.4lf  %.4lf   "  //27
   "%.6lf  %.5lf  %.4lf  %.6lf   %.5lf  %.5lf  %.5lf   %.7lf " ///35
   "%d  %d  %.1lf  %d  %.7lf  %.7lf   %.5lf  %.7lf  %.7lf   %.7lf  %.7lf  %.7lf   %.7lf  %.7lf  %.5lf  %.5lf\n", //49
   counter,s.lat, s.lon, s.right, s.decli, //5
   l.struc, l.Ml, l.Dl, l.vl, //9
   s.struc, s.cl, s.Ds, s.logT, s.logL, s.vs, s.Mab[2], s.Mab[6], s.Map[2], s.Map[6],//19
   s.mbs[0],s.mbs[1], s.fb[0], s.fb[1], s.Nblend[2], s.Nblend[6], s.Ai[2], s.Ai[6], //27
   l.tE, l.RE/AU, l.t0, l.mul, l.Vt, l.u0, s.opt*1.0e6,  l.tetE,//35
   flagf, flag_det, dchi, ndw, s.mus1, s.mus2, s.xi, s.mul1, s.mul2, l.piE, s.ampM, s.errM, s.ampA, s.errA, dchi2, dchia);//51
   fclose(fil1);
   cout<<"************* End of saving in the file ***********************"<<endl;
///####################################################################################
   
   Delta1[5]=0.05;///for mbs   
   Delta1[3]=0.1; //for xi in radian 
   if(l.t0>10.0) Delta1[0]=+10.0; 
   else          Delta1[0]=float(l.t0-2.0); 
   if(l.u0>0.1)  Delta1[1]=0.1; 
   else          Delta1[1]=float(l.u0*0.4); 
   if(l.tE>5.0)  Delta1[2]=5.0;          
   else          Delta1[2]=float(l.tE*0.4); 
   double bb= 0.5*(s.fb[0]+s.fb[1]);
   if(bb>0.1 and bb<0.9)   Delta1[4]=0.1;
   else if (bb<=0.1)       Delta1[4]=float(bb*0.5);
   else if (bb>=0.9)       Delta1[4]=float(1.0-bb)*0.5;
   Delta1[6]= float(0.25*l.piE);
   //if(l.piE>0.01) Delta1[6]=0.01;
   //else           Delta1[6]=float(l.piE*0.3); 
///######################################################3
   if(l.tetE>0.2) Delta2[0]=0.2;
   else           Delta2[0]=l.tetE*0.4;   
   Delta2[1]= s.mus1*0.25;  
   Delta2[2]= s.mus2*0.25; 
   Delta2[3]= float(0.25*l.piE); 
   //if(l.piE>0.01) Delta2[3]=0.01;
   //else           Delta2[3]=float(l.piE*0.3);  
///######################################################3   
   
   for(int j=0; j<Nx; ++j){
   for(int k=0; k<Nx; ++k){co.FishM[j][k]=0.0;}}
   for(int j=0; j<Ny; ++j){
   for(int k=0; k<Ny; ++k){co.FishA[j][k]=0.0;}}
 
///######################################################  
  
   for(int i=0; i<ndw; ++i){
   tt=int(tele[i]);
   
   for(int j=0; j<Nx; ++j){
   
   for(int h=0; h<2; ++h){
   if(j==0)  l.t0 += double(Delta1[0]*sig[h]);
   if(j==1)  l.u0 += double(Delta1[1]*sig[h]);
   if(j==2)  l.tE += double(Delta1[2]*sig[h]);
   if(j==3)  s.xi += double(Delta1[3]*sig[h]);
   if(j==4)  s.fb[tt] += double(Delta1[4]*sig[h]);
   if(j==5)  s.mbs[tt]+= double(Delta1[5]*sig[h]);
   if(j==6)  l.piE+= double(Delta1[6]*sig[h]);
      
   lightcurve(s,l,as,timn[i]);
   s.Astar= (s.ut*s.ut+2.0)/sqrt(s.ut*s.ut*(s.ut*s.ut+4.0));
   magw= s.mbs[tt] - 2.5*log10(s.Astar*s.fb[tt] + 1.0 - s.fb[tt]); 
   derm1[h]= float(magw-magn[i])/(Delta1[j]*sig[h]);
   
   if(j==0) l.t0 -= double(Delta1[0]*sig[h]);  
   if(j==1) l.u0 -= double(Delta1[1]*sig[h]);  
   if(j==2) l.tE -= double(Delta1[2]*sig[h]);  
   if(j==3) s.xi -= double(Delta1[3]*sig[h]);  
   if(j==4) s.fb[tt]  -= double(Delta1[4]*sig[h]);   
   if(j==5) s.mbs[tt] -= double(Delta1[5]*sig[h]); 
   if(j==6) l.piE-= double(Delta1[6]*sig[h]);}
   
   derm1f= double(derm1[0]+derm1[1])*0.5; 
   
   
   for(int k=0; k<=j; ++k){
   
   for(int h=0; h<2; ++h){
   if(k==0) l.t0+= double(Delta1[0]*sig[h]);
   if(k==1) l.u0+= double(Delta1[1]*sig[h]);
   if(k==2) l.tE+= double(Delta1[2]*sig[h]);
   if(k==3) s.xi+= double(Delta1[3]*sig[h]);
   if(k==4) s.fb[tt]+= double(Delta1[4]*sig[h]);
   if(k==5) s.mbs[tt]+=double(Delta1[5]*sig[h]);
   if(k==6) l.piE+=double(Delta1[6]*sig[h]);

   lightcurve(s,l,as,timn[i]);
   s.Astar= (s.ut*s.ut+2.0)/sqrt(s.ut*s.ut*(s.ut*s.ut+4.0));
   magw= s.mbs[tt] - 2.5*log10(s.Astar*s.fb[tt] + 1.0 - s.fb[tt]); 
   derm2[h]= float(magw-magn[i])/(Delta1[k]*sig[h]);
   
   if(k==0) l.t0 -= double(Delta1[0]*sig[h]);  
   if(k==1) l.u0 -= double(Delta1[1]*sig[h]);  
   if(k==2) l.tE -= double(Delta1[2]*sig[h]);  
   if(k==3) s.xi -= double(Delta1[3]*sig[h]);  
   if(k==4) s.fb[tt] -= double(Delta1[4]*sig[h]);   
   if(k==5) s.mbs[tt]-= double(Delta1[5]*sig[h]); 
   if(k==6) l.piE-= double(Delta1[6]*sig[h]);}
   
   derm2f= double(derm2[0]+derm2[1])*0.5;
   
   co.FishM[j][k] += derm1f*derm2f/(errm[i]*errm[i]);}
   
   }//end of for J
   }//end of data for
   
   for(int j=0; j<Nx; ++j){
   for(int k=(j+1);k<Nx; ++k) co.FishM[j][k]=co.FishM[k][j]; }
 
  
    fil2=fopen("./files/MONT1/IBH_MONT2c.dat","a+");
    fprintf(fil2,"%d  %.7lf  %.7lf %.7lf  %.7lf %.7lf   %.7lf %.7lf   %.7lf  %.7lf  %.8e   %.8e %.8e  %.8e %.8e %.8e   %.8e %.8e %.8e %.8e   %.8e %.8e  %.8e %.8e %.8e   %.8e %.8e %.8e %.8e %.8e %.8e  %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
    counter,l.t0, l.u0, l.tE, s.xi, s.fb[0], s.mbs[0],s.fb[1], s.mbs[1], l.piE, 
    co.FishM[0][0], 
    co.FishM[1][0], co.FishM[1][1], 
    co.FishM[2][0], co.FishM[2][1], co.FishM[2][2],
    co.FishM[3][0], co.FishM[3][1], co.FishM[3][2],co.FishM[3][3],
    co.FishM[4][0], co.FishM[4][1], co.FishM[4][2],co.FishM[4][3],co.FishM[4][4],
    co.FishM[5][0], co.FishM[5][1], co.FishM[5][2],co.FishM[5][3],co.FishM[5][4],co.FishM[5][5], 
    co.FishM[6][0], co.FishM[6][1], co.FishM[6][2],co.FishM[6][3],co.FishM[6][4],co.FishM[6][5], co.FishM[6][6]); //38
    fclose(fil2); 
   
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   

   
   for(int i=0; i<ndw; ++i){
   
   for (int j=0; j<Ny; ++j) {
   
   for(int h=0; h<2; ++h){
   if(j==0) l.tetE += double(Delta2[0]*sig[h]);
   if(j==1) s.mus1 += double(Delta2[1]*sig[h]);
   if(j==2) s.mus2 += double(Delta2[2]*sig[h]);
   if(j==3) l.piE  += double(Delta2[3]*sig[h]); 
   lightcurve(s,l,as,timn[i]);
   dera1[h]= float(s.pos1-soux[i])/(Delta2[j]*sig[h]);  
   derb1[h]= float(s.pos2-souy[i])/(Delta2[j]*sig[h]);
   if(j==0) l.tetE -= double(Delta2[0]*sig[h]);
   if(j==1) s.mus1 -= double(Delta2[1]*sig[h]); 
   if(j==2) s.mus2 -= double(Delta2[2]*sig[h]);
   if(j==3)  l.piE -= double(Delta2[3]*sig[h]);}
   dera1f=(dera1[0]+ dera1[1])*0.5;  
   derb1f=(derb1[0]+ derb1[1])*0.5;  
  
   for(int k=0; k<=j; ++k){   
   
   for(int h=0; h<2; ++h){
   if(k==0) l.tetE += double(Delta2[0]*sig[h]);
   if(k==1) s.mus1 += double(Delta2[1]*sig[h]);
   if(k==2) s.mus2 += double(Delta2[2]*sig[h]);
   if(k==3) l.piE  += double(Delta2[3]*sig[h]); 
   lightcurve(s,l,as,timn[i]);
   dera2[h]= float(s.pos1-soux[i])/(Delta2[k]*sig[h]);  
   derb2[h]= float(s.pos2-souy[i])/(Delta2[k]*sig[h]);  
   if(k==0) l.tetE -= double(Delta2[0]*sig[h]);
   if(k==1) s.mus1 -= double(Delta2[1]*sig[h]); 
   if(k==2) s.mus2 -= double(Delta2[2]*sig[h]);
   if(k==3)  l.piE -= double(Delta2[3]*sig[h]); }
   dera2f=(dera2[0]+ dera2[1])*0.5;  
   derb2f=(derb2[0]+ derb2[1])*0.5;  
   
   co.FishA[j][k] += (dera1f*dera2f + derb1f*derb2f)/(erra[i]*erra[i]*2.0);}
  
   }//end of loop J
   }//end of loop data
   
   for(int j=0; j<Ny; ++j){
   for(int k=(j+1);k<Ny; ++k) co.FishA[j][k]= co.FishA[k][j];}

  

   fil3=fopen("./files/MONT1/IBH_MONT3c.dat","a+");
   fprintf(fil3, "%d  %.7lf    %.7lf    %.7lf   %.9e    %.9e    %.9e   %.9e   %.9e   %.9e    %.9e   %.9e   %.9e   %.9e\n",
   counter,l.tetE, s.mus1, s.mus2,
   co.FishA[0][0], 
   co.FishA[1][0], co.FishA[1][1], 
   co.FishA[2][0], co.FishA[2][1], co.FishA[2][2],
   co.FishA[3][0], co.FishA[3][1], co.FishA[3][2], co.FishA[3][3]);
   fclose(fil3); 
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
/*
   cout<<"Inverse Matrix Magnification ****************************"<<endl;
   for(int i=0;i<Nx; ++i){
   for(int j=0;j<Nx; ++j)  cout << co.InvM[i][j] << "     ";
   cout << endl;}
   cout<<"Inverse Matrix Astrometry ***************** **************"<<endl;
   for(int i=0;i<Ny; ++i){
   for(int j=0;j<Ny; ++j)    cout << co.InvA[i][j] << "     ";
   cout << endl;}
   cout<<"************************************************"<<endl;*/
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   
   
   
   if(int(nml)%1==0){
   cout<<"============================================================="<<endl;
   cout<<"counter:     "<<counter<<"\t number_ML_events_detected:   "<<nml<<endl;
   cout<<"lat:  "<<s.lat<<"\t lon:  "<<s.lon<<endl;
   cout<<"********************** SOURCE **************************"<<endl;
   cout<<"Ds:  "<<s.Ds<<"\t nums:  "<<s.nums<<"\t strucs:  "<<s.struc<<endl;
   cout<<"cl:  "<<s.cl<<"\t logT:  "<<s.logT<<endl;
   cout<<"s.logL:  "<<s.logL<<"\t cl:  "<<s.cl<<endl;
   cout<<"vs:  "<<s.vs<<"\t Mag_I:  "<<s.Mab[1]<<"\t Mag_W149:  "<<s.Mab[4]<<endl;
   cout<<"********************** LENS **************************"<<endl;
   cout<<"Dl:  "<<l.Dl<<"\t numl:  "<<l.numl<<"\t strucl:  "<<l.struc<<endl;
   cout<<"mass:  "<<l.Ml<<"vl:  "<<l.vl<<endl;
   cout<<"*********************** LENSING ***********************"<<endl;
   cout<<"tE:  "<<l.tE<<"\t RE(AU):  "<<l.RE/AU<<"\t t0:  "<<l.t0<<endl;
   cout<<"Vt:  "<<l.Vt<<"\t mul(mas/days):  "<<l.mul<<"\t u0:  "<<l.u0<<endl;
   cout<<"************ WFIRST & OGLE ***************************"<<endl;
   cout<<"flag_det:  "<<flag_det<<"\t dchi:  "<<dchi<<"\t **************ndw:  "<<ndw<<endl;
   cout<<"t0:  "<<l.t0<<"\t u0:  "<<l.u0<<"\t tE:  "<<l.tE<<endl;
   cout<<"ksi:  "<<s.xi<<"\t fb:  "<<s.fb<<"\t mbs:  "<<s.mbs<<endl;
   cout<<"piE:   "<<l.piE<<"\t l.tetE:  "<<l.tetE <<endl;
   cout<<"mus1:  "<<s.mus1<<"\t mus2:  "<<s.mus2<<endl;
   cout<<"==============================================================="<<endl; }  
   }//if detection_wfirst>0
   }//if magnitude of baseline 
   }while(nml<700);//loop icon
   }//end of loop icon 
   delete [] magn, timn, soux, souy, errm, erra, tele; 
   return(0);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Light curves and Astrometry                //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void lightcurve(source & s, lens & l, astromet & as, double timh)
{ 
   double dvex, dvey, Ve_x, tt;
   double pis= double(1.0/s.Ds)/l.tetE;
   
   double int1=0.0, int2=0.0; 
   for (int ig=0; ig<2; ++ig){
   if(ig==0) tt=timh; 
   if(ig==1) tt=0.0;    
   dvex= +as.vearth* sin(as.omegae* tt + M_PI/2.0) / as.omegae;//km/s
   dvey= -as.vearth* cos(as.omegae* tt + M_PI/2.0) / as.omegae;//km/s
   as.Ve_n1= cos(as.tetp)*dvex*sin(s.deltao) - dvey*cos(s.deltao);
   Ve_x=    -cos(as.tetp)*dvex*cos(s.deltao) - dvey*sin(s.deltao);
   as.Ve_n2=-sin(s.FI)*Ve_x + cos(s.FI)*sin(as.tetp)*dvex;
   if(ig==0) {int1 = as.Ve_n1; int2 = as.Ve_n2; }
   if(ig==1) {int1-= as.Ve_n1; int2-= as.Ve_n2; }}
   as.ue_n1= int1; 
   as.ue_n2= int2; 
      
   s.ux=-l.u0*sin(s.xi) + (timh-l.t0)*cos(s.xi)/l.tE  + l.piE * as.ue_n1;
   s.uy= l.u0*cos(s.xi) + (timh-l.t0)*sin(s.xi)/l.tE  + l.piE * as.ue_n2;
   
   s.ut0=sqrt( (s.ux-l.piE*as.ue_n1)*(s.ux-l.piE*as.ue_n1) + (s.uy-l.piE*as.ue_n2)*(s.uy-l.piE*as.ue_n2) );
   s.ut= sqrt( s.ux*s.ux+ s.uy*s.uy);  
   
   if(s.ut==0.0)   s.ut=1.0e-50; 
   if(s.ut0==0.0)  s.ut0=1.0e-50; 
   
   
   s.def1a=  (s.ux-l.piE*as.ue_n1)*l.tetE/(s.ut0*s.ut0 + 2.0);//x-deflection, without parallax
   s.def2a=  (s.uy-l.piE*as.ue_n2)*l.tetE/(s.ut0*s.ut0 + 2.0);//y-deflection, without parallax 
   
   s.def1b=  (s.ux-l.piE*as.ue_n1)*l.tetE/(s.ut*s.ut + 2.0); //x-deflection, first sentence
   s.def2b=  (s.uy-l.piE*as.ue_n2)*l.tetE/(s.ut*s.ut + 2.0); //y-deflection, first sentence 
   
   s.def1c=  s.ux*l.tetE/(s.ut*s.ut + 2.0); //x-deflection
   s.def2c=  s.uy*l.tetE/(s.ut*s.ut + 2.0); //y-deflection 
   
   
   s.pos10= -l.u0*l.tetE*sin(s.xi) + s.mus1*(timh-l.t0) + s.def1a;//x-source trajectory, without parallax  
   s.pos20= +l.u0*l.tetE*cos(s.xi) + s.mus2*(timh-l.t0) + s.def2a;//y-source trajectory, without parallax  
  
   s.pos1= -l.u0*l.tetE*sin(s.xi) + s.mus1*(timh-l.t0) - as.ue_n1*pis  +s.def1c;//x-source trajectory
   s.pos2= +l.u0*l.tetE*cos(s.xi) + s.mus2*(timh-l.t0) - as.ue_n2*pis  +s.def2c;//y-source trajectory        
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double RandR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double errlsst2(double maga, int fi, double airm) 
{  
        double seei, mfive, x, Delta1;    
        
        seei=double(seeing[fi]+RandN(0.2,1));
        if(seei<=0.4) seei=0.4; 
        
        mfive=Cm[fi]+0.0+0.5*(msky[fi]-21.0)+2.5*log10(0.7/seei)+1.25*log10(texp/30.0)-km[fi]*(airm-1.0) ;
        x=pow(10.0 , 0.4*(maga-mfive));
        Delta1=sqrt(fabs((0.04-gama[fi])*x+gama[fi]*x*x)); 
        if(Delta1<0.0001)   Delta1=0.0001; 
        
        if(mfive<0.0 or mfive>40.0  or Delta1<0.00001 or Delta1>0.9) {
        cout<<"Error(LSST_phot_acuracy):  mfive:    "<<mfive<<"\t Delta1:  "<<Delta1<<"\t maga:  "<<maga<<endl;
        cout<<"msky[fi]:  "<<msky[fi]<<"\t seei:  "<<seei<<"\t airm:   "<<airm<<"\t  fi: "<<fi<<endl;
        int uue; cin>>uue; 
        }
        return( sqrt(delta2*delta2+Delta1*Delta1) ); 
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double errlsst(lsst & ls, double ghadr, int flag)
{
     double error=0.0, shib=0.0; 
    /*    
     if(flag==0){///calculating photometry errors LSST
     
     if(ghadr<ls.magl[0] or  ghadr==ls.magl[0])          error= ls.errl[0];
     
     else if(ghadr>ls.magl[Nm-1] or ghadr==ls.magl[Nm-1]) {
     shib=(ls.errl[Nm-1]-ls.errl[Nm-2])/(ls.magl[Nm-1]-ls.magl[Nm-2]); 
     error=ls.errl[Nm-1]+shib*(ghadr-ls.magl[Nm-1]);    }
     
     else{
     for(int i=1; i<Nm; ++i){
     if(double((ghadr-ls.magl[i])*(ghadr-ls.magl[i-1]))<0.0 or  ghadr==ls.magl[i-1]){
     shib=(ls.errl[i]-ls.errl[i-1])/(ls.magl[i]-ls.magl[i-1]); 
     error=ls.errl[i-1]+shib*(ghadr-ls.magl[i-1]); 
     break;}}}} */
     
     //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     
     if(flag==1){///calculating astrometry error LSST
     
     if(ghadr<ls.maga[0] or  ghadr==ls.maga[0])            error= ls.erra[0]; 
     
     else if(ghadr>ls.maga[Na-1] or ghadr==ls.maga[Na-1]){
     shib=(ls.erra[Na-1]-ls.erra[Na-2])/(ls.maga[Na-1]-ls.maga[Na-2]); 
     error=ls.erra[Na-1]+shib*(ghadr-ls.maga[Na-1]);   }   
       
     else{
     for(int i=1; i<Na; ++i){
     if(double((ghadr-ls.maga[i])*(ghadr-ls.maga[i-1]))<0.0 or ghadr==ls.maga[i-1]){
     shib=(ls.erra[i]-ls.erra[i-1])/(ls.maga[i]-ls.maga[i-1]);
     error=ls.erra[i-1] +  shib*(ghadr-ls.maga[i-1]); 
     break;}}}}
     
     return(error);     
}     
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double errELT(lsst & ls, double ghadr, int flag)
{
     double error=0.0, shib=0.0; 
     
     if(flag==0){//photometry-ELT
     if(ghadr<ls.mage[0] or ghadr==ls.mage[0])       error= ls.ere1[0];
     
     else if(ghadr>ls.mage[Nel-1] or ghadr==ls.mage[Nel-1]){
     shib=double(ls.ere1[Nel-1]-ls.ere1[Nel-2])/(ls.mage[Nel-1]-ls.mage[Nel-2]); 
     error=ls.ere1[Nel-1]+shib*(ghadr-ls.mage[Nel-1]);     }
     
     else{
     for(int i=1; i<Nel; ++i){
     if(double((ghadr-ls.mage[i])*(ghadr-ls.mage[i-1]))<0.0 or  ghadr==ls.mage[i-1]){
     shib=double(ls.ere1[i]-ls.ere1[i-1])/(ls.mage[i]-ls.mage[i-1]); 
     error=ls.ere1[i-1]+shib*(ghadr-ls.mage[i-1]); 
     break;}}} 
     }
     
     
     if(flag==1){//astrometry-ELT
     if(ghadr<ls.mage[0] or  ghadr==ls.mage[0])      error= ls.ere2[0];
      
     else if(ghadr>ls.mage[Nel-1] or ghadr==ls.mage[Nel-1]){
     shib=double(ls.ere2[Nel-1]-ls.ere2[Nel-2])/(ls.mage[Nel-1]-ls.mage[Nel-2]); 
     error=ls.ere2[Nel-1]+shib*(ghadr-ls.mage[Nel-1]); }   

     else{
     for(int i=1; i<Nel; ++i){
     if(double((ghadr-ls.mage[i])*(ghadr-ls.mage[i-1]))<0.0 or ghadr==ls.mage[i-1]){
     shib=double(ls.ere2[i]-ls.ere2[i-1])/(ls.mage[i]-ls.mage[i-1]);
     error=ls.ere2[i-1] +  shib*(ghadr-ls.mage[i-1]); 
     break;}}}
     }
 
     return(error);     
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///==============================================================//6
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*Msun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opt=0.0;
    s.od_dlmc=s.od_hlmc=s.od_blmc=0.0; 
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;
    s.od_dlmc +=  s.rho_dlmc[k]*x*(1.0-x)*dx*CC;
    s.od_blmc +=  s.rho_blmc[k]*x*(1.0-x)*dx*CC; 
    s.od_hlmc +=  s.rho_hlmc[k]*x*(1.0-x)*dx*CC;  }
    s.opt= fabs(s.od_disk+s.od_ThD+s.od_bulge+s.od_halo+ s.od_dlmc + s.od_blmc + s.od_hlmc );///total
}
///==============================================================//
///                                                              //
///                  func_source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD  &  cm)
{
    int num,struc,nums,yye;
    double rho, rf;
    double Ds, Ai[M+1], Av;
    double Mab[M+1], Map[M+1];
    double maxnb=0.0;


    for(int i=0; i<int(M+1); ++i){
    s.Nblend[i]=s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0);
    s.Nblend[i]=s.Nblend[i]+RandN(sqrt(s.Nblend[i]),1.5);
    if(s.Nblend[i]<=1.0) s.Nblend[i]=1.0;
    if(s.Nblend[i]>maxnb) maxnb=s.Nblend[i];}
    
    for(int i=0; i<int(M+1); ++i){
    s.Fluxb[i]=0.0;
    Ai[i]=0.0;  Mab[i]=0.0; Map[i]=0.0; 
    s.magb[i]=0.0; s.Ai[i]=0.0; s.Map[i]=0.0; s.Mab[i]=0.0;}



    for(int k=1; k<=int(maxnb); ++k){


    do{
    num=int((Num-15.00)*(double)rand()/(double)(RAND_MAX+1.) +10.00);
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    Ds=(double)(num*step);
    }while(rho>s.Rostari[num] or Ds<0.1 or Ds>MaxD);///distance larger than 50.
    
    
    if(Ds>MaxD or Ds<0.1){cout<<"ERROR (1): Ds: "<<Ds<<"\t MaxD: "<<MaxD<<"\t step: "<<step<<"\t num: "<<num<<endl; cin>>yye;}
    nums=num;
    if(k==1){s.Ds=Ds;  s.nums=nums; }



    rf= RandR(0.0, s.Rostar0[nums]); 
     if (rf<= s.rho_disk[nums])                    struc=0;///thin disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums])) struc=3;///halo
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums]+s.rho_hlmc[nums]))struc=4;///halo_lmc
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums]+s.rho_hlmc[nums]+s.rho_dlmc[nums]))struc=5;//disk_lmc
else struc=6;///bulge_lmc
     
      if(k==1) s.struc=struc; 

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(struc==0 or struc==5){//thin disk Milky Way and LMC
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-2.0));
    for(int i=0; i<M; ++i)   Mab[i]=cm.Mab_d[i][num]; 
    Mab[M]=cm.Mk_d[num];
    if(k==1){
    s.logT=cm.logT_d[num];
    s.logL= cm.logL_d[num];
    s.cl= cm.cl_d[num];}
    if(int(cm.cl_d[num])>6 or cm.logT_d[num]<0.0){
    cout<<"Error(thin disk) LogT: "<<cm.logT_d[num]<<"\t cl_d: "<<cm.cl_d[num]<<"\t counter: "<<num<<endl; cin>>yye;}}



    if(struc==1 or struc==6){///bulge  MilkyWay and LMC
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-2.0));
    for(int i=0; i<M; ++i) Mab[i]=cm.Mab_b[i][num]; 
    Mab[M]=cm.Mk_b[num];
    if(k==1){
    s.logT=cm.logT_b[num];
    s.logL=cm.logL_b[num];
    s.cl= cm.cl_b[num];}
    if(int(cm.cl_b[num])>6 or cm.logT_b[num]<0.0){
    cout<<"Error(bulge) LogT_b: "<<cm.logT_b[num]<<"\t cl_b: "<<cm.cl_b[num]<<"\t counter: "<<num<<endl; cin>>yye;}}




    if(struc==2){///thick disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-2.0));
    for(int i=0; i<M; ++i) Mab[i]=cm.Mab_t[i][num]; 
    Mab[M]=cm.Mk_t[num];
    if(k==1){
    s.logT=cm.logT_t[num];
    s.logL=cm.logL_t[num];
    s.cl= cm.cl_t[num];}
    if(int(cm.cl_t[num])>6 or cm.logT_t[num]<0.0){
    cout<<"Error(Thick disk) LogT_b: "<<cm.logT_t[num]<<"\t cl_b: "<<cm.cl_t[num]<<"\t counter: "<<num<<endl; cin>>yye;}}



    if(struc==3 or struc==4){///stellar halo
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-2.0));
    for(int i=0; i<M; ++i) Mab[i]=cm.Mab_h[i][num]; 
    Mab[M]=cm.Mk_h[num];
    if(k==1){
    s.logT=cm.logT_h[num];
    s.logL=cm.logL_h[num];
    s.cl= cm.cl_h[num]; }
    if(int(cm.cl_h[num])>6 or cm.logT_h[num]<0.0){
    cout<<"Error(Thick disk) LogT_b: "<<cm.logT_h[num]<<"\t cl_b: "<<cm.cl_h[num]<<"\t counter: "<<num<<endl; cin>>yye;}}
    
    
    
    Av=RandR(1.6, 2.1);///because we simulate only center of LMC    
    if(Av<0.0)    Av=0.0;
    for(int i=0;  i<int(M+1); ++i){
    Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i], 1.5); //extinction in other bands
    if(Ai[i]<0.0) Ai[i]=0.0;
    Map[i]=Mab[i] + 5.0*log10(Ds*100.0) + Ai[i];
    if(s.Nblend[i]>=k){ s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i]));}}
    if(k==1){
    for(int i=0; i<int(M+1);  ++i){s.Ai[i]=Ai[i]; s.Map[i]=Map[i]; s.Mab[i]= Mab[i];}}
    }///loop over the stars




    for(int i=0; i<int(M+1); ++i){
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    s.magb[i]= -2.5*log10( s.Fluxb[i] );    
    if(int(s.Nblend[i])<1.0 or (s.Nblend[i]==1.0 and s.blend[i]<1.0) or s.blend[i]>1.0 or s.blend[i]<0.0 or
    s.Fluxb[i]<-0.009543 or s.Ai[i]<0.0){
    cout<<"BIGG ERRROR nsbl: "<<s.Nblend[i]<<"\t Nblend[i]: "<<s.Nblend[i]<<"\t belnd: "<<s.blend[i]<<endl;
    cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<"\t No. filter:  "<<i<<"\t ext:  "<<s.Ai[i]<<endl; cin>>yye;} }
    
    
    s.fb[0]=s.blend[2];//r-LSST 
    s.fb[1]=s.blend[6];//K-band 
    s.mbs[0]=s.magb[2];//r-LSST      
    s.mbs[1]=s.magb[6];//K-band 
    
    //cout<<"************** End of func_source  ****************"<<endl;
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s)
{
    double f,test;
    double rholens[s.nums+2];
    l.rhomaxl=0.0;
    int il, yye;
    double mmin=2.0;
    double mmax=50.0;

    for(int k=1;k<=s.nums;++k){
    rholens[k]=0.0;
    l.Dl=k*step;
    l.xls=l.Dl/s.Ds;
    if(l.Dl>s.Ds){cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;  cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}


    do{
    l.numl=(int)((double)rand()*1000./((double)(RAND_MAX+1.)*1000.)*(s.nums-2.0)+1.0);
    test  =((double)rand()/(double)(RAND_MAX+1.)*l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){
    cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;  int ue; cin>>ue;}
    }while(test>rholens[l.numl]);
    l.Dl=double(l.numl*step);



double randflag=RandR(0.0,s.Rostar0[l.numl]);
if (randflag<=fabs(s.rho_disk[l.numl]) ) l.struc=0;//disk
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1;//bulge
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2;//thick
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]) ) l.struc=3;//halo
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]+s.rho_hlmc[l.numl]))l.struc=4;//halo_lmc
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]+s.rho_hlmc[l.numl]+s.rho_dlmc[l.numl])) l.struc=5;//disk_lmc
else if (randflag<=fabs( s.Rostar0[l.numl])) l.struc=6;//bar_lmc
else {  cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}

   
    
    
    l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(mmax - mmin)) + mmin;
    l.xls=l.Dl/s.Ds;
    l.RE=sqrt(4.0*G*l.Ml*Msun*s.Ds*KP)/velocity;
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    vrel(s,l);
    l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///in day
    l.mul=l.Vt*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days
    l.tetE= double(l.RE/AU/l.Dl);///marcs
    
    l.u0=RandR(0.0,1.0);
    l.t0=RandR(2.0 , Tobs-2.0 );// [2 days, 10years-2 days] in days
    l.pirel=double(1.0/l.Dl- 1.0/s.Ds);     

    if(l.tE<=0.0 or l.Dl>s.Ds or l.Vt<=0.0 or l.Ml<0.0 or l.t0>Tobs or l.t0<0.0 or l.Ml>mmax or l.Ml<mmin or l.pirel<0.0 or l.pirel==0.0){
    cout<<"Error (fun_lens)   Vt:  "<<l.Vt<<endl;
    cout<<"RE: "<<l.RE/AU<<"\t xls:  "<<l.xls<<"\t tE: "<<l.tE<<endl;
    cout<<"Ml:  "<<l.Ml<<"\t Dl:  "<<l.Dl<<"\t Ds:  "<<s.Ds<<"\t pirel:  "<<l.pirel<<endl;
    cout<<"numl:  "<<l.numl<<"\t Vt:  "<<l.Vt<<"\t mul:  "<<l.mul<<endl;  cin>>yye;}
    //cout<<"************** End of func_Lens  ****************"<<endl;
}
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double ylsst(yfilter & yf, double metal, double age, double magz)
{
    int h, g, k1,k2; 
    double magy;  
    //cout<<"metal_0:  "<<yf.Metal[0]<<"\t metal_1:  "<<yf.Metal[1]<<"\t metal_3:   "<<yf.Metal[3]<<endl;

    h=-1; 
    if(metal<=yf.Metal[0])          h=0; 
    else if(metal>=yf.Metal[met-1]) h=met-1;
    else{ 
    for(int i=0; i<int(met-1); ++i){
    if(float((metal-yf.Metal[i])*(metal-yf.Metal[i+1]))<0.0 or metal==yf.Metal[i]) {h=i; break;}}}
       
       
    g=-1;    
    k1=int(yf.count[h]);   
    k2=int(yf.count[h]+yf.number[h]);
    if(k1<0)           k1=0;
    if(k2>int(YZ-1))   k2=int(YZ-1);
    
    if(age<=yf.Age[k1])          g=int(k1); 
    else if(age>=yf.Age[k2-1])   g=int(k2-1); 
    
    else{
    for(int k=k1;  k<int(k2-1);  ++k){
    if(float((age-yf.Age[k])*(age-yf.Age[k+1]))<0.0  or  age==yf.Age[k]){g=k; break;}}}
    magy = double(yf.B[g] + yf.M[g]*magz);   
  
    if(g<0 or g>YZ or g<k1 or g>k2 or h<0 or h>met or k1>k2 or k2>YZ or fabs(age-yf.Age[g])>5.0 ){
    cout<<"ERROR2 g:"<<g<<"\t h:  "<<h<<"\t age:"<<age<<"\t metal: "<<metal<<"\t Age[g]:  "<<yf.Age[g]<<endl;
    cout<<"k1: "<<k1<<"\t k2:  "<<k2<<"\t metal[h]:  "<<yf.Metal[h]<<endl;  int uui;  cin>>uui;}
    return( magy ); 
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm, yfilter & yf)
{
//  mass  logT    Age  logL logG  Z  u_LSST  g_LSST  r_LSST i_LSST z_LSST  B  V  R    I   J    H    K    CL    Type//1402/04/06
//  0      1       2     3    4   5   6       7       8      9      10     11 12  13  14  15   16   17   18     19
    int yye, age;     
    double mass,gravity,metal,B,V,R,I,J, H,type; 
    char filename[40];
    FILE *fp2;
   
////=================================== THIN DISK ============================== 
    int j=0; 
    sprintf(filename,"./files/CMDBESANCON14020406/%c%c%c%c%c%c%d.dat",'C','M','D','T','i','b',3);///changed 1402/04/06
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTi.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %lf %d %lf\n",
    &mass,&cm.logT_d[j],&age,&cm.logL_d[j],&gravity,&metal,
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.Mab_d[3][j],&cm.Mab_d[4][j],&B,&V,&R,&I,&J,&H,&cm.Mk_d[j],&cm.cl_d[j],&type);
    //cout<<"metal:   "<<metal<<"\t age:  "<<age<<"\t z_LSST:   "<<cm.Mab_d[4][j]<<endl;
    cm.Mab_d[5][j]= ylsst(yf, metal, float(age*1.00), cm.Mab_d[4][j]);
    if(mass<0.0 or cm.logT_d[j]<0.0 or cm.Mab_d[2][j]>20.0 or metal>0.12 or age>10.0 or cm.cl_d[j]>7 or type>9.0){
    cout<<"ERROR(cmd thin_disk) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N1){cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}




////=================================== BULGE ================================== 
    j=0; 
    sprintf(filename,"./files/CMDBESANCON14020406/%c%c%c%c%c%d.dat",'C','M','D','B','b',3);
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDB.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %lf %d %lf\n",
    &mass,&cm.logT_b[j],&age,&cm.logL_b[j],&gravity,&metal,
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.Mab_b[3][j],&cm.Mab_b[4][j],&B,&V,&R,&I,&J,&H,&cm.Mk_b[j],&cm.cl_b[j],&type);
    cm.Mab_b[5][j]= ylsst( yf, metal, float(age*1.00), cm.Mab_b[4][j]);
    if(mass<0.0 or cm.logT_b[j]<0.0 or cm.Mab_b[2][j]>18.0 or metal>2.0 or age>10 or cm.cl_b[j]>7 or type>9.0){
    cout<<"ERROR(cmd bulge) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N2<<endl;  cin>>yye;}


////=================================== THICK DISK =============================
    j=0; 
    sprintf(filename,"./files/CMDBESANCON14020406/%c%c%c%c%c%c%d.dat",'C','M','D','T','k','b',3);
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTk.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %lf %d %lf\n",
    &mass,&cm.logT_t[j],&age,&cm.logL_t[j],&gravity,&metal,
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.Mab_t[3][j],&cm.Mab_t[4][j],&B,&V,&R,&I,&J,&H,&cm.Mk_t[j],&cm.cl_t[j],&type);
    cm.Mab_t[5][j]= ylsst( yf, metal, float(age*1.00), cm.Mab_t[4][j]);
    if(mass<0.0 or cm.logT_t[j]<0.0 or cm.Mab_t[2][j]>20.0 or metal>0.03 or age>8 or cm.cl_t[j]>7 or type>9.0){
    cout<<"ERROR(cmd Thick disk) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N3<<endl;  cin>>yye;}


////=================================== STELLAR HALO =========================== 
    j=0; 
    sprintf(filename,"./files/CMDBESANCON14020406/%c%c%c%c%c%d.dat",'C','M','D','H','b',3);
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDH.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %lf %d %lf\n",
    &mass,&cm.logT_h[j],&age,&cm.logL_h[j],&gravity,&metal,
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.Mab_h[3][j],&cm.Mab_h[4][j],&B,&V,&R,&I,&J,&H,&cm.Mk_h[j],&cm.cl_h[j],&type);
    cm.Mab_h[5][j]= ylsst( yf, metal, float(age*1.00), cm.Mab_h[4][j]);
    if(mass<0.0 or cm.logT_h[j]<0.0 or cm.Mab_h[2][j]>20.0 or metal>0.02 or age>9 or cm.cl_h[j]>7 or type>9.0){
    cout<<"ERROR(cmd Halo) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N4<<endl;  cin>>yye;}

   cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
   
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s, int yy)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nnf=0.4/0.8;
   double mBarre;//stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   s.Romins=10000000000.0;
   double fd=1.0;//see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;//0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;//No limitation 
   double Rdd=2.17;//2.53;///2.17;
   double Rhh=1.33;//1.32;//1.33;
   double alm= 2.0;//2KPC,  Mancini 2004
   double rlm_0=1.76*0.01;//solar mass/ Pc^-3
   double M_dlm=2.6;//[Msun]mass of disk LMC  from Mancini 2004 paper
   double M_blm=0.15*M_dlm;//[Msun] mass of bar LMC
   double xb0, yb0, zb0;//kPC
   double frac=0.05;//fraction of halo in the form of compact objects
   double qq=0.688; 
   double Rd0=1.54;
   double xol, yol, zol, x0, y0, z0; 
   double r0,  pos1, inc1, Rlm; 
   
   char filename[40];
   FILE *fil1;
   FILE *fill;
   
        
   int flagf=0;
   if(yy<10){
   flagf=1; 
   fil1=fopen("./files/density/right.txt", "a+"); 
   sprintf(filename,"./files/density/%c%.2lf%c%.2lf.dat",'D',s.lat,'_',s.lon);
   fill=fopen(filename,"w");}


for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;
   s.rho_hlmc[i]=s.rho_dlmc[i]=s.rho_blmc[i]=0.0; 


   x=double(i*step);
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = Dsun- x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);
   double rsun= sqrt(zb*zb+ yb*yb+ xb*xb); 


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///Msun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nnf)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*exp(nnf)*exp(-fabs(zb)/0.8)/(1.0+0.5*nnf);///Msun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=frac*(0.932*0.00001/867.067)*pow(0.5/Dsun,-2.44);
   else            s.rho_halo[i]=frac*(0.932*0.00001/867.067)*pow(rdi/Dsun,-2.44);///Msun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///Msun/pc^3
///==================================================================



///=====================  LMC density profiles ====================== 
   x0 = -x*cos(s.decli/RA)*sin(s.right/RA-RaLMC/RA);
   y0 =  x*sin(s.decli/RA)*cos(DecLMC/RA)-x*cos(s.decli/RA)*sin(DecLMC/RA)*cos(s.right/RA-RaLMC/RA);
   z0 = DLMC- x*cos(s.decli/RA)*cos(DecLMC/RA)*cos(s.right/RA-RaLMC/RA)-x*sin(s.decli/RA)*sin(DecLMC/RA);
   if(fabs(x-DLMC)<double(step*1.5)){s.xv=x0;  s.yv=y0;   s.zv=z0;}
   r0= sqrt(x0*x0+y0*y0+z0*z0); 
 
  
///=====================  LMC density HALO ======================
   if(r0<15.0) s.rho_hlmc[i]=fabs(rlm_0*frac/(1.0+r0*r0/alm/alm) );  ///Msun/PC^3
   else        s.rho_hlmc[i]=0.0;   
  
  
///=====================  LMC density Disk ======================
   pos1=(170.0-90.0)*M_PI/180.0; 
   inc1=34.7*M_PI/180.0;  
   xol= x0*cos(pos1)+ y0* sin(pos1);
   yol=-x0*sin(pos1)*cos(inc1) + y0*cos(pos1)*cos(inc1) -z0*sin(inc1);
   zol=-x0*sin(pos1)*sin(inc1) + y0*cos(pos1)*sin(inc1) +z0*cos(inc1);    
   Rlm= sqrt(xol*xol + yol*yol ); 
   double zd0=0.3;//KP Kim's paper
   double Rd1=1.8;///KPc  Kim's paper
   s.rho_dlmc[i]=fabs(M_dlm/(4.0*M_PI*zd0*Rd1*Rd1)*exp(-Rlm/Rd1)*exp(-fabs(zol/zd0))); ////kim [Msun/Pc^3]  Kim 2000
   

///=====================  LMC density Bulge ======================
   xb0= 1.2;  yb0=zb0= 0.44;
   pos1=(110.0-90.0)*M_PI/180.0; 
   inc1=0.0*M_PI/180.0;  
   xol= x0* cos(pos1) +  y0* sin(pos1);
   yol=-x0*sin(pos1)*cos(inc1) + y0*cos(pos1)*cos(inc1) -z0*sin(inc1);
   zol=-x0*sin(pos1)*sin(inc1) + y0*cos(pos1)*sin(inc1) +z0*cos(inc1);    
   s.rho_blmc[i]=M_blm*pow(2.0*M_PI,-1.5)/(xb0*yb0*zb0)*exp(-0.5*(xol*xol/xb0/xb0+yol*yol/yb0/yb0+zol*zol/zb0/zb0));///[Msun/Pc^3]  
///===============================================================================


s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]) +fabs(s.rho_hlmc[i])+fabs(s.rho_dlmc[i])+ fabs(s.rho_blmc[i]);///[Msun/pc^3]
s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Msun/deg^2]
s.Nstari[i]= binary_fraction*((s.rho_disk[i]+s.rho_dlmc[i])*fd/0.403445+s.rho_ThD[i]*fh/0.4542+(s.rho_halo[i]+s.rho_hlmc[i])*fh/0.4542+(s.rho_bulge[i]+s.rho_blmc[i])*fb/0.308571);////[Nt/pc^3] 
s.Nstari[i]=s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]


s.Nstart  +=  s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
if(s.Rostari[i]<s.Romins) s.Romins=s.Rostari[i];///source selection


  if(flagf>0)
  fprintf(fill,"%e   %e   %e   %e   %e  %e  %e   %e   %e   %e\n",
  x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.rho_dlmc[i],s.rho_blmc[i],s.rho_hlmc[i], s.Rostar0[i],s.Nstari[i]);  }

  if(flagf>0){
  fprintf(fil1,"%.5lf  %.5lf  %.5lf  %.5lf\n",s.right, s.decli, s.TET*RA, s.FI*RA);
  fclose(fill); 
  fclose(fil1);}

   //cout<<"xLMC:  "<<xLMC<<"\t yLMC:  "<<yLMC<<"\t zLMC:  "<<zLMC<<endl;
   //cout<<"LMC_distance:  "<<sqrt(xLMC*xLMC + yLMC*yLMC +  zLMC*zLMC)<<endl;
   //cout<<"xol:  "<<xol<<"\t yol:  "<<yol<<"\t  zol:  "<<zol<<endl;
   //cout<<"distance_LMC_Center:  "<<rlm<<"\t projected_distance:  "<<Rlm<<endl;
  // cout<<"Nstart [Nt/deg^2]: "<<s.Nstart<<"\t Ro_star [Mass/deg^2]: "<<s.Rostart<<endl;
   //cout<<">>>>>>>>>>>>>>>>>>>>>>>>> END OF DISK MODLE <<<<<<<<<<<<<<<<<<<<"<<endl;
   //exit(0);
}
///===========================================================================//
double RandN(double sigma, double nnd){
   double rr,f,frand;
   do{
   rr=double(((double)rand()/(double)(RAND_MAX+1.))*2.0-1.0)*sigma*nnd; ///[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/sigma/sigma);
   frand=fabs((double)rand()/((double)(RAND_MAX+1.))*1.0);
   }while(frand>f);
   return rr;
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
  if (l.Dl==0.0) l.Dl=0.00034735;
  double pi=M_PI;
  double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*l.Dl*cos(s.TET)*cos(s.FI));
  double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*s.Ds*cos(s.TET)*cos(s.FI));
  if(Rlc==0.0) Rlc=0.0000000000034346123;
  if(Rsc==0.0) Rsc=0.000000000004762654134; 
 
  double LVx, SVx;
  double SVT, SVR, SVZ, LVT, LVR, LVZ;
  double fv, testfv, test, age;
  double  VSunx, vls2, vls1;
  double betal, betas, deltal, deltas, tetd ;



  double NN=2.5;
  double sigma_R_Disk, sigma_T_Disk, sigma_Z_Disk;
  double sigma_R_DiskL,  sigma_T_DiskL, sigma_Z_DiskL;
  double sigma_R_DiskS,  sigma_T_DiskS, sigma_Z_DiskS;
  double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
  double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
  double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
  double Rho[8]={00.0}; double maxr=0.0;
  for(int i=0; i<8; ++i){ Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}

 
for (int i=0;i<2; ++i){
 test= ((double)rand()/(double)(RAND_MAX+1.))*maxr; ///total ages
     if(test<=Rho[0])                       {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0; age= 0.075;}
else if(test<=(Rho[0]+Rho[1]))              {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0; age=0.575; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]))       {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;age=1.5;  }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3])){sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2; age=2.5; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))
                                       {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8; age=4.0; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5]))
                                       {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4; age=6.0; }
else if(test<=maxr)                    {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5; age=8.5; }
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}
    if(i==0) {
    sigma_R_DiskS= sigma_R_Disk;
    sigma_T_DiskS= sigma_T_Disk;
    sigma_Z_DiskS= sigma_Z_Disk;}
    if(i==1){
    sigma_R_DiskL= sigma_R_Disk;
    sigma_T_DiskL= sigma_T_Disk;
    sigma_Z_DiskL= sigma_Z_Disk; }}
    
    
    double v_R_lmc=-57.0;
    double v_T_lmc=-226.0; 
    double v_Z_lmc= 221.0;
    double sigma_LMC=20.2; 
    double err_rlmc= 13.0; ///error of global velocity
    double err_tlmc= 15.0; 
    double err_zlmc= 19.0; 


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_DiskS, NN);
    SVT= RandN(sigma_T_DiskS, NN);
    SVZ= RandN(sigma_Z_DiskS, NN); }

    else if(s.struc==1){///Galactic bulge
    SVR= RandN(sigma_R_Bulge, NN);
    SVT= RandN(sigma_T_Bulge, NN);
    SVZ= RandN(sigma_Z_Bulge, NN); }

    else if(s.struc==2){///thick disk
    SVR= RandN(sigma_R_TDisk, NN);
    SVT= RandN(sigma_T_TDisk, NN);
    SVZ= RandN(sigma_Z_TDisk, NN); }

    else if(s.struc==3){///stellar halo
    SVR= RandN(sigma_R_halo, NN);
    SVT= RandN(sigma_T_halo, NN);
    SVZ= RandN(sigma_Z_halo, NN); }
    
    
    else if(s.struc>3){
    SVR= RandN(sigma_LMC, NN); 
    SVT= RandN(sigma_LMC, NN); 
    SVZ= RandN(sigma_LMC, NN); 
    SVZ +=  v_Z_lmc +RandR(-err_zlmc, err_zlmc);
    SVR +=  v_R_lmc +RandR(-err_rlmc, err_rlmc);
    SVT +=  v_T_lmc +RandR(-err_tlmc, err_tlmc);}
        
    if(s.struc==0 or s.struc==2)  SVT =SVT+ vro_sun*(1.00762*pow(Rsc/Dsun,0.0394) + 0.00712);
    
    s.vs=sqrt( SVR*SVR + SVT*SVT + SVZ*SVZ );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   if(l.struc==0){///Galactic disk
   LVR= RandN(sigma_R_DiskL, NN) ;
   LVT= RandN(sigma_T_DiskL, NN) ;
   LVZ= RandN(sigma_Z_DiskL, NN) ; }

   else if(l.struc==1){///Galactic bulge
   LVR= RandN(sigma_R_Bulge, NN) ;
   LVT= RandN(sigma_T_Bulge, NN) ;
   LVZ= RandN(sigma_Z_Bulge, NN) ; }

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN) ;
   LVT= RandN(sigma_T_TDisk, NN) ;
   LVZ= RandN(sigma_Z_TDisk, NN) ; }

   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN);
   LVT= RandN(sigma_T_halo, NN);
   LVZ= RandN(sigma_Z_halo, NN); }
   
   
   else if(l.struc>3){
   LVR= RandN(sigma_LMC, NN); 
   LVT= RandN(sigma_LMC, NN); 
   LVZ= RandN(sigma_LMC, NN); 
   LVZ +=  v_Z_lmc +RandR(-err_zlmc, err_zlmc);
   LVR +=  v_R_lmc +RandR(-err_rlmc, err_rlmc);
   LVT +=  v_T_lmc +RandR(-err_tlmc, err_tlmc); }
   
   
   if(l.struc==0 or l.struc==2)  LVT = LVT+vro_sun *(1.00762*pow(Rlc/Dsun,0.0394) + 0.00712);
   l.vl=sqrt( LVT*LVT + LVZ*LVZ + LVR*LVR );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  BETA  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    betal=betas=0.0;
    tetd= s.TET;
    test= double(l.Dl*cos(s.FI)*sin(tetd)/Rlc);
    if(fabs(test-1.0)<0.01)       betal= pi/2.0;
    else if(fabs(test+1.0)<0.01)  betal=-pi/2.0;
    else                          betal=asin(test);
    
    test= double(s.Ds*cos(s.FI)*sin(tetd)/Rsc); 
    if( fabs(test-1.0)<0.01)     betas=pi/2.0;
    else if(fabs(test+1.0)<0.01) betas=-pi/2.0;
    else                         betas=asin(test);
    
    if(Dsun < fabs(l.Dl*cos(s.FI)*cos(tetd))) betal= pi-betal; 
    if(Dsun < fabs(s.Ds*cos(s.FI)*cos(tetd))) betas= pi-betas; 

   if(fabs(l.Dl*cos(s.FI)*sin(tetd))>Rlc or fabs(test)>1.0){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<tetd<<"\t betal:  "<<betal<<endl;
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<"\t betas:   "<<betas<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(tetd)/Rlc<<"\t sin(s): "<<test<<endl;
   int ew; cin>>ew;}
///HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.TET>pi)  tetd=s.TET-2.0*pi; 
    deltal= pi - fabs(tetd) -fabs(betal);
    deltas= pi - fabs(tetd) -fabs(betas);  
    if(betal<0.0)  deltal= -1.0*deltal;
    if(betas<0.0)  deltas= -1.0*deltas;
    s.deltao= pi-fabs(tetd);
    if(tetd<0.0)  s.deltao=-1.0*s.deltao; 

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    s.SV_n1 =+SVR * sin(deltas)- SVT * cos(deltas);
    s.LV_n1 =+LVR * sin(deltal)- LVT * cos(deltal);
    s.VSun_n1=+VSunR*sin(s.deltao)-VSunT*cos(s.deltao);
    
    SVx= -SVR*cos(deltas)- SVT*sin(deltas);
    LVx= -LVR*cos(deltal)- LVT*sin(deltal);
    VSunx= -VSunR*cos(s.deltao) -VSunT*sin(s.deltao);
    
    s.SV_n2=-sin(s.FI)*(SVx) + cos(s.FI)*SVZ;
    s.LV_n2=-sin(s.FI)*(LVx) + cos(s.FI)*LVZ;
    s.VSun_n2=-sin(s.FI)*(VSunx)+cos(s.FI)*(VSunZ);
 
    
    vls1= l.xls*s.SV_n1 - s.LV_n1 +(1.0-l.xls)*s.VSun_n1;  ///Source - lens 
    vls2= l.xls*s.SV_n2 - s.LV_n2 +(1.0-l.xls)*s.VSun_n2;  /// Source -lens
    l.Vt=sqrt(fabs( vls1*vls1 + vls2*vls2 ) );
    

    if(vls1>0.0 and vls2>=0.0)         s.xi=atan(fabs(vls2)/fabs(vls1));//in radian
    else if(vls1<=0.0 and vls2>0.0)    s.xi=atan(fabs(vls1)/fabs(vls2)) + M_PI/2.0;
    else if(vls1<0.0 and vls2<=0.0)    s.xi=atan(fabs(vls2)/fabs(vls1)) + M_PI;
    else if(vls1>=0.0 and vls2<0.0)    s.xi=atan(fabs(vls1)/fabs(vls2)) + 3.0*M_PI/2.0;
    else if(s.xi > 2.0*M_PI)           s.xi-= 2.0*M_PI; 
    
    
    
    if(vls1==0.0 and vls2==0.0){
    cout<<"Big error both zeros:  "<<vls1<<"\t vls_b:  "<<vls2<<endl;
    int iee;  cin>>iee;}



    if (l.Vt<0.0 or l.Vt>1.0e6){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<s.vs<<endl;
    cout<<"source type:  "<<s.cl<<"\t s.struc:  "<<s.struc<<endl;
    int yee; cin>>yee;}
    //cout<<"(vrel_func):   s.deltao:  "<<s.deltao<<"FI: "<<s.FI<<endl;   
///*****************************************************
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
