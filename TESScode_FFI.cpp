#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include <cmath>
#include "VBBinaryLensingLibrary.h"

using std::cout;
using std::endl;
using std::cin;
///=============================================================================
const int    Num=40000;
const double MaxD=7.5;///kpc
const double RA=180.0/M_PI;
const double step=MaxD/(double)Num/1.0;///step in kpc
const double KP=3.08568025*pow(10.,19); // in meter.
const double G= 6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity=299792458.0;//velosity of light
const double Z_sun=float(0.0152);//#solar metalicity
const double M_sun=1.98892*pow(10.,30.0); //in [kg].
const double Mjupiter=1.898*pow(10,27.0); 
const double Mearth=  5.9722*pow(10.0,24.0);
const double mmin= double(13.0*Mjupiter/M_sun); //threshold for BD 
const double mmin0=double(0.000001*Mearth/M_sun);//in M_sun  low limit for planets

const double ml1=log10(mmin0); 
const double ml2=log10(2.0); 
const double dl1=0.0;
const double dl2=2.5;  

const double vro_sun=226.0;
const double AU=1.495978707*pow(10.0,11.0);
const double Rsun=6.9634*pow(10.0,8.0); ///solar radius [meter]

const double binary_fraction=double(2.0/3.0);
const double VSunR =11.1;
const double VSunT =vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ =7.25;
const double Avks=double(8.20922);
const int    Nco=int(48);  
const int    GG=200;///number of bins of tE distribution
const double ddeg=double(24.0/Nco/1.0);
///=============================================================================
const double Dsun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const int    N1=36224, N2=25000, N3=3818, N4=3500;///CMD_BESANCON, thinD, bulge, thickD, halo
const int    M=2;    ///No. of filter G, T
const double FWHM[2]= {0.059*5.0, 21.0*3.0};//20 Sep 

const double wave[4]={0.673,0.7865,0.532,0.797};//https://en.wikipedia.org/wiki/Photometric_system  G, T, BP, RP
const double AlAv[4]= {0.791986524645539, 0.617155245862836  ,1.0386670640894,0.601810722049874};
const double sigma[4]={0.017, 0.02,0.02, 0.02};// G, T, BP, RP  Table (2) Cardeli

//const int    Nes=int(26.4/(2.0/60.0/24.0)+5);

///=============================================================================
///SECTOR 12 https://tasoc.dk/docs/release_notes/tess_sector_12_drn17_v02.pdf
///SECTOR 39 https://tasoc.dk/docs/release_notes/tess_sector_39_drn56_v01.pdf
///SECTOR 4  https://tasoc.dk/docs/release_notes/tess_sector_04_drn05_v03.pdf
///SECTOR 31 https://tasoc.dk/docs/release_notes/tess_sector_31_drn47_v02.pdf
///SECTOR 13 https://tasoc.dk/docs/release_notes/tess_sector_13_drn18_v02.pdf

const int sector=1;  ///updated 6 Feb
const int camera=1; 
const int nmmax=159842;


///=============================================================================
struct source{
    int    nums,struc,cl;
    double Ds,TET,FI,ratios;
    double nsbl, blend, Fluxb, magb, Ai, Mab, Map;
    double lon, lat;
    double type, Tstar, logg, Rstar, mass, metal;
    double od_disk, od_ThD, od_bulge, od_halo,opd;
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs,nstart,nstarti;
    double nsdis[Num],nddis[Num];
    double ros,deltao;
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
};
struct lens{
    double Ml,Dl,ratiol,vl,vs,Vt,xls;
    double rhomaxl,tE,RE,t0, u0, tetE, mul;
    int    numl,struc;
    double piE, pirel;  
};
struct detection{
   int    det, numt;
   double dchi,t,tmin,tmax;
   double effs[GG+1],   effd[GG+1];  
   double effsC[GG+1],  effdC[GG+1]; 
   double Efte[GG+1][2], Efml[GG+1][2],Efdl[GG+1][2], Efmb[GG+1][2];   
   double timet;
   double weight; 
   double opta, optb, tE, fb, mbase; 
   double Dl, Ds, piE, pirel, ros, RE, vrel, u0, ml;  
   double fracB, fracP, fracR; 
   double Fb, Fp, Fr;
   double Nsim, Ndet;  
   double fbase, fside;//flag
   double fracS, fracBL;  
   double mst, emst, error; 
   double cadence, detect, satu, Tobs;  
   double gap[26], dur;
};
struct CMD{
    double Teff_d[N1+N3],logg_d[N1+N3],Mab_d[N1+N3],Rs_d[N1+N3],mass_d[N1+N3],type_d[N1+N3],metal_d[N1+N3];int cl_d[N1+N3];
    double Teff_b[N2],logg_b[N2],Mab_b[N2],Rs_b[N2],mass_b[N2],type_b[N2],metal_b[N2]; int cl_b[N2];  /// bulge
   // double Teff_t[N3],logg_t[N3],Mab_t[N3],Rs_t[N3],mass_t[N3],type_t[N3],metal[N3]; int cl_t[N3];  ///thick disk
    double Teff_h[N4],logg_h[N4],Mab_h[N4],Rs_h[N4],mass_h[N4],type_h[N4],metal_h[N4]; int cl_h[N4];  /// halo
};
struct extinc{
   double dis[100];///distance
   double Extks[100];///ks-band extinction
   double Aks;
   int flag; 
};
struct tessfile{
   int cam;  
   double Map, ra, dec, error, blend, radius;
   double Teff, logg, metal, elong, elat, glong, glat;  
};
struct coordinate{    
   double elo, ela;  
   double ra,  dec;  
   double l,   b; 
   double Nstar;
};
///=============================================================================
int  TEdet(detection & d, double mx, double mi, double mf);
int  Extinction(extinc & ex,source & s);
void read_cmd(CMD & cm);
void func_lens(lens & l, source & s);
void vrel(source & s,lens & l);
void Disk_model(source & s, int, int );
void optical_depth(source & s);
double ErrorTESS(detection & d, double maga); 
double Interpol(double ds, extinc & ex);
double RandN(double sigma, double);
double RandR(double down, double up);
time_t _timeNow;
unsigned int _randVal;
unsigned int _dummyVal;
FILE * _randStream;
///===========================================================================//
///                                                                           //
///                  Main program                                             //
///                                                                           //
///===========================================================================//	
int main()
{

   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   time( &_timeNow);
   printf("START time:   %s",ctime(&_timeNow));
   
   double nnt=double(2363.0+592.0+89.0+72.0+38.0+30.0);
   double frac[14]={0.0, 
   double(2363.0), 
   double(2363.0+592.0), 
   double(2363.0+592.0+89.0), 
   double(2363.0+592.0+89.0+9.0*1.0),
   double(2363.0+592.0+89.0+9.0*2.0), 
   double(2363.0+592.0+89.0+9.0*3.0), 
   double(2363.0+592.0+89.0+9.0*4.0),
   double(2363.0+592.0+89.0+9.0*5.0), 
   double(2363.0+592.0+89.0+9.0*6.0), 
   double(2363.0+592.0+89.0+9.0*7.0),
   double(2363.0+592.0+89.0+9.0*8.0), 
   double(2363.0+592.0+89.0+9.0*8.0+38.0), 
   double(2363.0+592.0+89.0+9.0*8.0+38.0+30.0)};
   
   VBBinaryLensing vbb;
   vbb.Tol=1.e-4;
   vbb.a1 =0.0;  
   vbb.LoadESPLTable("./ESPL.tbl");
  
   source s;
   lens l;
   detection d;
   CMD cm;
   extinc ex;
   coordinate co;
   tessfile t;
   read_cmd(cm);
     
     
   FILE* film;  
   FILE* fild; 
   FILE* result;  
   FILE* distr;
   FILE* tempo;    
   FILE* tefile;
   FILE* coord; 
   FILE* magdi;    

   
   if(sector==0){
   d.cadence=double(2.0/60.0/24.0); 
   d.detect=16.5; 
   d.satu= 3.5; }
   //d.Tobs= 13.7*13.0*2.0;
   //for(int i=0;i<26; ++i){
   //d.gap[i]= double(i+(i+1.0)*12.7);}
   //d.dur= 1.0;}
   
   
   if(sector==1){
   d.cadence=double(30.0/60.0/24.0); 
   d.detect=18.0;
   d.satu=  1.94;}
   
  
   if(sector==2){
   d.cadence=double(30.0/60.0/24.0);
   d.detect=18.0;
   d.satu=  -1.0;}
   //d.Tobs=  27.4132;
   //d.gap[0]=13.0521;
   //d.dur=   1.441;}
   
   if(sector==4){
   d.cadence=double(30.0/60.0/24.0);
   d.detect=18.0;
   d.satu=  0.0;
   d.Tobs=  25.95;}
   //d.gap[0]=7.6375; 
   //d.dur=   2.6743;
   //d.gap[1]=d.gap[0]+d.dur+2.2973;}
   
   
   if(sector==10){
   d.cadence=double(30.0/60.0/24.0); 
   d.detect= 18.0;
   d.satu=  -1.0;}
  // d.Tobs=   26.2486;
  // d.gap[0]=   12.3528;
  // d.dur=   0.9778;}
  
  
   if(sector==12){ 
   d.cadence=double(30.0/60.0/24.0);//in days 
   d.detect=18.0;
   d.satu=  1.0;}
   //d.Tobs=     27.942;
   //d.gap[0]=   14.046;
   //d.dur=     1.0375;}
   
   
   if(sector==13){
   d.cadence=double(30.0/60.0/24.0);
   d.detect=17.0;
   d.satu=  1.0;}
   //d.Tobs=   28.4417;
   //d.gap[0]= 13.775;
   //d.dur=    0.9292;}
  
  
   if(sector==31){
   d.cadence=double(10.0/60.0/24.0);
   d.detect=18.0;
   d.satu=  2.0;}
   //d.Tobs=  25.43473;
   //d.gap[0]=  12.94445;
   //d.dur=  1.4028;}    
      
      
   if(sector==39){
   d.cadence=double(10.0/60.0/24.0); 
   d.detect=18.0;
   d.satu=  1.0;}
   //d.Tobs=   27.9514;
   //d.gap[0]= 13.09;
   //d.dur= 1.00;}
    
     
   ///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   char   filenam0[40],filenam1[40], filenam2[40], filenam3[40], filenam4[40], filenam5[40], filenam6[40], filenam7[40];
   double tEmin=double(d.cadence*0.01);
   double tEmax=1.5*d.Tobs;      
   int    num, num2,nums, no; 
   int    nsim, ndet, los=0, countb, ffd, v0, v1, v2, v3, v4, detectf, fflag;
   double rf,testf, effe,Av, nwei, fte,ddo, nwea, meter,par0, par1, par2, par3, par4; 
   double lonn, countd, tei, tE, totn, wei,eff0, eff1, eff2, eff3, eff4, deltaA, module;   
   double gama, Nevent, Ntot=0.0;
   double maga, maga2, chi2, chi1,emt, u, As, Astar, timee,dt, dist, dmin, cdpp; 
   double flag0,flag1,flag2,  dista;
   int ID;
   double ra, dec; 
   for(int i=0; i<(GG+1); ++i){
   d.Efml[i][0]=0.0;   d.Efdl[i][0]=0.0;  d.Efmb[i][0]=0.0; d.Efte[i][0]=0.0;   
   d.Efml[i][1]=0.0;   d.Efdl[i][1]=0.0;  d.Efmb[i][1]=0.0; d.Efte[i][1]=0.0;   
   d.effdC[i]=0.0;     d.effsC[i]=0.0;} 
   d.Fb=0.0;     d.Fp=0.0;   d.Fr=0.0;
   d.Nsim=0.0;   d.Ndet=0.0;    
   d.fracS=0.0;  d.fracBL=0.0; 
   nsim=0;       ndet=0; 
   nwei=0.0;    nwea=0.0;  
   d.opta=0.0;  d.optb=0.0; 
   d.fracB=0.0; d.fracP=0.0; d.fracR=0.0; 
   d.tE=0.0;    d.fb=0.0;    d.mbase=0.0;
   d.mst=0.0;   d.emst=0.0; 
   d.Dl=0.0;    d.Ds=0.0; d.piE=0.0;  d.pirel=0.0;
   d.ros=0.0;   d.RE=0.0; d.vrel=0.0; d.u0=0.0; d.ml=0.0;
   for(int i=0; i<(GG+1);++i){d.effd[i]=0.0;   d.effs[i]=0.0;}
   for(int i=0; i<Num;  ++i){ s.nsdis[i]=0.0;  s.nddis[i]=0.0;}
 


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   sprintf(filenam0,"./files/light/lcF/%c%d%c%d.dat",'D',sector,'_',camera);
   distr=fopen(filenam0,"w");
   
   sprintf(filenam3,"./files/light/lcF/%c%c%c%c%d%c%d.dat",'T','E','S','S',sector,'_',camera);
   result=fopen(filenam3,"w");
   fclose(result);  
  
   sprintf(filenam4,"./files/light/lcF/%c%c%c%c%c%d%c%d.dat",'E','T','E','S','S',sector,'_',camera);
   tefile=fopen(filenam4,"w");
   fclose(tefile);  
 
 
   sprintf(filenam6,"./files/%c%c%c%d.dat",'T','E','m',sector);
   tempo=fopen(filenam6,"w");
   fclose(tempo);  
   
        
   for(int bn1=0; bn1<500000; ++bn1){
   
   num=int(RandR(0.0,nmmax-1.0));
   sprintf(filenam7,"./files/%c%c%c%d.dat",'F','F','I',sector);
   magdi=fopen(filenam7,"r"); 
   if(!magdi){cout<<"cannot read FFI1.dat "<<sector<<endl; exit(0);}
   for(int i=0; i<num; ++i){
   fscanf(magdi,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %d  %lf  %lf  %lf %lf\n",
   &t.Map,&t.ra,&t.dec,&t.error,&cdpp,&t.blend,&t.radius,&t.Teff,&t.logg,&t.metal,&t.cam,&t.elong,&t.elat,&t.glong,&t.glat);}
   fclose(magdi); 
   cout<<"file Information was read !!!!"<<endl; 
   cout<<"t.error:    "<<t.error<<endl;

   testf=RandR(0.0,nnt);  
   for(int i=0;i<13; ++i){
   if(double((testf-frac[i])*(testf-frac[i+1]))<0.0 or testf==frac[i]){no=i+1;  break;}}
   d.Tobs=double(13.7*no*2.0);
   for(int i=0;i<26;++i)   d.gap[i]=-1.0; 
   for(int i=1;i<=int(2.0*no); ++i) d.gap[i-1]=double(i-1.0+i*12.7);
   d.dur=1.0;
   cout<<"testf:  "<<testf<<"\t no:  "<<no<<"\t Tobs:  "<<d.Tobs<<endl;
   for(int i=0;i<int(2.0*no);++i)  cout<<"i:  "<<i<<"\t gap[i]:  "<<d.gap[i]<<endl;
   //int uue; cin>>uue;
   
   
   
   s.Map=  t.Map;   
   s.blend=t.blend;  
   s.Tstar=t.Teff;  
   s.Rstar=t.radius;
   s.logg= t.logg;
   if(t.metal>-9.5)s.metal= pow(10.0,t.metal)*Z_sun; 
   else            s.metal=-100.0; 
   s.lat=double(t.glat);  
   lonn= double(t.glong);
   if(lonn<=0.0) s.lon=360.0+lonn;
   else          s.lon=lonn; 
   s.TET=double(360.0-s.lon)/RA;///radian
   s.FI =double(s.lat/RA);///radian
   Disk_model(s , 1, 1);
   ex.flag=0;
   ex.flag=Extinction(ex,s);
   cout<<"ex.flag:  "<<ex.flag<<endl;
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   num=-1; 
   if(s.Tstar>cm.Teff_d[N1+N3-1])  num=int(N1+N3-1); 
   else if(s.Tstar<cm.Teff_d[0] )  num=0;  
   else{
   for(int j=1; j<int(N1+N3); ++j){
   if(float((s.Tstar-cm.Teff_d[j])*(s.Tstar-cm.Teff_d[j-1]))<=0.0){
   num=j; break; j=int(N1+N3);}}}
   if(num<0){cout<<"There is an error num<0:  "<<num<<"\t Tstar:  "<<s.Tstar<<endl;  int yye;  cin>>yye;  }
   cout<<"Tstar:  "<<s.Tstar<<"\t Teff_d[num]:  "<<cm.Teff_d[num]<<"\t Teff_d[num-1]:  "<<cm.Teff_d[num-1]<<endl;
    
    
    
   dmin=100000; dist=0.0;
   for(int j=-15; j<15; ++j){
   num2=int(num+j);
   if(num2>=0 and num2<int(N1+N3)){
   if(s.metal>-9.5 and s.logg>0.0)
                      dist=sqrt(pow(s.logg-cm.logg_d[num2],2.0)+pow(s.Rstar-cm.Rs_d[num2],2.0)+pow(s.metal-cm.metal_d[num2],2.0));
   else if(s.logg>0.0)dist=sqrt(pow(s.logg-cm.logg_d[num2],2.0)+pow(s.Rstar-cm.Rs_d[num2],2.0));
   else               dist=fabs(s.Rstar-cm.Rs_d[num2]);
   if(dist<dmin){dmin=dist;  nums=num2;}}}
   cout<<"nums:  "<<nums<<endl;
   cout<<"num2:    "<<nums<<"\t    logg  "  <<s.logg<<"\t logg_cmd:  "<<cm.logg_d[nums]<<endl;
   cout<<"metal:  "<<s.metal<<"\t metal:  "<<cm.metal_d[nums]<<endl;
   cout<<"Rstar:  "<<s.Rstar<<"\t Rs_d:  "  <<cm.Rs_d[nums]<<"\t dist:  "<<dist<<endl;
   s.mass= cm.mass_d[nums];
   s.cl=   cm.cl_d[nums]; 
   s.type= cm.type_d[nums];
   s.Mab=  cm.Mab_d[nums];
   dista=pow(10.0,0.2*fabs(s.Map-s.Mab))/100.0;//kpc 
   if(dista<step){cout<<"dista/step:  "<<double(dista/step)<<"\t dista:  "<<dista<<endl; int iie;  cin>>iie;}
   dmin=1000.0; dist=0.0; s.Ds=0;  
   for(double Ds=step; Ds<=dista; Ds=Ds+step){
   module=5.0*log10(Ds*100.0);
   if(ex.flag>0)Av=double(Interpol(Ds,ex)*Avks);
   else         Av=double(0.7*Ds); 
   if(Av<0.0)   Av=0.0;
   s.Ai=fabs(Av*AlAv[1])+RandN(sigma[1],1.5);
   if(s.Ai<0.0) s.Ai=0.0;   
   dist=fabs(s.Map - module- s.Ai - s.Mab);  
   if(dist<dmin){dmin=dist;  s.Ds=Ds;}}
   cout<<"Ds:  "<<s.Ds<<"\t dmin:  "<<dmin<<"\t Map:  "<<s.Map<<"\t Mab:  "<<s.Mab<<endl;
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

   s.nums=int(s.Ds/step); 
   rf=RandR(0.0,  s.Rostar0[s.nums]);
        if (rf<= s.rho_disk[s.nums]) s.struc=0;///thin disk
   else if (rf<=(s.rho_disk[s.nums]+s.rho_bulge[s.nums])) s.struc=1;/// bulge structure
   else if (rf<=(s.rho_disk[s.nums]+s.rho_bulge[s.nums]+s.rho_ThD[s.nums])) s.struc=2;///thick disk
   else if (  rf<=s.Rostar0[s.nums]) s.struc=3;///halo
   s.Fluxb=fabs(pow(10.0,-0.4*s.Map)/s.blend);//T-band  baseline flux
   s.magb =s.Map+2.5*log10(s.blend);//T-band baseline magnitude
   s.nsbl =double(1.0/s.blend);
   if(ex.flag>0) Av=double(Interpol(s.Ds,ex)*Avks); 
   else          Av=double(0.7*s.Ds);
   if(Av<0.0)    Av=0.0;
   s.Ai=fabs(Av*AlAv[1])+RandN(sigma[1],1.5);
   if(s.Ai<0.0) s.Ai=0.0;  
   s.Mab=s.Map-s.Ai-5.0*log10(s.Ds*100.0);  
   func_lens(l,s);

  
   l.t0=RandR(0.0,d.Tobs);
   optical_depth(s);
   s.nsdis[s.nums]+=1.0;
   d.weight=1.0; 
   nsim+=1;
   s.nddis[s.nums]+=1.0*d.weight;
   d.opta +=s.opd*1.0e8*d.weight;  
   nwei   +=1.0*d.weight; 
   d.Nsim +=1.0*d.weight;

    
   dt=double(d.cadence/3.0); 
   if(dt>double(l.tE*0.05))  dt=double(l.tE*0.05);  
   d.det=0;      d.numt=0;   countb=0;
   chi1=0.0;     chi2=0.0;   d.dchi=0.0;
   timee=0.0;    flag0=0.0;  d.error=0.0; 
   flag1=0.0;    flag2=0.0;
   d.fbase=0.0;  d.fside=0.0;  
   detectf=0;    ffd=0;  
  
  
   if(nsim<=50){
   ffd=1;
   sprintf(filenam1,"./files/light/lcF/%c%d%c%d%c%d%c%d.dat",'L',sector,'_',camera,'_',1,'_',nsim);
   fild=fopen(filenam1,"w");
   sprintf(filenam2,"./files/light/lcF/%c%d%c%d%c%d%c%d.dat",'M',sector,'_',camera,'_',1,'_',nsim);
   film=fopen(filenam2,"w");}
 
  
   for(d.t=-5.0; d.t<=d.Tobs+5.0; d.t=d.t+dt){
   u=sqrt(l.u0*l.u0+(d.t-l.t0)*(d.t-l.t0)/l.tE/l.tE);
   if(s.ros<99.9999) As=vbb.ESPLMag2(u, s.ros);
   else              As=1.0+2.0/s.ros/s.ros;  
   Astar=double(As*s.blend+1.0-s.blend);   
   maga=s.magb-2.5*log10(Astar);
   if(ffd>0)    fprintf(film,"%.4lf    %.4lf   %.4lf\n",d.t, maga, Astar);
   if(d.t>0 and d.t<=d.Tobs and 
        float((d.t-d.gap[0])*(d.t-d.gap[0]-d.dur))>=0.0 and float((d.t-d.gap[1])*(d.t-d.gap[1]-d.dur))>=0.0 
   and  float((d.t-d.gap[2])*(d.t-d.gap[2]-d.dur))>=0.0 and float((d.t-d.gap[3])*(d.t-d.gap[3]-d.dur))>=0.0
   and  float((d.t-d.gap[4])*(d.t-d.gap[4]-d.dur))>=0.0 and float((d.t-d.gap[5])*(d.t-d.gap[5]-d.dur))>=0.0
   and  float((d.t-d.gap[6])*(d.t-d.gap[6]-d.dur))>=0.0 and float((d.t-d.gap[7])*(d.t-d.gap[7]-d.dur))>=0.0
   and  float((d.t-d.gap[8])*(d.t-d.gap[8]-d.dur))>=0.0 and float((d.t-d.gap[9])*(d.t-d.gap[9]-d.dur))>=0.0
   and float((d.t-d.gap[10])*(d.t-d.gap[10]-d.dur))>=0.0 and float((d.t-d.gap[11])*(d.t-d.gap[11]-d.dur))>=0.0
   and float((d.t-d.gap[12])*(d.t-d.gap[12]-d.dur))>=0.0 and float((d.t-d.gap[13])*(d.t-d.gap[13]-d.dur))>=0.0
   and float((d.t-d.gap[14])*(d.t-d.gap[14]-d.dur))>=0.0 and float((d.t-d.gap[15])*(d.t-d.gap[15]-d.dur))>=0.0
   and float((d.t-d.gap[16])*(d.t-d.gap[16]-d.dur))>=0.0 and float((d.t-d.gap[17])*(d.t-d.gap[17]-d.dur))>=0.0
   and float((d.t-d.gap[18])*(d.t-d.gap[18]-d.dur))>=0.0 and float((d.t-d.gap[19])*(d.t-d.gap[19]-d.dur))>=0.0
   and float((d.t-d.gap[20])*(d.t-d.gap[20]-d.dur))>=0.0 and float((d.t-d.gap[21])*(d.t-d.gap[21]-d.dur))>=0.0
   and float((d.t-d.gap[22])*(d.t-d.gap[22]-d.dur))>=0.0 and float((d.t-d.gap[23])*(d.t-d.gap[23]-d.dur))>=0.0
   and float((d.t-d.gap[24])*(d.t-d.gap[24]-d.dur))>=0.0 and float((d.t-d.gap[25])*(d.t-d.gap[25]-d.dur))>=0.0){
   timee+=dt; 
   if(timee>=d.cadence){
   timee-=d.cadence; 
   if(maga>=d.satu and maga<=d.detect){
   emt= fabs(ErrorTESS(d,maga));
   deltaA=fabs(pow(10.0,-0.4*emt)-1.0)*Astar; 
   maga2=maga+RandN(emt, 3.0);
   d.error+= emt; 
   if(Astar<1.1) countb+=1;
   chi1 += pow((maga2-s.magb)/emt,2.0);
   chi2 += pow((maga2-     maga)/emt,2.0);
   if(d.det==0){
   flag2=0.0;
   if( float(s.magb-maga2)>float(5.0*emt))     flag2=1.0;
   if(d.numt>1 and float(flag2+flag1+flag0)>2.0)  d.det=1;
   flag0=flag1;
   flag1=flag2;}
   if(ffd>0) fprintf(fild,"%.4lf  %.5lf  %.7lf   %.5lf   %.7lf\n",d.t, maga2, emt, Astar+RandN(deltaA,3.0), deltaA);
   d.numt+=1;
   //if(d.numt>Nes){cout<<"Error  numt:  "<<d.numt<<"\t Nes:  "<<Nes<<endl;  exit(0);}
   }}}}//time loop 
   if(countb<5) d.fbase=0.0;
   else          d.fbase=1.0;  
   if(fabs(l.t0)<float(2.0*d.cadence) or fabs(d.Tobs-l.t0)<float(2.0*d.cadence))  d.fside=0.0;
   else            d.fside=1.0;    
   d.error=fabs(d.error/(d.numt+0.00000078)); 
   d.dchi=fabs(chi2-chi1);
   if(ffd>0){fclose(film);  fclose(fild);}
   //meter= metric(d); 
   
   

   if(d.dchi>=500.0 and d.det>0 and d.fbase>0.0 and d.fside>0.0){//Lensing detection
   detectf=  1;  
   ndet  +=  1;
   nwea  +=  1.0*d.weight; 
   d.Ndet+=  1.0*d.weight;      
   d.optb+=  s.opd*1.0e8*d.weight;  
   d.tE  +=  l.tE*d.weight; 
   d.fb  +=  s.blend*d.weight; 
   d.mbase+= s.magb*d.weight; 
   d.Dl  +=  l.Dl*d.weight; 
   d.Ds  +=  s.Ds*d.weight; 
   d.piE +=  l.piE*d.weight; 
   d.pirel+= l.pirel*d.weight; 
   d.ros +=  s.ros*d.weight; 
   d.RE  +=  l.RE/AU *d.weight; 
   d.vrel+=  l.Vt*d.weight; 
   d.u0  +=  l.u0*d.weight;
   d.ml  +=  l.Ml*d.weight;  
   d.mst +=  s.Map*d.weight; 
   d.emst+=  d.error*d.weight; 
   if(l.Ml<0.08 and l.Ml>mmin)      {d.fracB+=1.0*d.weight; d.Fb+=d.weight;}//fraction of BD
   else if(l.Ml<mmin and l.Ml>mmin0){d.fracP+=1.0*d.weight; d.Fp+=d.weight;}//fraction of FFPs
   if(s.cl<5)                       {d.fracR+=1.0*d.weight; d.Fr+=d.weight;}//fraction of giant stars 
   if(d.fside>0.0)                   d.fracS+=1.0*d.weight; 
   if(d.fbase>0.0)                   d.fracBL+=1.0*d.weight;}
   
 
   v0=int(TEdet(d,l.tE, tEmin, tEmax) );  
   v1=int(TEdet(d, log10(l.tE), log10(tEmin), log10(tEmax) ) ); 
   v2=int(TEdet(d,l.Dl, dl1, dl2) ); 
   v3=int(TEdet(d,log10(l.Ml), ml1, ml2));
   v4=int(TEdet(d,s.magb, d.satu, d.detect) );
   
   
   d.effs[v0] +=1.0;  
   d.effsC[v0]+=1.0;
   d.Efte[v1][0]+=1.0;  
   d.Efdl[v2][0]+=1.0;
   d.Efml[v3][0]+=1.0;
   d.Efmb[v4][0]+=1.0;
   
   if(detectf>0){
   d.effd[v0] +=1.0;  
   d.effdC[v0]+=1.0;  
   d.Efte[v1][1]+=1.0;
   d.Efdl[v2][1]+=1.0;
   d.Efml[v3][1]+=1.0; 
   d.Efmb[v4][1]+=1.0;}   
   

    
   if(detectf>0){ 
   sprintf(filenam6,"./files/%c%c%c%d.dat",'T','E','m',sector);
   tempo=fopen(filenam6,"a+");
   fprintf(tempo, "%.8lf   %.2lf   %d\n",l.tE, d.weight, v0);  
   fclose(tempo);}
  
 
   fprintf(distr,
   "%d  %d   %.5lf  %.5lf  "///4
   "%d  %.8lf   %.6lf   %.5lf  "///8
   "%d   %d  %.5lf   %.5lf  %.4lf  %.4lf  %.4lf  %.2lf  %.3lf  %.4lf  %.4lf  "   //19
   "%.5lf  %.8lf    %.2lf   %.4lf   "  //23
   "%.6lf  %.5lf  %.5lf  %.7lf  %.5lf  %.5lf  %.10lf  %.7lf  %.9lf  " ///32
   "%d  %.1lf  %d   %d  %.5lf  %.5lf   %.5lf %d %d\n", //42
   los, nsim, s.lat, lonn, //4
   l.struc, log10(l.Ml), l.Dl, l.vl, //8
   s.struc, s.nums, s.logg, s.Ds, s.Tstar, s.Rstar,s.metal, s.type, l.vs, s.Mab, s.Map,//19
   s.magb, s.blend, s.nsbl, s.Ai, //23
   l.tE, l.RE/AU, l.t0, l.mul, l.Vt, l.u0, s.opd*1.0e8, s.ros, l.tetE,//32
   detectf,d.dchi,d.det,d.numt, l.piE, meter,s.mass, s.cl,no);//42
   }
 



///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   
   sprintf(filenam4,"./files/light/lcF/%c%c%c%c%c%d%c%d.dat",'E','T','E','S','S',sector,'_',camera);
   tefile=fopen(filenam4,"w");   
   for(int j=0; j<(GG+1); ++j){
   d.effd[j]=double(d.effd[j]/(d.effs[j]+0.01));
   par0=double(tEmin+ (tEmax-tEmin)*j/GG/1.00);  
   par1=double(log10(tEmin)+(log10(tEmax)-log10(tEmin))*j/GG/1.00);  
   par2=double(dl1+ (dl2-dl1)*j/GG/1.00);  
   par3=double(ml1+ (ml2-ml1)*j/GG/1.00);  
   par4=double(d.satu+( d.detect-d.satu )*j/GG/1.00); 
   eff0=double(d.effdC[j]*100.000/(d.effsC[j]  +0.01));
   eff1=double(d.Efte[j][1]*100.0/(d.Efte[j][0]+0.01));
   eff2=double(d.Efdl[j][1]*100.0/(d.Efdl[j][0]+0.01));
   eff3=double(d.Efml[j][1]*100.0/(d.Efml[j][0]+0.01)); 
   eff4=double(d.Efmb[j][1]*100.0/(d.Efmb[j][0]+0.01)); 
   fprintf(tefile,"%.6lf %.10lf  %.6lf  %.10lf  %.6lf  %.10lf  %.6lf  %.10lf  %.6lf  %.10lf\n",       
   par0,eff0, par1,eff1, par2,eff2, par3,eff3, par4,eff4);}
   fclose(tefile);  
    
  
   countd=0.0;   tei=0.0; 
   sprintf(filenam6,"./files/%c%c%c%d.dat",'T','E','m',sector);
   tempo=fopen(filenam6,"r"); 
   for(int i=0; i<ndet;++i){
   fscanf(tempo,"%lf   %lf  %d\n", &tE, &wei, &v1); 
   countd+=wei;  
   tei   +=fabs(d.effd[v1]/(tE+0.000001))*wei;}
   tei= double(tei/countd); 
   fclose(tempo); 
   
   
   
   if(tei==0.0 or tei<0.0 or fabs(countd-ndet)>0.5){
   cout<<"Error tei: "<<tei<<"\t ndet:  "<<ndet<<endl;
   cout<<"last low:  "<<tE<<wei<<fflag<<v1<<"\t  countd:  "<<countd<<endl;  
   int yye; cin>>yye;}
   
    
   s.nstart=0.0;
   for(int i=1; i<Num; ++i){
   effe= double(s.nddis[i]/(s.nsdis[i]+0.01));
   s.nstart+= s.Rostari[i]*(s.Nstart/s.Rostart)*effe;}// per deg^{-2}
   d.opta=double(d.opta/nwei);  
   gama=2.0*d.opta*1.0e-8*tei/M_PI;//[N_event/days
   Nevent=gama*1.0*d.Tobs;  
   Ntot += Nevent*ddeg*ddeg;//comulative
   if(s.nstart<=0.0 or d.opta<=0.0 or nwei<=0.0 or gama<=0.0 or Nevent<0.0 or d.Tobs<=0.0 or ddeg<=0.0 or 
   nwea<=0.0 or s.Rostart<=0.0 or s.Nstart<=0.0){
   cout<<"Error nstart:  "<<s.nstart<<"\t opta:  "<<d.opta<<"\t nwei:  "<<nwei<<endl;
   cout<<"gama:  "<<gama<<"\t Nevent:  "<<Nevent<<"\t Tobs:  "<<d.Tobs<<"\t nwea:  "<<nwea<<endl;
   cout<<"ddeg:  "<<ddeg<<"\t Ntot:  "<<Ntot<<"\t co.Nstar[los]: "<<co.Nstar<<endl;  
   cout<<"Nstart:  "<<s.Nstart<<"\t Rostart:   "<<s.Rostart<<endl;  int yye;  cin>>yye;}
   
  
   d.tE =   double(d.tE/nwea);
   d.fb=    double(d.fb/nwea); 
   d.mbase= double(d.mbase/nwea); 
   d.optb=  double(d.optb/nwea); 
   d.Dl =   double(d.Dl/nwea);  
   d.Ds  =  double(d.Ds/nwea); 
   d.piE =  double(d.piE/nwea); 
   d.pirel= double(d.pirel/nwea); 
   d.ros=   double(d.ros/nwea); 
   d.RE =   double(d.RE/nwea); 
   d.vrel=  double(d.vrel/nwea); 
   d.u0=    double(d.u0/nwea); 
   d.ml=    double(d.ml/nwea); 
   d.mst=   double(d.mst/nwea);  
   d.emst=  double(d.emst/nwea);    
   d.fracB= double(d.fracB*100.0/nwea);  
   d.fracP= double(d.fracP*100.0/nwea);  
   d.fracR= double(d.fracR*100.0/nwea);  
   sprintf(filenam3,"./files/light/lcF/%c%c%c%c%d%c%d.dat",'T','E','S','S',sector,'_', camera);
   result=fopen(filenam3,"a+");
   fprintf(result,"%d  %d  %.5lf   %.5lf   %.5lf   %.5lf  %.5lf   %.5lf   %d  %d  "//10
   "%.6lf  %.6lf  %.6lf   %.6lf "  //14
   "%.7lf  %.7lf  %.8lf  %.8lf  %.8lf  " //19
   "%.5lf  %.5lf  %.6lf   %.6lf  %.6lf   %.8lf   %.8lf   %.8lf   %.6lf  %.5lf   %.6lf   %.8lf   %.6lf  %.6lf  %.6lf  "//34
   "%.7lf   %.7lf  %.7lf   %.7lf   %.7lf  %.4lf  %.4lf   %.5lf   %.9lf\n",//46  
   camera,los, co.ela, co.elo, co.b, co.l, co.dec, co.ra, nsim, ndet,//10
   log10(s.Nstart),log10(s.Rostart), log10(tei), log10(s.nstart), //14
   log10(d.opta), log10(d.optb), log10(gama), Nevent, Ntot, //19
   d.tE, d.mbase, d.fb, d.Dl, d.Ds, log10(d.piE), log10(d.pirel), log10(d.ros), d.RE, d.vrel, d.u0, log10(d.ml),
   d.fracB,d.fracP, d.fracR, //34
   double(d.Fb*100.0/d.Ndet), double(d.Fp*100.0/d.Ndet), double(d.Fr*100.0/d.Ndet), //37
   double(d.fracS*100.0/d.Ndet), double(d.fracBL*100.0/d.Ndet), d.Nsim, d.Ndet,d.mst,log10(d.emst));//46
   fclose(result);
    cout<<"=============================================================="<<endl;
    cout<<"camera:  "<<camera<<endl;
    cout<<"latitude: "<<s.lat<<"\t longtitude: "<<s.lon<<"\t los:  "<<los<<endl;
    cout<<"nsim:  "<<nsim<<"\t ndet:  "<<ndet<<endl;
    cout<<"nsim:      "<<nsim<<"\t ndet:  "<<ndet<<endl;
    cout<<"tei:     "<<tei<<"\t nstart:  "<<s.nstart<<endl;
    cout<<"gama:  "<<gama<<"\t Nevent:  "<<Nevent<<endl;
    cout<<"tE:  "<<d.tE<<"\t mbase:  "<<d.mbase<<"\t fb:  "<<d.fb<<endl;
    cout<<"DL:   "<<d.Dl<<"\t Ds:  "<<d.Ds<<"\t pirel:  "<<d.pirel<<endl;
    cout<<"Ro_star:  "<<d.ros<<"\t RE:  "<<d.RE<<"\t u0:  "<<d.u0<<endl; 
    cout<<"vrel:  "<<d.vrel<<"\t piE:  "<<d.piE<<"\t ML:  "<<d.ml<<endl;
    cout<<"No.  total ML events:  "<<Ntot<<endl;
    cout<<"mst:  "<<d.mst<<"\t emst:  "<<d.emst<<endl;
    cout<<"fracB:  "<<d.fracB<<"\t fracP:  "<<d.fracP<<"\t fracR:  "<<d.fracR<<endl;
    cout<<"fracB_Com:  "<<double(d.Fb*100.0/d.Ndet)<<endl;
    cout<<"fracP_Com:  "<<double(d.Fp*100.0/d.Ndet)<<endl;
    cout<<"fracR_Com:  "<<double(d.Fr*100.0/d.Ndet)<<endl;
    cout<<"frac detecting one side:  "<<double(d.fracS*100.0/d.Ndet)<<endl;
    cout<<"frac NOT detecting baseline:  "<<double(d.fracBL*100.0/d.Ndet)<<endl;
    cout<<"Total simulated events:  "<<d.Nsim<<"\t total detected events:  "<<d.Ndet<<endl;
    cout<<"=============================================================="<<endl;
    fclose(distr);
    fclose(_randStream);
    return(0);
}
///==========================================================================//
///                                                                          //
///                   metric                                                 //
///                                                                          //
///==========================================================================//
/*double metric( detection & d)
{
    double avem, median, mean, cutoff,STD, meter;
    double temp=0.0; int cob=0;  median=0.0;  
    double base[d.numt]={0.0};  
    int numn;  
    
    for(int i=0;   i<d.numt; ++i){
    for(int j=i+1; j<d.numt; ++j){ 
    if(d.magn[i]>d.magn[j]){
    temp = d.magn[i];
    d.magn[i]=d.magn[j];
    d.magn[j] = temp;}}} 
    cob= int(d.numt/2.0);
    if(d.numt%2==1)  median=d.magn[cob];  
    else             median=double(d.magn[cob]+d.magn[cob-1])*0.5;  
    
    
    
    cutoff=0.0;  
    for(int i=0; i<d.numt; ++i)    
    d.magn[i]=double(d.magn[i]-median);//REL_MAG  
    if(d.numt%2==1) cutoff=d.magn[cob];  
    else            cutoff=(d.magn[cob]+d.magn[cob-1])*0.5;  



    numn=0;  
    mean=0.0; 
    for(int i=0; i<d.numt; ++i){
    if(d.magn[i]>cutoff or d.magn[i]==cutoff){//baseline
    base[numn]=d.magn[i]; 
    mean+= d.magn[i]; 
    numn+=1;}}
    mean=double(mean/(numn+1.0e-10)/1.0);  


    STD=0.0;  
    meter=0.0;  
    for(int i=0; i<numn; ++i){
    STD+=double(base[i]-mean)*(base[i]-mean);}  
    STD=sqrt(STD/(numn+1.0e-10)/1.0); 
    if(STD==0.0) STD=1.0e-10;  
    meter=fabs(d.magn[d.numt-1]-d.magn[0])/(STD+1.0e-10);


    if(meter<0 or d.numt==0 or fabs(cutoff)>0.01 or mean==0.0 or STD==0.0){
    cout<<"Error metric: "<<meter<<"\t num:  "<<d.numt<<"\t STD:  "<<STD<<endl;
    cout<<"mean:  "<<mean<<"\t error:  "<<d.error<<"\t dchi:  "<<d.dchi<<endl;
    exit(0);}
    return(meter);    
}
*/
///#############################################################################
double ErrorTESS(detection & d, double maga){  //checked 6 Feb  2024
   double emt=-1.0, m;     
   
   if(sector==12){
   if(maga<8.2) emt=double(0.195*maga-5.6  ); 
   else         emt=double(0.26 *maga-6.133);}
   
   if(sector==39){
   if(maga<8.2) emt=double(0.198*maga -5.38); 
   else         emt=double(0.265*maga -5.93);}
   
   if(sector==4){
   if(maga<8.0)                emt=double(0.20*maga -5.64); 
   if(maga>=8.0 and maga<11.0) emt=double(0.22*maga -5.8); 
   if(maga>=11.0)              emt=double(0.32*maga -6.9);}
   
   if(sector==31){
   if(maga<8.0)                emt=double(0.2*maga  -5.408); 
   if(maga>=8.0 and maga<11.0) emt=double(0.224*maga-5.6 );  
   if(maga>=11.0)              emt=double(0.32*maga -6.656);}
   
   if(sector==13 or sector==10){
   if(maga<8.0)                emt=double(0.2*maga -5.64); 
   if(maga>=8.0 and maga<11.5) emt=double(0.242*maga-5.976);
   if(maga>=11.5)              emt=double(0.3*maga-6.643); }
   
   if(sector==0){///CTL
   if(maga<8.2)    emt=double(0.198*maga-5.61); 
   if(maga>=8.2)   emt=double(0.250*maga-6.03); }

   if(sector==1){
   if(maga<10.3)    emt=double(0.198*maga-5.61); 
   if(maga>=10.3)   emt=double(0.299*maga-6.65); }     
     
   emt=emt+RandN(0.1,3.0);
   if(emt<-5.0)   emt=-5.0;  
   emt=pow(10.0,emt);

   if(emt<0.00001 or emt>0.5 or maga<0.0){
   cout<<"Error emt:  "<<emt<<"\t maga:  "<<maga<<endl;}
   return(emt); 
}
///#############################################################################
int TEdet(detection & d, double mx, double mi, double mf){
    int vf=-1; 
    double sstep= double(mf-mi)/GG/1.00;  
    if(mx<mi or mx==mi)     vf=0;  
    else if(mx>mf or mx==mf)vf=GG;  
    else                    vf=int(double(mx-mi)/sstep );  
    if(vf<0 or vf>GG){
    cout<<"Error vf: "<<vf<<"\t mx:  "<<mx<<"\t mi:  "<<mi<<"\t mf:  "<<mf<<endl;
    int rre; cin>>rre;} 
    return(vf); 
}
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex)
{
  double F=-1.0;
  if(ds<ex.dis[0])        F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{ 
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] and ds<ex.dis[i+1]){
  F = ex.Extks[i]+(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0 or F<0.0){cout<<"ERROR big Extinction(ds): "<<F<<"\t ds: "<<ds<<endl; exit(0); }
  return(F);
}
///==================================================================
double RandN(double sigma, double nn){
   double rr,f,frand;
   do{
   rr=RandR(-sigma*nn , sigma*nn); ///[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/(sigma*sigma));
   frand=RandR(0.0 , 1.0);
   }while(frand>f);
   return(rr);
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm)
{
    //mass, Teff, Age, logL,  log(g),  Z,  Rs,  MB, MV, MI, MK, Cl, type (13)
    int yye;
    double age, lumi, MB, MV, MK, MI, Gg, T, GBP, GRP;
    char filename[40];
    FILE *fp2;
////=================================== THIN DISK ==============================
    int j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c%c%c.dat",'C','M','D','T','i','W','S','B');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTiW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.Teff_d[j],&age,&lumi,&cm.logg_d[j],&cm.metal_d[j],&cm.Rs_d[j],&MB,&MV, &MI, &MK, &cm.cl_d[j], &cm.type_d[j]);
    
    Gg=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428- 0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                   +RandN(0.04474, 1.5);  
    
    if(double(GBP-GRP)<=6.0 and double(GBP-GRP)>=-1.0) 
    T= Gg-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473+RandN(0.006, 1.5); 
    else T=Gg-0.430 + RandN(0.6,1.5);
    cm.Mab_d[j]= T;
    
    if(cm.mass_d[j]<=0.0 or cm.Teff_d[j]<0.0 or cm.metal_d[j]>0.2 or int(cm.cl_d[j])==6 or cm.type_d[j]>8.0 or cm.type_d[j]<1.0){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<"\t cm.metal[j]:  "<<cm.metal_d[j]<<"\t age:  "<<age<<endl;
    cout<<"Teff_d[j]:  "<<cm.Teff_d[j]<<"\t Teff_d[j-1]:  "<<cm.Teff_d[j-1]<<endl;
    cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=(N1+N3)){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1+N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;





////=================================== BULGE ==================================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','b','W','S');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDbW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.Teff_b[j],&age,&lumi,&cm.logg_b[j],&cm.metal_b[j],&cm.Rs_b[j],&MB,&MV,&MI,&MK,&cm.cl_b[j],&cm.type_b[j]);
    
    Gg=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428-0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                    +RandN(0.04474, 1.5);  
    if(double(GBP-GRP)<=6.0 and double(GBP-GRP)>=-1.0) 
    T= Gg-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473 +RandN(0.006, 1.5); 
    else T=Gg-0.430;
    
    //cm.Mab_b[0][j]= G;  
    cm.Mab_b[j]= T;
    if(cm.mass_b[j]<=0.0 or cm.Teff_b[j]<0.0 or age>10 or cm.metal_b[j]>0.9 or cm.cl_b[j]==6 or cm.type_b[j]>=8.0 or 
    (j>0 and cm.Teff_b[j]<cm.Teff_b[j-1])  ){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<"\t age:  "<<age<<"\t cm.metal[j]:"<<cm.metal_b[j]<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<"\t Teff:  "<<cm.Teff_b[j]<<endl; 
    cout<<"Teff_b[j]:  "<<cm.Teff_b[j]<<"\t Teff_b[j-1]:  "<<cm.Teff_b[j-1]<<endl; cin>>yye;}
    j++;} 
    fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
    cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;


////=================================== THICK DISK =============================
    /*
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c%c.dat",'C','M','D','T','k','W','S');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTkW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.Teff_t[j],&age,&lumi,&cm.logg_t[j],&cm.metal_t[j],&cm.Rs_t[j],&MB,&MV,&MI,&MK,&cm.cl_t[j],&cm.type_t[j]);
    G=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428-0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                    +RandN(0.04474, 1.5);  
    if(GBP-GRP<6.0 and GBP-GRP>-1.0) 
    T= G-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473 +RandN(0.006, 1.5); 
    else T=G-0.430;
    //cm.Mab_t[0][j]= G;  
    cm.Mab_t[j]= T;
    
    if(cm.mass_t[j]<=0.0 or cm.Teff_t[j]<0.0 or cm.metal_t[j]>0.2 or cm.cl_t[j]==6 or cm.type_t[j]>=8.0 or
      (j>0 and cm.Teff_t[j]<cm.Teff_t[j-1]) ){
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; 
    cout<<"Teff_t[j]:  "<<cm.Teff_t[j]<<"\t Teff_t[j-1]:  "<<cm.Teff_t[j-1]<<endl;   cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;
*/



////=================================== STELLAR HALO ===========================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','h','W','S');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDhW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.Teff_h[j],&age,&lumi,&cm.logg_h[j],&cm.metal_h[j],&cm.Rs_h[j],&MB,&MV,&MI,&MK,&cm.cl_h[j],&cm.type_h[j]);
    
    Gg=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428-0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                    +RandN(0.04474, 1.5);  
    if(double(GBP-GRP)<=6.0 and double(GBP-GRP)>=-1.0) 
    T= Gg-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473 +RandN(0.006, 1.5); 
    else T=Gg-0.430;
    cm.Mab_h[j]= T;
    
    if(cm.mass_h[j]<=0.0 or age<0 or cm.cl_h[j]==6 or cm.Teff_h[j]<0.0 or cm.metal_h[j]>0.1 or cm.cl_h[j]>7 or 
    cm.type_h[j]>9 or (cm.cl_h[j]<5 and int(cm.type_h[j])==9) or (j>0 and cm.Teff_h[j]<cm.Teff_h[j-1]) ){
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"age:  "<<age<<"\t cm.metalh[j]:  "<<cm.metal_h[j]<<"\t mass:  "<<cm.mass_h[j]<<endl; 
    cout<<"Teff_h[j]:  "<<cm.Teff_h[j]<<"\t Teff_h[j-1]:  "<<cm.Teff_h[j-1]<<endl;cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
   cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
   cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
}
///==============================================================//
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*M_sun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opd=0.0;
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;}
    s.opd= s.od_disk+s.od_ThD+s.od_bulge+s.od_halo;///total
    //cout<<"total_opticalD: "<<s.opd<<"\t od_disk: "<<s.od_disk<<endl;
   // cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///==============================================================//
///                                                              //
///                  Func source   Initial amounts               //
///                                                              //
///==============================================================//
/*
void func_source(source & s, CMD & cm, extinc & ex, tessfile & t, double elong, double elat, int c)
{
    int    struc,nums,num, num2=0;
    double rf, Ds,Av, dist, dmin, module, dista;
    s.Fluxb=0.0; s.magb=0.0; s.Ai=0.0; s.Map=0.0; s.Mab=0.0;
    s.nsbl=0.0;  s.blend=0.0;  
   
    do{
    num=int(RandR(0.0,nm[c]-1.1));
    num2+=1; 
    if(num2>=nm[c]) break;
    }while(fabs(t.elong[num]-elong)>3.0  or fabs(t.elat[num]-elat)>3.0  );


    s.Map  = t.Map[num]; 
    s.blend= t.blend[num];  
    s.Tstar= t.teff[num];  
    s.Rstar= t.radius[num]; 
    s.logg = t.logg[num];  
    if(t.mh[num]>-9.5) s.metal= pow(10.0,t.mh[num])*Z_sun; 
    else               s.metal=-10.0;  
    
    if(s.Tstar<=0.0  or s.logg<0.0 or s.blend<=0.0 or  s.blend>1.0 or s.Rstar<=0.0 or s.Map<=0.0){
    cout<<"Error Tstar:  "<<s.Tstar<<"\tlogg:  "<<s.logg<<"\t blend:  "<<s.blend<<endl;
    cout<<"Rstar:  "<<s.Rstar<<"\t elong:  "<<elong<<"\t elat:  "<<elat<<endl; 
    cout<<"elong[num]:  "<<t.elong[num]<<"\t elat[num]:  "<<t.elat[num]<<endl;
    cout<<"num:  "<<num<<"\t  nm[c]:  "<<nm[c]<<"\t camera:  "<<c<<endl;  
    int eei; cin>>eei;}
    
    
    num=-1; 
    if(s.Tstar>cm.Teff_d[N1+N3-1])  num=int(N1+N3-1); 
    else if(s.Tstar<cm.Teff_d[0] )  num=0;  
    else{
    for(int j=1; j<int(N1+N3); ++j){
    if(float((s.Tstar-cm.Teff_d[j])*(s.Tstar-cm.Teff_d[j-1]))<=0.0){
    num=j; break; j=int(N1+N3);}}}
    if(num<0){cout<<"There is an error num<0:  "<<num<<"\t Tstar:  "<<s.Tstar<<endl;  int yye;  cin>>yye;  }
    //cout<<"Tstar:  "<<s.Tstar<<"\t Teff_d[num]:  "<<cm.Teff_d[num]<<"\t Teff_d[num-1]:  "<<cm.Teff_d[num-1]<<endl;
    
    
    dmin=100000; dist=0.0;
    for(int j=-15; j<15; ++j){
    num2=int(num+j);
    if(num2>=0 and num2<int(N1+N3) ){
    if(s.metal>-9.5) dist=sqrt(pow(s.logg-cm.logg_d[num2],2.0)+pow(s.Rstar-cm.Rs_d[num2],2.0)+pow(s.metal-cm.metal_d[num2],2.0));
    else             dist=sqrt(pow(s.logg-cm.logg_d[num2],2.0)+pow(s.Rstar-cm.Rs_d[num2],2.0));
    //cout<<"num2:    "<<num2<<"\t    logg  "  <<s.logg<<"\t logg_cmd:  "<<cm.logg_d[num2]<<endl;
    //cout<<"metal:  "<<s.metal<<"\t metal:  "<<cm.metal_d[num2]<<endl;
    //cout<<"Rstar:  "<<s.Rstar<<"\t Rs_d:  "  <<cm.Rs_d[num2]<<"\t dist:  "<<dist<<endl;
    if(dist<dmin){dmin=dist;  nums=num2;}}}
    //cout<<"nums:  "<<nums<<endl;
    
    
    s.mass= cm.mass_d[nums];
    s.cl=   cm.cl_d[nums]; 
    s.type= cm.type_d[nums];
    s.Mab=  cm.Mab_d[nums];
    dista=pow(10.0,0.2*fabs(s.Map-s.Mab))/100.0;//KPC  
    if(dista<step){cout<<"dista/step:  "<<double(dista/step)<<"\t dista:  "<<dista<<endl; int iie;  cin>>iie;}

    
    dmin=1000.0; dist=0.0; s.Ds=0;  
    for(Ds=0.01; Ds<=dista; Ds=Ds+step){
    module=5.0*log10(Ds*100.0); 
    ex.Aks=Interpol(Ds,ex);
    Av=ex.Aks*Avks;
    if(Av<0.0) Av=0.0;
    s.Ai=fabs(Av*AlAv[1])+RandN(sigma[1],1.5);
    if(s.Ai<0.0) s.Ai=0.0;   
    dist= fabs(s.Map - module- s.Ai - s.Mab);  
    if(dist<dmin){dmin=dist; s.Ds=Ds; }}
    //cout<<"Ds:  "<<s.Ds<<"\t dmin:  "<<dmin<<"\t Map:  "<<s.Map<<"\t Mab:  "<<s.Mab<<endl;
    
    
    ex.Aks=Interpol(s.Ds,ex);
    Av=ex.Aks*Avks;
    if(Av<0.0)  Av=0.0;
    s.Ai=fabs(Av*AlAv[1])+RandN(sigma[1],1.5);
    if(s.Ai<0.0) s.Ai=0.0;  
    //cout<<"Extinction:  "<<s.Ai<<endl; 
    
    
   
    nums=int(s.Ds/step); 
    rf=RandR(0.0, s.Rostar0[nums]);
         if (rf<= s.rho_disk[nums]) struc=0;///thin disk
    else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
    else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
    else if (  rf<=s.Rostar0[nums]) struc=3;///halo
    s.struc=struc;
    //cout<<"struc:  "<<s.struc<<endl;
    
    
    s.Fluxb=fabs(pow(10.0,-0.4*s.Map)/s.blend);//T-band  baseline flux
    s.magb =s.Map+2.5*log10(s.blend);//T-band baseline magnitude
    s.nsbl =double(1.0/s.blend);
    s.nums=int(s.Ds/step);
   
    //cout<<"*******************************************"<<endl;
    //cout<<"Ds:     "<<s.Ds<<"\t   Map:  "<<s.Map<<endl; 
    //cout<<"Mab:    "<<s.Mab<<"\t  Ai:  "<<s.Ai<<endl;
    //cout<<"Fluxb:  "<<s.Fluxb<<"\tmagb:  "<<s.magb<<endl; 
   // cout<<"Tstar:  "<<s.Tstar<<"\tnsbl:  "<<s.nsbl<<endl;
   // cout<<"Mass:   "<<s.mass<<"\t Rstar:  "<<s.Rstar<<endl; 
   // cout<<"type:   "<<s.type<<"\t cl:  "<<s.cl<<endl; 
   // cout<<"logg:   "<<s.logg<<"\t metal:  "<<s.metal<<endl;
   // cout<<"nums:  "<<s.nums<<endl;
   // cout<<" end of  func_source!!!! "<<endl;
}
*/
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
     double sig,Lon,Lat;
     int uue, flag=0;
     if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
     else sig=1.0;
     double delt=fabs(s.lon)-floor(fabs(s.lon));
     if(delt>1.0 or delt<0.0){cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
     else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
     else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
     else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
     else               Lon=(floor(fabs(s.lon))+0.75)*sig;
     if(fabs(s.lon)<0.24999999)     Lon=360.00;
     if(fabs(s.lon-360.0)<0.2499999)  Lon=360.00;
     //cout<<"s.lat:  "<<s.lon<<"\t s.lat:  "<<s.lat<<endl;
     //cout<<"Lon:    "<<Lon<<"\t     Lat:  "<<Lat<<endl;


     if(s.lat<0.0) sig=-1.0;
     else sig=1.0;
     delt=fabs(s.lat)-floor(fabs(s.lat));
     if(delt>1.0 or delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
     else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
     else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
     else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
     else                Lat=(floor(fabs(s.lat))+0.75)*sig;
     
     
     
     if(fabs(Lon)<0.2499999) Lon=360.00; 
     if(Lat==-0.00)  Lat=0.00;
    // if(Lat>10.0)    Lat=10.0; 
   //  if(Lat<-10.0)   Lat=-10.0;
    // if(Lon>100.0 and Lon<260.0) Lon=100.0; 
     if(fabs(Lon)<0.24999999)    Lon=360.00;
     cout<<"Lon:    "<<Lon<<"\t     Lat:  "<<Lat<<endl;
    // if(Lon>360.000 or Lon<0.25 or fabs(Lat)>10.0 or (Lon>100 and Lon<260)){
    // cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}
     //cout<<"(2) Lon:    "<<Lon<<"\t     Lat:  "<<Lat<<endl;

     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");


     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     //FILE *SD;
    // SD=fopen("./files/Ext/saved_direction.txt","r");
    // for(int i=0; i<64881; ++i) {
    // fscanf(SD,"%lf %lf \n",&latit,&lonti);
    // if(fabs(Lat-latit)<0.1 and fabs(Lon-lonti)<0.1){
    // cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
    // cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
     if(ex.dis[i]<0.2  or ex.dis[i]>50.0 or ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0; }}
     fclose(fpd);}
     //cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<endl;
     
     return(flag);
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s)
{
    double f,test,t1,t2,test1;
    double rholens[s.nums+2];
    l.rhomaxl=0.0;
    
    for(int k=1; k<s.nums; ++k){
    rholens[k]=0.0;
    l.Dl =k*step;
    l.xls=l.Dl/s.Ds;
    //if(l.Dl>s.Ds){cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl; l.Dl=s.Ds*0.9;}//  int yye; cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}
   

    do{
    l.numl=int(RandR(1.0,s.nums-1.0));
    test =     RandR(0.0,l.rhomaxl);
    }while(test>rholens[l.numl]);
    l.Dl= double(l.numl*step);


   double randflag=RandR(0.0,s.Rostar0[l.numl]);
       if (randflag<=s.rho_disk[l.numl])    l.struc=0;///thin disk
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1; // bulge structure
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2; //thick disk
  else if (randflag<= s.Rostar0[l.numl]) l.struc=3;//halo
  else    {cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}

///#############################################################################
  test1=RandR(0.0,1.0);
  if(test1<0.2){//Mroz et al., 2017 Nature
  do{
  l.Ml=RandR(mmin0, mmin);
  f=double(1.0/l.Ml); 
  test=RandR(1.0/mmin , 1.0/mmin0);//log uniform
  }while(test>f);}
  
  
  else{
  if(l.struc==0){//thin disk
  t1=pow(2.0,-3.0);
  t2=pow(mmin,-0.7)*pow(0.08,0.7-1.6);
  do{
  l.Ml=RandR(mmin,2.0);
  test=RandR(t1  , t2);
  if(l.Ml<0.08)                     f=pow(l.Ml,-0.7)*pow(0.08,0.7-1.6);   
  else if(l.Ml<1.0 and l.Ml>=0.08)  f=pow(l.Ml,-1.6);
  else if(l.Ml>=1.0)                f=pow(l.Ml,-3.0);
  if(t1>t2 or f<t1 or f>t2){cout<<"Error_1 t1:  "<<t1<<"\t t2: "<<t2<<"\t f:  "<<f<<"\t Ml: "<<l.Ml<<endl;  exit(0);}
  }while(test>f);}
  //=========================================================
  if(l.struc==1){///Galactic bulge
  t1=pow(1.4,-2.35)*pow(0.08,2.35-0.7); 
  t2=pow(mmin,-0.7);
  do{
  l.Ml=RandR(mmin,1.4);
  test=RandR(t1  , t2);
  if(l.Ml<0.08)         f=pow(l.Ml, -0.7); 
  else if(l.Ml>=0.08)   f=pow(l.Ml,-2.35)*pow(0.08,2.35-0.7);
  if(t1>t2 or f<t1 or f>t2 or l.Ml>1.4){cout<<"Error_2 t1:  "<<t1<<"\t t2: "<<t2<<"\t f:  "<<f<<"\t Ml: "<<l.Ml<<endl;  exit(0);}
  }while(test>f);}
  //========================================================= 
  if(l.struc==2  or l.struc==3){///thick disk  &&  Stellar Halo
  t1=pow(1.4,-0.5)*pow(0.08,0.5-0.7); 
  t2=pow(mmin,-0.7);
  do{
  l.Ml=RandR(mmin, 1.4); 
  test=RandR(t1,t2); 
  if(l.Ml<0.08)          f=pow(l.Ml,-0.7);   
  else if(l.Ml>=0.08)    f=pow(l.Ml,-0.5)*pow(0.08,0.5-0.7);
  if(t1>t2 or f<t1 or f>t2 or l.Ml>1.4){cout<<"Error_3 t1:  "<<t1<<"\t t2: "<<t2<<"\t f:  "<<f<<"\t Ml: "<<l.Ml<<endl;  exit(0);}
  }while(test>f);}
  }
///#############################################################################

  l.xls=l.Dl/s.Ds;
  l.RE=sqrt(4.0*G*l.Ml*M_sun*s.Ds*KP)/velocity;
  l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
  vrel(s,l);
  l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///in day
  s.ros=double(s.Rstar*Rsun*l.xls/l.RE);
  l.mul=l.Vt*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days
  l.u0=fabs(RandR(0.0,1.0));
  if(l.u0==0.0) l.u0=1.0e-50; 
  l.tetE=  double(l.RE/AU/l.Dl);///marcs
  l.pirel= double(1.0/l.Dl - 1.0/s.Ds);//mas
  l.piE=double(l.pirel/l.tetE); 


  if(l.tE<0.0 or l.tE==0.0 or l.RE<0.0 or  l.Dl>s.Ds  or l.xls>=1.0 or s.ros<0.0 or l.Ml<0.0 or l.t0<0.0 or 
  l.Dl<=0.0 or s.Ds<=0.0 or l.xls>1.0){ 
  cout<<"BIG ERROR te: "<<l.tE<<"\t RE: "<<l.RE<<"\t V_rel: "<<l.Vt<<endl; 
  cout<<"xls:  "<<l.xls<<"\t ros:  "<<s.ros<<"\t t0:  "<<l.t0<<endl;  
  cout<<"Dl:  "<<l.Dl<<"\t Ds:  "<<s.Ds<<"\t xls:  "<<l.xls<<endl; }  
  ///int rre; cin>>rre; }   
  
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
     if(test<=Rho[0])                           {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0; age= 0.075;}
else if(test<=(Rho[0]+Rho[1]))                  {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0; age=0.575; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]))           {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;age=1.5;  }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]))    {sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2; age=2.5; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4])){sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8; age=4.0; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5])){sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4; age=6.0; }
else if(test<=maxr)                                       {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5; age=8.5; }
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}
    if(i==0) {
    sigma_R_DiskS= sigma_R_Disk;
    sigma_T_DiskS= sigma_T_Disk;
    sigma_Z_DiskS= sigma_Z_Disk;}
    if(i==1){
    sigma_R_DiskL= sigma_R_Disk;
    sigma_T_DiskL= sigma_T_Disk;
    sigma_Z_DiskL= sigma_Z_Disk; }}


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
    if(s.struc==0 or s.struc==2)  SVT =SVT+ vro_sun*(1.00762*pow(Rsc/Dsun,0.0394) + 0.00712);
    l.vs=sqrt( SVR*SVR + SVT*SVT + SVZ*SVZ );
///======================================================================================
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
   
    if (l.Vt<0.0 or l.Vt>1.0e6 or l.Vt==0.0){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;   int yee; cin>>yee;}
    //cout<<"Vt: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;
}
///==================================================================
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double RandR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s, int counter, int camera)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nn=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   double fd=1.0; ///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;///No limitation 
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;
   
   int ffd=0; 
   char filename[40];
   FILE *fill;
   if(counter%50==0){ 
   ffd=1; 
   sprintf(filename,"./files/light/lcF/%c%c%c%d%c%d%c%d.dat",'D','E','N',sector,'_',camera,'_',counter);
   fill=fopen(filename,"w");
   if(!fill){cout<<"cannot open file (density) : "<<endl;  exit(0);}}



for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;

   x=i*step;
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = Dsun-x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///M_sun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nn)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*exp(nn)*exp(-fabs(zb)/0.8)/(1.0+0.5*nn);///M_sun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(0.5/Dsun,-2.44);
   else            s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(rdi/Dsun,-2.44);///M_sun/pc^3
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
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///M_sun/pc^3
///=================================================
///در اینجا اینکه تعداد ستاره ه
///ا به این بستگی دارد که ما  نمودار قدر رنگ مطلق ستاره ها را چگونه درست کرده باشیم.اگر هیچ گونه محدودیتی برای
///درست کردن آن  در نظر نگرفته ایم،  پس تعداد کل ستاره ها را نظر  میگیریم.
/// ولی بهتر است که ما رابطه بین قدر مطلق و جرم را تعیین کنیم. در این صورت می توانیم  ورودی قدر رنگ
///وارد شده به کد را خودمان محدود به ستاره های روشن کنیم تا سرعت اجرای برنامه بالارود.
///averaged mass are the same as the previous work!!! because we did not change the besancon model


s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]);///[M_sun/pc^3]
s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[M_sun/deg^2]
s.Nstari[i]=binary_fraction*(s.rho_disk[i]*fd/0.403445+s.rho_ThD[i]*fh/0.4542+s.rho_halo[i]*fh/0.4542+s.rho_bulge[i]*fb/0.308571);////[Nt/pc^3] 
s.Nstari[i]= s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]
s.Nstart  +=s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection

 if(ffd>0)
 fprintf(fill,"%e %e %e %e %e  %e  %e\n",x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.Rostar0[i],s.Nstari[i]);
 }
   if(ffd>0)   fclose(fill);
 // cout<<"Nstart [Nt/deg^2]: "<<s.Nstart<<"\t Ro_star [Mass/deg^2]: "<<s.Rostart<<endl;
 //cout<<">>>>>>>>>>>>>>>>>>>>>>>>> END OF DISK MODLE <<<<<<<<<<<<<<<<<<<<"<<endl;
}
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
