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
const double RA=180.0/M_PI;
const double KP=3.08568025*pow(10.,19); // in meter.
const double G= 6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity=299792458.0;//velosity of light
const double Z_sun=float(0.0152);//#solar metalicity
const double Msun=1.98892*pow(10.,30.0); //in [kg].
const double Mjupiter=1.898*pow(10,27.0); 
const double Mearth=  5.9722*pow(10.0,24.0);
const double mmin= double(13.0*Mjupiter/Msun); //threshold for BD 
const double mmin0=double(0.01*Mearth/Msun);//in Msun  low limit for planets
const double AU=1.495978707*pow(10.0,11.0);
const double Rsun=6.9634*pow(10.0,8.0); ///solar radius [meter]
const double ddeg=0.5;
const double Avks=double(8.20922);
const int    Nco=int(24.0/ddeg);  

///=============================================================================
const int    N1=36224, N2=25000, N3=3818, N4=3500;///CMD_BESANCON, thinD, bulge, thickD, hal
const double FWHM[2]= {0.059*5.0, 21.0*3.0};//20 Sep 



const double wave[4]={0.673,0.7865,0.532,0.797};//https://en.wikipedia.org/wiki/Photometric_system  G, T, BP, RP
const double AlAv[4]= {0.791986524645539, 0.617155245862836  ,1.0386670640894,0.601810722049874};
const double sigma[4]={0.017, 0.02,0.02, 0.02};// G, T, BP, RP  Table (2) Cardeli
//x=1.0/wave    ///https://heasarc.gsfc.nasa.gov/docs/tess/the-tess-space-telescope.html
//y=x-1.82
//a=1.0+0.17699*y-0.50447*y*y-0.02427*y*y*y+0.72085*y*y*y*y+ 0.01979*y*y*y*y*y-0.7753*y*y*y*y*y*y+0.32999*y*y*y*y*y*y*y
//b=1.41338*y+2.28305*y*y+1.07233*y*y*y-5.38434*y*y*y*y-0.62251*y*y*y*y*y+5.3026*y*y*y*y*y*y-2.09002*y*y*y*y*y*y*y
//AlAv=a+b/3.1

const double thre=0.001; 
const int    NB=int(1000); 
const int nbh=int(50);  
const int Nc=int(9586);
///=============================================================================
///SECTOR 12 https://tasoc.dk/docs/release_notes/tess_sector_12_drn17_v02.pdf
///SECTOR 39 https://tasoc.dk/docs/release_notes/tess_sector_39_drn56_v01.pdf
///SECTOR 4  https://tasoc.dk/docs/release_notes/tess_sector_04_drn05_v03.pdf
///SECTOR 31 https://tasoc.dk/docs/release_notes/tess_sector_31_drn47_v02.pdf
///SECTOR 13 https://tasoc.dk/docs/release_notes/tess_sector_13_drn18_v02.pdf

const int sector=12;  ///updated 6 Feb
const int nm=24282;//  39:9293/  4: 11682/   31:9341/    13:11264

const double cadence=double(30.0/60.0/24.0);//in days 
const double detect=18.0;
const double satu=  1.0;
const double Tobs=  27.942;//days updated 3 Feb  [21/5/2019,   18/6/2019]
const double gap1=  14.046;//[21/5/2019,  04/6/2019]  OBSERVATION
const double dur1=  1.0375;//[4/6/2019,  5/6/2019]  GAP
const double gap2=  20.0;  
const double dur2=  0.0;
const int nd= int(3.0+double(Tobs-dur1-dur2)/cadence);         
const int nw=1772;///rows in WDCat.dat
///=============================================================================
struct source{
    int    nums,struc,cl;
    double Ds,TET,FI;
    double nsbl, blend, magb, Ai, Mab, Map;
    double lon, lat;
    double type, Tstar, logg, Rstar, mass, metal;
    double ros, limb;
    double Logg[Nc], Z[Nc], Tef[Nc], Limb[Nc];  
    double cdpp;  
};
struct lens{
    int num;  
    double ecen, inc, tet, tp, period, a; 
    double phi, RE, stepb;   
    double magG, magBP, magRP;  
    double ratio, q, MBH, RBH, ext, Map, Mab;
    double mass[nw], radius[nw], logg[nw], teff[nw], G[nw], BP[nw], RP[nw], dist[nw];
};
struct CMD{
    double Teff_d[N1+N3],logg_d[N1+N3],Mab_d[N1+N3],Rs_d[N1+N3],mass_d[N1+N3],type_d[N1+N3],metal_d[N1+N3];int cl_d[N1+N3];
    double Teff_b[N2],logg_b[N2],Mab_b[N2],Rs_b[N2],mass_b[N2],type_b[N2],metal_b[N2]; int cl_b[N2];  /// bulge
    double Teff_h[N4],logg_h[N4],Mab_h[N4],Rs_h[N4],mass_h[N4],type_h[N4],metal_h[N4]; int cl_h[N4];  /// halo
};
struct detection{
    int    numt; 
    double dchi;
    double magn[nd];     
};
struct extinc{
   double dis[100];///distance
   double Extks[100];///ks-band extinction
   double Aks;
};
struct tessfile{
   double Map[nm],blend[nm],teff[nm],radius[nm],mh[nm],logg[nm],ra[nm],dec[nm];
   double cdpp[nm];   
};
struct coordinate{
   double elong[Nco*Nco], elat[Nco*Nco];  
   double    ra[Nco*Nco],  dec[Nco*Nco];  
   double     l[Nco*Nco],    b[Nco*Nco]; 
   double Nstar[Nco*Nco];
};
///=============================================================================
int  Extinction(extinc & ex,source & s);
void read_cmd(CMD & cm);
void func_source(source & s,CMD & cm, extinc & ex, tessfile & t, double , double);
void func_lens(lens & l, source & s,extinc & ex);
double metric(detection & d); 
double ErrorTESS(double maga); 
double Interpol(double ds, extinc & ex);
double RandN(double sigma,double);
double RandR(double down, double up);
double Fluxlimb(double limb, double rstar);  
double Kepler(double phi, double ecen); 
double Bessel(int n,double x); 
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
      
   
   VBBinaryLensing vbb;
   vbb.Tol=1.e-4;
   vbb.a1 =0.0;  
   vbb.LoadESPLTable("./ESPL.tbl");
  
  
   source s;
   lens l;
   CMD cm;
   extinc ex;
   detection d; 
   coordinate co;
   tessfile t;
   read_cmd(cm);
     
     
   FILE* film;  
   FILE* fild; 
   FILE* limbt;  
   FILE* distr;
   FILE* tefile;
   FILE* coord; 
   FILE* magdi; 
   FILE* wdcat;    
  
    /*  
   if(sector==39){
   d.cadence=double(10.0/60.0/24.0);//in days 
   d.detect=18.0;
   d.satu=  1.0;
   d.Tobs=     27.9514;//days   
   d.gap1=     13.09;//[27/5/2021, 9/6/2021] OBSERVATION
   d.dur1=     1.00; //[9/6/2021, 10/6/2021] GAP
   d.gap2=     20.0;   
   d.dur2=     0.0;}    
        
   if(sector==4){
   d.cadence=double(30.0/60.0/24.0);//in days 
   d.detect=18.0;
   d.satu=  0.0;
   d.Tobs=     25.95;//days
   d.gap1=     7.6375;//[19/10/2018 - 27/10/2018]          OBSER 
   d.dur1=     2.6743; //[27/10/2018- 29/10/2018]   GAP
   d.gap2=     d.gap1+ d.dur1 + 2.2973;//[29/10/2018,  1/11/2018]     OBSER
   d.dur2=     1.0403;}//[1/11/2018,   2/11/2018]  G
  
   if(sector==31){
   d.cadence=double(10.0/60.0/24.0);//in days 
   d.detect=18.0;
   d.satu=  2.0;
   d.Tobs=     25.43473;//days [22/10/202,   16/11/2020]
   d.gap1=     12.94445;//[22/10/2020,   03/11/2020]
   d.dur1=     1.4028; //[03/11/2020,  05/11/2020]
   d.gap2=     20.0;   
   d.dur2=     0.0;}
       
   if(sector==13){
   d.cadence=double(30.0/60.0/24.0);//in days 
   d.detect=17.0;
   d.satu=  1.0;
   d.Tobs=     28.4417;//days
   d.gap1=     13.775;//
   d.dur1=     0.9292;
   d.gap2=     20.0;   
   d.dur2=     0.0;}    */
  
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   int    num1, num2, numm, ddo, ntran; 
   char   filenam0[40],filenam1[40], filenam2[40], filenam3[40], filenam4[40], filenam5[40], filenam7[40];
   int    nsim, ndet, los=0, countb, ffd, detectf;
   double fte, meter; 
   double lonn, deltaA;   
   double chi2,chi1, emt, timee; 
   double ksi, x0, y0, x1, y1, z1, dis, yc, zc; 
   double phase, RE, u, us, Astar, frac, frac0, flux, maga, Occul;   
   double yb,zb, zlim, rstar, maga2;  
   double dt=double(cadence/3.0);  
   double Amin, Amax=0.0;  
   double  depth, snr; 
   double error;
   double Nsim=0.0, Ndet=0.0; 
   double paral;  
   long int id; 
///=============================================================================
    
   sprintf(filenam5,"./files/%c%c%c%c%c%c%d.dat",'C','o','o','r','S','B',sector);
   coord=fopen(filenam5,"r");
   if(!coord){ cout<<"cannot read CoordSB.dat "<<sector<<endl; exit(0);}
   for(int i=0; i<int(Nco*Nco); ++i){
   fscanf(coord,"%d %d %lf %lf %lf %lf %lf %lf %lf\n",
   &num1,&num2,&co.elat[i],&co.elong[i],&co.b[i],&co.l[i],&co.dec[i],&co.ra[i],&co.Nstar[i]);}
   fclose(coord);
   
       
   wdcat=fopen("./files/WDCat.dat","r");
   if(!wdcat){ cout<<"cannot read wdcat.dat "<<sector<<endl; exit(0);}
   for(int i=0; i<nw; ++i){// ID, Teff, logg, mass, parallax, G, BP,  RP, radius 
   fscanf(wdcat,"%d  %lf  %lf %lf %lf %lf %lf %lf %lf\n",
   &id, &l.teff[i], &l.logg[i], &l.mass[i], &paral, &l.G[i], &l.BP[i], &l.RP[i], &l.radius[i]);  
   l.G[i] -= 5.0*log10(100.0/paral);  
   l.BP[i] -= 5.0*log10(100.0/paral);  
   l.RP[i] -= 5.0*log10(100.0/paral);
   l.dist[i]=double(1.0/paral); //kpc 
   }
   fclose(wdcat); 
   cout<<"**** File WDcat.dat was read ****"<<endl;    
      
    
      

   sprintf(filenam7,"./files/%c%c%c%d%c%c%c%c%c.txt",'s','e','c',sector,'_','p','a','r','B');
   magdi=fopen(filenam7,"r"); 
   if(!magdi){cout<<"cannot read sec12_parB.txt "<<sector<<endl; exit(0);}
   for(int i=0; i<nm; ++i){
   fscanf(magdi,"%d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
   &numm,&t.Map[i],&t.blend[i],&t.teff[i],&t.radius[i],&t.mh[i],&t.logg[i],&t.ra[i],&t.dec[i],&t.cdpp[i]);}
   fclose(magdi); 
   cout<<"file Information was read !!!!"<<endl; 
         
         
         
  
   limbt=fopen("./files/table24bClaret.dat","r");
   for(int i=0; i<Nc; ++i){
   fscanf(limbt,"%lf  %lf  %lf %lf\n",&s.Tef[i], &s.Logg[i], &s.Z[i], &s.Limb[i]);}
   fclose(limbt); 
   cout<<"**** File table24bClaret.dat was read ****"<<endl;    
   
     
               
      
///=============================================================================
   for(int bn1=0; bn1<int(Nco); ++bn1){  
   for(int bn2=0; bn2<int(Nco); ++bn2){  
   s.lat=double(co.b[los]);  
   lonn= double(co.l[los]);
   if(lonn<=0.0) s.lon=360.0+lonn;
   else          s.lon=lonn;
   s.TET=double(360.0-s.lon)/RA;///radian
   s.FI =double(s.lat/RA);///radian
   if(Extinction(ex,s)==1){
   
   
    
   sprintf(filenam0,"./files/light/lcF2/%c%d%c%d.dat",'D',sector,'_',los);
   distr=fopen(filenam0,"w"); 
   nsim=0; ndet=0; 

 

   do{
   func_source(s,cm,ex, t, double(co.ra[los]), double(co.dec[los]) );
   func_lens(l, s, ex);   
   nsim+=1;
   Nsim+=1.0; 
   
   if(double(Tobs/l.period)>1.0){
   cout<<">>>>Tobs:  "<<Tobs<<"\t period:  "<<l.period<<"\t Tobs/period:  "<<Tobs/l.period<<endl;


   ffd=1;
   sprintf(filenam1,"./files/light/lcF2/%c%d%c%d%c%d.dat",'L',sector,'_',los,'_',nsim);
   fild=fopen(filenam1,"w");
   sprintf(filenam2,"./files/light/lcF2/%c%d%c%d%c%d.dat",'M',sector,'_',los,'_',nsim);
   film=fopen(filenam2,"w");
   l.inc=double(l.inc/RA);
   l.tet=double(l.tet/RA);
  

   Amin=1000.0;  Amax=0.0;     error=0.0;      
   d.numt=0;  chi1=0.0; chi2=0.0;  timee=0.0;  
   for(double tim=0.0; tim<=Tobs; tim=tim+dt){
   l.phi=double((tim-l.tp)*2.0*M_PI/l.period); 
   if(l.ecen<0.01) ksi=l.phi;
   else            ksi=Kepler(l.phi , l.ecen);
   x0=l.a*(cos(ksi)-l.ecen);
   y0=l.a*sin(ksi)*sqrt(1.0-l.ecen*l.ecen); 
   y1=              y0*cos(l.tet)+x0*sin(l.tet);
   x1= cos(l.inc)*(-y0*sin(l.tet)+x0*cos(l.tet));
   z1=-sin(l.inc)*(-y0*sin(l.tet)+x0*cos(l.tet));
   dis=sqrt(x1*x1+y1*y1+z1*z1)+1.0e-50; 
   phase=acos(-x1/dis);
   //cout<<"tim:  "<<tim<<"\t phi:  "<<l.phi*RA<<"\t ksi:  "<<ksi*RA<<endl;
   //cout<<"Orbital componenets:  x0/a:  "<<x0/l.a<<"\t y0/a:  "<<y0/l.a<<endl;
   //cout<<"x1/a:  "<<x1/l.a<<"\t y1/a:  "<<y1/l.a<<"\t z1/a:  "<<z1/l.a<<endl;
   //cout<<"distance/a:  "<<dis/l.a<<"\t phase(deg):  "<<phase*RA<<endl;
  
   l.RE=sqrt(4.0*G*Msun*l.MBH)*sqrt(fabs(x1)*s.Ds*KP/(s.Ds*KP+fabs(x1)))/velocity+1.0e-50;
   s.ros=fabs(s.Rstar*Rsun*s.Ds*KP/(s.Ds*KP+fabs(x1))/l.RE); 
   u= sqrt(y1*y1+z1*z1)/l.RE;
   if(sqrt(y1*y1+z1*z1)<double(l.RBH)) us=0.0; 
   else us=fabs(sqrt(y1*y1 + z1*z1)-l.RBH)/(s.Rstar*Rsun);  
   if((x1<0.0 and fabs(phase)>M_PI/2.0) or l.RE<=0.0 or s.ros<=0.0 or dis<=0.0 or x1>dis){ 
   cout<<"x1:  "<<x1<<"\t phase: "<<phase<<"\t RE:  "<<RE<<"\t limb:  "<<s.limb<<endl;int uue; cin>>uue; }
   //cout<<"RE/AU:  "<<l.RE/AU<<"\t ros:  "<<s.ros<<"\t u:  "<<u<<"\t us:  "<<us<<endl;

   

   vbb.a1=s.limb; 
   Astar=1.0;     
   Occul=1.0;   
   
   
   if(x1<-0.001){//Self-lensing Check 
   if(s.ros>100.0){  
   if(u<s.ros) Astar=double(1.0+2.0/s.ros/s.ros);
   else        Astar=double(u*u+2.0)/sqrt(u*u*(u*u+4.0));}
   else        Astar=vbb.ESPLMag2(u, s.ros);}
   if(x1>0.0 and l.MBH<2.17){//occultation of the WD brightness  
   if(us<=1.0){
   if(fabs(sqrt(y1*y1 + z1*z1)+l.RBH)<double(s.Rstar*Rsun) ){frac=frac0=1.0;}
   else{   
   frac=0.0;  frac0=0.0;    
   for(int i=0; i<nbh; ++i){
   for(int j=0; j<nbh; ++j){
   yb=double(-l.RBH + i*l.stepb);   
   zb=double(-l.RBH + j*l.stepb); 
   zlim=sqrt(l.RBH*l.RBH - yb*yb); 
   if(fabs(zb)<=zlim){
   yc=y1 + yb;    
   zc=z1 + zb;  
   rstar=sqrt(yc*yc+zc*zc)/(s.Rstar*Rsun); 
   frac0+=1.0; 
   if(rstar<=1.0) frac+=1.0; //fabs(Fluxlimb(s.limb,rstar));
   }}}}   
   //frac0= M_PI*s.Rstar*Rsun*s.Rstar*Rsun*(1.0-s.limb/3.0)/(l.stepb*l.stepb); 
   Occul=double(1.0-frac/frac0);}}
   

   
   
   Astar=(Astar+l.ratio*Occul)/(1.0+l.ratio);  
   Astar=double(Astar*s.blend+1.0-s.blend);
   if(Astar<Amin)  Amin=Astar;   
   if(Astar>Amax)  Amax=Astar;
   maga=s.magb-2.5*log10(Astar);
   if(ffd>0) fprintf(film,"%.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf %.4lf  %.4lf\n",
   tim, maga, Astar,double(x1/l.a),double(y1/l.a),double(z1/l.a),double(phase*RA),double(l.RE/Rsun),s.ros,u,us);
   
   
   
   if(tim>0 and tim<=Tobs and float((tim-gap1)*(tim-gap1-dur1))>=0.0 and float((tim-gap2)*(tim-gap2-dur2))>=0.0){
   timee+=dt; 
   if(timee>=cadence){
   timee-=cadence; 
   if(maga>=satu and maga<=detect){
   emt= fabs(ErrorTESS(maga));
   deltaA=fabs(pow(10.0,-0.4*emt)-1.0)*Astar; 
   maga2=maga+RandN(emt, 1.5);
   error += emt;  
   d.magn[d.numt]=maga2; 
   chi1+= pow((maga2-s.magb)/emt,2.0);
   chi2+= pow((maga2-  maga)/emt,2.0);
   if(ffd>0) fprintf(fild,"%.4lf %.5lf %.7lf\n", tim, maga2, emt);//, Astar+RandN(deltaA,1.5), deltaA);
   d.numt+=1;
   if(d.numt>=nd){cout<<"Error  numt:  "<<d.numt<<"\t Nes:  "<<nd<<endl;  exit(0);}
   }}}
   }//time loop 
   d.dchi= fabs(chi2-chi1)/d.numt;
   if(fabs(Amax-1.0)>fabs(1.0-Amin)) depth=fabs(Amax-1.0); 
   else                              depth=fabs(1.0-Amin);   
   ntran= int(Tobs/l.period);  
   snr=  double(sqrt(ntran*1.0) *depth*1000000.0 / s.cdpp);    
   error=double(error/d.numt);
   if(snr>5.0){ Ndet+=1.0;   ndet+=1; }  
   if(ffd>0){fclose(film); fclose(fild);}
   meter= metric(d); 
   fprintf(distr,
   "%d   %d   %.5lf  %.5lf  "///4
   "%.8lf  %.7lf   %.5lf   %.5lf   %.6lf  %.5lf  %.6lf  %.5lf  "//12
   "%d   %.2lf  %.5lf  %.6lf  %.2lf  %.5lf  %.5lf  %.5lf  %.5lf  "//21
   "%.5lf  %.5lf  %.5lf  %.5lf   %.2lf   %.4lf  "  //27
   "%.1lf  %d  %.5lf   %.7lf  %.7lf  %.7lf  %.7lf   %.4lf   %.4lf  %.8lf  %.4lf  %.2lf  %d  "//40
   "%.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf\n", //46
   los, nsim, s.lat, lonn, //4
   l.MBH, l.RBH/Rsun, l.inc*RA, l.tet*RA, l.ecen, l.period,  log10(l.a/(s.Rstar*Rsun)), l.tp,//12
   s.cl, s.type, s.Ds, s.mass, s.Tstar, s.Rstar, s.logg, s.metal, s.limb, //21
   s.Mab, s.Map, s.magb, s.blend, s.nsbl, s.Ai,//27
   d.dchi, d.numt, meter, error, snr, depth, s.cdpp, l.Map, l.Mab, l.ratio, l.logg[l.num], l.teff[l.num], l.num,//40
   l.magG, l.magBP, l.magRP,l.G[l.num], l.BP[l.num], l.RP[l.num] );//46
   
   cout<<"=============================================================="<<endl;
   cout<<"*** l.Mab:  "<<l.Mab<<"\t l.Map:  "<<l.Map<<"\t ratio:  "<<l.ratio<<endl;
   cout<<"latitude: "<<s.lat<<"\t longtitude: "<<s.lon<<"\t los:  "<<los<<endl;
   cout<<"nsim:      "<<nsim<<"\t l.MBH: "<<l.MBH<<"\t RBH: "<<l.RBH<<endl;
   cout<<"Map:  "<<s.Map<<"\t magb:  "<<s.magb<<"\t l.period:  "<<l.period<<endl;
   cout<<"dchi:  "<<d.dchi<<"\t d.numt: "<<d.numt<<"\t metric:  "<<meter<<endl;
   cout<<"SNR:  "<<snr<<"\t N_transit:   "<<ntran<<"\t error:  "<<error<<endl;
   cout<<"Amax:  "<<Amax<<"\t Amin:  "<<Amin<<"\t depth :    "<<depth<<endl;
   cout<<"cdpp:  "<<s.cdpp<<endl;
   cout<<"Nsim:  "<<Nsim<<"\t Ndet:  "<<Ndet<<"\t fraction:  "<<Ndet*100.0/Nsim<<endl;
   cout<<"=============================================================="<<endl;
   }
   }while(ndet<10);
   fclose(distr);}//end of extinction
   else{cout<<"There is no Extinction file !!!"<<endl;  exit(0);}

    los+=1;}}
    fclose(_randStream);
    return(0);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double Kepler(double phi, double ecen){  
    double ksi=0;   
    double term, term0;  
    phi=double(phi*RA); 
    while(phi>360.0) phi=phi-360.0; 
    while(phi<0.0)   phi=phi+360.0;      
    if(phi>180)      phi=double(phi-360.0);
    if(phi<-181.0 or phi>181.0){ 
    cout<<"Error :  Phi:  "<<phi<<"\t ecent:  "<<ecen<<endl;   int yye;  cin>>yye;}
    phi=double(phi/RA);
    ksi=phi; 
    for(int i=1; i<NB; ++i){
    term= Bessel(i,i*ecen)*sin(i*phi)*2.0/i;  
    ksi+=term; 
    if(i==1)  term0=fabs(term); 
    if(fabs(term)<double(thre*term0) and i>5)  break;}        
    return(ksi); 
}    
///==========================================================================//
///                                                                          //
///                   metric                                                 //
///                                                                          //
///==========================================================================//
double metric( detection & d)
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
    cout<<"mean:  "<<mean<<"\t dchi:  "<<d.dchi<<endl;
    exit(0);}
    cout<<"meter:  "<<meter<<endl;
    return(meter);    
}
///#############################################################################
double ErrorTESS(double maga){//checked 6 Feb  2024
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
   
   if(sector==13){
   if(maga<8.0)                emt=double(0.2*maga -5.64); 
   if(maga>=8.0 and maga<11.5) emt=double(0.242*maga-5.976);
   if(maga>=11.5)              emt=double(0.3*maga-6.643); }
   
   emt=emt+RandN(0.1,3.0);
   if(emt<-5.0)   emt=-5.0;  
   emt=pow(10.0,emt);

   if(emt<0.00001 or emt>0.5 or maga<0.0){
   cout<<"Error emt:  "<<emt<<"\t maga:  "<<maga<<endl;}
   return(emt); 
}
///#############################################################################
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
///#############################################################################
double RandN(double sigma, double nn){
   double rr,f,frand;
   do{
   rr=RandR(-sigma*nn , sigma*nn); ///[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/(sigma*sigma));
   frand=RandR(0.0 , 1.0);
   }while(frand>f);
   return(rr);
}
///#############################################################################
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm)
{
    int yye;
    double age, lumi, MB, MV, MK, MI, G, T, GBP, GRP;
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
    
    G=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428- 0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                   +RandN(0.04474, 1.5);  
    if(GBP-GRP<6.0 and GBP-GRP>-1.0) 
    T= G-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473+RandN(0.006, 1.5); 
    else T=G-0.430 + RandN(0.6,1.5);
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
    
    G=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428-0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                    +RandN(0.04474, 1.5);  
    if(GBP-GRP<6.0 and GBP-GRP>-1.0) 
    T= G-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473 +RandN(0.006, 1.5); 
    else T=G-0.430;
    
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
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;  */



////=================================== STELLAR HALO ===========================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','h','W','S');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDhW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.Teff_h[j],&age,&lumi,&cm.logg_h[j],&cm.metal_h[j],&cm.Rs_h[j],&MB,&MV,&MI,&MK,&cm.cl_h[j],&cm.type_h[j]);
    
    G=  MV - 0.01746+0.008092*(MV-MI)-0.2810*pow(MV-MI,2.0)+0.03655*pow(MV-MI,3.0) +RandN(0.04670, 1.5);  
    GBP=MV - 0.05204+  0.4830*(MV-MI)-0.2001*pow(MV-MI,2.0)+0.02186*pow(MV-MI,3.0) +RandN(0.04483, 1.5);  
    GRP=MV +0.0002428-0.8675*(MV-MI)-0.02866*pow(MV-MI,2.0)+0.0                    +RandN(0.04474, 1.5);  
    if(GBP-GRP<6.0 and GBP-GRP>-1.0) 
    T= G-0.00522555*pow(GBP-GRP,3.0)+0.0891337*pow(GBP-GRP,2.0)-0.633923*(GBP-GRP)+0.0324473 +RandN(0.006, 1.5); 
    else T=G-0.430;
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
///#############################################################################
///==============================================================//
///                                                              //
///                  Func source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD & cm, extinc & ex, tessfile & t, double ral, double decl)
{
    int    struc,nums,num, num2=0;
    double rf,Ds,Av, dist, dmin, module, dista;
    s.magb=0.0; s.Ai=0.0; s.Map=0.0; s.Mab=0.0;
    s.nsbl=0.0;  s.blend=0.0;  
    //cout<<"func_source star "<<endl; 
 
 
   
    do{
    num=int(RandR(0.0,nm-1.0));
    num2+=1; 
    if(num2>nm) break;  
    }while(fabs(t.ra[num]-ral)>3.0  or fabs(t.dec[num]-decl)>3.0);
    //cout<<"ral:  "<<ral<<"\t decl:      "<<decl<<"\t counter: "<<num<<endl; 
    //cout<<"tra:  "<<t.ra[num]<<"\t tdec:    "<<t.dec[num]<<endl;
    
    s.Map  = t.Map[num]; 
    s.blend= t.blend[num];  
    s.Tstar= t.teff[num];  
    s.Rstar= t.radius[num]; 
    s.logg = t.logg[num];  
    s.cdpp=  t.cdpp[num];  
    if(t.mh[num]>-9.5) s.metal= pow(10.0,t.mh[num])*Z_sun; 
    else               s.metal=-10.0;  

    
    if(s.Tstar<0.0  or s.logg<0.0 or s.blend<0.0 or  s.blend>1.0 or s.Rstar<0.0){
    cout<<"Tstar:  "<<s.Tstar<<"\tlogg:  "<<s.logg<<"\t blend:  "<<s.blend<<endl;
    cout<<"Rstar:  "<<s.Rstar<<"\t ral:  "<<ral<<"\t decl:  "<<decl<<endl; 
    cout<<"ra[num]:  "<<t.ra[num]<<"\t dec[num]:  "<<t.dec[num]<<endl;  int eei; cin>>eei;}
    
    
    num=-1; 
    if(s.Tstar>cm.Teff_d[N1+N3-1])  num=int(N1+N3-1); 
    else if(s.Tstar<cm.Teff_d[0] )  num=0;  
    else{
    for(int j=1; j<int(N1+N3); ++j){
    if(float((s.Tstar-cm.Teff_d[j])*(s.Tstar-cm.Teff_d[j-1]))<=0.0){
    num=j-1; break; j=int(N1+N3);}}}
    
    if(num<0 or fabs(s.Tstar-cm.Teff_d[num])>2000.0){
    cout<<"There is an error num<0:  "<<num<<"\t Tstar:  "<<s.Tstar<<"\t teffcm:  "<<cm.Teff_d[num]<<endl;  
    int yye;  cin>>yye;  }

    
    dmin=100000; dist=0.0;nums=-1;  
    for(int j=-20; j<20; ++j){
    num2=int(num+j);
    if(num2>=0 and num2<int(N1+N3) ){
    if(s.metal>-9.5) dist=sqrt(pow(s.logg-cm.logg_d[num2],2.0)+pow(s.Rstar-cm.Rs_d[num2],2.0)+pow(s.metal-cm.metal_d[num2],2.0));
    else             dist=sqrt(pow(s.logg-cm.logg_d[num2],2.0)+pow(s.Rstar-cm.Rs_d[num2],2.0));
    if(dist<dmin){dmin=dist;  nums=num2;}}}
    s.mass= cm.mass_d[nums];
    s.cl=   cm.cl_d[nums]; 
    s.type= cm.type_d[nums];
    s.Mab=  cm.Mab_d[nums];
    dista=pow(10.0,0.2*fabs(s.Map-s.Mab))/100.0; 
    //cout<<"Tstar:  "<<s.Tstar<<"\t Teff_d[num]:  "<<cm.Teff_d[num]<<endl;
    //cout<<"logg:  "<<s.logg<<"\t logg_d:  "<<cm.logg_d[nums]<<endl;
    //cout<<"Rstar:  "<<s.Rstar<<"\t Rs_d:  "<<cm.Rs_d[nums]<<endl;

//##############################################################################
    num=-1; nums=-1;  
    if(s.Tstar>s.Tef[Nc-1])  num=int(Nc-1); 
    else if(s.Tstar<s.Tef[0])  num=0;  
    else{
    for(int j=1; j<int(Nc); ++j){
    if(float((s.Tstar-s.Tef[j])*(s.Tstar-s.Tef[j-1]))<=0.0){
    num=j-1; break; j=int(Nc);}}}
    if(num<0){cout<<"There is an error num<0:  "<<num<<"\t Tstar:  "<<s.Tstar<<endl;  int yye;  cin>>yye;  }
    //cout<<"Tstar:  "<<s.Tstar<<"\t s.Tef[num]:  "<<s.Tef[num]<<endl;
    
    
    
    dmin=100000; dist=0.0;
    for(int j=-10; j<10; ++j){
    num2=int(num+j);
    if(num2>=0 and num2<int(Nc)){
    if(s.metal>-9.5) dist=sqrt(pow(s.logg-s.Logg[num2],2.0)+pow(s.metal-s.Z[num2],2.0));
    else             dist=fabs(s.logg-s.Logg[num2]);
    if(dist<dmin){dmin=dist; nums=num2;}}}
    //cout<<"s.logg:  "<<s.logg<<"\t logg[num2]:  "<<s.Logg[num2]<<endl;
    s.limb=s.Limb[nums];  
    
    
    

//##############################################################################
    dmin=1000.0; dist=0.0; s.Ds=0;  
    for(Ds=0.01; Ds<=dista; Ds=Ds+0.005){
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
   
   

    
    //s.Fluxb=fabs(pow(10.0,-0.4*s.Map)/s.blend);//T-band  baseline flux
    s.magb =s.Map+2.5*log10(s.blend);//T-band baseline magnitude
    s.nsbl =double(1.0/s.blend);
   
    
    //cout<<"*******************************************"<<endl;
    cout<<"Ds:     "<<s.Ds<<"\t   Map:  "<<s.Map<<endl; 
    cout<<"Mab:    "<<s.Mab<<"\t  Ai:  "<<s.Ai<<endl;
    cout<<"Tstar:  "<<s.Tstar<<"\tnsbl:  "<<s.nsbl<<endl;
    cout<<"Mass:   "<<s.mass<<"\t Rstar:  "<<s.Rstar<<endl; 
    cout<<"type:   "<<s.type<<"\t cl:  "<<s.cl<<endl; 
    cout<<"logg:   "<<s.logg<<"\t metal:  "<<s.metal<<endl;
    cout<<"Tstar:  "<<s.Tstar<<"\t logg:  "<<s.logg<<"\t metal:  "<<s.metal<<endl;
    cout<<"limb:  "<<s.limb<<"\t mass:  "<<s.mass<<"\t type:  "<<s.type<<endl;
    cout<<"num:  "<<num<<"\t nums:  "<<nums<<"\t cl:  "<<s.cl<<endl;
    cout<<"Mab:  "<<s.Mab<<"\t Map:  "<<s.Map<<"\t dista:  "<<dista<<endl;
    //cout<<"*******************************************"<<endl; 
    cout<<" end of  func_source!!!! "<<endl;
    //int yye;    cin>>yye;  
}
///#############################################################################
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
     if(Lat>10.0)    Lat=10.0; 
     if(Lat<-10.0)   Lat=-10.0;
     if(Lon>100.0 and Lon<260.0) Lon=100.0; 
     if(fabs(Lon)<0.24999999)    Lon=360.00;
     cout<<"Lon:    "<<Lon<<"\t     Lat:  "<<Lat<<endl;
     if(Lon>360.000 or Lon<0.25 or fabs(Lat)>10.0 or (Lon>100 and Lon<260)){
     cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}
     //cout<<"(2) Lon:    "<<Lon<<"\t     Lat:  "<<Lat<<endl;

     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");


     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     FILE *SD;
     SD=fopen("./files/Ext/saved_direction.txt","r");
     for(int i=0; i<64881; ++i) {
     fscanf(SD,"%lf %lf \n",&latit,&lonti);
     if(fabs(Lat-latit)<0.1 and fabs(Lon-lonti)<0.1){
     cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
     cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
     if(ex.dis[i]<0.2  or ex.dis[i]>50.0 or ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0; }
     }}
     //cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<endl;
     fclose(fpd);
     return(flag);
}
///==============================================================//
///                                                              //
///                  func lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s, extinc & ex)
{  
  
  double emax, Roche, color; 
  double Av, extG, extBP, extRP;  
  l.inc=RandR(0.0, 2.0);
  l.tet=RandR(0.1,359.9);
  
 
  do{
  l.num=int(RandR(0.0,nw-1.0));
  l.MBH=l.mass[l.num];
  l.a=RandR(log10(3.0*s.Rstar),log10(1000000.0*s.Rstar));
  l.a=pow(10.0,l.a)*Rsun;//meter 
  l.period=sqrt(4.0*M_PI*M_PI/(G*Msun*(s.mass+l.MBH)))*pow(l.a,1.5)/(3600.0*24.0);//days
  l.q=double(l.MBH/s.mass);
  Roche=l.a*0.49*pow(l.q,2.0/3.0)/(0.6*pow(l.q,2.0/3.0)+log(1.0+pow(l.q,1.0/3.0)));
  cout<<"Roche:  "<<Roche/s.Rstar/Rsun<<"\t q:  "<<l.q<<endl; 
  }while(Roche<double(s.Rstar*Rsun) );  
  
 
 
  emax=double(0.8-8.0*exp(-pow(6.0*l.period,0.35)));//Fig 7
  if(emax<0.0)  emax=0.0; 
  if(emax>1.0)  emax=1.0;
  if(emax<0.0 or emax>1.0){cout<<"Error emax:  "<<emax<<"\t period:  "<<l.period<<endl; int yyw; cin>> yyw; } 
  l.ecen=RandR(0.0, emax); 
  
  
  l.tp= RandR(0.0,l.period); 
  l.RBH=l.radius[l.num]*Rsun;//meter 
  l.stepb=double(2.0*l.RBH/nbh/1.0); 


  ex.Aks=Interpol(l.dist[l.num],ex);
  Av=double(ex.Aks*Avks);
  if(Av<0.0) Av=0.0;
  extG= fabs(Av*AlAv[0])+RandN(sigma[0],1.5);
  extBP=fabs(Av*AlAv[2])+RandN(sigma[2],1.5);
  extRP=fabs(Av*AlAv[3])+RandN(sigma[3],1.5);
  if(extG<0.0)  extG=0.0;
  if(extBP<0.0) extBP=0.0;
  if(extRP<0.0) extRP=0.0;
  l.magG=  l.G[l.num] -extG;  
  l.magBP= l.BP[l.num]-extBP;  
  l.magRP= l.RP[l.num]-extRP;  
  

  
  color=double(l.magBP-l.magRP); 
  if(color<=6.0 and color>=-1.0)
       l.Mab=l.magG-0.00522555*pow(color,3.0)+0.0891337*pow(color,2.0)-0.633923*color+0.0324473+RandN(0.006,1.5); 
  else l.Mab=l.magG-0.430 + RandN(0.6,1.5);
 
  l.Map=l.Mab + s.Ai + 5.0*log10(s.Ds*100.0);  
  l.ratio=pow(10.0,-0.4*fabs(l.Map-s.Map)); 
  
   
  cout<<"extG:  "<<extG<<"\t extBP:  "<<extBP<<"\t extRP:  "<<extRP<<endl; 
  cout<<"Map(WD):  "<<l.Map<<"\t Map(source):  "<<s.Map<<"\t l.ratio:  "<<l.ratio<<endl;
  cout<<"inc:  "<<l.inc<<"\t tet:  "<<l.tet<<"\t MBH:  "<<l.MBH<<endl; 
  cout<<"period:  "<<l.period<<"\t semi:  "<<l.a/Rsun/s.Rstar<<"\t ecen_max:  "<<emax<<endl;
  cout<<"ecen:  "<<l.ecen<<"\t RBH(km):  "<<l.RBH/1000<<"\t tp:   "<<l.tp<<endl;
  cout<<"The end of func_lens **"<<endl;
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double RandR(double down, double up){
   double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
   return(p*(up-down)+down);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double Fluxlimb(double limb, double rstar){
    return ( double(1.0-limb*(1.0-sqrt(fabs(1.0-rstar*rstar)))) );
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double Bessel(int n,double x)
{
    double j1=0.00000001,tet;
    int kmax=10000;
    for(int k=0; k<kmax; ++k){
    tet=double(k*M_PI/kmax);
    j1+=double(M_PI/kmax/1.0)*cos(n*tet-x*sin(tet)); }
    return(j1/M_PI);
}    
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
