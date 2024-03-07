from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
import numpy as np
import seaborn as sns
import pandas as pd 
import matplotlib as mpl
from scipy.stats import iqr
import pylab
rcParams["font.size"] = 18
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"


AU=1.495978707*pow(10.0,11.0);
M_sun=1.98892*pow(10.,30.0); #//in [kg].
Mjupiter=1.898*pow(10,27.0); 
Mearth=  5.9722*pow(10.0,24.0);
mmin= float(13.0*Mjupiter/M_sun); #//threshold for BD 
mmin0=float(0.000001*Mearth/M_sun);#//in M_sun  low limit for planets
KP=3.08568025*pow(10.,19); #// in meter.
G= 6.67384*pow(10.,-11.0);#// in [m^3/s^2*kg].
velocity=299792458.0;#//velosity of light

'''
vrel=40.0*1000.0
xx=0.5
Ds=0.2
Dl= xx*Ds
Ml=0.1
RE=np.sqrt(4.0*G*Ml*M_sun*Ds*KP)*np.sqrt(xx*(1.0-xx))/velocity
tetE=float(RE/AU/Dl) 
print("tE:  ",   RE/vrel/(3600.0*24.0))
print("piE:  ", (1.0/Dl-1.0/Ds)/tetE , tetE, np.sqrt(1.0/Dl-1.0/Ds)) 
input("Enter a number ")
'''


nam0=[r"$\log_{10}[t_{\rm{E}}(\rm{days})]$", r"$\log_{10}[\rm{Priority}]$", r"$D_{\rm{l}}/D_{\rm{s}}$", r"$\rm{Galactic}~\rm{Longitude} (\rm{deg})$", 
r"$\rm{Lens}~\rm{Structure}$", r"$\log_{10}[M_{\rm{l}}(M_{\odot})]$", r"$D_{\rm{l}}(\rm{kpc})$", r"$v_{\rm{l}}(\rm{km}/s)$",
r"$\rm{source}~\rm{Structure}$",r"$\rm{nums}~\rm{source}$",r"$M_{\rm{l}}(M_{\odot})$",r"$D_{\rm{s}}(\rm{kpc})$",r"$T_{\star}(\rm{k})$",r"$R_{\star}(R_{\odot})$", r"$\rm{Priority}$", r"$\rm{source}~\rm{type}$", r"$v_{\star}(\rm{km}/s)$", 
r"$M_{\star,~\rm{T}}(\rm{mag})$", r"$m_{\star,~\rm{T}}(\rm{mag})$", r"$m_{\rm{base},~\rm{T}}(\rm{mag})$", 
r"$f_{\rm{b}}$", r"$N_{\rm{b}}$", r"$\rm{Extinction}$", 
r"$t_{\rm{E}}(\rm{days})$", r"$R_{\rm{E}}(\rm{AU})$", r"$t_{0}(\rm{days})$", r"$\mu_{\rm{rel}}(\rm{mas}/yrs)$", 
r"$V_{\rm{rel}}(\rm{km}/s)$", r"$u_{0}$", r"$\tau \times 10^{8}$", r"$\rho_{\star}$", r"$\theta_{\rm{E}}(\rm{mas})$", 
r"$Flag_{\rm{detect}}$", r"$\Delta \chi^{2}$", r"$\rm{Detect}$", r"$\log_{10}[\pi_{\rm{E}}]$", r"$\pi_{\rm{E}}$", r"$\log_{10}[\mathcal{M}]$"]

nam4=[r"$\log_{10}[\rm{Priority}]$", r"$D_{\rm{s}}(\rm{kpc})$", r"$D_{\rm{l}}(\rm{kpc})$", r"$\log_{10}[M_{\rm{l}}(M_{\odot})]$", r"$m_{\rm{TESS}}(\rm{mag})$", r"$R_{\star}[R_{\odot}]$", r"$t_{\rm{E}}(\rm{days})$", r"$\log_{10}[\tau \times 10^{6}]$", r"$\log_{10}[\rho_{\star}]$", r"$\log_{10}[\pi_{\rm{E}}]$"]

##==============================================================================
nx=13
nc=38
sector=int(0)
Tobs=[27.4,2.0*27.4,3.0*27.4,4*27.4,5*27.4,6*27.4,7*27.4,8*27.4,9*27.4,10*27.4,11*27.4,12*27.4,13.0*27.4]##13
cam= [1,2,3,4,5,6,7,8,9,10,11,12,13]#13
Nstar=[125000.0,40000.0,15000.0,3125.0,3125.0,3125.0, 3125.0,3125.0, 3125.0, 3125.0, 3125.0,6000.0,5000.0]#

n=np.zeros((nx))
for i in range(nx):
    f1=open("./files/D{0:d}_{1:d}.dat".format(sector,cam[i]),"r")
    n[i]=int(sum(1 for line in f1)) 

nf=int(np.max(n))
par= np.zeros((nf,nc,nx))
pard=np.zeros((nf,nc,nx))
for i in range(nx):
    par[:int(n[i]),:,i]= np.loadtxt("./files/D{0:d}_{1:d}.dat".format( sector, int(cam[i]) ) ) 
    par[:int(n[i]),0,i]= np.log10( par[:int(n[i]),23,i] )##logtE
    par[:int(n[i]),1,i]= np.log10( par[:int(n[i]),14,i] )##logP
    par[:int(n[i]),2,i]= par[:int(n[i]),6,i]/par[:int(n[i]),11,i]##Dl/Ds 
    par[:int(n[i]),10,i]= np.power(10.0,par[:int(n[i]),5,i])##ML 
    par[:int(n[i]),35,i]= np.log10(par[:int(n[i]),36,i])##log piE
    print("***********************step:  ", i) 


ndd=15
ntot=int(np.sum(n))
optd=np.zeros(( ntot , ndd ))
diss=np.zeros(( ntot , ndd ))
num= np.zeros((ndd))
prio=np.zeros((ndd))
pmin=np.min(np.concatenate((par[:int(n[0]),1,0],par[:int(n[5]),1,5],par[:int(n[8]),1,8],par[:int(n[12]),1,12]),axis=0))*1.0000001
pmax=np.max(np.concatenate((par[:int(n[0]),1,0],par[:int(n[5]),1,5],par[:int(n[8]),1,8],par[:int(n[12]),1,12]),axis=0))
for i in range(ndd):  
    prio[i]=float(pmin+i*(pmax-pmin)/ndd/1.00)
    
################################################################################
for k in range(nx): 
    df= pd.read_csv("./files/D{0:d}_{1:d}.dat".format(sector, cam[k]),sep=" ",skipinitialspace=True,header=None,usecols=[14,11,2,5,18,13,23,29,30,35],names=nam4)
    corrM = df.corr()
    fig, ax = plt.subplots(figsize=(10, 8))
    corrM.style.background_gradient(cmap='coolwarm').set_precision(2)
    ax=sns.heatmap(corrM, annot=True, xticklabels=nam4, yticklabels=nam4, annot_kws={"size":13.5},square=True,linewidth=1.0,cbar_kws={"shrink":.99}, linecolor="k",fmt=".1f", cbar=True, vmax=1, vmin=-1, center=0.0, robust=True)
    cbar=ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=17, pad=0.0)
    plt.xticks(rotation=35,horizontalalignment='right',fontweight='light', fontsize=15)
    plt.yticks(rotation=0, horizontalalignment='right',fontweight='light', fontsize=15)
    plt.title(r"$\rm{Correlation}~\rm{Matrix}$", fontsize=17)
    fig.tight_layout()
    plt.savefig("./Figs/Histo/corr{0:d}.jpg".format(cam[k]),dpi=200)
print("**** Correlation matrix was calculated ******** ")
################################################################################
nk=np.zeros((nx))
for k in range(nx):
    ms=0.0; bd=0.0;  ffp=0.0;
    nk[k]=int(0)
    for i in range(int(n[k])):
        prior,detflag, opt, lml, Ds =par[i,1,k],  par[i,32,k], par[i,29,k], par[i,5,k], par[i,11,k]
        f=int((prior-pmin)/((pmax-pmin)/ndd/1.00))
        if(f>=0 and f<ndd): 
            optd[int(num[f]),f]= opt*10.0
            diss[int(num[f]),f]= Ds
            num[f]+=1

        if(detflag>0): 
            pard[int(nk[k]),:,k]=par[i,:,k]  
            nk[k]+=1
            if(par[i,10,k]>0.08):    ms+=1.0
            elif(par[i,10,k]>mmin):  bd+=1.0
            else:                    ffp+=1.0
    print("MS, BD, FFPs[\%]: ",k,  round(ms*100.0/nk[k],1),"\t", round(bd*100.0/nk[k],1),"\t" ,round(ffp*100.0/nk[k],1) )        
################################################################################         
Effi=np.zeros((ndd,nc,3,nx))
eft= np.zeros((nx))
min1=np.zeros(nc) 
max1=np.zeros(nc)
step=np.zeros(nc)

for j in range(nc):
    min1[j]=np.min(np.concatenate((par[:int(n[0]),j,0],par[:int(n[5]),j,5],par[:int(n[8]),j,8],par[:int(n[12]),j,12]),axis=0))
    max1[j]=np.max(np.concatenate((par[:int(n[0]),j,0],par[:int(n[5]),j,5],par[:int(n[8]),j,8],par[:int(n[12]),j,12]),axis=0))
    if(j==30): max1[j]= 150.0
    if(j==1):  max1[j]=-1.75;
    if(j==6):  max1[j]= 0.5;
    if(j==11): max1[j]= 0.45
    if(j==27): max1[j]= 170.0
    if(j==35): max1[j]= 5.0
    step[j]=float(max1[j]-min1[j])/ndd/1.0+0.0000000000003

for i in range(nx):      
    for k in range(nc):         
        for j in range(int(n[i])): 
            nn=int((par[j,k,i]-min1[k])/step[k])
            if(nn>=ndd): nn=ndd 
            if(nn>=0 and nn<ndd):  
                Effi[nn,k,1,i]+=1.0
                if(par[j,32,i]>0): Effi[nn,k,2,i]+=1.0
    for j in range (ndd):
        for k in range(nc):  
            Effi[j,k,0,i]=float(min1[k]+j*step[k])
            Effi[j,k,2,i]=float(Effi[j,k,2,i]/(Effi[j,k,1,i]+0.00001))
    for j in range(int(nk[i])):
        nn=int((pard[j,23,i]-min1[23])/step[23])
        eft[i]+=float(Effi[nn,23,2,i]/pard[j,23,i]/1.0)     
    eft[i]=float(eft[i]/nk[i]/1.0)                       
                
                
                
for i in range(nc):
    if(i==5 or i==6 or i==11 or i==14 or i==0 or i==30 or i==1 or i==2 or i==24 or i==27 or i==35 or i==18): #log Ml, Dl,Ds, Priority,log tE , rho*
        plt.clf()
        plt.cla()
        fig=plt.figure(figsize=(8,6))
        plt.step(Effi[:,i,0,0],Effi[:,i,2,0]*100.0,where='mid',c='k',ls=':',lw=2.,label=r"$T_{\rm{obs}}=27.4~\rm{days}$")
        plt.step(Effi[:,i,0,1],Effi[:,i,2,1]*100.0,where='mid',c='r',ls='-', lw=2.,label=r"$T_{\rm{obs}}=54.8~\rm{days}$")
        plt.step(Effi[:,i,0,2],Effi[:,i,2,2]*100.0,where='mid',c='b',ls='--',lw=2.,label=r"$T_{\rm{obs}}=82.2~\rm{days}$")
        #plt.step(Effi[:,i,0,3],Effi[:,i,2,3]*100.0,where='mid',c='m',ls='-.',lw=2.,label=r"$T_{\rm{obs}}=109.6~\rm{days}$")
        plt.step(Effi[:,i,0,11],Effi[:,i,2,11]*100.0,where='mid',c='g',ls='-.',lw=2.,label=r"$T_{\rm{obs}}=328.8~\rm{days}$")
        #plt.step(Effi[:,i,0,12],Effi[:,i,2,12]*100.0,where='mid',c='g',ls='-.',lw=2.,label=r"$T_{\rm{obs}}=356.2~\rm{days}$")
        plt.xlabel(str(nam0[i]),fontsize=19)
        plt.ylabel(r"$\rm{Detection}~\rm{Efficiency}[\%]$", fontsize=19)
        plt.xticks(fontsize=19, rotation=0)
        plt.yticks(fontsize=19, rotation=0)
        #plt.ylim([-0.1,102.0])
        plt.grid("True")
        plt.grid(linestyle='dashed')
        if(i==5): 
            plt.legend()
            plt.legend(loc='best',fancybox=True, shadow=True)
            plt.legend(prop={"size":20})
        fig=plt.gcf()
        fig.tight_layout()
        fig.savefig("./Figs/Histo/Effi{0:d}.jpg".format(i),dpi=200)
print("Efficieny is plotted ********************* ")   
################################################################################

for i in range(nx):
    ns=int(n[i])
    nd=int(nk[i])
    Gama=  2.0*np.mean(par[:ns,29,i])*1.0e-8*eft[i]/np.pi 
    Neven= Gama*Tobs[i]*Nstar[i] 
    print("i:  ", i)   
    print("Dl: ", np.mean(pard[:nd,6,i]) , "\t", np.round(np.percentile(pard[:nd,6 ,i],[25,75]),2))   
    print("Ds: ", np.mean(pard[:nd,11,i]), "\t", np.round(np.percentile(pard[:nd,11,i],[25,75]),2))     
    print("ML: ", np.mean(pard[:nd,10,i]), "\t", np.round(np.percentile(pard[:nd,10,i],[25,75]),2))
    print("tE: ", np.mean(pard[:nd,23,i]), "\t", np.round(np.percentile(pard[:nd,23,i],[25,75]),2))         
    print("RE: ", np.mean(pard[:nd,24,i]), "\t", np.round(np.percentile(pard[:nd,24,i],[25,75]),2))       
    print("Vt: ", np.mean(pard[:nd,27,i]), "\t", np.round(np.percentile(pard[:nd,27,i],[25,75]),2)) 
    print("Ros:", np.mean(pard[:nd,30,i]), "\t", np.round(np.percentile(pard[:nd,30,i],[25,75]),2))
    print("piE: ",np.mean(pard[:nd,36,i]), "\t", np.round(np.percentile(pard[:nd,36,i],[25,75]),2)) 
    print("Map: ",np.mean(pard[:nd,18,i]), "\t", np.round(np.percentile(pard[:nd,18,i],[25,75]),2))  
    print("mbg: ",np.mean(pard[:nd,19,i]), "\t", np.round(np.percentile(pard[:nd,19,i],[25,75]),2))
    print("fb:  ",np.mean(pard[:nd,20,i]), "\t", np.round(np.percentile(pard[:nd,20,i],[25,75]),2))
    print("u0: ", np.mean(pard[:nd,28,i]), "\t", np.round(np.percentile(pard[:nd,28,i],[25,75]),2)) 
    print("<tau*10^8>_th:  ", np.mean(par[:ns,29,i]) , "\pm" ,  np.std(par[:ns,29,i]))      
    print("<tau*10^8>_ob:  ", np.mean(pard[:nd,29,i]) ,"\pm" , np.std(pard[:nd,29,i]))     
    print("Event rate(per :  ", Gama )  
    print("<Eff(tE)/tE>:     ", round(eft[i],2)   )
    print("Number of BG stars: ", Nstar[i] )                        
    print("Number of events: ",   round(Neven,3) )
    print("Efficiency[\%]: ",     round(nd*100.0/ns,2) )
    print("**********************************************************************")
#input("enter a number ") 
################################################################################
            
plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
for i in range(ndd): 
    plt.errorbar(prio[i],np.mean(optd[:int(num[i]),i]), yerr= np.std(optd[:int(num[i]),i]),fmt='o',markersize=8.0,color='g', ecolor='lime',elinewidth=1.2,capsize=0)
plt.xlabel(r"$\log_{10}[\rm{Priority}]$",     fontsize=18)
plt.ylabel(r"$\overline{\tau} \times 10^{9}$",fontsize=18)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./Figs/Histo/OptP.jpg" , dpi=200)

################################################################################ 
   
plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
for i in range(ndd): 
    plt.errorbar(prio[i],np.mean(diss[:int(num[i]),i]), yerr= np.std(diss[:int(num[i]),i]),fmt='o',markersize=8.0,color='g', ecolor='lime',elinewidth=1.2,capsize=0)
plt.xlabel(r"$\log_{10}[\rm{Priority}]$",       fontsize=18)
plt.ylabel(r"$\overline{D_{\rm{s}}}(\rm{kpc})$",fontsize=18)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./Figs/Histo/DisP.jpg",dpi=200)


################################################################################ 
ns=int(np.min(n))
print(ns)   
for i in range(nx):
    nd=int(nk[i])
    nb=int(n[i])
    par[:nb,30,i]=np.log10(par[:nb,30,i])
    par[:nb,36,i]=np.log10(par[:nb,36,i])
    par[:nb,29,i]=np.log10(par[:nb,29,i])
    pard[:nd,30,i]=np.log10(pard[:nd,30,i])
    pard[:nd,36,i]=np.log10(pard[:nd,36,i])
    pard[:nd,29,i]=np.log10(pard[:nd,29,i])
nam0[30]=r"$\log_{10}[\rho_{\star}]$"
nam0[36]=r"$\log_{10}[\pi_{\rm{E}}]$"    
nam0[29]=r"$\log_{10}[\tau \times 10^{8}]$"    
for i in range(nc):
    plt.clf()
    plt.cla()
    fig=plt.figure(figsize=(8,6))
    ax= plt.gca()              
    plt.hist(np.concatenate((par[:int(n[0]),i,0],par[:int(n[1]),i,1],par[:int(n[2]),i,2],par[:int(n[4]),i,4],par[:int(n[8]),i,8],par[:int(n[12]),i,12]),axis=0), density=True,bins=40,histtype='bar',ec='grey',facecolor='gray',alpha=0.65,rwidth=1.0, label=r"$\rm{All}~\rm{simulated}~\rm{events}$")
    plt.hist(pard[:int(nk[0]),i,0],density=True,bins=40,histtype='step',color='k',lw=1.9,ls=':',label=r"$\rm{Detected}~\rm{events,}~T_{\rm{obs}}=27.4~\rm{days}$")
    plt.hist(pard[:int(nk[1]),i,1],density=True,bins=40,histtype='step',color='r',lw=1.9,ls='-',label=r"$\rm{Detected}~\rm{events,}~T_{\rm{obs}}=54.8~\rm{days}$")
    plt.hist(pard[:int(nk[2]),i,2],density=True,bins=40,histtype='step',color='b',lw=1.9,ls='--',label=r"$\rm{Detected}~\rm{events,}~T_{\rm{obs}}=82.2~\rm{days}$")
    plt.hist(pard[:int(nk[11]),i,11],density=True,bins=40,histtype='step',color='g',lw=1.9,ls='-.',label=r"$\rm{Detected}~\rm{events,}~T_{\rm{obs}}=328.8~\rm{days}$")
    #y_vals =ax.get_yticks()
    #ax.set_yticks(y_vals)
    #ax.set_yticklabels(['{:.2f}'.format(float(1.0*x*(1.0/ns))) for x in y_vals]) 
    #y_vals = ax.get_yticks()
    #plt.ylim([np.min(y_vals), np.max(y_vals)])
    if(i==0):  plt.xlim([-4.5,2.5])
    if(i==5):  plt.xlim([-11.8,0.4])
    if(i==30): plt.xlim([-4.0,4.0])
    if(i==36): plt.xlim([-1.0,6.5])
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=18,labelpad=0.1)
    ax.set_xlabel(str(nam0[i]),fontsize=18,labelpad=0.1)
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    #plt.grid("True")
    #plt.grid(linestyle='dashed')
    if(i==5): 
        plt.legend()
        plt.legend(loc='best',fancybox=True, shadow=True)
        plt.legend(prop={"size":18})
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./Figs/Histo/HCTL{0:d}.jpg".format(i),dpi=200)
print("****  All histo_simulated_parameters are plotted ********************")

################################################################################



plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
plt.scatter(pard[:int(nk[0]),0,0],pard[:int(nk[0]),2,0],marker='o',c='b', s=7.0,alpha=0.3, label=r"$\rm{Camera:~1}$")
plt.scatter(par[:int(n[0]),0,0],par[:int(n[0]),2,0],marker='*',c='r', s=3.0,alpha=0.3,label=r"$\rm{Camera:~1}$")
plt.xlabel(r"$log tE$", fontsize=18)
plt.ylabel(r"$x$", fontsize=18)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
plt.legend(prop={"size":16.5}, loc='best')
fig=plt.gcf()
fig.savefig("./Figs/Histo/scatter.jpg",dpi=200)



























































