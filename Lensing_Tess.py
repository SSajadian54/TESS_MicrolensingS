from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
import numpy as np
import matplotlib as mpl
import pylab
rcParams["font.size"] = 18
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
##==============================================================================
sector=int(13)
if(sector==12): Tobs=27.942
if(sector==39): Tobs=27.9514
if(sector==4):  Tobs=25.95
if(sector==31): Tobs=25.43473
if(sector==13): Tobs=28.4417

ddeg=0.5
num=int(24.0/ddeg) 
Nlos=4###int(num*num)
nc= 38

###########################################################
N0=100000
no= int(N0/Nlos)
pars=np.zeros((N0,nc))
pard=np.zeros((N0,nc))
k1=0;  k2=0
############################################################


for los0 in range(Nlos): 
    los=los0+0
    f1=open("./files/D{0:d}_{1:d}.dat".format(sector,los),"r")
    nf= sum(1 for line in f1)  
    par=np.zeros((nf,nc))
    par=np.loadtxt("./files/D{0:d}_{1:d}.dat".format(sector, los)) 
    for i in range(nf):
        los,nsim,lat,lon                        =par[i,0],par[i,1],par[i,2],par[i,3]
        strucl, lml, Dl, vl, strucs, cl, mass,Ds=par[i,4],par[i,5],par[i,6],par[i,7],par[i,8],par[i,9],par[i,10],par[i,11]
        Tstar,Rstar, logg, typs, vs             =par[i,12],par[i,13],par[i,14],par[i,15],par[i,16]
        Mab, Map, magb, fb, nsbl, ext           =par[i,17],par[i,18],par[i,19],par[i,20],par[i,21],par[i,22]
        tE, RE, t0, mul, Vt, u0, opt            =par[i,23],par[i,24],par[i,25],par[i,26],par[i,27],par[i,28],par[i,29]
        ros,tetE, detflag, dchi, det            =par[i,30],par[i,31],par[i,32],par[i,33],par[i,34]
        numt, timet,meter                       =par[i,35],par[i,36],par[i,37]     
        los=int(los);  nsim=int(nsim) 
        Ml= pow(10.0,lml)  
        pirel=float(1.0/Dl-1.0/Ds)  
        par[i,36]= np.log10(pirel/tetE)  
        dchin=dchi/numt
        ############################################################
        if(i<no and k1<int(N0)):  
            pars[k1,:]=par[i,:];  
            k1+=1
            if(detflag>0): 
                pard[k2,:]=par[i,:]  
                k2+=1
        ########################################################################         
        if(nsim<50 and los%1==0): 
            print("****************************************************")
            print("Counter:  ", los,   nsim,    dchin,     meter)
            print("tE, ros, pirel:       ",  tE, ros, pirel)
            print("photometry_parameters:      ", Mab, Map, magb, fb, nsbl)
            print("****************************************************")   
            nd=-1;  nm=-1;  
            try: 
                f1=open('./files/L{0:d}_{1:d}_{2:d}.dat'.format(sector,los,nsim) )
                nd=int(len(f1.readlines()))
                print(nd)
            except: 
                print("file does not exist",  nsim)    
            try:
                f2=open('./files/M{0:d}_{1:d}_{2:d}.dat'.format(sector,los,nsim) )
                nm=int(len(f2.readlines()))  
                print (nm)
            except: 
                print("file does not exist",  nsim)        
            print("nsim,   nd,  nm:  ", nsim,   nd,  nm)   
            if(nd>1 and nm>1): 
                dat=np.zeros((nd,5))
                dat=np.loadtxt('./files/L{0:d}_{1:d}_{2:d}.dat'.format(sector,los,nsim)) 
                mod=np.zeros((nm,3))
                mod=np.loadtxt('./files/M{0:d}_{1:d}_{2:d}.dat'.format(sector,los,nsim)) 
                stat=[]
                if(detflag>0):
                    stat=r"$\rm{Detected}$";   col='g'
                else:
                    stat=r"$\rm{Not}-\rm{detected}$";   col='r'
                ymax1, ymin1= np.max(mod[:,1]), np.min(mod[:,1])
                ymax2, ymin2= np.max(dat[:,1]), np.min(dat[:,1])
                ymax, ymin= np.max(np.array([ymax1,ymax2])) , np.min(np.array([ymin1, ymin2]))
                plt.clf()
                plt.cla()
                fig=plt.figure(figsize=(8,6))
                ax1=fig.add_subplot(111)
                plt.plot(mod[:,0],mod[:,1],'k--',label=r"$\rm{Model}~\rm{Light}~\rm{Curve}$", lw=2.1,alpha=1.0)
                plt.errorbar(dat[:,0],dat[:,1],yerr=dat[:,2],fmt=".",markersize=6.8,color='m',ecolor='gray',elinewidth=0.2, capsize=0,alpha=0.7,label=r"$\rm{TESS}~\rm{Data}$")
                #plt.plot([],[],' ', color=col, label= str(stat))
                plt.ylabel(r"$\rm{TESS}-\rm{magnitude}$",fontsize=18)
                plt.xlabel(r"$\rm{time}(\rm{days})$",fontsize=18)
                plt.title(
                r"$t_{\rm{E}}(\rm{days})=$"+'{0:.1f}'.format(tE)+
                r"$,~\rho_{\star}=$"+'{0:.3f}'.format(ros)+  
                r"$,~f_{\rm{b}}=$"+'{0:.2f}'.format(fb)+
                r"$,~\log_{10}[M_{\rm{l}}(M_{\odot})]=$"+'{0:.2f}'.format(lml)+
                r"$,~\pi_{\rm{rel}}(\rm{mas})=$"+'{0:.2f}'.format(pirel)+
                r"$,~\mathcal{M}=$"+'{0:.1f}'.format(meter),fontsize=13.0,color='k')
                pylab.ylim([ymin, ymax])
                plt.xlim([0.0-1.0,Tobs+1.0])
                plt.xticks(fontsize=17,rotation=0)
                plt.yticks(fontsize=17,rotation=0)
                plt.gca().invert_yaxis()
                ax1.grid("True")
                ax1.grid(linestyle='dashed')
                ax1.legend(title=stat, prop={"size":14.})
                fig=plt.gcf()
                fig.tight_layout()
                fig.savefig("./Figs/Light{0:d}_{1:d}_{2:d}.jpg".format(sector,los,nsim),dpi=200)
             
################################################################################ 
nam0=['los', 'nsim', r"$\rm{Galactic}~\rm{Latitude}~(\rm{deg})$", r"$\rm{Galactic}~\rm{Longitude} (\rm{deg})$", 
r"$\rm{Lens}~\rm{Structure}$", r"$\log_{10}[M_{\rm{l}}(M_{\odot})]$", r"$D_{\rm{l}}(\rm{kpc})$", r"$v_{\rm{l}}(\rm{km}/s)$",
r"$\rm{source}~\rm{Structure}$",r"$\rm{CL}~\rm{source}$",r"$M_{\star}(M_{\odot})$",r"$D_{\rm{s}}(\rm{kpc})$",r"$T_{\star}(\rm{k})$",
r"$R_{\star}(R_{\odot})$", r"$\log_{10}[g]$", r"$\rm{source}~\rm{type}$", r"$v_{\star}(\rm{km}/s)$", 
r"$M_{\star,~\rm{T}}(\rm{mag})$", r"$m_{\star,~\rm{T}}(\rm{mag})$", r"$m_{\rm{base},~\rm{T}}(\rm{mag})$", 
r"$f_{\rm{b}}$", r"$N_{\rm{b}}$", r"$\rm{Extinction}$", 
r"$t_{\rm{E}}(\rm{days})$", r"$R_{\rm{E}}(\rm{AU})$", r"$t_{0}(\rm{days})$", r"$\mu_{\rm{rel}}(\rm{mas}/yrs)$", 
r"$V_{\rm{rel}}(\rm{km}/s)$", r"$u_{0}$", r"$\tau \times 10^{6}$", r"$\log_{10}[\rho_{\star}]$", r"$\theta_{\rm{E}}(\rm{mas})$", 
r"$Flag_{\rm{detect}}$", r"$\Delta \chi^{2}$", r"$\rm{Detect}$", r"$No.$", r"$\log_{10}[\pi_{\rm{E}}]$", r"$\log_{10}[\mathcal{M}]$"]
  
 
pars[:k1,30]=np.log10(pars[:k1,30])
pard[:k2,30]=np.log10(pard[:k2,30])
pars[:k1,36]=np.log10(pars[:k1,36])
pard[:k2,36]=np.log10(pard[:k2,36])


for i in range(nc):
    plt.clf()
    plt.cla()
    fig=plt.figure(figsize=(8,6))
    ax= plt.gca()              
    plt.hist(pars[:k1,i],30,histtype='bar',ec='darkgreen',facecolor='green', alpha=0.4,rwidth=1.0)
    plt.hist(pard[:k2,i],30,histtype='step',color='k', alpha=1.0, lw=1.5,ls='--')
    y_vals =ax.get_yticks()
    ax.set_yticks(y_vals)
    ax.set_yticklabels(['{:.2f}'.format(float(1.0*x*(1.0/k1))) for x in y_vals]) 
    y_vals = ax.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=18,labelpad=0.1)
    ax.set_xlabel(str(nam0[i]),fontsize=18,labelpad=0.1)
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./Figs/Histo/Histo{0:d}_{1:d}.jpg".format(sector,i),dpi=200)
    print("****  All histo_simulated_parameters are plotted ********************")

################################################################################

































































