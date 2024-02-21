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
sector=int(12)
if(sector==12): Tobs=27.942
if(sector==39): Tobs=27.9514
if(sector==4):  Tobs=25.95
if(sector==31): Tobs=25.43473
if(sector==13): Tobs=28.4417



ddeg=0.5
num=int(24.0/ddeg) 
Nlos=3##int(num*num)
nc= 46

###########################################################
N0=100000
no= int(N0/Nlos)
pars=np.zeros((N0,nc))
pard=np.zeros((N0,nc))
k1=0;  k2=0
###########################################################
for los0 in range(Nlos): 
    los=los0+0
    f1=open("./files2/D{0:d}_{1:d}.dat".format(sector,los),"r")
    nf= sum(1 for line in f1)  
    par=np.zeros((nf,nc))
    par=np.loadtxt("./files2/D{0:d}_{1:d}.dat".format(sector, los)) 
    for i in range(nf):
        los,nsim,lat,lon   =par[i,0],par[i,1],par[i,2],par[i,3]
        Ml, rml, inc, tet, ecen, period, semi,tp=par[i,4],par[i,5],par[i,6],par[i,7],par[i,8],par[i,9],par[i,10],par[i,11]
        scl, typ, Ds, mass, Tstar=par[i,12],par[i,13],par[i,14],par[i,15],par[i,16]
        Rstar, logg, metal, limb, Mab, Map=par[i,17],par[i,18],par[i,19],par[i,20],par[i,21],par[i,22]
        magb, fb,Nbl, Ext, dchi,numt,meter=par[i,23],par[i,24],par[i,25],par[i,26],par[i,27],par[i,28],par[i,29]
        error, snr, depth, cdpp= par[i,30],par[i,31],par[i,32],par[i,33]
        lMap,lMab, ratio, loggl, teffl=  par[i,34],par[i,35],par[i,36],par[i,37],par[i,38]
        numl, MagG, MagBP, MagRP = par[i,39], par[i,40], par[i,41], par[i,42]
        magG, magBP, magRP       = par[i,43], par[i,44], par[i,45]
        
        los=int(los);  nsim=int(nsim) 
        ntrans= Tobs/period
        #par[i,4]=pow(10.0,par[i,4])
        #Ml= pow(10.0,lml)      
        ############################################################
        if(1>0):##i<no and k1<int(N0)):  
            pars[k1,:]=par[i,:];  
            k1+=1
            if(snr>7.0 and ntrans>2.0): 
                pard[k2,:]=par[i,:]  
                k2+=1
        ########################################################################         
        if(nsim>=0 and los%1==0): 
            print("****************************************************")
            print("Counter:  ", los,   nsim,    dchi,     meter, snr, error, depth, cdpp )
            print("tp, ecen, lml, rml:       ",  tp, ecen, Ml, rml)
            print("photometry_parameters:      ", Mab, Map, magb, fb, Nbl)
            print("****************************************************")   
            nd=-1;  nm=-1;  
            try: 
                f1=open('./files2/L{0:d}_{1:d}_{2:d}.dat'.format(sector,los,nsim) )
                nd=int(len(f1.readlines()))
                print(nd)
            except: 
                print("file does not exist",  nsim)    
            try:
                f2=open('./files2/M{0:d}_{1:d}_{2:d}.dat'.format(sector,los,nsim) )
                nm=int(len(f2.readlines()))  
                print (nm)
            except: 
                print("file does not exist",  nsim)        
            print("nsim,   nd,  nm:  ", nsim,   nd,  nm)   
            if(nd>1 and nm>1): 
                dat=np.zeros((nd,3))
                dat=np.loadtxt('./files2/L{0:d}_{1:d}_{2:d}.dat'.format(sector,los,nsim)) 
                mod=np.zeros((nm,11))
                mod=np.loadtxt('./files2/M{0:d}_{1:d}_{2:d}.dat'.format(sector,los,nsim)) 
                stat=[]
                if(snr>7.0 and ntrans>2.0):
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
                plt.plot(mod[:,0],mod[:,1],'k--',label=r"$\rm{Model}~\rm{Light}~\rm{Curve}$", lw=1.3,alpha=1.0)
                plt.errorbar(dat[:,0],dat[:,1],yerr=dat[:,2],fmt=".",markersize=4.8,color='m',ecolor='gray',elinewidth=0.1, capsize=0,alpha=0.7,label=r"$\rm{TESS}~\rm{Data}$")
                #plt.plot([],[],' ', color=col, label= str(stat))
                plt.ylabel(r"$\rm{TESS}-\rm{magnitude}$",fontsize=18)
                plt.xlabel(r"$\rm{time}(\rm{days})$",fontsize=18)
                plt.title(
                r"$M_{\rm{WD}}(M_{\odot})=$"+'{0:.1f}'.format(Ml)+
                r"$,~\epsilon=$"+'{0:.2f}'.format(ecen)+  
                r"$,~\rm{Period}(\rm{days})=$"+'{0:.2f}'.format(period)+
                r"$,~\log_{10}[\mathcal{R}]=$"+'{0:.2f}'.format(np.log10(ratio))+
                r"$,~\rm{SNR}=$"+'{0:.3f}'.format(snr),fontsize=13.0,color='k')
                pylab.ylim([ymin, ymax])
                plt.xlim([0.0,Tobs])
                plt.xticks(fontsize=17, rotation=0)
                plt.yticks(fontsize=17, rotation=0)
                plt.gca().invert_yaxis()
                ax1.grid("True")
                ax1.grid(linestyle='dashed')
                ax1.legend(title=stat, prop={"size":13.})
                fig=plt.gcf()
                fig.tight_layout()
                fig.savefig("./Figs2/Light{0:d}_{1:d}_{2:d}.jpg".format(sector,los,nsim),dpi=200)
                ##################################################################3
                '''
                plt.cla()
                plt.clf()
                plt.figure(figsize=(8, 6))
                plt.plot(mod[:,0]/period,mod[:,3], color="r",label="x1(LoS)",  lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,4], color="b",label="y1(RIGHT)",lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,5], color="g",label="z1(UP)",   lw=1.2, alpha=0.95)
                plt.xticks(fontsize=17, rotation=0)
                plt.yticks(fontsize=17, rotation=0)
                plt.xlabel(r"$time$", fontsize=18)
                plt.grid("True")
                plt.legend()
                fig=plt.gcf()
                fig.savefig("./Figs2/xyz{0:d}_{1:d}_{2:d}.jpg".format(sector,los,nsim), dpi=200)
                ##################################################################
                plt.cla()
                plt.clf()
                plt.figure(figsize=(8, 6))
                plt.plot(mod[:,0]/period,mod[:,2], color="r",label="Magnification",  lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,7]/np.max(mod[:,7]), color="b",label="RE",lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,8]/np.max(mod[:,8]), color="g",label="ros",   lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,9]/np.max(mod[:,9]), color="m",label="u",   lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,10]/np.max(mod[:,10]), color="k",label="us",   lw=1.2, alpha=0.95)
                plt.xticks(fontsize=17, rotation=0)
                plt.yticks(fontsize=17, rotation=0)
                plt.xlabel(r"$time$", fontsize=18)
                plt.grid("True")
                plt.legend()
                fig=plt.gcf()
                fig.savefig("./Figs2/Astar{0:d}_{1:d}_{2:d}.jpg".format(sector,los,nsim), dpi=200)
                '''
print ("fraction of detection:  ", k2*100.0/k1)                
################################################################################ 
nam0=['los', 'nsim', 's.lat', 'lonn','MBH', 'RBH/Rsun', 'inc', 'tet', 'ecen','period',  
'semi/Rsun', 'tp','cl', 'type', 'Ds','mass', 'Tstar', 'Rstar', 'logg', 'metal', 'limb', 
'Mab', 'Map', 'magb', 'blend', 'nsbl', 'Ai','dchi', 'numt', 'meter', 'error', 'snr', 'depth','cdpp',
'l.Map', 'l.Mab', 'l.ratio', 'l.logg', 'l.teff', 'l.num', 'l.MagG', 'l.MagBP', 'l.MagRP', 'l.magG', 'l.magBP', 'l.magRP']

#pars[:k1,30]=np.log10(pars[:k1,30])
#pard[:k2,30]=np.log10(pard[:k2,30])
#pars[:k1,36]=np.log10(pars[:k1,36])
#pard[:k2,36]=np.log10(pard[:k2,36])


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
    fig.savefig("./Figs2/Hist/Histo{0:d}_{1:d}.jpg".format(sector,i),dpi=200)
    print("****  All histo_simulated_parameters are plotted ********************")

################################################################################
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8,6))
plt.plot( pars[:k1,41]-pars[:k1,42], pars[:k1,40], "go",  ms=3.5, label=r"$\rm{Absolute}~\rm{Magnitude}$")
plt.xlabel(r"$\rm{M}_{\rm{G}}(\rm{mag})$",  fontsize=18)
plt.ylabel(r"$G_{\rm{BP}}-G_{\rm{RP}}(\rm{mag})$",fontsize=18)
#plt.xlim([2.5, 17.5])
#plt.ylim([-5.2,-1.0])
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
plt.gca().invert_yaxis()
plt.legend(prop={"size":15}, loc='best')
fig.tight_layout()
fig=plt.gcf()
fig.savefig("./Figs2/Hist/CM_WD1.jpg")


################################################################################
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8,6))
plt.plot( pars[:k1,44]-pars[:k1,45], pars[:k1,43], "go",  ms=3.5, label=r"$\rm{Apparent}~\rm{Magnitude}$")
plt.xlabel(r"$\rm{m}_{\rm{G}}(\rm{mag})$",  fontsize=18)
plt.ylabel(r"$g_{\rm{BP}}-g_{\rm{RP}}(\rm{mag})$",fontsize=18)
#plt.xlim([2.5, 17.5])
#plt.ylim([-5.2,-1.0])
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
plt.gca().invert_yaxis()
plt.legend(prop={"size":15}, loc='best')
fig.tight_layout()
fig=plt.gcf()
fig.savefig("./Figs2/Hist/CM_WD2.jpg")





################################################################################














































