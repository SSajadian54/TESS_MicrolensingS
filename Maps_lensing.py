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
if(sector==12): 
    elati=float(-18.0)
    elong=float(253.56) 

    
if(sector==39): 
    elati=float(-18.0)
    elong=float(259.2)
    
if(sector==4): 
    elati=float(-18.0)
    elong=float(38.48)
    
if(sector==31): 
    elati=float(-18.0)
    elong=float(42.4)
    
if(sector==13): 
    elati=float(-18.0)
    elong=float(281.12)

lef, rig, bot, top=  elong-12.0, elong+12.0, elati-12.0, elati+12.0
##==============================================================================

if(sector==12): Tobs=27.942
if(sector==39): Tobs=27.9514
if(sector==4):  Tobs=25.95
if(sector==31): Tobs=25.43473
if(sector==13): Tobs=28.4417
nam=[r"$\rm{Los}$", r"$\rm{Ecliptic}~\rm{Latitude}(\rm{deg})$", r"$\rm{Ecliptic}~\rm{Longitude}(\rm{deg})$", 
r"$\rm{Galactic}~\rm{Latitude}(\rm{deg})$", r"$\rm{Galactic}~\rm{Longitude}(\rm{deg})$", r"$\rm{Right}~\rm{Ascention}(\rm{deg})$", 
r"$\rm{Declination}(\rm{deg})$", r"$\rm{Nsim}$" , r"$\rm{Ndet}$", r"$\log_{10}[N_{\star}(\rm{deg}^{-2})]$", 
r"$\log_{10}[\rho_{\star}(M_{\odot}/\rm{deg}^{2})]$", r"$\log_{10}[\overline{\epsilon(t_{\rm{E}})/t_{\rm{E}}}]$",r"$\log_{10}[N_{\star, \rm{TESS}}(\rm{deg}^{-2})]$", r"$\log_{10}[\tau_{\rm{t}}\times 10^{6}]$",r"$\log_{10}[\tau\times 10^{6}]$", 
r"$\log_{10}[\Gamma(\rm{star}^{-1}\rm{days}^{-1})]$", r"$\log_{10}[N_{\rm{e}}(\rm{deg}^{-2})]$", r"$N_{\rm{com}}$", 
r"$\overline{t_{\rm{E}}}(\rm{days})$", r"$\overline{m_{\rm{base}}}(\rm{mag})$", r"$\overline{f_{\rm{b}}}$", 
r"$\overline{D_{\rm{l}}}(\rm{kpc})$", r"$\overline{D_{\rm{s}}}(\rm{kpc})$", r"$\log_{10}[\overline{\pi_{\rm{E}}}]$", 
r"$\log_{10}[\overline{\pi_{\rm{rel}}} (\rm{mas})]$", r"$\log_{10}[\overline{\rho_{\star}}]$" , r"$\overline{R_{\rm{E}}}(\rm{au})$", 
r"$\overline{v_{\rm{rel}}}(\rm{km}/s)$", r"$\overline{u_{0}}$", r"$\log_{10}[\overline{M_{\rm{l}}(M_{\odot})}]$", 
r"$f_{\rm{BD}}$", r"$f_{\rm{FFPs}}$", r"$f_{\rm{giant}}$", r"$f_{\rm{BD}, \rm{Com}}$", r"$f_{\rm{FFPs}, \rm{Com}}$", 
r"$f_{\rm{giants}, \rm{Com}}$", r"$f_{\rm{side}, \rm{Com}}$", r"$f_{\rm{BL}, \rm{Com}}$" , r"$N_{\rm{sim}, \rm{Com}}$", r"$N_{\rm{det}, \rm{Com}}$", r"$\overline{m_{\star,\rm{T}}} (\rm{mag})$", r"$\log_{10}[\overline{\sigma_{\rm{m}}}(\rm{mag})]$"]

################################################################################        
ddeg=0.5
num=int(24.0/ddeg) 
Nlos=int(num*num)
nc= 42
f1=open("./files/TESS{0:d}.dat".format(sector),"r")
dat=np.zeros((Nlos, nc))
dat=np.loadtxt("./files/TESS{0:d}.dat".format(sector)) 
nu=0
mapp=np.zeros((nc,num,num))
for i in range(num):
    for j in range(num):
        mapp[:,i,num-j-1]=dat[nu,:]
        nu+=1
################################################################################        
tt=int(9)
v=np.zeros((tt))
for i in range(nc):
    plt.cla()
    plt.clf()
    plt.figure(figsize=(8,8))
    plt.imshow(mapp[i,:,:],cmap='viridis',extent=(lef,rig,bot,top),interpolation='nearest',aspect='auto', origin='lower')
    plt.clim()
    plt.title(str(nam[i]), fontsize=18)
    minn=np.min(mapp[i,:,:])
    maxx=np.max(mapp[i,:,:])
    step=float((maxx-minn)/(tt-1.0));
    for m in range(tt):
        v[m]=round(float(minn+m*step),1)
    cbar=plt.colorbar(orientation='horizontal',shrink=0.85,pad=0.08,ticks=v)
    cbar.ax.tick_params(labelsize=17)
    plt.clim(v[0]-0.005*step,v[tt-1]+0.005*step)
    plt.xticks(np.arange(elong-11.5, elong+12.0, 4))
    plt.yticks(np.arange(elati-11.5, elati+12.0, 4))
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    ax=plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.tick_params(direction='out',pad=5,top=False, right=False)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    plt.xlabel(r"$\rm{Ecliptic}~\rm{longitude}~(\rm{deg})$",fontsize=18,labelpad=0.05)
    plt.ylabel(r"$\rm{Ecliptic}~\rm{latitude}~(\rm{deg})$", fontsize=18,labelpad=0.05)
    fig=plt.gcf()
    fig.tight_layout(pad=0.3)
    fig.savefig("./Figs/Histo/MapB{0:d}_{1:d}.jpg".format(sector,i),dpi=200)
    print ("map is plotted  i , sector:    ",  sector,  i)

################################################################################        


   





























