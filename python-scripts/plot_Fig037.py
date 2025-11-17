import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt
import sys
from matplotlib import rcParams
rcParams.update({'font.size': 14})  # set default font size

#--- Star Type ---
star_type = sys.argv[1]  # read from command line
star_types = ['M0','M2']
if star_type not in star_types:
    print("Please enter M0 or M2 for the desired star type.")
    sys.exit(0)

#--- Constants ---
m_earth = 3.0034e-6  # Earth mass in solar masses
r_earth = 4.26352e-5  # Earth radius in AU
fs = 'x-large'  # axis font size

def plot_col(mp, ap, color, ax_col, xmin, xmax, sublbl):
    ax_list_index = ax_list.index(ax_col)
    for k in range(0,2):  # tidal data first, secular second
        if k == 0:
            data = np.genfromtxt("../data/%s/reb_tide[%1.2f,%1.3f]_698.txt" % (figname,mp,ap),delimiter=',',comments='#')
        else:
            data = np.genfromtxt("../data/%s/m%s_sec_%1.1f_%1.2f.txt" % (figname,star_type[1],mp,ap),delimiter=', ')
        plotlabel = "%1.2f au, %1.1f $M_\\oplus$" % (ap,mp)
        a_p = data[:,3]
        R_H = a_p*(mp*m_earth/(3*Mcen))**(1/3.)  # Hill radius
        time, a_m, e_m = data[:,0], data[:,1], data[:,2]
        for i in range(0,3):
            if k == 1 and i==2:  # skip secular eccentricity
                continue
            ax = ax_col[i]
            if i == 0:  # top panel
                if k == 0:
                    ax.plot(time/T_p,a_m*r_earth/R_H,'.',ms=2,color=color,label=plotlabel)
                else:
                    ax.plot(time/T_p,a_m/R_H,'--',color='k')
                ax.set_ylim(0.1,0.4)
                if ax_list_index == 0: ax.set_ylabel("$a_\\mathrm{m} (R_\\mathrm{H})$",fontsize=20)
                else: ax.set_yticklabels([])
                ax.set_yticks(np.arange(0.2,0.5,0.1))
                ax.legend(loc='lower right',markerscale=5,fontsize=10)
                ax_top = ax.twiny()
                ax_top.minorticks_on()
                ax_top.tick_params(which='major',direction='out',length=8.0,width=4.0,labelsize=fs)
                ax_top.tick_params(which='minor',direction='out',length=4.0,width=4.0,labelsize=fs)
                ax_top.set_xlim(xmin*T_p,xmax*T_p)
                ax_top.set_xscale('log')
            elif i == 1:  # middle panel
                T_p_nonconstant = np.sqrt(a_p**3/Mcen)  # evolving planet period
                if k == 0:
                    T_m = np.sqrt((a_m*r_earth)**3/(mp*m_earth))
                    ax.plot(time/T_p,T_p_nonconstant/T_m,'.',ms=2,color=color,label=plotlabel)
                else:
                    T_m = np.sqrt(a_m**3/(mp*m_earth))
                    ax.plot(time/T_p,T_p_nonconstant/T_m,'--',color='k')
                ax.set_ylim(5,30)
                if ax_list_index == 0: ax.set_ylabel("Period Ratio $P_\\mathrm{p}/P_\\mathrm{m}$",fontsize=16)
                else: ax.set_yticklabels([])
            elif i == 2:  # bottom panel
                time_med, ecc_med = [], []
                left_bin,right_bin = 1e6,2e6
                while right_bin < data[-1,0]/T_p:
                    time_med.append(0.5*(right_bin+left_bin))
                    temp_idx = np.where(np.logical_and(left_bin<=time/T_p, right_bin>=time/T_p))[0]
                    ecc_med.append(np.median(e_m[temp_idx]))
                    left_bin,right_bin = right_bin,right_bin+1e6
                ax.plot(time/T_p,e_m,'.',ms=5,color=color,label=plotlabel)
                ax.plot(time_med,ecc_med,'.',alpha=1,ms=8,color='k')
                ax.set_ylim(0,0.2)
                if ax_list_index == 0: ax.set_ylabel("$e_\\mathrm{m}$",fontsize=20)
                else: ax.set_yticklabels([])
            ax.text(0.05,0.95,sublbl[i],transform=ax.transAxes,fontsize=20,fontweight='bold',va='top')
            ax.grid(True,color='gray',alpha=0.7,lw=1.5)
            ax.set_xlim(xmin,xmax)
            ax.set_xscale('log')
            ax.minorticks_on()
            ax.tick_params(which='major',direction='out',length=8.0,width=4.0,labelsize=fs)
            ax.tick_params(which='minor',direction='out',length=4.0,width=4.0,labelsize=fs)
            if i < 2: ax.set_xticklabels([])
    return

# Stellar parameters and plot settings by type
if star_type == 'M2':
    Mcen = 0.44
    ap_in = 0.17
    ap_out = 0.34
    xmin_list = [1e6,1e6,1e6]
    xmax_list = [2.3e7,9.9e7,5.3e8]
    figname = "Fig03"
else:
    Mcen = 0.57
    ap_in = 0.27
    ap_out = 0.52
    xmin_list = [1e6,1e6,1e7]
    xmax_list = [4.7e7,7.5e8,1.2e9]
    figname = "Fig07"

T_p = np.sqrt(ap_in**3/Mcen)  # reference orbital period

fig = plt.figure(figsize=(11,7),dpi=150)
ax_list = []
for i in range(3):  # 3 columns
    ax_col = []
    for j in range(3):  # 3 rows
        ax_col.append(fig.add_subplot(3,3,3*j+i+1))
    ax_list.append(ax_col)

# plot three planet configs
plot_col(2,ap_in,'deepskyblue',ax_list[0],xmin_list[0],xmax_list[0],['a','d','g'])
plot_col(1,ap_out,'red',ax_list[1],xmin_list[1],xmax_list[1],['b','e','h'])
plot_col(2,ap_out,'red',ax_list[2],xmin_list[2],xmax_list[2],['c','f','i'])

fig.subplots_adjust(hspace=0.1,wspace=0.1)
fig.text(0.5,0.01,"Normalized Planetary Orbits",ha='center',fontsize=18)
fig.text(0.5,0.97,"Time (yr)",fontsize=22,ha='center')

fig.savefig("../Figs/%s.png" % figname,bbox_inches='tight',dpi=300)
