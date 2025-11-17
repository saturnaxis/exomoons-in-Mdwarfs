import numpy as np
import matplotlib 
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt

def plot_fit(x,y,t0,order,ax):
    fit_idx = np.where(x>t0)[0]
    fit_t = x[fit_idx]
    fit_y = y[fit_idx]
    coeff = np.polyfit(fit_t,fit_y,order)
    fit_t_plot = np.arange(t0,9.5,0.01)
    ax.plot(10**fit_t_plot,np.poly1d(coeff)(fit_t_plot),'r--',lw=1,alpha=0.05)
    return coeff

# --- physical + file params ---
mp = 2.0
ap0 = 0.52
fname = "../data/Fig08/reb_tide[%1.2f,%1.3f]_698.txt" % (mp,ap0)
time,a_m,e_m,a_p,e_p=np.genfromtxt(fname,delimiter=',',comments='#',usecols=(0,1,2,3,4),skip_header=2,unpack=True)

# -- Constants --- 
M_star = 0.57      # Msun
M_E = 3.003e-6     # Msun
R_E = 4.26352e-5   # AU

R_H = a_p*(mp*M_E/(3*M_star))**(1./3.)
a_m *= R_E/R_H
stab_limit=0.4031*(1-1.123*e_p-0.1862*e_m) #from Rosario-Franco et al. (2020)

fs = 'x-large'
ms = 5
# --- plot ---
fig = plt.figure(figsize=(10,4),dpi=150)
ax = fig.add_subplot(111)

ax.plot(time,a_m,'k.',ms=ms)
ax.plot(time,stab_limit,'--',lw=2)

# --- random crossing estimates ---
t_rng=np.random.uniform(7.3,8,1000)
cross_x=np.zeros(len(t_rng))

for i,t_c in enumerate(t_rng):
    polys = []
    for y,order in [(a_m,1),(stab_limit,1)]:
        coeff = plot_fit(np.log10(time),y,t_c,order,ax)
        polys.append(np.poly1d(coeff))
    f1,f2 = polys
    fit_time_plot = np.arange(t_c,9.5,0.0001)
    diff = np.abs(f1(fit_time_plot)-f2(fit_time_plot))
    x_idx = np.where(diff<1e-5)[0]
    if len(x_idx)==0: continue
    inter_x = fit_time_plot[x_idx[0]]
    cross_x[i] = inter_x
    ax.axvline(10**inter_x,color='m',linestyle='-',lw=1,alpha=0.05)

x_mean,x_std = np.mean(cross_x),np.std(cross_x)

ax.set_xlim(10**6.5,10**9.5)
ax.set_ylim(0.1,0.4)
ax.set_xlabel("Time (yr)",fontsize=24)
ax.set_ylabel("$a_\\mathrm{m}\ (R_\\mathrm{H})$",fontsize=24)
ax.minorticks_on()
ax.tick_params(which='major',direction='out',length=8.0,width=4.0,labelsize=fs)
ax.tick_params(which='minor',direction='out',length=4.0,width=4.0,labelsize=fs)
ax.set_xscale('log')
ax.grid(True,color='gray',alpha=0.7,lw=1.5)

fig.savefig("../Figs/Fig08.png",bbox_inches='tight',dpi=300)
