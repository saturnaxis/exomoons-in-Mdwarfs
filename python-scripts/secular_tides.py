from scipy.integrate import solve_ivp,odeint
from scipy.optimize import minimize,fmin
import numpy as np
import threading
import os
import multiprocessing as mp
import sys

home = "/home/u.bq252623/Shaan_Tides/"

#constants
G = 4.*np.pi**2
M_sun = 1.989e33 #Mass of Sun in g
AU = 149597870700. #AU in m
AU_cm = AU*100.

M_E = 3.0035e-6 #M_sun
M_J = 9.54e-4 #M_sun
R_E = 4.2635e-5 #AU
# M_Moon = M_E/81.
M_Moon = 3.6943e-8
rho_E = 5.515 # g/cm^3
rho_Mars = 3.93
rho_Moon = 3.34
rho_N = 1.64 # g/cm^3
rho_J = 1. # g/cm^3
rho_moon = 3.34 #g/cm^3

psi_E = np.radians(0)


def write2file(fname,f_str,head):
    lock.acquire() # thread blocks at this line until it can obtain lock
    if head == "head":
        f = open(fname,'w')
    else:
        f = open(fname,'a')
    f.write(f_str)
    f.close()
    lock.release()

def f_ecc(f_idx,ecc):
    if f_idx == 1:
        return 1. + 31./2.*ecc**2 + 255./8.*ecc**4 + 185./16.*ecc**6 + 25./64.*ecc**8
    elif f_idx == 2:
        return 1. + 15./2.*ecc**2 + 45./8.*ecc**4 + 5./16.*ecc**6
    elif f_idx == 3:
        return 1. + 15./4.*ecc**2 + 15./8.*ecc**4 + 5./64.*ecc**6
    elif f_idx == 4:
        return 1. + 3./2.*ecc**2 + 1./8.*ecc**4
    elif f_idx == 5:
        return 1. + 3.*ecc**2 + 3./8.*ecc**4

def deriv_t_Hut(t,y,kp,tau_p,M_p,R_p,M_m,M_star):
    #y = [a_m,e_m,a_p,e_p,O_p]
    #secular tidal evolution (Hut 1981)

    #calculate factors
    T_Hut = R_p**3/(G*M_p*tau_p)
    q_pm = M_m/M_p
    q_ps = M_p/M_star
    dam_fact = -6*kp/T_Hut*q_pm*(1+q_pm)*(R_p/y[0])**8*y[0]
    dem_fact = -27*kp/T_Hut*q_pm*(1+q_pm)*(R_p/y[0])**8*y[1]
    em_beta = np.sqrt(1.-y[1]**2)

    dap_fact = -6*kp/T_Hut*q_ps*(1+q_ps)*(R_p/y[2])**8*y[2]
    dep_fact = -27*kp/T_Hut*q_ps*(1+q_ps)*(R_p/y[2])**8*y[3]
    ep_beta = np.sqrt(1.-y[3]**2)

    r_g = np.sqrt(C_p)
    n_m = np.sqrt(G*(M_p+M_m)/y[0]**3)
    n_p = np.sqrt(G*(M_star+M_p)/y[2]**3)
    dOp_pm = 3*kp/T_Hut*(q_pm/r_g)**2*(R_p/y[0])**6*n_m
    dOp_ps = 3*kp/T_Hut*(q_ps/r_g)**2*(R_p/y[2])**6*n_p

    #moon semimajor axis derivative
    damdt = dam_fact*(f_ecc(1,y[1])/em_beta**15-f_ecc(2,y[1])/em_beta**12*y[4]/n_m)
    #moon eccentricity derivative
    demdt = dem_fact*(f_ecc(3,y[1])/em_beta**13-(11./18)*f_ecc(4,y[1])/em_beta**10*y[4]/n_m)

    #planet semimajor axis derivative
    dapdt = dap_fact*(f_ecc(1,y[3])/em_beta**15-f_ecc(2,y[3])/ep_beta**12*y[4]/n_p)
    #planet eccentricity derivative
    depdt = dep_fact*(f_ecc(3,y[3])/ep_beta**13-(11./18)*f_ecc(4,y[3])/ep_beta**10*y[4]/n_p)


    dOpdt  = dOp_pm*(f_ecc(2,y[1])/em_beta**12-f_ecc(5,y[1])/em_beta**9*y[4]/n_m)
    dOpdt += dOp_ps*(f_ecc(2,y[3])/ep_beta**12-f_ecc(5,y[3])/ep_beta**9*y[4]/n_p)

    return [damdt,demdt,dapdt,depdt,dOpdt]

def deriv_t_Barnes(t,y,k2p,tau_p,M_p,R_p,M_m,M_star):
    #y = [a_m,e_m,a_p,e_p,O_p]
    #secular tidal evolution (Barnes 2017)

    #calculate factors
    dam_fact = 2*y[0]**2/(G*M_p*M_m)
    dem_fact = 11*y[0]*y[1]/(2*G*M_p*M_m) #
    em_beta = np.sqrt(1.-y[1]**2)
    Z_pm = 3*G**2*k2p*M_m**2*(M_p+M_m)*R_p**5/y[0]**9*tau_p

    dap_fact = 2*y[2]**2/(G*M_star*M_p)
    dep_fact = 11*y[2]*y[3]/(2*G*M_star*M_p)
    ep_beta = np.sqrt(1.-y[3]**2)
    Z_ps = 3*G**2*k2p*M_star**2*(M_p+M_star)*R_p**5/y[2]**9*tau_p

    r_g = np.sqrt(0.3308)
    n_m = np.sqrt(G*(M_p+M_m)/y[0]**3)
    n_p = np.sqrt(G*(M_star+M_p)/y[2]**3)
    dOp_pm = Z_pm/(2*M_p*r_g**2*R_p**2*n_m)
    dOp_ps = Z_ps/(2*M_p*r_g**2*R_p**2*n_p)

    #moon semimajor axis derivative
    damdt = dam_fact*Z_pm*(np.cos(psi_E)*f_ecc(2,y[1])/em_beta**12*y[4]/n_m-f_ecc(1,y[1])/em_beta**15)
    #moon eccentricity derivative
    demdt = dem_fact*Z_pm*(np.cos(psi_E)*f_ecc(4,y[1])/em_beta**10*y[4]/n_m-f_ecc(3,y[1])/em_beta**13)

    #planet semimajor axis derivative
    dapdt = dap_fact*Z_ps*(np.cos(psi_E)*f_ecc(2,y[3])/ep_beta**12*y[4]/n_p-f_ecc(1,y[3])/ep_beta**15)
    #planet eccentricity derivative
    depdt = dep_fact*Z_ps*(np.cos(psi_E)*f_ecc(4,y[3])/ep_beta**10*y[4]/n_p-f_ecc(4,y[3])/ep_beta**12)


    dOpdt  = dOp_pm*(2*np.cos(psi_E)*f_ecc(2,y[1])/em_beta**12-(1+np.cos(psi_E)**2)*f_ecc(5,y[1])/em_beta**9*y[4]/n_m)
    dOpdt += dOp_ps*(2*np.cos(psi_E)*f_ecc(2,y[3])/ep_beta**12-(1+np.cos(psi_E)**2)*f_ecc(5,y[3])/ep_beta**9*y[4]/n_p)

    return [damdt,demdt,dapdt,depdt,dOpdt]

def calc_Tide_evolution(par_idx):
    m_earth_mul, a_p = params[par_idx]
    # times = np.concatenate((np.arange(0,1000,100),np.arange(1000,t_fin+1000,1000)))
    times = np.arange(0,tscale+10000,10000)
    nsteps = len(times)
    M_p = m_earth_mul * M_E
    
    
    R_H = a_p*((M_p+M_sat)/(3.*M_star))**(1./3.)
    a_roche = 2.44*R_p * (5.515/3.34) ** (1./3.)
    a_m = 3 * a_roche
    e_m = 1e-4
    e_p = e_p0

    n_m = np.sqrt(G*(M_p+M_sat)/a_m**3) #mean motion of moon in rad/yr
    a_crit = 0.4031*(1.-1.123*e_p)*R_H
    n_p = np.sqrt(G*(M_star+M_p+M_sat)/a_p**3) #mean motion of planet in rad/yr
    T_p = 5./24./365.25 # rotation period in yr (Kokubo & Ida 2007)
    O_p = 2.*np.pi/T_p  #rotation rate of planet in rad/yr

    tau_p = tau/3600./24./365.25  #Hut-->550; Barnes-->600 ; Reboundx -->300
    moon_escape = False

    e_max = 0
    for i in range(1,nsteps):
        t_i = times[i-1]
        t_n = times[i]

        delta_t = [t_i,t_n]
        temp = solve_ivp(deriv_t_Hut,delta_t,[a_m,e_m,a_p,e_p,O_p],method='RK45',atol=1e-12,rtol=1e-12,args=(kp,tau_p,M_p,R_p,M_sat,M_star))

        if np.abs(temp.y[0,-1]-n_m)<1e-8:
            moon_escape = True
            break
        a_m,e_m,a_p,e_p,O_p = temp.y[0,-1],temp.y[1,-1],temp.y[2,-1],temp.y[3,-1],temp.y[4,-1]

        if e_m > e_max:
            e_max = e_m
        if a_m*(1.-e_m) < R_p: #check if q < R_p
            moon_escape = True
            break
        if a_m*(1.+e_m) > a_crit: #check if Q > a_crit
            moon_escape = True
            break
    P_n = 2.*np.pi/n_p
    numorbits = t_n/P_n
    outstg = "%1.2f, %1.3f, %1.3e, %1.3f, %1.3e, %1.3e\n" % (m_earth_mul,a_p,e_max,a_m/R_H,numorbits,t_n)
    write2file(fname,outstg,'foot')


lock = threading.Lock()

st_idx = int(sys.argv[1])
end_idx = st_idx + 5000
tau = int(sys.argv[2]) #in seconds
st_type = sys.argv[3]

fname = home + '%s_sec_%i.txt' % (st_type, tau)
if not os.path.exists(fname):
    outstg = "#M_p (M_E), a_p (AU), max e_m, final a_m (R_H), num of planetary orbits, lifetime (yr)\n"
    write2file(fname,outstg,'head')

e_p0 = 0.03
if st_type == 'M0':
    M_star = 0.57 #Mass in Msun
    a_i, a_o = 0.27, 0.521 #a_p in AU
elif st_type == 'M2':
    M_star = 0.44 #Mass in Msun
    a_i, a_o = 0.17, 0.351 #a_p in AU

R_p = R_E # R_E in AU
M_sat = M_Moon #Mass in Msun (of moon)
kp = 0.298/2 #k of host planet = k2p/2
C_p = 0.3308 #Norm MOI constant
tscale = 5e9 #max time of simulation

params = []
for m in np.arange(0.8,2.01,0.01):
    for a in np.arange(a_i, a_o,0.001):
        params.append((round(m, 6), round(a, 6)))
        
np.random.seed(42)
np.random.shuffle(params)

params = params[st_idx:end_idx]

#Load already processed parameters
datafile = np.genfromtxt(fname, delimiter=',', comments='#')
evaluated = set((round(row[0], 6), round(row[1], 6)) for row in datafile)
#Filter remaining params
params = [p for p in params if p not in evaluated]

print(st_idx,tau,st_type,"  ",len(params), "runs left.")
#sys.exit(0)

pool = mp.Pool(processes=60)
pool.map(calc_Tide_evolution, range(0,len(params)))
pool.close()