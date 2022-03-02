
import os
import numpy as np
from scipy.interpolate import interp2d
from copy import deepcopy
import pickle
import random
from datetime import datetime
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def nedderimp(t1,params):
# gives inner radius as a function of time for a cylinder implosion
# input t1: spits out inner radius times at each of the times in vector t1.
# the 6-vector params = [R1 lam s rho mrat u0]
#
# function that computes the time vs. y (inner radius) trace of the implosion
# according to Neddermeyer's formula given in his '43 paper.
# The basic form of the equation is of the form:
#
# y' = [2/(R1^2*f(y)^2)*{v0^2/2 - s/(2*rho)*g(y)}]^{-.5}
# 
# where
# f(y) = (y^2/(1-lam^2))*ln((y^2+1-lam^2)/y^2)   and
# g(y) = 1/(1-lam^2)*[y^2+1-lam^2]*ln(y^2+1-lam^2)-lam^2*ln(lam^2)
#
# we're simulating experiment 10 in Neddermeyer's '43 paper.
# 
# lam = 2/3   ratio of outer cylinder radius to the inner radius
# s = 3*10^5*6.84*10^4  yield stress in cgs units (300,000 lb/in^2)
# rho = 7.5 g/cm^3  -- specific density (relative to water)
# R1 = 1.5in  -- initial outer radius of steel cylinder being imploded
# R2 = 1.0in  -- initial inner radius of steel cylinder being imploded
# v0 = .3*10^5  cm/s -- initial velicity imparted on outer radius of
#                       cylinder from the HE
# mrat = mass ratio of HE to cylinder = rho_HE*h_HE/(rho_cyl*h_cyl) = .32
#        in expt 10
# u0 = energy per gram of exploded gas from the HE
# critical velocity given by sqrt(g(0)*s/rho) = 5.816*10^4 cm/s
# final scaled inner radius given by solving
#  v0^2/2 - s/(2*rho)*g(y)=0  => resluting y = .505
# y = r2/R1; x = r1/R1; x^2 - y^2 = lam^2 (conservation of mass)
#
# R1 = initial outer radius;  r1 = outer radius as a function of time...
# R2 = initial inner radius;  r2 = inner radius as function of time...
# lam = R2/R1 = 1.0in/1.5in = 2/3.

    R1 = params[0]; lam = params[1]; R2 = R1*lam; 
    s = params[2]; rho = params[3]; mrat = params[4]; u0 = params[5];
    if (2*u0/(1+mrat)) < 0: print('sqrt of',2*u0/(1+mrat))
    v0 = mrat*np.sqrt(2*u0/(1+mrat));
    # the 6-vector params = [R1 lam s rho mrat u0]
    nt = len(t1);
    dt = 1e-8;
    ninc = 8000;
    yout = np.zeros((ninc, 1));
    yout[0] = lam;
    tout = np.arange(0,dt*ninc-dt,dt)
    for i in range(1,ninc):
        yout[i] = yout[i-1] + (v(yout[i-1],params)*dt);
    
    ineg = (yout < 0);
    yout[ineg] = 0;
    
    #out = interp1(tout,yout,t1,'linear');
    return(np.interp(t1.squeeze(),tout,yout.squeeze()))

def f(y,params): # computes Neddermeyer's f(y)
    R1 = params[0]; lam = params[1]; R2 = R1*lam; 
    s = params[2]; rho = params[3]; mrat = params[4]; u0 = params[5];
    if (2*u0/(1+mrat))<0: print('sqrt of',2*u0/(1+mrat))
    v0 = mrat*np.sqrt(2*u0/(1+mrat));
    fout = y**2/(1-lam**2)*np.log((y**2 + 1 - lam**2)/y**2);
    return(fout)

def g(y,params): # computes Neddermeyer's g(y)
    R1 = params[0]; lam = params[1]; R2 = R1*lam; 
    s = params[2]; rho = params[3]; mrat = params[4]; u0 = params[5];
    if (2*u0/(1+mrat)) < 0: print('sqrt of',2*u0/(1+mrat))
    v0 = mrat*np.sqrt(2*u0/(1+mrat));
    gout = 1/(1-lam**2)*(y**2*np.log(y**2)-(y**2+1-lam**2)*np.log(y**2+1-lam**2)-lam**2*np.log(lam**2));
    return(gout)

def v(y,params): # computes the velocity as a function of y
    R1 = params[0]; lam = params[1]; R2 = R1*lam; 
    s = params[2]; rho = params[3]; mrat = params[4]; u0 = params[5];
    if (2*u0/(1+mrat)) < 0: print('sqrt of',2*u0/(1+mrat))
    v0 = mrat*np.sqrt(2*u0/(1+mrat));
    x1 = v0**2/2 - s*g(y,params)/(2*rho);
    ineg = x1 < 0;
    x1[ineg] = 0;
    if (2./(R1**2*f(y,params))*x1) < 0: print('sqrt of',2./(R1**2*f(y,params))*x1)
    vout = -np.sqrt(2./(R1**2*f(y,params))*x1);
    return(vout)

def dist2pi(x1,x2):
    # computes the distance assuming periodicity: 2pi=0
    # x1 and x2 are vectors with common length and values
    # between 0 and 2pi
    d = abs(x1-x2)
    iwrap = d > np.pi
    d[iwrap] = 2*np.pi - d[iwrap]
    return(d)
def dnorm(x,mu,scale):  
    # normal density in 1-d. 
    # It is scaled so that the 1-d integral is 1
    # mu and scale are scalars, x is an array...
    out=np.zeros(len(x))
    u=abs(x-mu)/scale
    out = (1.0/(np.sqrt(2*np.pi)*scale)) * np.exp(-.5 * u**2)
    return(out)
  
design = np.array([
  [0.7714,    0.4286,    0.0286],
  #[0.76923077,    .5,    0.5],
  [0.3714,    0.1143,    0.7143],
  [0.1714,    0.4571,    0.8857],
  [0.3429,    0.6000,    0.8000],
  [0.8000,    0.6286,    0.4000],
  [0.7429,    0.5429,         0],
  [0.6571,    1.0000,    0.6286],
  [0.2857,         0,    0.4571],
  [0.5143,    0.9429,    0.2286],
  [0.6857,    0.3143,    0.6571],
  [0.8286,    0.2000,    1.0000],
  [0.9714,    0.3429,    0.6000],
  [0.4000,    0.8000,    0.2000],
  [0.5429,    0.2857,    0.2857],
  [0.9143,    0.8857,    0.2571],
  [0.0571,    0.0286,    0.0857],
  [0.1143,    0.5714,    0.7429],
  [0.2000,    0.2286,    0.3714],
  [0.4571,    0.9143,    0.3429],
  [0.6286,    0.7143,    0.6857],
  [     0,    0.8286,    0.9429],
  [0.8857,    0.0857,    0.9714],
  [0.2286,    0.0571,    0.5714],
  [0.7143,    0.1714,    0.8571],
  [0.2571,    0.4857,    0.1429],
  [0.5714,    0.4000,    0.8286],
  [0.9429,    0.6857,    0.4857],
  [0.4857,    0.1429,    0.1143],
  [1.0000,    0.8571,    0.9143],
  [0.6000,    0.6571,    0.5143],
  [0.1429,    0.7429,    0.5429],
  [0.8571,    0.2571,    0.0571],
  [0.3143,    0.3714,    0.4286],
  [0.4286,    0.7714,    0.7714],
  [0.0286,    0.9714,    0.3143],
  [0.0857,    0.5143,    0.1714]])

# number of experiements and simulations
n = 4; m = design.shape[0]

# these parameter values simulate expt 10 in Neddermeyer '43
# params =             [R1      lam     s         rho   mratio   u0]
params10 = np.array([1.5*2.54, 2/3,  3e5*6.84e4,  7.5,   .32,   1.65e10])
paramslhs = np.zeros((m, 6))
for i in range(m):
  paramslhs[i,:] = params10*np.array([1, 1, design[i,1]*.2+.9, 1, design[i,0]*.65+.5, design[i,2]*.2+.9])
# the simulation runs will vary mratio from .32*[.5 to 1.15]
#                                s      from s0*[.9 to 1.1]
#                                u0     from u0*[.9 to 1.1]

nt = 22; nphi = 26
time = np.c_[np.linspace(0,5.0e-5,nt,endpoint=True)]
phi = np.linspace(0,1,nphi,endpoint=True) * 2*np.pi;
rinner = nedderimp(time,params10);
lam = params10[1]; R1 = params10[0];
router = np.sqrt(rinner**2 + 1 - lam**2);
xycirc = np.array([np.cos(phi),np.sin(phi)]).T
r = nedderimp(time,params10);

# fig, axs = plt.subplots(2,5,figsize=[10,4])
# tframes = np.arange(1,nt,2)
# time_text = ['0','1','2','3','4',r'$\times 10^{-5} s$']
# time_x_loc = np.arange(-4,2,1)
# time_y_loc = -3.9
# for j,ax in enumerate(axs.flatten()):
#     it = tframes[j]
#     outer_circle = xycirc*router[it]*R1
#     ax.fill(outer_circle[:,0],outer_circle[:,1],c='gray')
#     inner_circle = xycirc*rinner[it]*R1
#     ax.fill(inner_circle[:,0],inner_circle[:,1],c='white')
#     ax.set_xlim([-4.5,4.5]); ax.set_ylim([-4.5,4.5])
#     for i in range(len(time_text)):
#         ax.text(time_x_loc[i],time_y_loc,time_text[i],fontsize=9)
#     ax.plot((time[0:it]/(6e-5)*6)-4,np.repeat(-4,time[0:it].shape[0]),'lightgreen',linewidth=4)
#     ax.get_xaxis().set_ticks([])
#     if j == 0 or j == 5: ax.set_ylabel('cm')
#     else: ax.get_yaxis().set_visible(False)
#         
#plt.subplots_adjust(wspace=0,hspace=0.05)
#plt.savefig('plots/nedderdata.png',dpi=300,bbox_inches='tight')
#plt.show()

# fig = plt.figure(figsize=[6,12])
# ax1 = fig.add_subplot(211)
yr = np.zeros((m,nt))
for i in range(m):
     params = paramslhs[i,:]
     yr[i,:] = params[0]*nedderimp(time,params)
     
#     ax1.plot(time,yr[i,:],c='lightgreen')
# ax1.set_xlabel('time (s)'); ax1.set_ylabel('inner radius (cm)')
# ax1.set_title('Simulation runs')

y_sim = np.tile(yr,nphi)
# ax2 = fig.add_subplot(212, projection='3d')
# example_runs = [0,1,28]
# colors = ['tomato','deepskyblue','lightgreen']
# for i in range(len(example_runs)):
#     ax2.plot_wireframe(time,phi,y_sim[example_runs[i],:].reshape(nt,nphi,order='F'),color=colors[i],alpha=.75)
# ax2.set_xlabel('time (s)'); ax2.set_ylabel('phi (angle, radians)'); ax2.set_zlabel('inner radius (cm)')

#plt.savefig('plots/neddersims.png',dpi=300,bbox_inches='tight')
#plt.show()
# indices
y_sim_ind_time_phi = np.zeros((22*26,2))
y_sim_ind_time_phi[:,0] = np.repeat(time,26)
y_sim_ind_time_phi[:,1] = np.tile(phi,22)
x_sim = design[:,0].reshape(m,1)
t_sim = design[:,1:3]

phi_obs = np.arange(0,(2*np.pi-.1),(2*np.pi/16))
n_phi_obs = phi_obs.shape[0]
#time_obs = [np.array([1.5e-5, 2.7e-5, 4.5e-5]),np.array([4.5e-5]),np.array([2.5e-5, 4.5e-5])]
time_obs = [np.array([1.5e-5, 2.7e-5, 4.5e-5]),np.array([1.5e-5, 2.7e-5, 4.5e-5]),np.array([1.5e-5, 2.7e-5, 4.5e-5]),np.array([1.5e-5, 2.7e-5, 4.5e-5])]
n_time_obs = [tmp.shape[0] for tmp in time_obs]

phiknots = np.arange(0,2*np.pi-.1,2*np.pi/8)
dknots = np.expand_dims(np.array([.04, -.03, .03, -.03, .02, -.03, .03, -.03]),1)*2.5
pphiknots = len(phiknots)
Ddelt = np.zeros((phi_obs.shape[0], pphiknots));
datadelt = np.matmul(Ddelt,dknots)

# observations
ssq = 0
r_obs = [None]*n; y_obs = [None]*n
for i in range(n):
  obs_params = deepcopy(params10)
  if i==1: obs_params[4]=.17
  elif i==2: obs_params[4]=.36
  elif i==3: obs_params[4]=.25
  r_obs[i] = np.atleast_2d(obs_params[0]*nedderimp(time_obs[i],obs_params))
  y_obs[i] = np.tile(r_obs[i].T,phi_obs.shape[0]).reshape(n_phi_obs,n_time_obs[i])
  y_obs[i] += np.tile(datadelt,n_time_obs[i])
  y_obs[i] = (y_obs[i] + ssq*np.random.normal(size=y_obs[i].shape)).flatten()

# indices of observations
x_obs = ((np.array([params10[4], .17, .36, .25])/.32-.5)/.65).reshape(n,1)
# create y_ind_obs where each row is a (time, angle) pair. Pairs are grouped by time.
# experiment 1
y_ind_obs_1 = np.column_stack( ( np.concatenate((np.ones(phi_obs.shape[0])*time_obs[0][0],\
                                 np.ones(phi_obs.shape[0])*time_obs[0][1],\
                                     np.ones(phi_obs.shape[0])*time_obs[0][2])), np.tile(phi_obs,3).T ) )
# experiment 2

#y_ind_obs_2 = np.column_stack( ( (np.ones(phi_obs.shape[0])*time_obs[1]).reshape(16,1), phi_obs.T ) )
y_ind_obs_2 = y_ind_obs_1
# experiment 3
#y_ind_obs_3 = np.column_stack( ( np.concatenate((np.ones(phi_obs.shape[0])*time_obs[2][0],\
#                                   np.ones(phi_obs.shape[0])*time_obs[2][1])), np.tile(phi_obs,2).T ) )
y_ind_obs_3 = y_ind_obs_1
y_ind_obs_4 = y_ind_obs_1
# Store in a list containing all 3 experiments
y_ind_obs = [y_ind_obs_1, y_ind_obs_2, y_ind_obs_3, y_ind_obs_4]
del y_ind_obs_1, y_ind_obs_2, y_ind_obs_3, y_ind_obs_4

data_dict = dict([('Xsim',x_sim),('XTsim',design),\
              ('Tsim',t_sim),('Ysim',y_sim),\
              ('YindSim',y_sim_ind_time_phi),
              ('Xobs',x_obs),('Yobs',y_obs),('YindObs',y_ind_obs)])
