from __future__ import print_function
from __future__ import absolute_import
import extract
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
import yaml
from plots_dop_paper import PlotDOSPreview,PlotDOSAll, PlotDOSIP

# specific parameters
figname='fig1'
E_fermi=-5.22

# general 
fid= open('jobs/r_0/dop_pap.yaml','r')
settings = yaml.load(fid)
fid.close()
disorder = settings['layers'][0]['molecule_species'][0]['electron_affinity']['disorder']
dop_frac = settings['layers'][0]['molecule_species'][1]['concentration']

#hars coded: TODO: txt
NE	=	1000
Emin	=	-7.5
Emax	=	-0.0
Energy = np.linspace(Emin,Emax,NE)

nr = int(np.loadtxt('s4s.txt')[2])

J = np.zeros(nr)
JJ = 0.0

#for DOS
IP_0_0 = np.zeros([nr,NE])
IPIP_0_0 = np.zeros(NE)

EA_plus_0_0 = np.zeros([nr,NE])
EAEA_plus_0_0 = np.zeros(NE)

EA_1_0 = np.zeros([nr,NE])
EAEA_1_0 = np.zeros(NE)

IP_minus_1_0 = np.zeros([nr,NE])
IPIP_minus_1_0 = np.zeros(NE)


for i_r in range(0,nr):
    fn = 'jobs/r_{}/output_job_0'.format(i_r)
    print(fn)
    J[i_r] = extract.extract_float(fn, "Ionized dopant fraction: ")
    IP_0_0[i_r,:] = np.loadtxt("jobs/r_{}/experiments/tmp/IP_0_0.dat".format(i_r))
    EA_1_0[i_r,:] = np.loadtxt("jobs/r_{}/experiments/tmp/EA_1_0.dat".format(i_r))
    EA_plus_0_0[i_r,:] = np.loadtxt("jobs/r_{}/experiments/tmp/EA_plus_0_0.dat".format(i_r))
    IP_minus_1_0[i_r,:] = np.loadtxt("jobs/r_{}/experiments/tmp/IP_minus_1_0.dat".format(i_r))
    JJ 		= np.mean(J[:])
    IPIP_0_0[:]= np.mean(IP_0_0[:,:], 0)
    EAEA_1_0[:]= np.mean(EA_1_0[:,:], 0)
    EAEA_plus_0_0[:]= np.mean(EA_plus_0_0[:,:], 0)
    IPIP_minus_1_0[:]= np.mean(IP_minus_1_0[:,:], 0)
print("Av. ionized dopant fraction: {}".format(JJ))


#-----------------------------------
#    Generate template for Fig.
#-----------------------------------
ylim=[-6.1, -4.4]
resc_IP_h=1
resc_EA_h=1
resc_IP_d=1
resc_EA_d=1
rescalex=4
color_IP_h='blue'
color_EA_h='red'
color_IP_d='green'
color_EA_d='orange'
IP_h,EA_h,IP_d,EA_d = IPIP_0_0,EAEA_plus_0_0,IPIP_minus_1_0,EAEA_1_0
linewidth=2

fig,axes = plt.subplots(nrows=1,ncols=2,figsize = (4.5,7),sharey=True)
fig.subplots_adjust(wspace=0)
ax_h,ax_d=axes[0],axes[1]
max_DOS_h=np.max([EA_h,IP_h])
# Host
# plot curves
ax_h.plot(resc_IP_h*IP_h/max_DOS_h,         Energy, label = 'IP$_{h}$', color=color_IP_h,alpha=0.5,lw=linewidth)# IP host 
ax_h.fill_between(resc_IP_h*IP_h/max_DOS_h, Energy, facecolor=color_IP_h,alpha=0.5)
ax_h.plot(resc_EA_h*EA_h/max_DOS_h,    Energy, label = 'EA$^+_{h}$', color=color_EA_h,alpha=1,lw=linewidth)# EA host
ax_h.set_ylabel('Energy [eV]')
ax_h.set_xlabel('Host DOS')
ax_h.set_ylim(ylim)
ax_h.set_xlim(0,rescalex)
ax_h.spines['right'].set_visible(False)
ax_h.set_xticks(np.linspace(0,int(rescalex),3))
# Dopant
max_DOS_d=np.max([EA_d,IP_d])
# plot curves
ax_d.plot(IP_d*resc_IP_d/max_DOS_d, Energy, label = 'IP$^-_{d}$',color=color_IP_d,alpha=0.5,lw=linewidth)# IP dop
ax_d.fill_between(IP_d*resc_IP_d/max_DOS_d,Energy, facecolor=color_IP_d, alpha=0.5)
ax_d.plot(EA_d*resc_EA_d/max_DOS_d,       Energy, label = 'EA$_{d}$',  color=color_EA_d, alpha=1,lw=linewidth)# EA dop
ax_d.set_ylim(ylim)
ax_d.set_xlim(rescalex,0)
ax_d.spines['left'].set_visible(False)
ax_d.get_yaxis().set_visible(False)
ax_d.set_xticks(np.linspace(0,int(rescalex),3))
ax_d.set_xlabel('Dopant DOS')

fig.legend(loc=(0.91,0.11))

plt.savefig('{}_curves.svg'.format(figname))
plt.savefig('{}_curves.png'.format(figname))
plt.clf()

print('Max EA_h: {} eV'.format(np.max(EAEA_plus_0_0)))
print('Max IP_h: {} eV'.format(np.max(IPIP_0_0)))
print('Max EA_d: {} eV'.format(np.max(EAEA_1_0)))
print('Max IP_d: {} eV'.format(np.max(IPIP_minus_1_0)))