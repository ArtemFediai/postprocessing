
from __future__ import print_function
from __future__ import absolute_import
import extract
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
import yaml
from plots_dop_paper import PlotDOSPreview,PlotDOSAll, PlotDOSIP, PlotDOSDopant

# specific parameters
figname='fig4b'


# general 
fid= open('jobs/r_0/dop_pap.yaml','r')
settings = yaml.load(fid)
fid.close()
disorder = settings['layers'][0]['molecule_species'][0]['electron_affinity']['disorder']
dop_frac = settings['layers'][0]['molecule_species'][1]['concentration']
offset   = settings['layers'][0]['molecule_species'][0]['ionization_potential']['mean']-settings['layers'][0]['molecule_species'][1]['electron_affinity']['mean']

#hard coded: TODO: txt
NE	=	1000
Emin	=	-7.5
Emax	=	-0.0

if not os.path.exists('DOS_{}.txt'.format(figname)):
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

    # save DOS data
    with open('DOS_{}.txt'.format(figname),'w') as f:
        f.write("# Av. ionized dopant fraction: {}\n".format(JJ))
        f.write('# Offset: {} eV\n'.format(offset))
        f.write('# Max EA_h: {} eV\n'.format(np.max(EAEA_plus_0_0)))
        f.write('# Max IP_h: {} eV\n'.format(np.max(IPIP_0_0)))
        f.write('# Max EA_d: {} eV\n'.format(np.max(EAEA_1_0)))
        f.write('# Max IP_d: {} eV\n'.format(np.max(IPIP_minus_1_0)))
        f.write('# Av. DOS: IP_h   EA_h   IP_d   EA_d\n')
        for i in range(len(IPIP_0_0)):
            f.write('{} {} {} {}\n'.format(IPIP_0_0[i],EAEA_plus_0_0[i],IPIP_minus_1_0[i],EAEA_1_0[i]))
    f.close()

DOS_data=np.loadtxt('DOS_{}.txt'.format(figname))
IPIP_0_0       = DOS_data[:,0]
EAEA_plus_0_0  = DOS_data[:,1]
IPIP_minus_1_0 = DOS_data[:,2]
EAEA_1_0       = DOS_data[:,3]


Energy = np.linspace(Emin,Emax,NE)

# Plot DOS vs. Energy to see if smooth enough
PlotDOSPreview(IPIP_0_0,EAEA_plus_0_0,IPIP_minus_1_0,EAEA_1_0,Energy,figname)

#-----------------------------------
#    Generate template for Fig.4
#-----------------------------------
IP_h_mean = settings['layers'][0]['molecule_species'][0]['ionization_potential']['mean']
EA_d_mean = settings['layers'][0]['molecule_species'][1]['electron_affinity']['mean']

ylim_DOS_all=[-6.5, -3]
resc_IP_d=1
resc_EA_d=10
rescalex=1.2
PlotDOSDopant(IPIP_minus_1_0,EAEA_1_0,Energy, IP_h_mean, EA_d_mean, offset, figname,
           ylim_DOS_all,rescalex,resc_IP_d,resc_EA_d,
           labels=True)

print('Offset: {} eV'.format(offset))
print('Max EA_h: {} eV'.format(np.max(EAEA_plus_0_0)))
print('Max IP_h: {} eV'.format(np.max(IPIP_0_0)))
print('Max EA_d: {} eV'.format(np.max(EAEA_1_0)))
print('Max IP_d: {} eV'.format(np.max(IPIP_minus_1_0)))

