from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
import yaml


def PlotDOSPreview(IP_h,EA_h,IP_d,EA_d,Energy,figname,
                color_IP_h='blue',color_EA_h='red',
                color_IP_d='green',color_EA_d='orange'):
    plt.plot(Energy, IP_h, label = 'IP$_{h}$', color=color_IP_h,alpha=0.5,lw=2)# IP host 
    plt.plot(Energy, EA_h, label = 'EA$^+_{h}$', color=color_EA_h)# EA host
    plt.plot(Energy, IP_d, label = 'IP$^-_{d}$', color=color_IP_d)# IP dop
    plt.plot(Energy, EA_d, label = 'EA$_{d}$', color=color_EA_d)# EA dop
    plt.xlim([-7.5, -0])
    plt.legend()
    plt.ylabel('DOS')
    plt.xlabel('Energy [eV]')
    plt.savefig('{}_DOS_preview.pdf'.format(figname))
    plt.yscale('log')
    plt.savefig('{}_DOS_log_preview.pdf'.format(figname))
    plt.clf()

def PlotDOSAll(IP_h,EA_h,IP_d,EA_d,Energy,
               E_fermi,figname,ylim,rescalex=1.5,
               resc_IP_h=1,resc_EA_h=1,resc_IP_d=1,resc_EA_d=1,
               color_IP_h='blue',color_EA_h='red',color_IP_d='green',color_EA_d='orange',
               labels=True,linewidth=2):
    # Generate template for Fig.
    matplotlib.rcParams.update({'font.size': 16})
    fig,axes = plt.subplots(nrows=1,ncols=2,figsize = (4.5,7),sharey=True)
    fig.subplots_adjust(wspace=0)
    ax_h,ax_d=axes[0],axes[1]
    max_DOS_h=np.max([EA_h,IP_h])
    # Host
    # plot curves
    ax_h.plot(resc_IP_h*IP_h/max_DOS_h,         Energy, label = 'IP$_{h}$', color=color_IP_h,alpha=0.5,lw=linewidth)# IP host 
    ax_h.fill_between(resc_IP_h*IP_h/max_DOS_h, Energy, facecolor=color_IP_h,alpha=0.5)
    ax_h.plot(resc_EA_h*EA_h/max_DOS_h,    Energy, label = 'EA$^+_{h}$', color=color_EA_h,alpha=1,lw=linewidth)# EA host
    if labels:
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
    if labels:
        ax_d.set_xlabel('Dopant DOS')
    # determine and plot Fermi level
    ax_d.plot(np.linspace(0,rescalex*max_DOS_h,10), E_fermi*np.ones(10), label = 'E$_{\mathrm{F}}$',color='k',linestyle='-.',lw=linewidth)
    ax_h.plot(np.linspace(0,rescalex*max_DOS_h,10), E_fermi*np.ones(10), label = 'E$_{\mathrm{F}}$',color='k',linestyle='-.',lw=linewidth)
    if labels:
        fig.legend(loc=(0.91,0.11))
    
    plt.savefig('{}_DOS_all.svg'.format(figname))
    plt.savefig('{}_DOS_all.png'.format(figname))
    plt.clf()

def PlotDOSIP(IP_h,IP_d,Energy,E_fermi,disorder, dop_frac,IP_h_mean,EA_d_mean,
              figname,rescalex,ylim,
              color_IP_h='blue',color_IP_d='orange',labels=True,linewidth=2):
    matplotlib.rcParams.update({'font.size': 16})
    plt.figure(figsize = (4.5,4.5))
    max_IP_d=np.max(IP_d)
    max_IP_h=np.max(IP_h)
    # Rescale
    resc_IP_h=1/max_IP_h
    resc_IP_d=1/max_IP_d
    # plot curves
    # plt.plot(IP_d*resc_IP_d, Energy, label = 'IP$^-_{d}$',color=color_IP_d,alpha=0.5,lw=linewidth)# IP dop
    plt.fill_between(IP_d*resc_IP_d,Energy, facecolor=color_IP_d, alpha=0.5, label = 'HOMO$^-_{d}$')
    # plt.plot(resc_IP_h*IP_h,         Energy, label = 'IP$_{h}$', color=color_IP_h,alpha=0.5,lw=linewidth)# IP host 
    plt.fill_between(resc_IP_h*IP_h, Energy, facecolor=color_IP_h,alpha=0.5,label = 'HOMO$_{h}$')

    # plot mean
    # mean_IP_h = Energy[np.argmax(IP_h)]
    # mean_IP_d = Energy[np.argmax(IP_d)]
    # plt.plot(np.linspace(0,rescalex,10), mean_IP_h*np.ones(10), label = 'mean IP$_{h}$',color='k',linestyle=':')
    # plt.plot(np.linspace(0,rescalex,10), mean_IP_d*np.ones(10), label = 'mean IP$^-_{d}$',color='k',linestyle=':')
    # plot IP_h and EA_d-0.36
    plt.plot(np.linspace(0,rescalex,10), IP_h_mean*np.ones(10), label = 'HOMO$_{h}$(0)',color=color_IP_h,linestyle='--',lw=linewidth)
    plt.plot(np.linspace(0,rescalex,10), (EA_d_mean-0.36)*np.ones(10), label = 'LUMO$_{d}$ - |V$_{C}$|',color=color_IP_d,linestyle='--',lw=linewidth)
    # plot Fermi level
    plt.plot(np.linspace(0,rescalex,10), E_fermi*np.ones(10), label = 'E$_{F}$',color='k',linestyle='-.',lw=linewidth)
    if labels:
        plt.ylabel('Energy [eV]')
        plt.xlabel('DOS[a.u.]')
        plt.title('Disorder = {} eV, doping fraction = {} %'.format(disorder,dop_frac*100))
        plt.legend(framealpha=1)
    plt.ylim(ylim)
    plt.xlim(0,rescalex*1) # DOS normed to 1
    plt.xticks(np.linspace(0,1.0,6))
    
    plt.savefig('{}_DOS_IP.svg'.format(figname))
    plt.savefig('{}_DOS_IP.png'.format(figname))
    plt.clf()

def PlotDOSDopant(IP_d,EA_d,Energy,IP_h_mean,EA_d_mean,offset,
               figname,ylim,rescalex=1.5,
               resc_IP_d=1,resc_EA_d=1,
               color_IP_d='green',color_EA_d='orange',
               labels=True,linewidth=2):
    # Generate template for Fig.
    matplotlib.rcParams.update({'font.size': 16})
    plt.figure(figsize = (5,5))
    plt.subplots_adjust(left=0.2,bottom=0.15)
    # Dopant
    max_DOS_d=np.max([EA_d,IP_d])
    # plot curves
    plt.plot(IP_d*resc_IP_d/max_DOS_d, Energy, label = 'HOMO$^-_{d}$',color=color_IP_d,alpha=0.5,lw=linewidth)# IP dop
    plt.fill_between(IP_d*resc_IP_d/max_DOS_d,Energy, facecolor=color_IP_d, alpha=0.5)
    plt.plot(EA_d*resc_EA_d/max_DOS_d,       Energy, label = 'LUMO$_{d}$',  color=color_EA_d, alpha=1,lw=linewidth)# EA dop
    plt.plot(np.linspace(0,rescalex,10), (EA_d_mean-0.36)*np.ones(10), label = 'LUMO$_{d}^{(0)}$ - |V$_{C}$|',color=color_IP_d,linestyle='--',lw=linewidth)
    plt.plot(np.linspace(0,rescalex,10), EA_d_mean*np.ones(10), label = 'LUMO$_{d}^{(0)}$',color=color_EA_d,linestyle='--',lw=linewidth)
    plt.plot(np.linspace(0,rescalex,10), IP_h_mean*np.ones(10), label = 'HOMO$_{h}^{(0)}$',color='k',linestyle='--',lw=linewidth)
    plt.ylim(ylim)
    plt.xlim(0,rescalex)
    plt.xticks(np.linspace(0,1.0,6))
    if labels:
        plt.ylabel('Energy [eV]')
        plt.xlabel('Dopant DOS')
        plt.legend(framealpha=1)
    
    plt.savefig('{}_DOS_dopant.svg'.format(figname))
    plt.savefig('{}_DOS_dopant.png'.format(figname))
    plt.clf()