import numpy as np
import extract
import os
import shutil
import glob
import scipy.stats.mstats as stat
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
from sklearn.linear_model import LinearRegression
from scipy.constants import Boltzmann,e

class CurrTempSimulation:
    '''
    Analyse simulations for J(T), i.e. for several Temp. values and one 
    value for field, disorder, DMR.
    The directory must contain subdirectories in the format temp_*/r_*
    and a sim_data_file (default= 'sim_param.txt') containing the relevant 
    simulation parameters.
    '''
    # TODO: use yaml file for sim_param.txt
    # TODO: rewrite the code in a way that false format is recognised
    # TODO: add mode for writting the lambda (just a single value) for 
    #       Marcus rates in all plots and data files
    def __init__(self,source_dir='',dest_dir='analysis',sim_data_file='sim_param.txt'):
        sim_data = np.genfromtxt(source_dir+sim_data_file,dtype=float)
        self.field    = sim_data[1]
        self.dis      = sim_data[2]
        self.DMR      = sim_data[3]
        self.sys_size = sim_data[4]
        self.temp     = sim_data[5:]
        self.n_r      = int(sim_data[0])

        self.source_dir = source_dir
        self.dest_dir   = dest_dir

        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        self.current     = None
        self.std_current = None
        self.conduct     = None
        self.std_conduct = None
        self.act_energy  = None

    def get_current(self):
        '''
        Extract final current density from outputfile 'output_job_0' from 
        converged jobs and average for over all replicas for each temperature.
        The outputfile 'current.txt' contains (in columns): 
        - temperature values (K)
        - average current density values (A/m^2)
        - standard error values of current density (A/m^2)
        - logarithmic std. error of current density (A/m^2) 
          (to use for errorbars in logarithmic plots)
        
        '''
        n_temp      = len(self.temp)
        av_current  = np.zeros(n_temp)
        std_current = np.zeros(n_temp)
        log_std_current = np.zeros(n_temp)
        for i_t in range(n_temp):
            current=[]
            for i_r in range(self.n_r):
                outputdir=self.source_dir+'temp_{}/r_{}'.format(i_t, i_r)
                if os.path.exists(outputdir):
                    if 'final current density: ' in open(outputdir+'/output_job_0').read():
                        current_dens=extract.extract_e(outputdir+'/output_job_0', "final current density: ")
                        if current_dens<0:
                            print("Negative current density encountered at {}: {} A/m^2.\n Not counted into average current density.".format(outputdir,current_dens))
                        else:
                            current.append(current_dens)

            if len(current)>0:
                av_current[i_t]  = stat.gmean(current)
                std_current[i_t] = np.std(current)/np.sqrt(len(current))
                log_std_current[i_t] = np.std(np.log(current))/np.sqrt(len(current))

        self.current = av_current
        self.std_current = [std_current,log_std_current]
        # Write current, temperature, DMR to txt file
        if not os.path.isdir(self.dest_dir+'/data'):
            os.makedirs(self.dest_dir+'/data')
        with open(self.dest_dir+'/data/current.txt', 'w') as f:
                f.write('# Simulation for field = {} eV, disorder = {} eV, DMR = {}, system size = {}\n'.format(self.field,self.dis,self.DMR,self.sys_size))
                f.write('# Temperature(K)   Av. Current density(A/m$^2$)   Normal std. error(A/m$^2$)   Log. std. error(A/m$^2$)\n')
                for i_t in range(n_temp):
                    f.write('{0:<16}   {1:<26}   {2:<24}   {3:<30}\n'.format(self.temp[i_t], self.current[i_t], self.std_current[0][i_t], self.std_current[1][i_t]))
        f.close()


    def get_conductivity(self):
        '''
        Extract conductivity( = av. current density / field ) from outputfile 'current.txt'.
        The outputfile 'current.txt' contains (in columns): 
        - Temperature values (K) 
        - Average conductivity values (Av/eV/m^2)
        - Standard error values of current density (Av/eV/m^2)
        - Logarithmic std. error of current density (Av/eV/m^2)
          (to use for errorbars in logarithmic plots) 
        
        '''
        # check if CurrTempSimulation already has current etc. attribute, if not check if data is already generated and if yes read it in
        if not os.path.exists(self.dest_dir+'/data/current.txt'):
            self.get_current()   
        else:
            if np.any(self.current) == None:
                current_data = np.loadtxt(self.dest_dir+'/data/current.txt',comments='#',unpack=True)
                self.current, self.std_current = current_data[1], [current_data[2],current_data[3]] 
        self.conduct = self.current*1e-09/self.field
        self.std_conduct     = [self.std_current[0]*1e-09/self.field,self.std_current[1]*1e-09/self.field]
    
        with open(self.dest_dir+'/data/conductivity.txt', 'w') as f:
            f.write('# Simulation for field = {} eV, disorder = {} eV, DMR = {}, system size = {}\n'.format(self.field,self.dis,self.DMR,self.sys_size))
            f.write('# Temperature(K)   Av. Conductivity(A/eV/m$^2$)   Normal std. error(A/m$^2$)   Log. std. error(A/m$^2$)\n')
            n_temp      = len(self.temp)
            for i_t in range(n_temp):
                f.write('{0:<16}   {1:<26}   {2:<24}   {3:<30}\n'.format(self.temp[i_t], self.conduct[i_t], self.std_conduct[0][i_t], self.std_conduct[1][i_t]))    
        f.close()

    def get_conv_analysis(self):
        '''
        Analyse convergence parameters for each temperature value.
        The outputfile 'conv_analysis.txt' contains (in columns):
        - Temperature values (K)
        - Number of run, non-run, converged and non-converged jobs
        - Convergence, i.e. #conv. jobs/#run jobs (%)
        - Average, maximum, minimum and std. deviation of
          simulation time per replica (h)
        '''
        n_temp      = len(self.temp)
        # count run, non run, conv. and non conv. jobs and save to file    
        run=0
        conv=0
        tot_nonrun_jobs=[]
        tot_nonconv_jobs=[]
        # determine conv., perc. of run jobs and average convergence time for each system
        convergence   = np.zeros((n_temp))
        run_jobs_perc = np.zeros((n_temp))
        av_conv_time  = np.zeros((n_temp))
        min_conv_time = np.zeros((n_temp))
        max_conv_time = np.zeros((n_temp))
        std_conv_time = np.zeros((n_temp))
        analysis_data=[]
        for i_t in range(n_temp):
                run_jobs=[]
                nonrun_jobs=[]
                conv_jobs=[]
                nonconv_jobs=[]
                conv_times=[]
                for i_r in range(self.n_r):
                    job=str(i_r)
                    outputdir=self.source_dir+'temp_{}/r_{}'.format(i_t, i_r)
                    if os.path.exists(outputdir):
                        run_jobs.append(job)
                        run+=1
                        if 'final current density:' in open(outputdir+'/output_job_0').read():
                            conv+=1
                            conv_jobs.append(job)
                            if type(extract.extract_float(outputdir+'/output_job_0', "Time per simulation: ")) is float:
                                conv_times.append(extract.extract_float(outputdir+'/output_job_0', "Time per simulation: ")) # add convergence time
                        else:
                            nonconv_jobs.append(job)
                            tot_nonconv_jobs.append(job)
                    else:
                        nonrun_jobs.append(job)
                        tot_nonrun_jobs.append(job)      
                # write info about each system to txt file
                conv_times=np.array(conv_times)/3600                                                                                   
                run_jobs_perc[i_t]=100.0*len(run_jobs)/self.n_r   
                if len(run_jobs)>0:
                    convergence[i_t]=100.0*len(conv_jobs)/len(run_jobs)
                    if len(conv_jobs)>0:
                        if len(conv_times) != 0:
                            av_conv_time[i_t]=np.mean(conv_times) # take average over converged jobs
                            max_conv_time[i_t]=np.max(conv_times)
                            min_conv_time[i_t]=np.min(conv_times)
                            std_conv_time[i_t]=np.std(conv_times)
                    else:
                        av_conv_time[i_t]=0.0   # set av.conv. time to 0 as there is no convergence
                        max_conv_time[i_t]=0.0
                        min_conv_time[i_t]=0.0
                        std_conv_time[i_t]=0.0
                else:
                    convergence[i_t]=0         # set convergence to 0 if no job is run                        
                analysis_data.append(('{0:<3}   {1:<2}   {2:<2}   {3:<2}   {4:<2}   {5:<6}   {6:<20}   {7:<20}   {8:<20}   {9:<20}\n'.format( self.temp[i_t],
                                                                                                    len(run_jobs),
                                                                                                    len(nonrun_jobs),
                                                                                                    len(conv_jobs), 
                                                                                                    len(nonconv_jobs),
                                                                                                    convergence[i_t],
                                                                                                    av_conv_time[i_t],
                                                                                                    max_conv_time[i_t],
                                                                                                    min_conv_time[i_t],
                                                                                                    std_conv_time[i_t]
                                                                                                 )))
        if not os.path.isdir(self.dest_dir+'/data'):
            os.makedirs(self.dest_dir+'/data')        
        with open(self.dest_dir+'/data/conv_analysis.txt', 'w') as f:
            f.write('# Total number of jobs: {}\n'.format(n_temp*self.n_r))
            f.write('# Field: {} eV, Disorder: {} eV, Dop.fract.: {} %, System size: {} nm\n'.format(self.field,
                                                                                                   self.dis, 
                                                                                                   100.0*self.DMR,
                                                                                                   self.sys_size))
            f.write('# Number of replicas for each settings configuration: {} \n'.format(self.n_r))             
            f.write('# Overview:   {} % of jobs started to run ({}).\n'.format(100.0*run/n_temp/self.n_r,run))
            if run>0:
                f.write('#             {} % of run jobs have converged ({}). \n \n'.format(round(100.0*conv/run,1),conv))                
            f.write('# Temperature(K)     Number of jobs: run    non-run     conv.    non-conv.     Conv.(%)     Conv. time: Average(h)     Maximum(h)     Minimum(h)     Standarddev.(h)\n')
            for i_t in range(n_temp): f.write(analysis_data[i_t])           
        f.close()

    def get_act_energy(self,tlim_low=50,tlim_high=600,save_to_file=True):
        '''
        Fit log(J) vs. 1000/T with linear fit: y = a*x + b.
        Use specific temperature interval from tlim_low to tlim_high.
        The outputfile 'activation_energy.txt' contains the activation energy 
        (eV) for this temperature set (within the spec. limits).
        The plot 'current_lin_fit_dmr*.png' contains the log J vs. 1000/T
        with the curve from the linear regression.
        '''
        if not os.path.exists(self.dest_dir+'/data/current.txt'):
            self.get_current() 
        else:
            if np.any(self.current) == None:
                current_data = np.loadtxt(self.dest_dir+'/data/current.txt',comments='#',unpack=True)
                self.current, self.std_current = current_data[1], [current_data[2],current_data[3]]           
            if np.all(self.current == 0.0):
                print("Only zero values for current found, linear regretion on log(J) vs. 1/T not possible.")
            else:
                # check if CurrTempSimulation already has current etc. attribute, if not check if data is already generated and if yes read it in                      
                x_conv = np.where(self.current[:,]>0.0)
                current      = self.current[x_conv[0][0]:x_conv[0][-1]+1]
                std_current  = [self.std_current[0][x_conv[0][0]:x_conv[0][-1]+1],self.std_current[1][x_conv[0][0]:x_conv[0][-1]+1]]
                temperatures = self.temp[x_conv[0][0]:x_conv[0][-1]+1]  
                x = 1000/temperatures                  
                y = np.log(current)
                idx_low  = np.argwhere( temperatures >= tlim_low )[0,0]
                idx_high = np.argwhere( temperatures <= tlim_high )[-1,0]
                x_fit = (x[idx_low:idx_high+1]).reshape((-1,1))
                y_fit = y[idx_low:idx_high+1]
                if x_fit == [] or y_fit == []:
                    print("Data points to analyse don't lie within the specified temperature range, linear regretion on log(J) vs. 1/T not possible.")
                else:
                    linear_model = LinearRegression().fit(x_fit,y_fit)
                    a, b, r2 = linear_model.coef_[0], linear_model.intercept_, linear_model.score(x_fit,y_fit)
                    print('Linear fit: slope     = ', a)
                    print('            intercept = ',b)
                    print('                  R^2 = ',r2)
                    self.act_energy = -1000*1000*a*Boltzmann/e
                    print('Activation energy: E_A = {} meV'.format(self.act_energy))
                    if save_to_file:
                        if not os.path.isdir(self.dest_dir+'/data/act_energy'):
                            os.makedirs(self.dest_dir+'/data/act_energy')
                        with open(self.dest_dir+'/data/act_energy/activation_energy.txt','w') as f:
                            f.write('# Activation Energy (meV)\n{}'.format(self.act_energy))
                        f.close()

                    # Plot linear regression
                    plt.errorbar(x, y,yerr=std_current[1],linestyle = '',marker='o',color='r', label='Simulated points')
                    plt.plot(x,a*x + b, color='k',label='Linear fit: $y = a*x + b$\n$a = {0:.3}, b =  {1:.3}$\n$R^2 = {2:.3}$'.format(a,b,r2))
                    plt.legend()
                    plt.xlabel('1/T (1000/K)')
                    plt.ylabel('log J (A/m$^2$)')
                    plt.title('Simulation for field = {} eV, disorder = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.DMR,self.sys_size),fontsize=10)
                    if not os.path.isdir(self.dest_dir+'/act_energy'):
                        os.makedirs(self.dest_dir+'/act_energy')                    
                    plt.savefig(self.dest_dir+'/act_energy/current_lin_fit_dmr{}.png'.format(self.DMR))
                    if not os.path.isdir('analysis/act_energy'):
                        os.makedirs('analysis/act_energy')
                    plt.savefig('analysis/act_energy/current_lin_fit_dmr{}.png'.format(self.DMR))
                    plt.close()
    

    def plot_current(self,plot_log=True,errorbar=False, only_conv=True,curve_color='b', save_to_file=True):
        '''
        Plot current vs. 1000/T. 
        plot_log:     plot natural logarithm of current.
        errorbar:     display errorbars of each simulation points.
        only_conv:    plot only points for which at least one replica has 
                      converged.
        save_to_file: save the generated plot directly to 'current.png'.

        '''
        # check if CurrTempSimulation already has current etc. attribute, if not check if data is already generated and if yes read it in
        if not os.path.exists(self.dest_dir+'/data/current.txt'):
            self.get_current()   
        else:
            if np.any(self.current) == None:
                current_data = np.loadtxt(self.dest_dir+'/data/current.txt',comments='#',unpack=True)
                self.current, self.std_current = current_data[1], [current_data[2],current_data[3]]                 
            # only plot points with at least one converged job
            if only_conv:
                x = np.where(self.current[:,]>0.0)
                if len(x[0])>0:
                    current      = self.current[x[0][0]:x[0][-1]+1]
                    std_current  = [self.std_current[0][x[0][0]:x[0][-1]+1],self.std_current[1][x[0][0]:x[0][-1]+1]]
                    temperatures = self.temp[x[0][0]:x[0][-1]+1]
                else:
                    if plot_log:
                        current, std_current = np.ones(len(self.temp)), [np.ones(len(self.temp)),np.ones(len(self.temp))]
                    else:
                        current, std_current = self.current, self.std_current
                    temperatures = self.temp
            else:
                current      = self.current
                std_current  = self.std_current
                temperatures = self.temp
            # make logarithmic plot with appropiate std error
            if plot_log:
                if errorbar:
                    plt.errorbar(1000/temperatures,np.log(current),yerr=std_current[1],label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                else:
                    plt.plot(1000/temperatures,np.log(current),'-+',label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                ylabel = 'log J (A/m$^2$)'
            else:
                if errorbar:
                    plt.errorbar(1000/temperatures,current,yerr=std_current[0],label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                else:
                    plt.plot(1000/temperatures,current,'-+',label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                ylabel = 'J (A/m$^2$)'
            if save_to_file:
                plt.xlabel('1/T (1000/K)')
                plt.ylabel(ylabel)
                plt.title('Simulation for field = {} eV, disorder = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.DMR,self.sys_size),fontsize=10)
                if not os.path.isdir(self.dest_dir+'/plots'):
                    os.makedirs(self.dest_dir+'/plots')
                plt.savefig(self.dest_dir+'/plots/current.png')
                plt.close()

    def plot_conductivity(self, plot_log=True, errorbar=False, only_conv=True, curve_color='b', save_to_file=True):
        '''
        Plot conductivity vs. 1000/T. 
        plot_log:     plot natural logarithm of current.
        errorbar:     display errorbars of each simulation points.
        only_conv:    plot only points for which at least one replica has 
                      converged.
        save_to_file: save the generated plot directly to 'current.png'.

        '''
        # check if CurrTempSimulation already has current etc. attribute, if not check if data is already generated and if yes read it in
        if not os.path.exists(self.dest_dir+'/data/conductivity.txt'):
            print('Error: "{}" not found. Please generate data first.'.format(self.dest_dir+'/data/conductivity.txt'))   
        else:
            if np.any(self.conduct) == None:
                conduct_data = np.loadtxt(self.dest_dir+'/data/conductivity.txt',comments='#',unpack=True)
                self.conduct, self.std_conduct = conduct_data[1], [conduct_data[2],conduct_data[3]]
            # only plot points with at least one converged job
            if only_conv:
                x = np.where(self.conduct[:,]>0.0)
                if len(x[0])>0:
                    conduct      = self.conduct[x[0][0]:x[0][-1]+1]
                    std_conduct  = [self.std_conduct[0][x[0][0]:x[0][-1]+1],self.std_conduct[1][x[0][0]:x[0][-1]+1]]
                    temperatures = self.temp[x[0][0]:x[0][-1]+1]
                else:
                    if plot_log:
                        conduct, std_conduct = np.ones(len(self.temp)), [np.ones(len(self.temp)),np.ones(len(self.temp))]
                    else:
                        conduct, std_conduct = self.conduct, self.std_conduct
                    temperatures = self.temp                
            else:
                conduct      = self.conduct
                std_conduct  = self.std_conduct
                temperatures = self.temp
            # make logarithmic plot with appropiate std error
            if plot_log:
                if errorbar:
                    plt.errorbar(1000/temperatures,np.log(conduct),yerr=std_conduct[1],label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                else:
                    plt.plot(1000/temperatures,np.log(conduct),'-+',label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                ylabel = 'log $\sigma$ (A/eV/m)'
            else:
                if errorbar:
                    plt.errorbar(1000/temperatures,conduct,yerr=std_conduct[0],label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                else:
                    plt.plot(1000/temperatures,conduct,'-+',label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                ylabel = '$\sigma$ (A/eV/m)'
            if save_to_file:
                plt.xlabel('1/T (1000/K)')
                plt.ylabel(ylabel)
                plt.title('Simulation for field = {} eV, disorder = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.DMR,self.sys_size),fontsize=10)
                if not os.path.isdir(self.dest_dir+'/plots'):
                    os.makedirs(self.dest_dir+'/plots')
                plt.savefig(self.dest_dir+'/plots/conductivity.png')
                plt.close()

    def plot_conv_analysis(self,x_loc_scale=1.25,bar_loc=0.0,bar_width=0.25, 
                           bar_color='b',save_to_file=True,figure=None,axes=None):
        '''
        Plot convergence parameters for set of temperatures.
        The output plot 'conv_analysis.png' contains 6 subplots where
        - Convergence and number of run jobs
        - Average, maximum, minimum and std. dev. of simulation time
        are shown as bars for each temperature point.
        '''
        if not os.path.exists(self.dest_dir+'/data/conv_analysis.txt'):
            self.get_act_energy()  
        else:
            n_temp      = len(self.temp)
            conv_data = np.loadtxt(self.dest_dir+'/data/conv_analysis.txt', comments='#', unpack=True)
            temperatures                                              = conv_data[0]
            convergence                                               = conv_data[5]
            run_jobs_perc                                             = 100.0*conv_data[1]/self.n_r
            av_conv_time, max_conv_time, min_conv_time, std_conv_time = conv_data[6:]

            conv_quan = [[convergence, run_jobs_perc],
                        [av_conv_time, max_conv_time],
                        [min_conv_time, std_conv_time]]
            conv_desc = [['Convergence of run jobs (%)','Number of run jobs (%)'],
                        ['Average convergence time (h)','Maximum convergence time (h)'],
                        ['Minimum convergence time (h)','Standarddeviation convergence time (h)']]
            n_col, n_row = len(conv_quan[:]), len(conv_quan[0])
            # Visualize
            ind = np.arange(n_temp)*x_loc_scale  # the x locations for the groups
            if save_to_file:
                fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(2.5*n_temp*8/10, 10))
            else:
                fig, axes = figure, axes

            for i in range(n_col):
                for j in range(n_row):
                    if not save_to_file:
                        dmr_label = 'Dop.fract. {} %'.format(self.DMR*100)
                    else:
                        dmr_label = None
                    axes[j,i].bar(ind + bar_loc*bar_width/2, conv_quan[i][j], bar_width,color=bar_color,label=dmr_label)
                    axes[j,i].set_xlabel('Temperature (K)')
                    axes[j,i].set_ylabel(conv_desc[i][j])
                    axes[j,i].set_xticks(ind)
                    axes[j,i].set_xticklabels(temperatures)
                    axes[j,i].grid(axis="y")
                    axes[j,i].set_axisbelow(True)
                    if not save_to_file:
                        axes[j,i].legend(loc='lower left')
            if save_to_file:
                fig.suptitle('Convergence for field={} V/nm, disorder={} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.DMR,self.sys_size),size=10)
                if not os.path.isdir(self.dest_dir+'/plots'):
                    os.makedirs(self.dest_dir+'/plots')
                plt.savefig(self.dest_dir+'/plots/conv_analysis.png')
                plt.close()
            else:
                fig.suptitle('Convergence for field = {} V/nm, disorder = {} eV.'.format(self.field,self.dis),size=14)

    def get_nonconv_joblist(self, new_joblist_name = 'joblist_non_conv'):
        '''
        If simulation has not converged completely write a new joblist 
        with all temperatures that have not converged in order to restart 
        the simulation with this new joblist file.
        '''
        n_cpus = 1 # hardcoded
        n_temp = len(self.temp)
        n_replicas = self.n_r

        # find non converged and non run jobs and save to new joblist  
        jobfile = open(new_joblist_name+"_run_"+str(n_run), 'w')

        # count_nonconv_jobs = 0
        for i_t in reversed(range(n_temp)):
            for i_r in range(n_replicas):
                outputdir='{}temp_{}/r_{}'.format(self.source_dir,i_t, i_r)
                if os.path.exists(outputdir):
                    if not 'final current density:' in open(outputdir+'/output_job_0').read():
                        for i in range(n_replicas):
                            jobfile.write("{} JV.run_kmc {} {}\n".format(n_cpus, i_t, i))
                        break
                else:
                    for i in range(n_replicas):
                        jobfile.write("{} JV.run_kmc {} {}\n".format(n_cpus, i_t, i))
                    break
        print("Created new joblist: {}\n".format(new_joblist_name))

        
   

class CurrFieldSimulation:
    '''
    Analyse simulations for J(E), i.e. for several field values and one 
    value for temperature, disorder, DMR.
    The directory must contain subdirectories in the format f_*/r_*
    and a sim_data_file (default= 'sim_param.txt') containing the relevant 
    simulation parameters.

    Works analogously to CurrTempSimulation but without computation of
    activation energy.

    '''
    def __init__(self,source_dir='',dest_dir='analysis',sim_data_file='sim_param.txt'):
        sim_data = np.genfromtxt(source_dir+sim_data_file,dtype=float)
        self.n_r      = int(sim_data[0])
        self.temp     = sim_data[1]
        self.dis      = sim_data[2]
        self.DMR      = sim_data[3]
        self.sys_size = sim_data[4]
        self.field    = sim_data[5:]

        self.source_dir = source_dir
        self.dest_dir   = dest_dir

        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)

        self.current  = None
        self.std_current     = None
        self.conduct  = None
        self.std_conduct = None
        self.act_energy = None

    def get_current(self):
        '''
        Extract final current density from outputfiles from converged
        jobs and average for over all replicas for each temperature.
        This is done for one set of temperatures and only one field,
        disorder and DMR.
        The outputfile contains: 
        - average current density values(1st column) 
        - standard error values of current density (2nd column)
        - temperature values (3rd column)
        '''
        n_field      = len(self.field)
        av_current  = np.zeros(n_field)
        std_current = np.zeros(n_field)
        log_std_current = np.zeros(n_field)
        for i_f in range(n_field):
            current=[]
            for i_r in range(self.n_r):
                outputdir=self.source_dir+'f_{}/r_{}'.format(i_f, i_r)
                if os.path.exists(outputdir):
                    if 'final current density: ' in open(outputdir+'/output_job_0').read():
                        current_dens=extract.extract_e(outputdir+'/output_job_0', "final current density: ")
                        if current_dens<0:
                            print("Negative current density encountered at {}: {} A/m^2.\n Not counted into average current density.".format(outputdir,current_dens))
                        else:
                            current.append(current_dens)

            if len(current)>0:
                av_current[i_f]  = stat.gmean(current)
                std_current[i_f] = np.std(current)/np.sqrt(len(current))
                log_std_current[i_f] = np.std(np.log(current))/np.sqrt(len(current))

        self.current = av_current
        self.std_current = [std_current,log_std_current]
        # Write current, field, DMR to txt file
        if not os.path.isdir(self.dest_dir+'/data'):
            os.makedirs(self.dest_dir+'/data')
        with open(self.dest_dir+'/data/current.txt', 'w') as f:
                f.write('# Simulation for field = {} V/nm, disorder = {} eV, DMR = {}, system size = {}\n'.format(self.field,self.dis,self.DMR,self.sys_size))
                f.write('# Field(V/nm)   Av. Current density(A/m$^2$)   Normal std. error(A/m$^2$)   Log. std. error(A/m$^2$)\n')
                for i_f in range(n_field):
                    f.write('{0:<16}   {1:<26}   {2:<24}   {3:<30}\n'.format(self.field[i_f], self.current[i_f], self.std_current[0][i_f], self.std_current[1][i_f]))
        f.close()

    def get_conductivity(self):
        # TODO
        pass
    
    def get_conv_analysis(self):
        '''
        Analyse convergence(%); number of run jobs; average, maximum and minimum 
        convergence time (h); standarddeviation of convergence time (h) for each 
        temperature value for a given field, disorder and doping MR.
        '''
        n_field      = len(self.field)
        # count run, non run, conv. and non conv. jobs and save to file    
        run=0
        conv=0
        tot_nonrun_jobs=[]
        tot_nonconv_jobs=[]
        # determine conv., perc. of run jobs and average convergence time for each system
        convergence   = np.zeros((n_field))
        run_jobs_perc = np.zeros((n_field))
        av_conv_time  = np.zeros((n_field))
        min_conv_time = np.zeros((n_field))
        max_conv_time = np.zeros((n_field))
        std_conv_time = np.zeros((n_field))
        analysis_data=[]
        for i_f in range(n_field):
                run_jobs=[]
                nonrun_jobs=[]
                conv_jobs=[]
                nonconv_jobs=[]
                conv_times=[]
                for i_r in range(self.n_r):
                    job=str(i_r)
                    outputdir=self.source_dir+'f_{}/r_{}'.format(i_f, i_r)
                    if os.path.exists(outputdir):
                        run_jobs.append(job)
                        run+=1
                        if 'final current density:' in open(outputdir+'/output_job_0').read():
                            conv+=1
                            conv_jobs.append(job)
                            if type(extract.extract_float(outputdir+'/output_job_0', "Time per simulation: ")) is float:
                                conv_times.append(extract.extract_float(outputdir+'/output_job_0', "Time per simulation: ")) # add convergence time
                        else:
                            nonconv_jobs.append(job)
                            tot_nonconv_jobs.append(job)
                    else:
                        nonrun_jobs.append(job)
                        tot_nonrun_jobs.append(job)      
                # write info about each system to txt file
                conv_times=np.array(conv_times)/3600                                                                                   
                run_jobs_perc[i_f]=100.0*len(run_jobs)/self.n_r   
                if len(run_jobs)>0:
                    convergence[i_f]=100.0*len(conv_jobs)/len(run_jobs)
                    if len(conv_jobs)>0:
                        if len(conv_times) != 0:
                            av_conv_time[i_f]=np.mean(conv_times) # take average over converged jobs
                            max_conv_time[i_f]=np.max(conv_times)
                            min_conv_time[i_f]=np.min(conv_times)
                            std_conv_time[i_f]=np.std(conv_times)
                    else:
                        av_conv_time[i_f]=0.0   # set av.conv. time to 0 as there is no convergence
                        max_conv_time[i_f]=0.0
                        min_conv_time[i_f]=0.0
                        std_conv_time[i_f]=0.0
                else:
                    convergence[i_f]=0         # set convergence to 0 if no job is run                        
                analysis_data.append(('{0:<3}   {1:<2}   {2:<2}   {3:<2}   {4:<2}   {5:<6}   {6:<20}   {7:<20}   {8:<20}   {9:<20}\n'.format( self.field[i_f],
                                                                                                    len(run_jobs),
                                                                                                    len(nonrun_jobs),
                                                                                                    len(conv_jobs), 
                                                                                                    len(nonconv_jobs),
                                                                                                    convergence[i_f],
                                                                                                    av_conv_time[i_f],
                                                                                                    max_conv_time[i_f],
                                                                                                    min_conv_time[i_f],
                                                                                                    std_conv_time[i_f]
                                                                                                 )))
        if not os.path.isdir(self.dest_dir+'/data'):
            os.makedirs(self.dest_dir+'/data')
        with open(self.dest_dir+'/data/conv_analysis.txt', 'w') as f:
            f.write('# Total number of jobs: {}\n'.format(n_field*self.n_r))
            f.write('# Field: {} V/nm, Disorder: {} eV, Dop.fract.: {} %, System size: {} nm\n'.format(self.field,
                                                                                                   self.dis, 
                                                                                                   100.0*self.DMR,
                                                                                                   self.sys_size))
            f.write('# Number of replicas for each settings configuration: {} \n'.format(self.n_r))             
            f.write('# Overview:   {} % of jobs started to run ({}).\n'.format(100.0*run/n_field/self.n_r,run))
            if run>0:
                f.write('#             {} % of run jobs have converged ({}). \n \n'.format(round(100.0*conv/run,1),conv))                
            f.write('# Field(V/nm)     Number of jobs: run    non-run     conv.    non-conv.     Conv.(%)     Conv. time: Average(h)     Maximum(h)     Minimum(h)     Standarddev.(h)\n')
            for i_t in range(n_field): f.write(analysis_data[i_t])           
        f.close()

    def plot_current(self,plot_log=True,errorbar=False, only_conv=True,curve_color='g', save_to_file=True):
        # check if CurrTempSimulation already has current etc. attribute, if not check if data is already generated and if yes read it in
        if not os.path.exists(self.dest_dir+'/data/current.txt'):
            self.get_current()   
        else:
            if np.any(self.current) == None:
                current_data = np.loadtxt(self.dest_dir+'/data/current.txt',comments='#',unpack=True)
                self.current, self.std_current = current_data[1], [current_data[2],current_data[3]]                 
            # only plot points with at least one converged job
            if only_conv:
                x = np.where(self.current[:,]>0.0)
                if len(x[0])>0:
                    current      = self.current[x[0][0]:x[0][-1]+1]
                    std_current  = [self.std_current[0][x[0][0]:x[0][-1]+1],self.std_current[1][x[0][0]:x[0][-1]+1]]
                    fields = self.field[x[0][0]:x[0][-1]+1]
                else:
                    if plot_log:
                        current, std_current = np.ones(len(self.field)), [np.ones(len(self.field)),np.ones(len(self.field))]
                    else:
                        current, std_current = self.current, self.std_current
                    fields = self.field
            else:
                current      = self.current
                std_current  = self.std_current
                fields = self.field
            # make logarithmic plot with appropiate std error
            if plot_log:
                if errorbar:
                    plt.errorbar(fields,np.log(current),yerr=std_current[1],label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                else:
                    plt.plot(fields,np.log(current),'-+',label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                ylabel = 'log J (A/m$^2$)'
            else:
                if errorbar:
                    plt.errorbar(fields,current,yerr=std_current[0],label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                else:
                    plt.plot(fields,current,'-+',label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                ylabel = 'J (A/m$^2$)'
            if save_to_file:
                plt.xlabel('E (V/nm)')
                plt.ylabel(ylabel)
                plt.title('Sim. for temperature = {} K, disorder = {} eV, DMR = {}, system size = {} nm.'.format(self.temp,self.dis,self.DMR,self.sys_size),fontsize=10)
                if not os.path.isdir(self.dest_dir+'/plots'):
                    os.makedirs(self.dest_dir+'/plots')
                plt.savefig(self.dest_dir+'/plots/current.png')
                plt.close()

    def plot_conductivity(self):
        # TODO
        pass
    
    def plot_conv_analysis(self,x_loc_scale=1.25,bar_loc=0.0,bar_width=0.25, 
                           bar_color='g',save_to_file=True,figure=None,axes=None):
        
        if not os.path.exists(self.dest_dir+'/data/conv_analysis.txt'):
            print('Error: "{}" not found. Please generate data first.'.format(self.dest_dir+'/data/conv_analysis.txt'))   
        else:
            n_field      = len(self.field)
            conv_data = np.loadtxt(self.dest_dir+'/data/conv_analysis.txt', comments='#', unpack=True)

            temperatures                                              = conv_data[0]
            convergence                                               = conv_data[5]
            run_jobs_perc                                             = 100.0*conv_data[1]/self.n_r
            av_conv_time, max_conv_time, min_conv_time, std_conv_time = conv_data[6:]

            conv_quan = [[convergence, run_jobs_perc],
                        [av_conv_time, max_conv_time],
                        [min_conv_time, std_conv_time]]
            conv_desc = [['Convergence of run jobs (%)','Number of run jobs (%)'],
                        ['Average convergence time (h)','Maximum convergence time (h)'],
                        ['Minimum convergence time (h)','Standarddeviation convergence time (h)']]
            n_col, n_row = len(conv_quan[:]), len(conv_quan[0])
            # Visualize
            ind = np.arange(n_field)*x_loc_scale  # the x locations for the groups
            if save_to_file:
                fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(2.5*n_field*8/10, 10))
            else:
                fig, axes = figure, axes

            for i in range(n_col):
                for j in range(n_row):
                    if not save_to_file:
                        dmr_label = 'Dop.fract. {} %'.format(self.DMR*100)
                    else:
                        dmr_label = None
                    axes[j,i].bar(ind + bar_loc*bar_width/2, conv_quan[i][j], bar_width,color=bar_color,label=dmr_label)
                    axes[j,i].set_xlabel('Field (V/nm)')
                    axes[j,i].set_ylabel(conv_desc[i][j])
                    axes[j,i].set_xticks(ind)
                    axes[j,i].set_xticklabels(temperatures)
                    axes[j,i].grid(axis="y")
                    axes[j,i].set_axisbelow(True)
                    if not save_to_file:
                        axes[j,i].legend(loc='lower left')
            if save_to_file:
                fig.suptitle('Convergence for temperature={} K, disorder={} eV, DMR = {}, system size = {} nm.'.format(self.temp,self.dis,self.DMR,self.sys_size),size=10)
                if not os.path.isdir(self.dest_dir+'/plots'):
                    os.makedirs(self.dest_dir+'/plots')
                plt.savefig(self.dest_dir+'/plots/conv_analysis.png')
                plt.close()
            else:
                fig.suptitle('Convergence for temperature = {} K, disorder = {} eV.'.format(self.temp,self.dis),size=14)

    def get_nonconv_joblist(self, new_joblist_name = 'joblist_non_conv'):
        # TODO
        pass


class CurrSysSizeSimulation:
    '''
    Analyse simulations for J(E), i.e. for several system size values and
    one value for field, temperature, disorder, DMR.
    The directory must contain subdirectories in the format sys_*/r_*
    and a sim_data_file (default= 'sim_param.txt') containing the relevant 
    simulation parameters.

    Works analogously to CurrTempSimulation but without computation of
    activation energy.
    '''
    def __init__(self,source_dir='',dest_dir='analysis',sim_data_file='sim_param.txt'):
        sim_data = np.genfromtxt(sim_data_file,dtype=float)

        self.n_r      = int(sim_data[0])
        self.field    = sim_data[1]
        self.temp     = sim_data[2]
        self.dis      = sim_data[3]
        self.DMR      = sim_data[4]
        self.sys_size = sim_data[5:]

        self.source_dir = source_dir
        self.dest_dir   = dest_dir

        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        shutil.copyfile(sim_data_file,dest_dir+'/'+sim_data_file)

        self.current  = None
        self.std_current     = None


    def get_current(self):
        n_sys      = len(self.sys_size)
        av_current  = np.zeros(n_sys)
        std_current = np.zeros(n_sys)
        log_std_current = np.zeros(n_sys)
        for i_s in range(n_sys):
            current=[]
            for i_r in range(self.n_r):
                outputdir=self.source_dir+'sys_{}/r_{}'.format(i_s, i_r)
                if os.path.exists(outputdir):
                    if 'final current density: ' in open(outputdir+'/output_job_0').read():
                        current_dens=extract.extract_e(outputdir+'/output_job_0', "final current density: ")
                        if current_dens<0:
                            print("Negative current density encountered at {}: {} A/m^2.\n Not counted into average current density.".format(outputdir,current_dens))
                        else:
                            current.append(current_dens)

            if len(current)>0:
                av_current[i_s]  = stat.gmean(current)
                std_current[i_s] = np.std(current)/np.sqrt(len(current))
                log_std_current[i_s] = np.std(np.log(current))/np.sqrt(len(current))

        self.current = av_current
        self.std_current = [std_current,log_std_current]
        # Write current, temperature, DMR to txt file
        if not os.path.isdir(self.dest_dir+'/data'):
            os.makedirs(self.dest_dir+'/data')
        with open(self.dest_dir+'/data/current.txt', 'w') as f:
                f.write('# Simulation for field = {} eV, disorder = {} eV, DMR = {}, temperature = {} K\n'.format(self.field,self.dis,self.DMR,self.temp))
                f.write('# System size(nm)   Av. Current density(A/m$^2$)   Normal std. error(A/m$^2$)   Log. std. error(A/m$^2$)\n')
                for i_s in range(n_sys):
                    f.write('{0:<16}   {1:<26}   {2:<24}   {3:<30}\n'.format(self.sys_size[i_s], self.current[i_s], self.std_current[0][i_s], self.std_current[1][i_s]))
        f.close()
    
    def get_conv_analysis(self):
        n_sys      = len(self.sys_size)
        # count run, non run, conv. and non conv. jobs and save to file    
        run=0
        conv=0
        tot_nonrun_jobs=[]
        tot_nonconv_jobs=[]
        # determine conv., perc. of run jobs and average convergence time for each system
        convergence   = np.zeros((n_sys))
        run_jobs_perc = np.zeros((n_sys))
        av_conv_time  = np.zeros((n_sys))
        min_conv_time = np.zeros((n_sys))
        max_conv_time = np.zeros((n_sys))
        std_conv_time = np.zeros((n_sys))
        analysis_data=[]
        for i_s in range(n_sys):
                run_jobs=[]
                nonrun_jobs=[]
                conv_jobs=[]
                nonconv_jobs=[]
                conv_times=[]
                for i_r in range(self.n_r):
                    job=str(i_r)
                    outputdir=self.source_dir+'sys_{}/r_{}'.format(i_s, i_r)
                    if os.path.exists(outputdir):
                        run_jobs.append(job)
                        run+=1
                        if 'final current density:' in open(outputdir+'/output_job_0').read():
                            conv+=1
                            conv_jobs.append(job)
                            if type(extract.extract_float(outputdir+'/output_job_0', "Time per simulation: ")) is float:
                                conv_times.append(extract.extract_float(outputdir+'/output_job_0', "Time per simulation: ")) # add convergence time
                        else:
                            nonconv_jobs.append(job)
                            tot_nonconv_jobs.append(job)
                    else:
                        nonrun_jobs.append(job)
                        tot_nonrun_jobs.append(job)      
                # write info about each system to txt file
                conv_times=np.array(conv_times)/3600                                                                                   
                run_jobs_perc[i_s]=100.0*len(run_jobs)/self.n_r   
                if len(run_jobs)>0:
                    convergence[i_s]=100.0*len(conv_jobs)/len(run_jobs)
                    if len(conv_jobs)>0:
                        if len(conv_times) != 0:
                            av_conv_time[i_s]=np.mean(conv_times) # take average over converged jobs
                            max_conv_time[i_s]=np.max(conv_times)
                            min_conv_time[i_s]=np.min(conv_times)
                            std_conv_time[i_s]=np.std(conv_times)
                    else:
                        av_conv_time[i_s]=0.0   # set av.conv. time to 0 as there is no convergence
                        max_conv_time[i_s]=0.0
                        min_conv_time[i_s]=0.0
                        std_conv_time[i_s]=0.0
                else:
                    convergence[i_s]=0         # set convergence to 0 if no job is run                        
                analysis_data.append(('{0:<3}   {1:<2}   {2:<2}   {3:<2}   {4:<2}   {5:<6}   {6:<20}   {7:<20}   {8:<20}   {9:<20}\n'.format( self.sys_size[i_s],
                                                                                                    len(run_jobs),
                                                                                                    len(nonrun_jobs),
                                                                                                    len(conv_jobs), 
                                                                                                    len(nonconv_jobs),
                                                                                                    convergence[i_s],
                                                                                                    av_conv_time[i_s],
                                                                                                    max_conv_time[i_s],
                                                                                                    min_conv_time[i_s],
                                                                                                    std_conv_time[i_s]
                                                                                                 )))
        if not os.path.isdir(self.dest_dir+'/data'):
            os.makedirs(self.dest_dir+'/data')
        with open(self.dest_dir+'/data/conv_analysis.txt', 'w') as f:
            f.write('# Total number of jobs: {}\n'.format(n_sys*self.n_r))
            f.write('# Field: {} eV, Temperature: {} K, Disorder: {} eV, Dop.fract.: {} %\n'.format(self.field,
                                                                                                    self.temp,
                                                                                                    self.dis, 
                                                                                                    100.0*self.DMR))
            f.write('# Number of replicas for each settings configuration: {} \n'.format(self.n_r))             
            f.write('# Overview:   {} % of jobs started to run ({}).\n'.format(100.0*run/n_sys/self.n_r,run))
            if run>0:
                f.write('#             {} % of run jobs have converged ({}). \n \n'.format(round(100.0*conv/run,1),conv))                
            f.write('# System size(nm)     Number of jobs: run    non-run     conv.    non-conv.     Conv.(%)     Conv. time: Average(h)     Maximum(h)     Minimum(h)     Standarddev.(h)\n')
            for i_s in range(n_sys): f.write(analysis_data[i_s])           
        f.close()

    def plot_current(self,y_logmin=-15,y_logmax=15,plot_log=True,errorbar=False, only_conv=True,curve_color='b', save_to_file=True):
        # check if CurrTempSimulation already has current etc. attribute, if not check if data is already generated and if yes read it in
        if not os.path.exists(self.dest_dir+'/data/current.txt'):
            self.get_current()   
        else:
            if np.any(self.current) == None:
                current_data = np.loadtxt(self.dest_dir+'current.txt',comments='#',unpack=True)
                self.current, self.std_current = current_data[1], [current_data[2],current_data[3]]                 
            # only plot points with at least one converged job
            if only_conv:
                x = np.where(self.current[:,]>0.0)
                if len(x[0])>0:
                    current      = self.current[x[0][0]:x[0][-1]+1]
                    std_current  = [self.std_current[0][x[0][0]:x[0][-1]+1],self.std_current[1][x[0][0]:x[0][-1]+1]]
                    sys_size = self.sys_size[x[0][0]:x[0][-1]+1]
                else:
                    if plot_log:
                        current, std_current = np.ones(len(self.sys_size)), [np.ones(len(self.sys_size)),np.ones(len(self.sys_size))]
                    else:
                        current, std_current = self.current, self.std_current
                    sys_size = self.sys_size
            else:
                current      = self.current
                std_current  = self.std_current
                sys_size     = self.sys_size
            # make logarithmic plot with appropiate std error
            if plot_log:
                if errorbar:
                    plt.errorbar(sys_size,np.log(current),yerr=std_current[1],label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                else:
                    plt.plot(sys_size,np.log(current),'-+',label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                ylabel = 'log J (A/m$^2$)'
            else:
                if errorbar:
                    plt.errorbar(sys_size,current,yerr=std_current[0],label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                else:
                    plt.plot(sys_size,current,'-+',label='Dop.fract. {} %'.format(self.DMR*100),color=curve_color)
                ylabel = 'J (A/m$^2$)'
            if save_to_file:
                if not os.path.isdir(self.dest_dir+'/plots'):
                    os.makedirs(self.dest_dir+'/plots')
                plt.xlabel('system size (nm)')
                plt.ylim(y_logmin,y_logmax)
                plt.ylabel(ylabel)
                plt.title('Simulation for field = {} eV, temperature = {} K, disorder = {} eV, DMR = {}.'.format(self.field,self.temp,self.dis,self.DMR,self.sys_size),fontsize=10)
                plt.savefig(self.dest_dir+'/plots/current.png')
                plt.close()

    def plot_conv_analysis(self,x_loc_scale=1.25,bar_loc=0.0,bar_width=0.25, 
                           bar_color='b',dest_dir='',save_to_file=True,figure=None,axes=None):
        n_sys      = len(self.sys_size)
        if not os.path.exists(self.dest_dir+'/data/conv_analysis.txt'):
            self.get_current()   
        else:
            conv_data = np.loadtxt(self.dest_dir+'/data/conv_analysis.txt', comments='#', unpack=True)

        sys_sizes                                                 = conv_data[0]
        convergence                                               = conv_data[5]
        run_jobs_perc                                             = 100.0*conv_data[1]/self.n_r
        av_conv_time, max_conv_time, min_conv_time, std_conv_time = conv_data[6:]

        conv_quan = [[convergence, run_jobs_perc],
                     [av_conv_time, max_conv_time],
                     [min_conv_time, std_conv_time]]
        conv_desc = [['Convergence of run jobs (%)','Number of run jobs (%)'],
                     ['Average convergence time (h)','Maximum convergence time (h)'],
                     ['Minimum convergence time (h)','Standarddeviation convergence time (h)']]
        n_col, n_row = len(conv_quan[:]), len(conv_quan[0])
        # Visualize
        ind = np.arange(n_sys)*x_loc_scale  # the x locations for the groups
        if save_to_file:
            fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(2.5*n_sys*8/10, 10))
        else:
            fig, axes = figure, axes

        for i in range(n_col):
            for j in range(n_row):
                if not save_to_file:
                    dmr_label = 'Dop.fract. {} %'.format(self.DMR*100)
                else:
                    dmr_label = None
                axes[j,i].bar(ind + bar_loc*bar_width/2, conv_quan[i][j], bar_width,color=bar_color,label=dmr_label)
                axes[j,i].set_xlabel('System size (nm)')
                axes[j,i].set_ylabel(conv_desc[i][j])
                axes[j,i].set_xticks(ind)
                axes[j,i].set_xticklabels(sys_sizes)
                axes[j,i].grid(axis="y")
                axes[j,i].set_axisbelow(True)
                if not save_to_file:
                    axes[j,i].legend(loc='lower left')
        if save_to_file:
            fig.suptitle('Convergence for field={} V/nm, temperature = {} K, disorder = {} eV, DMR = {}.'.format(self.field,self.temp,self.dis,self.DMR),size=10)
            if not os.path.isdir(self.dest_dir+'/plots'):
                os.makedirs(self.dest_dir+'/plots')
            plt.savefig(self.dest_dir+'/plots/conv_analysis.png')
            plt.close()
        else:
            fig.suptitle('Convergence for field = {} V/nm, temperature = {} K, disorder = {} eV.'.format(self.field,self.temp,self.dis),size=14)

class CurrTempDMRset:
    '''
    Analyse J(T)-sets for different doping molar ratios (DMR)
    and one value for field and disorder.
    The directory must contain subdirectories in the format DMR_*/temp_*/r_*
    and the J(T)-sets have to consist of the same sim. parameters.
    '''
    def __init__(self,source_dir='',dest_dir='analysis'):
        '''
        Initialize attributes of DMRset: number of replicas, field, 
        disorder, doping molar ratios (DMRs), system sizes, temperatures.
        Initialize colorset for plots.
        '''
        self.source_dir = source_dir
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)        
        self.dest_dir   = dest_dir+'/'
        self.list_DMR_dirs = sorted(glob.glob('DMR_*'))
        if self.list_DMR_dirs == []:
            print("SimAnalysisError: No compatible Simulations found.")
        else:
            print(self.list_DMR_dirs)
            for i_dmr,dir in enumerate(self.list_DMR_dirs):
                dmr_sim = CurrTempSimulation(source_dir=dir+'/',dest_dir='analysis/'+dir)
                if i_dmr==0:
                    self.n_r      = dmr_sim.n_r
                    self.field    = dmr_sim.field
                    self.dis      = dmr_sim.dis
                    self.DMR      = [dmr_sim.DMR]
                    self.sys_size = [dmr_sim.sys_size]
                    self.temp     = dmr_sim.temp
                else:
                    if dmr_sim.n_r != self.n_r  or dmr_sim.field != self.field or dmr_sim.dis != self.dis or np.any(dmr_sim.temp != self.temp):
                        print("SimAnalysisError: default sim. parameters of CurrTemp simulations do not agree, DMR set cannot be analysed.")
                        break
                    else:
                        self.DMR.append(dmr_sim.DMR)
                        self.sys_size.append(dmr_sim.sys_size)

        self.act_energy = None 
        self.colorset   = np.array(['salmon','orangered','darkred',
                                       'violet','blueviolet','darkblue',
                                       'royalblue','dodgerblue','mediumseagreen',
                                       'green'])

        sim_data = open(self.dest_dir+'/sim_param_DMR_set.txt','w')
        sim_data.write("# n_replicas\n{}\n".format(self.n_r))
        sim_data.write("# field\n{}\n".format(self.field))
        sim_data.write("# disorder\n{}\n".format(self.dis))
        sim_data.write("# doping ratios\n")
        for dmr in self.DMR: sim_data.write("{}\n".format(dmr))
        sim_data.write("# system size\n")
        for sys in self.sys_size: sim_data.write("{}\n".format(sys))
        sim_data.write("# temperatures\n")
        for t in self.temp: sim_data.write("{}\n".format(t))



    def plot_current(self,errorbar=False,plot_log=True, only_conv=True):
        DMR_dir = self.list_DMR_dirs
        curr_dmr_colors  = self.colorset
        for i_dmr,dir in enumerate(DMR_dir):
            dmr_sim = CurrTempSimulation(source_dir=dir+'/',dest_dir='analysis/'+dir)
            if not os.path.exists(self.dest_dir+dir+"/data/current.txt"):    
                dmr_sim.get_current()
            dmr_sim.plot_current( plot_log, errorbar, only_conv, curr_dmr_colors[i_dmr], False)
        if plot_log:
            ylabel = 'log J (A/m$^2$)'
        else:
            ylabel = 'J (A/m$^2$)'
        plt.xlabel('1/T (1000/K)')
        plt.title('Simulation for field = {} V/nm, disorder = {} eV.'.format(self.field,self.dis),fontsize=14)
        plt.ylabel(ylabel)
        plt.legend()    
        plt.savefig(self.dest_dir+'/current_DMR_set.png')
        plt.close()

    def plot_conductivity(self,errorbar=False,plot_log=True, only_conv=True,log10=False):
        DMR_dir = self.list_DMR_dirs
        curr_dmr_colors  = self.colorset
        for i_dmr,dir in enumerate(DMR_dir):
            dmr_sim = CurrTempSimulation(source_dir=dir+'/',dest_dir='analysis/'+dir)
            if not os.path.exists(self.dest_dir+dir+"/data/conductivity.txt"):    
                dmr_sim.get_conductivity()
            dmr_sim.plot_conductivity( plot_log, errorbar, only_conv, curr_dmr_colors[i_dmr], False)
        if plot_log:
            ylabel = 'log $\sigma$ (A/eV/m)'
        else:
            ylabel = '$\sigma$ (A/eV/m)'
        if log10:
            plt.yscale("log")
        plt.xlabel('1/T (1000/K)')
        plt.ylabel(ylabel)
        plt.legend()
        plt.title('Simulation for field = {} V/nm, disorder = {} eV.'.format(dmr_sim.field,dmr_sim.dis),fontsize=14)
        plt.savefig(self.dest_dir+'/conductivity_DMR_set.png')
        plt.close()


    def plot_conv_analysis(self):
        DMR_dir = self.list_DMR_dirs
        n_DMR    = len(self.DMR)
        # General bar plot settings
        x_loc_scale = n_DMR*0.4
        width   = round(0.25*n_DMR/4,2)  # the width of the bars
        bar_loc = np.linspace(-3,3,n_DMR)
        curr_dmr_colors  = self.colorset
        
        for i_dmr, dir in enumerate(DMR_dir):
            dmr_sim = CurrTempSimulation(source_dir=dir+'/',dest_dir='analysis/'+dir)
            if not os.path.exists(self.dest_dir+dir+"/data/conv_analysis.txt"):    
                dmr_sim.get_conv_analysis()
            if i_dmr == 0:
                n_temp = len(self.temp)
                fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(2.5*n_temp*8/10, 10))
            dmr_sim.plot_conv_analysis(x_loc_scale=x_loc_scale,bar_loc=bar_loc[i_dmr],bar_width=width,
                                    bar_color=curr_dmr_colors[i_dmr],
                                    save_to_file=False,figure=fig,axes=axes)

        plt.savefig(self.dest_dir+'/conv_analysis_DMR_set.png')
        plt.close()  

    def get_act_energy(self,tlim_low=250,tlim_high=350):
        DMR_dir = sorted(glob.glob('DMR_*'))
        n_DMR    = len(DMR_dir)  
        act_energy = np.zeros(n_DMR)
        dmr        = self.DMR * 100
        for i_dmr, dir in enumerate(DMR_dir):
            dmr_sim = CurrTempSimulation(source_dir=dir+'/',dest_dir='analysis/'+dir)
            dmr_sim.get_act_energy(tlim_low=tlim_low,tlim_high=tlim_high)
            act_energy[i_dmr] = np.loadtxt(self.dest_dir+dir+"/data/act_energy/activation_energy.txt",unpack=True)
        if not os.path.isdir(self.dest_dir+'/act_energy'):
            os.makedirs(self.dest_dir+'/act_energy')       
        with open(self.dest_dir+"act_energy/act_energy_DMR_set.txt",'w') as f:
            f.write('# Simulation for field = {} V/nm, disorder = {} eV\n'.format(dmr_sim.field,dmr_sim.dis))
            f.write('# DMR(%)   Activation energy(meV)\n')
            for i in range(n_DMR):f.write('{}   {}\n'.format(dmr[i],act_energy[i]))
        f.close()
        self.act_energy = act_energy

    def plot_act_energy(self):
        if not os.path.exists(self.dest_dir+"act_energy/act_energy_DMR_set.txt"):
            self.get_act_energy()     
        act_energy = (np.loadtxt(self.dest_dir+"act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[1]
        dmr        = self.DMR
        plt.plot(dmr,act_energy)
        plt.xlabel('DMR (%)')
        plt.xscale("log")
        plt.ylabel('$E_A$ (meV)')
        plt.title('Simulation for field = {} V/nm, disorder = {} eV'.format(self.field,self.dis))
        plt.savefig(self.dest_dir+'act_energy/act_energy_DMR_set.png')
        plt.close()        

class CurrFieldDMRset:
    '''
    Analyse J(E)-sets for different doping molar ratios (DMR)
    and one value for field and disorder.
    The directory must contain subdirectories in the format DMR_*/temp_*/r_*
    and the J(E)-sets have to consist of the same sim. parameters.
    '''
    def __init__(self,source_dir='',dest_dir='analysis'):
        '''
        Initialize attributes of DMRset: number of replicas, field, 
        disorder, doping molar ratios (DMRs), system sizes, temperatures.
        Initialize colorset for plots.
        '''
        self.source_dir = source_dir
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)        
        self.dest_dir   = dest_dir+'/'
        self.list_DMR_dirs = sorted(glob.glob('DMR_*'))
        if self.list_DMR_dirs == []:
            print("SimAnalysisError: No compatible Simulations found.")
        else:
            print(self.list_DMR_dirs)
            for i_dmr,dir in enumerate(self.list_DMR_dirs):
                dmr_sim = CurrFieldSimulation(source_dir=dir+'/',dest_dir='analysis/'+dir)
                if i_dmr==0:
                    self.n_r      = dmr_sim.n_r
                    self.temp     = dmr_sim.temp
                    self.dis      = dmr_sim.dis
                    self.DMR      = [dmr_sim.DMR]
                    self.sys_size = [dmr_sim.sys_size]
                    self.field    = dmr_sim.field
                else:
                    if dmr_sim.n_r != self.n_r  or dmr_sim.temp != self.temp or dmr_sim.dis != self.dis or np.any(dmr_sim.field != self.field):
                        print("SimAnalysisError: default sim. parameters of CurrField simulations do not agree, DMR set cannot be analysed.")
                        break
                    else:
                        self.DMR.append(dmr_sim.DMR)
                        self.sys_size.append(dmr_sim.sys_size)

        self.act_energy = None 
        self.colorset   = np.array(['salmon','orangered','darkred',
                                       'violet','blueviolet','darkblue',
                                       'royalblue','dodgerblue','mediumseagreen',
                                       'green'])

        sim_data = open(self.dest_dir+'/sim_param_DMR_set.txt','w')
        sim_data.write("# n_replicas\n{}\n".format(self.n_r))
        sim_data.write("# temperature\n{}\n".format(self.temp))
        sim_data.write("# disorder\n{}\n".format(self.dis))
        sim_data.write("# doping ratios\n")
        for dmr in self.DMR: sim_data.write("{}\n".format(dmr))
        sim_data.write("# system size\n")
        for sys in self.sys_size: sim_data.write("{}\n".format(sys))
        sim_data.write("# fields\n")
        for f in self.field: sim_data.write("{}\n".format(f))



    def plot_current(self,errorbar=False,plot_log=True, only_conv=True):
        DMR_dir = self.list_DMR_dirs
        curr_dmr_colors  = self.colorset
        for i_dmr,dir in enumerate(DMR_dir):
            dmr_sim = CurrFieldSimulation(source_dir=dir+'/',dest_dir='analysis/'+dir)
            if not os.path.exists(self.dest_dir+dir+"/data/current.txt"):    
                dmr_sim.get_current()
            dmr_sim.plot_current( plot_log, errorbar, only_conv, curr_dmr_colors[i_dmr], False)
        if plot_log:
            ylabel = 'log J (A/m$^2$)'
        else:
            ylabel = 'J (A/m$^2$)'
        plt.xlabel('E (V/nm)')
        plt.title('Simulation for temperature = {} K, disorder = {} eV.'.format(self.temp,self.dis),fontsize=14)
        plt.ylabel(ylabel)
        plt.legend()    
        plt.savefig(self.dest_dir+'/current_DMR_set.png')
        plt.close()

    def plot_conductivity(self,errorbar=False,plot_log=True, only_conv=True,log10=False):
        # TODO: CurrFieldSimulation.get_conductivity() (if necessary)
        pass


    def plot_conv_analysis(self):
        DMR_dir = self.list_DMR_dirs
        n_DMR    = len(self.DMR)
        # General bar plot settings
        x_loc_scale = n_DMR*0.4
        width   = round(0.25*n_DMR/4,2)  # the width of the bars
        bar_loc = np.linspace(-3,3,n_DMR)
        curr_dmr_colors  = self.colorset
        
        for i_dmr, dir in enumerate(DMR_dir):
            dmr_sim = CurrFieldSimulation(source_dir=dir+'/',dest_dir='analysis/'+dir)
            if not os.path.exists(self.dest_dir+dir+"/data/conv_analysis.txt"):    
                dmr_sim.get_conv_analysis()
            if i_dmr == 0:
                n_field = len(self.field)
                fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(2.5*n_field*8/10, 10))
            dmr_sim.plot_conv_analysis(x_loc_scale=x_loc_scale,bar_loc=bar_loc[i_dmr],bar_width=width,
                                    bar_color=curr_dmr_colors[i_dmr],
                                    save_to_file=False,figure=fig,axes=axes)

        plt.savefig(self.dest_dir+'/conv_analysis_DMR_set.png')
        plt.close()  

        