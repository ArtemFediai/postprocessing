import numpy as np
import extract
import os
import shutil
import glob
import scipy.stats.mstats as stat
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
plt.switch_backend('Agg')
from sklearn.linear_model import LinearRegression
from scipy.constants import Boltzmann,e

# TODO: use try except statements instead of if else.
# TODO: use yaml file for sim_param.txt
# TODO: rewrite the code in a way that false format is recognised
# TODO: add analyse_Ecoul in case of miller rates

class CurrTempSimulation:
    '''
    Analyse simulations for J(T), i.e. for several Temp. values and one 
    value for field, disorder, DMR.
    The directory must contain subdirectories in the format temp_*/r_*
    and a sim_data_file (default= 'sim_param.txt') containing the relevant 
    simulation parameters.
    '''
    def __init__(self,rates="Miller",analyse_Ecoul=False, source_dir='',dest_dir='analysis',sim_data_file='sim_param.txt'):
        sim_data = np.genfromtxt(source_dir+sim_data_file,dtype=float)
        self.rates         = rates
        self.analyse_Ecoul = analyse_Ecoul
        self.n_r      = int(sim_data[0])
        self.field    = sim_data[1]
        self.dis      = sim_data[2]
        if self.rates=="Miller":
            self.DMR      = sim_data[3]
            self.sys_size = sim_data[4]
            self.temp     = sim_data[5:]
        elif self.rates=="Marcus":
            self.lam      = sim_data[3]
            if analyse_Ecoul:
                self.Ecoul    = sim_data[4]
                self.DMR      = sim_data[5]
                self.sys_size = sim_data[6]
                self.temp     = sim_data[7:]
            else:               
                self.DMR      = sim_data[4]
                self.sys_size = sim_data[5]
                self.temp     = sim_data[6:]
        else:
            print("SimAnalysisError: False rates specified.")
            exit()

        self.source_dir = source_dir
        self.dest_dir   = dest_dir+"/"

        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        self.current     = None
        self.std_current = None
        self.conduct     = None
        self.std_conduct = None
        self.act_energy  = None

    def collect_current_data(self):
        '''
        Extract final current density from outputfile 'output_job_0' from 
        converged jobs for every replica. If no final current is found, 
        current is set to 0.
        The outputfile 'J_temp_*.txt' contains (in columns): 
        - replica number
        - final current density values (A/m^2)
        
        '''
        n_temp      = len(self.temp)
        current_data = np.zeros((n_temp,self.n_r))
        if not os.path.isdir(self.dest_dir+'current_data'):
            os.makedirs(self.dest_dir+'current_data')
        for i_t in range(n_temp):
            for i_r in range(self.n_r):
                outputdir=self.source_dir+'temp_{}/r_{}'.format(i_t, i_r)
                if os.path.exists(outputdir+'/output_job_0'):
                    if 'final current density: ' in open(outputdir+'/output_job_0').read():
                        current_dens=extract.extract_e(outputdir+'/output_job_0', "final current density: ")
                        if current_dens<0:
                            print("Negative current density encountered at {}: {} A/m^2.\n Not counted to current densities.".format(outputdir,current_dens))
                        else:
                            current_data[i_t,i_r] = current_dens
            with open(self.dest_dir+'current_data/curr_temp_{}.txt'.format(i_t), 'w') as f: 
                    f.write("# Raw current data for temp_{} (T={} K)\n".format(i_t,self.temp[i_t]))
                    if self.rates == "Marcus":
                        if self.analyse_Ecoul:  
                            f.write('# Field = {} V/nm, disorder = {} eV, lambda = {} eV, Ecoul = {} eV, DMR = {}, system size = {} nm\n'.format(self.field,self.dis,self.lam, self.Ecoul,self.DMR,self.sys_size))    
                        else:
                            f.write('# Field = {} V/nm, disorder = {} eV, lambda = {} eV, DMR = {}, system size = {} nm\n'.format(self.field,self.dis,self.lam,self.DMR,self.sys_size))
                    else:
                        f.write('# Field = {} V/nm, disorder = {} eV, DMR = {}, system size = {}\n'.format(self.field,self.dis,self.DMR,self.sys_size))
                    f.write('#\n# Replica Nr.     Current density(A/m$^2$)\n')
                    for i_r in range(self.n_r):
                        f.write('  {0:<13}   {1:<1}\n'.format(i_r, current_data[i_t,i_r]))
            f.close()

    def get_av_current(self):
        '''
        Average final current density over all replicas for each temperature.
        Get data from '/current_data'.
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
            if not os.path.exists(self.dest_dir+'current_data/curr_temp_{}.txt'.format(i_t)):
                self.collect_current_data()
            current = (np.loadtxt(self.dest_dir+'current_data/curr_temp_{}.txt'.format(i_t)))[:,1]
            current = current[np.nonzero(current)]
            # Calculate average current
            if len(current)>0:
                av_current[i_t]  = stat.gmean(current)
                std_current[i_t] = np.std(current)/np.sqrt(np.count_nonzero(current))
                log_std_current[i_t] = np.std(np.log(current))/np.sqrt(np.count_nonzero(current))

        self.current = av_current
        self.std_current = [std_current,log_std_current]
        # Write current, temperature, DMR to txt file
        with open(self.dest_dir+'av_current.txt', 'w') as f:
            f.write("# Average current and related values.\n")
            if self.rates == "Marcus":
                if self.analyse_Ecoul:
                    f.write('# Field = {} eV, disorder = {} eV, lambda = {} eV, Ecoul = {} eV, DMR = {}, system size = {} nm\n'.format(self.field,self.dis,self.lam,self.Ecoul,self.DMR,self.sys_size))
                else:
                    f.write('# Field = {} eV, disorder = {} eV, lambda = {} eV, DMR = {}, system size = {} nm\n'.format(self.field,self.dis,self.lam,self.DMR,self.sys_size))
            else:
                f.write('# Field = {} eV, disorder = {} eV, DMR = {}, system size = {} nm\n'.format(self.field,self.dis,self.DMR,self.sys_size))
            f.write('#\n# Temperature(K)   Av. Current density(A/m$^2$)   Normal std. error(A/m$^2$)   Log. std. error(A/m$^2$)\n')
            for i_t in range(n_temp):
                f.write('  {0:<14}   {1:<28}   {2:<26}   {3:<1}\n'.format(self.temp[i_t], self.current[i_t], self.std_current[0][i_t], self.std_current[1][i_t]))
        f.close()

    def get_av_conductivity(self):
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
        if not os.path.exists(self.dest_dir+'av_current.txt'):
            self.get_av_current()   
        else:
            if np.any(self.current) == None:
                current_data = np.loadtxt(self.dest_dir+'av_current.txt',comments='#',unpack=True)
                self.current, self.std_current = current_data[1], [current_data[2],current_data[3]] 
        self.conduct = self.current*1e-09/self.field
        self.std_conduct     = [self.std_current[0]*1e-09/self.field,self.std_current[1]*1e-09/self.field]
    
        with open(self.dest_dir+'av_conductivity.txt', 'w') as f:
            f.write("# Average conductivity and related values.\n")
            if self.rates == "Marcus":
                if self.analyse_Ecoul:
                    f.write('# Field = {} eV, disorder = {} eV, lambda = {} eV, Ecoul = {} eV, DMR = {}, system size = {} nm\n'.format(self.field,self.dis,self.lam,self.Ecoul,self.DMR,self.sys_size))
                else:
                    f.write('# Field = {} eV, disorder = {} eV, lambda = {} eV, DMR = {}, system size = {} nm\n'.format(self.field,self.dis,self.lam,self.DMR,self.sys_size))
            else:
                f.write('# Field = {} eV, disorder = {} eV, DMR = {}, system size = {} nm\n'.format(self.field,self.dis,self.DMR,self.sys_size))
            f.write('#\n# Temperature(K)   Av. Conductivity(A/eV/m$^2$)   Normal std. error(A/m$^2$)   Log. std. error(A/m$^2$)\n')
            n_temp      = len(self.temp)
            for i_t in range(n_temp):
                f.write('  {0:<14}   {1:<28}   {2:<26}   {3:<1}\n'.format(self.temp[i_t], self.conduct[i_t], self.std_conduct[0][i_t], self.std_conduct[1][i_t]))    
        f.close()

    def get_conv_analysis(self):
        # TODO: collect here also sim.time and convergence for each replica
        #       as for current
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
                    if os.path.exists(outputdir+'/output_job_0'):
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
                            av_conv_time[i_t]=np.round(np.mean(conv_times),4) # take average over converged jobs
                            max_conv_time[i_t]=np.round(np.max(conv_times),4)
                            min_conv_time[i_t]=np.round(np.min(conv_times),4)
                            std_conv_time[i_t]=np.round(np.std(conv_times),4)
                    else:
                        av_conv_time[i_t]=0.0   # set av.conv. time to 0 as there is no convergence
                        max_conv_time[i_t]=0.0
                        min_conv_time[i_t]=0.0
                        std_conv_time[i_t]=0.0
                else:
                    convergence[i_t]=0         # set convergence to 0 if no job is run                        
                analysis_data.append(('  {0:<14}   {1:<4}   {2:<7}   {3:<6}   {4:<11}   {5:<12}   {6:<19}   {7:<20}   {8:<20}   {9:<20}\n'.format( self.temp[i_t],
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
        with open(self.dest_dir+'conv_analysis.txt', 'w') as f:
            f.write('# Convergence analysis.\n')
            f.write('# Total number of jobs: {}\n'.format(n_temp*self.n_r))
            if self.rates == "Marcus":
                if self.analyse_Ecoul:
                    f.write('# Field: {} eV, Disorder: {} eV, Lambda: {} eV, Ecoul: {} eV, Dop.fract.: {} %, System size: {} nm\n'.format(self.field,self.dis,self.lam,self.Ecoul,100.0*self.DMR,self.sys_size))
                else:
                    f.write('# Field: {} eV, Disorder: {} eV, Lambda: {} eV, Dop.fract.: {} %, System size: {} nm\n'.format(self.field,self.dis,self.lam,100.0*self.DMR,self.sys_size))
            else:
                f.write('# Field: {} eV, Disorder: {} eV, Dop.fract.: {} %, System size: {} nm\n'.format(self.field,
                                                                                                         self.dis, 
                                                                                                         100.0*self.DMR,
                                                                                                         self.sys_size))             
            f.write('# Overview:   {} % of jobs started to run ({}).\n'.format(100.0*run/n_temp/self.n_r,run))
            if run>0:
                f.write('#             {} % of run jobs have converged ({}). \n \n'.format(round(100.0*conv/run,1),conv))                
            f.write('#\n# Temperature(K)   run    non-run   conv.    non-conv.     Conv.(%)  Conv. time: Average(h)     Maximum(h)     Minimum(h)     Standarddev.(h)\n')
            for i_t in range(n_temp): f.write(analysis_data[i_t])           
        f.close()

    def get_act_energy(self,Tlim_low=50,Tlim_high=600,save_to_file=True):
        '''
        Fit log(J) vs. 1000/T with linear fit: y = a*x + b.
        Use specific temperature interval from tlim_low to tlim_high.
        The outputfile 'lin_fit_data.txt' contains the activation energy 
        (eV) for this temperature set (within the spec. limits).
        The plot 'current_lin_fit_dmr*.png' contains the log J vs. 1000/T
        with the curve from the linear regression.
        '''
        if not os.path.exists(self.dest_dir+'av_current.txt'):
            self.get_av_current() 
        else:
            if np.any(self.current) == None:
                current_data = np.loadtxt(self.dest_dir+'av_current.txt',comments='#',unpack=True)
                self.current, self.std_current = current_data[1], [current_data[2],current_data[3]]           
            if np.all(self.current == 0.0):
                print("Only zero values for current found, linear regretion on log(J) vs. 1/T not possible.")
                self.act_energy = 0.0
                a, b, r2 = 0.0, 0.0, 0.0
            else:
                # check if CurrTempSimulation already has current etc. attribute, if not check if data is already generated and if yes read it in                      
                x_conv = np.where(self.current[:,]>0.0)
                current      = self.current[x_conv[0][0]:x_conv[0][-1]+1]
                std_current  = [self.std_current[0][x_conv[0][0]:x_conv[0][-1]+1],self.std_current[1][x_conv[0][0]:x_conv[0][-1]+1]]
                temperatures = self.temp[x_conv[0][0]:x_conv[0][-1]+1]  
                # Only use non zero current for linear regression
                nonzero_ids = np.nonzero(current)
                current = current[nonzero_ids]
                std_current = [std_current[0][nonzero_ids],std_current[1][nonzero_ids]]
                temperatures = temperatures[nonzero_ids]
                # define x and y for linear regression
                x = 1000/temperatures                  
                y = np.log(current)
                if np.min(temperatures)>Tlim_high or np.max(temperatures)<Tlim_low:
                    print("Data points to analyse don't lie within the specified temperature range, linear regretion on log(J) vs. 1/T not possible.")
                    self.act_energy = 0.0
                    a, b, r2 = 0.0, 0.0, 0.0
                else:
                    idx_low  = np.argwhere( temperatures >= Tlim_low )[0,0]
                    idx_high = np.argwhere( temperatures <= Tlim_high )[-1,0]
                    x_fit = (x[idx_low:idx_high+1]).reshape((-1,1))
                    y_fit = y[idx_low:idx_high+1]
                    linear_model = LinearRegression().fit(x_fit,y_fit)
                    a, b, r2 = linear_model.coef_[0], linear_model.intercept_, linear_model.score(x_fit,y_fit)
                    print('Linear fit: slope     = ', a)
                    print('            intercept = ',b)
                    print('                  R^2 = ',r2)
                    self.act_energy = -1000*1000*a*Boltzmann/e
                    print('Activation energy: E_A = {} meV'.format(self.act_energy))
                    # Plot linear regression
                    plt.errorbar(x, y,yerr=std_current[1],linestyle = '',marker='o',color='r', label='Simulated points')
                    plt.plot(x,a*x + b, color='k',label='Linear fit: $y = a*x + b$\n$a = {0:.3}, b =  {1:.3}$\n$R^2 = {2:.3}$'.format(a,b,r2))
                    plt.legend()
                    plt.xlabel('1/T (1000/K)')
                    plt.ylabel('log J (A/m$^2$)')
                    if self.rates == "Marcus":
                        if self.analyse_Ecoul:
                            plt.title('Field = {} eV, disorder = {} eV, $\lambda$ = {} eV, Ecoul= {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.lam,self.Ecoul, self.DMR,self.sys_size),fontsize=10)
                        else:
                            plt.title('Field = {} eV, disorder = {} eV, $\lambda$ = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.lam,self.DMR,self.sys_size),fontsize=10)
                    else:
                        plt.title('Field = {} eV, disorder = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.DMR,self.sys_size),fontsize=10)
                    if not os.path.isdir(self.dest_dir+'act_energy'):
                        os.makedirs(self.dest_dir+'act_energy')
                    plt.savefig(self.dest_dir+'act_energy/current_lin_fit_dmr{}.png'.format(self.DMR))
                    plt.close()
                if save_to_file:
                    if not os.path.isdir(self.dest_dir+'act_energy'):
                        os.makedirs(self.dest_dir+'act_energy')
                    with open(self.dest_dir+'act_energy/lin_fit_data.txt','w') as f:
                        f.write('# Linear regression information\n')
                        f.write('# Slope (AK/1000/m^2) \n{}\n'.format(a))
                        f.write('# Intercept (A/m^2)\n{}\n'.format(b))
                        f.write('# Coefficient of determination(R^2)\n{}\n'.format(r2))
                        f.write('# Lower T limit (K)\n{}\n'.format(Tlim_low))
                        f.write('# Upper T limit (K)\n{}\n'.format(Tlim_high))                            
                        f.write('# Activation Energy (meV)\n{}'.format(self.act_energy))
                    f.close()
        

    def plot_av_current(self,Tlim_low=None,Tlim_high=None,Jlim_low=None,Jlim_high=None,plot_log=True,errorbar=False, only_conv=True,
                     curve_color='b', save_to_file=True
                     ):
        '''
        Plot current vs. 1000/T. 
        plot_log:     plot natural logarithm of current.
        errorbar:     display errorbars of each simulation points.
        only_conv:    plot only points for which at least one replica has 
                      converged.
        save_to_file: save the generated plot directly to 'av_current.png'.

        '''
        # check if CurrTempSimulation already has current etc. attribute, 
        # if not check if data is already generated and if yes read it in
        if not os.path.exists(self.dest_dir+'av_current.txt'):
            self.get_av_current()   
        else:
            if np.any(self.current) == None:
                current_data = np.loadtxt(self.dest_dir+'av_current.txt',comments='#',unpack=True)
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
            # adjust plot limits
            if Tlim_low!=None:
                x = np.where(temperatures[:]>=Tlim_low)
                if len(x[0])>0:
                    current      = current[x[0][0]:x[0][-1]+1]
                    std_current  = [std_current[0][x[0][0]:x[0][-1]+1],std_current[1][x[0][0]:x[0][-1]+1]]
                    temperatures = temperatures[x[0][0]:x[0][-1]+1]
            if Tlim_high!=None:
                x = np.where(temperatures[:]<=Tlim_high)
                if len(x[0])>0:
                    current      = current[x[0][0]:x[0][-1]+1]
                    std_current  = [std_current[0][x[0][0]:x[0][-1]+1],std_current[1][x[0][0]:x[0][-1]+1]]
                    temperatures = temperatures[x[0][0]:x[0][-1]+1]
            if Jlim_low!=None:
                x = np.where(current[:]>=Jlim_low)
                if len(x[0])>0:
                    current      = current[x[0][0]:x[0][-1]+1]
                    std_current  = [std_current[0][x[0][0]:x[0][-1]+1],std_current[1][x[0][0]:x[0][-1]+1]]
                    temperatures = temperatures[x[0][0]:x[0][-1]+1]
            if Jlim_high!=None:
                x = np.where(current[:]<=Jlim_high)
                if len(x[0])>0:
                    current      = current[x[0][0]:x[0][-1]+1]
                    std_current  = [std_current[0][x[0][0]:x[0][-1]+1],std_current[1][x[0][0]:x[0][-1]+1]]
                    temperatures = temperatures[x[0][0]:x[0][-1]+1]
                
            # make logarithmic plot with appropiate std error
            if plot_log:
                if errorbar:
                    plt.errorbar(1000/temperatures,np.log(current),yerr=std_current[1],label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                else:
                    plt.plot(1000/temperatures,np.log(current),'-+',label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                ylabel = 'log J (A/m$^2$)'
            else:
                if errorbar:
                    plt.errorbar(1000/temperatures,current,yerr=std_current[0],label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                else:
                    plt.plot(1000/temperatures,current,'-+',label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                ylabel = 'J (A/m$^2$)'
            if save_to_file:
                plt.xlabel('1/T (1000/K)')
                plt.ylabel(ylabel)
                if self.rates == "Marcus":
                    if self.analyse_Ecoul:
                        plt.title('Field = {} eV, disorder = {} eV, $\lambda$ = {} eV, Ecoul = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.lam,self.Ecoul,self.DMR,self.sys_size),fontsize=10)
                    else:
                        plt.title('Field = {} eV, disorder = {} eV, $\lambda$ = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.lam,self.DMR,self.sys_size),fontsize=10)
                else:
                    plt.title('Field = {} eV, disorder = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.DMR,self.sys_size),fontsize=10)
                plt.savefig(self.dest_dir+'av_current.png')
                plt.close()

    def plot_av_conductivity(self, plot_log=True, errorbar=False, only_conv=True, curve_color='b', save_to_file=True):
        '''
        Plot conductivity vs. 1000/T. 
        plot_log:     plot natural logarithm of current.
        errorbar:     display errorbars of each simulation points.
        only_conv:    plot only points for which at least one replica has 
                      converged.
        save_to_file: save the generated plot directly to 'av_current.png'.

        '''
        # check if CurrTempSimulation already has current etc. attribute, if not check if data is already generated and if yes read it in
        if not os.path.exists(self.dest_dir+'av_conductivity.txt'):
            print('Error: "{}" not found. Please generate data first.'.format(self.dest_dir+'av_conductivity.txt'))   
        else:
            if np.any(self.conduct) == None:
                conduct_data = np.loadtxt(self.dest_dir+'av_conductivity.txt',comments='#',unpack=True)
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
                    plt.errorbar(1000/temperatures,np.log(conduct),yerr=std_conduct[1],label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                else:
                    plt.plot(1000/temperatures,np.log(conduct),'-+',label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                ylabel = 'log $\sigma$ (A/eV/m)'
            else:
                if errorbar:
                    plt.errorbar(1000/temperatures,conduct,yerr=std_conduct[0],label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                else:
                    plt.plot(1000/temperatures,conduct,'-+',label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                ylabel = '$\sigma$ (A/eV/m)'
            if save_to_file:
                plt.xlabel('1/T (1000/K)')
                plt.ylabel(ylabel)
                if self.rates == "Marcus":
                    if self.analyse_Ecoul:
                        plt.title('Field = {} eV, disorder = {} eV, $\lambda$ = {} eV, Ecoul = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.lam,self.Ecoul,self.DMR,self.sys_size),fontsize=10)
                    else:
                        plt.title('Field = {} eV, disorder = {} eV, $\lambda$ = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.lam,self.DMR,self.sys_size),fontsize=10)
                else:
                    plt.title('Field = {} eV, disorder = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.DMR,self.sys_size),fontsize=10)
                plt.savefig(self.dest_dir+'av_conductivity.png')
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
        if not os.path.exists(self.dest_dir+'conv_analysis.txt'):
            self.get_act_energy()  
        else:
            n_temp      = len(self.temp)
            conv_data = np.loadtxt(self.dest_dir+'conv_analysis.txt', comments='#', unpack=True)
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
                        dmr_label = 'Dop.fract. {} %'.format(round(self.DMR*100,2))
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
                if self.rates == "Marcus":
                    if self.analyse_Ecoul:
                        fig.suptitle('Field={} V/nm, disorder={} eV, $\lambda$ = {} eV, Ecoul = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.lam,self.Ecoul,self.DMR,self.sys_size),size=10)
                    else:
                        fig.suptitle('Field={} V/nm, disorder={} eV, $\lambda$ = {} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.lam,self.DMR,self.sys_size),size=10)
                else:
                    fig.suptitle('Field={} V/nm, disorder={} eV, DMR = {}, system size = {} nm.'.format(self.field,self.dis,self.DMR,self.sys_size),size=10)
                plt.savefig(self.dest_dir+'conv_analysis.png')
                plt.close()
            else:
                if self.rates == "Marcus":
                    if self.analyse_Ecoul:
                        fig.suptitle('Field={} V/nm, disorder={} eV, $\lambda$ = {} eV, Ecoul = {} eV.'.format(self.field,self.dis,self.lam,self.Ecoul),size=10)
                    else:
                        fig.suptitle('Field = {} V/nm, disorder = {} eV, $\lambda$ = {} eV.'.format(self.field,self.dis,self.lam),size=14)
                else:
                    fig.suptitle('Field = {} V/nm, disorder = {} eV.'.format(self.field,self.dis),size=14)



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
    # TODO: if needed add "Marcus" rates mode with lambda
    # TODO: as in CurrTempSimulation, collect current data seperately
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

    def get_av_current(self):
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
                if os.path.exists(outputdir+'/output_job_0'):
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
        if not os.path.isdir(self.dest_dir+'data'):
            os.makedirs(self.dest_dir+'data')
        with open(self.dest_dir+'av_current.txt', 'w') as f:
                f.write('# Field = {} V/nm, disorder = {} eV, DMR = {}, system size = {}\n'.format(self.field,self.dis,self.DMR,self.sys_size))
                f.write('# Field(V/nm)   Av. Current density(A/m$^2$)   Normal std. error(A/m$^2$)   Log. std. error(A/m$^2$)\n')
                for i_f in range(n_field):
                    f.write('{0:<16}   {1:<26}   {2:<24}   {3:<30}\n'.format(self.field[i_f], self.current[i_f], self.std_current[0][i_f], self.std_current[1][i_f]))
        f.close()

    def get_av_conductivity(self):
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
                    if os.path.exists(outputdir+'/output_job_0'):
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
        if not os.path.isdir(self.dest_dir+'data'):
            os.makedirs(self.dest_dir+'data')
        with open(self.dest_dir+'conv_analysis.txt', 'w') as f:
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

    def plot_av_current(self,plot_log=True,errorbar=False, only_conv=True,curve_color='g', save_to_file=True):
        # check if CurrTempSimulation already has current etc. attribute, if not check if data is already generated and if yes read it in
        if not os.path.exists(self.dest_dir+'av_current.txt'):
            self.get_av_current()   
        else:
            if np.any(self.current) == None:
                current_data = np.loadtxt(self.dest_dir+'av_current.txt',comments='#',unpack=True)
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
                    plt.errorbar(fields,np.log(current),yerr=std_current[1],label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                else:
                    plt.plot(fields,np.log(current),'-+',label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                ylabel = 'log J (A/m$^2$)'
            else:
                if errorbar:
                    plt.errorbar(fields,current,yerr=std_current[0],label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                else:
                    plt.plot(fields,current,'-+',label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                ylabel = 'J (A/m$^2$)'
            if save_to_file:
                plt.xlabel('E (V/nm)')
                plt.ylabel(ylabel)
                plt.title('Sim. for temperature = {} K, disorder = {} eV, DMR = {}, system size = {} nm.'.format(self.temp,self.dis,self.DMR,self.sys_size),fontsize=10)
                plt.savefig(self.dest_dir+'av_current.png')
                plt.close()

    def plot_av_conductivity(self):
        # TODO
        pass
    
    def plot_conv_analysis(self,x_loc_scale=1.25,bar_loc=0.0,bar_width=0.25, 
                           bar_color='g',save_to_file=True,figure=None,axes=None):
        
        if not os.path.exists(self.dest_dir+'conv_analysis.txt'):
            print('Error: "{}" not found. Please generate data first.'.format(self.dest_dir+'conv_analysis.txt'))   
        else:
            n_field      = len(self.field)
            conv_data = np.loadtxt(self.dest_dir+'conv_analysis.txt', comments='#', unpack=True)

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
                        dmr_label = 'Dop.fract. {} %'.format(round(self.DMR*100,2))
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
                plt.savefig(self.dest_dir+'conv_analysis.png')
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
    # TODO: if needed add "Marcus" rates mode with lambda
    def __init__(self,rates="Miller",source_dir='',dest_dir='analysis',sim_data_file='sim_param.txt'):
        sim_data = np.genfromtxt(sim_data_file,dtype=float)

        self.rates    = rates
        self.n_r      = int(sim_data[0])
        self.field    = sim_data[1]
        self.temp     = sim_data[2]
        self.dis      = sim_data[3]
        if self.rates=="Marcus":
            self.lam = sim_data[4]
            self.DMR = sim_data[5]
            sys_size = sim_data[6:] 
        else:
            self.DMR      = sim_data[4]
            sys_size = sim_data[5:]   
        # count only sys_sizes not n_cpus
        self.sys_size = []
        for i in sys_size:
            if i > 10.0:
                    (self.sys_size).append(i)
        print(self.sys_size)

        self.source_dir = source_dir
        self.dest_dir   = dest_dir+"/"

        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        shutil.copyfile(sim_data_file,self.dest_dir+sim_data_file)

        self.current  = None
        self.std_current     = None

    def collect_current_data(self):
        '''
        Extract final current density from outputfile 'output_job_0' from 
        converged jobs for every replica. If no final current is found, 
        current is set to 0.
        The outputfile 'J_temp_*.txt' contains (in columns): 
        - replica number
        - final current density values (A/m^2)
        
        '''
        n_sys      = len(self.sys_size)
        current_data = np.zeros((n_sys,self.n_r))
        if not os.path.isdir(self.dest_dir+'current_data'):
            os.makedirs(self.dest_dir+'current_data')
        for i_s in range(n_sys):
            for i_r in range(self.n_r):
                outputdir=self.source_dir+'sys_{}/r_{}'.format(i_s, i_r)
                if os.path.exists(outputdir+'/output_job_0'):
                    if 'final current density: ' in open(outputdir+'/output_job_0').read():
                        current_dens=extract.extract_e(outputdir+'/output_job_0', "final current density: ")
                        if current_dens<0:
                            print("Negative current density encountered at {}: {} A/m^2.\n Not counted to current densities.".format(outputdir,current_dens))
                        else:
                            current_data[i_s,i_r] = current_dens
            with open(self.dest_dir+'current_data/curr_sys_{}.txt'.format(i_s), 'w') as f: 
                    f.write("# Raw current data for sys_{} (T={} K)\n".format(i_s,self.sys_size[i_s]))
                    if self.rates == "Marcus":
                        f.write('# Temp. = {} K, Field = {} eV, disorder = {} eV, lambda = {} eV, DMR = {}\n'.format(self.temp, self.field,self.dis,self.lam, self.DMR))
                    else:
                        f.write('# Temp. = {} K, Field = {} eV, disorder = {} eV, DMR = {}\n'.format(self.temp, self.field,self.dis,self.DMR))
                    f.write('#\n# Replica Nr.     Current density(A/m$^2$)\n')
                    for i_r in range(self.n_r):
                        f.write('  {0:<13}   {1:<1}\n'.format(i_r, current_data[i_s,i_r]))
            f.close()

    def get_av_current(self):
        '''
        Average final current density over all replicas for each system size.
        Get data from '/current_data'.
        The outputfile 'current.txt' contains (in columns): 
        - system sizes (nm)
        - average current density values (A/m^2)
        - standard error values of current density (A/m^2)
        - logarithmic std. error of current density (A/m^2) 
          (to use for errorbars in logarithmic plots)
        
        '''
        
        n_sys      = len(self.sys_size)
        av_current  = np.zeros(n_sys)
        std_current = np.zeros(n_sys)
        log_std_current = np.zeros(n_sys)
        for i_s in range(n_sys):
            if not os.path.exists(self.dest_dir+'current_data/curr_sys_{}.txt'.format(i_s)):
                self.collect_current_data()
            current = (np.loadtxt(self.dest_dir+'current_data/curr_sys_{}.txt'.format(i_s)))[:,1]
            current = current[np.nonzero(current)]
            # Calculate average current
            if len(current)>0:
                av_current[i_s]  = stat.gmean(current)
                std_current[i_s] = np.std(current)/np.sqrt(np.count_nonzero(current))
                log_std_current[i_s] = np.std(np.log(current))/np.sqrt(np.count_nonzero(current))
        
        self.current = av_current
        self.std_current = [std_current,log_std_current]
        # Write current, temperature, DMR to txt file
        with open(self.dest_dir+'av_current.txt', 'w') as f:
            f.write("# Average current and related values.\n")
            if self.rates == "Marcus":
                f.write('# Temp. = {} K, Field = {} eV, disorder = {} eV, lambda = {} eV, DMR = {}\n'.format(self.temp, self.field,self.dis,self.lam, self.DMR))
            else:
                f.write('# Temp. = {} K, Field = {} eV, disorder = {} eV, DMR = {}\n'.format(self.temp,self.field,self.dis,self.DMR))
            f.write('# System size(nm)   Av. Current density(A/m$^2$)   Normal std. error(A/m$^2$)   Log. std. error(A/m$^2$)\n')
            for i_s in range(n_sys):
                f.write('  {0:<14}   {1:<28}   {2:<26}   {3:<1}\n'.format(self.sys_size[i_s], self.current[i_s], self.std_current[0][i_s], self.std_current[1][i_s]))
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
                    if os.path.exists(outputdir+'/output_job_0'):
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
                analysis_data.append(('  {0:<14}   {1:<4}   {2:<7}   {3:<6}   {4:<11}   {5:<12}   {6:<19}   {7:<20}   {8:<20}   {9:<20}\n'.format( self.sys_size[i_s],
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
        with open(self.dest_dir+'conv_analysis.txt', 'w') as f:
            f.write('# Total number of jobs: {}\n'.format(n_sys*self.n_r))
            if self.rates == "Marcus":
                f.write('# Temp.: {} K, Field: {} eV,  Disorder: {} eV, Lambda = {} eV, Dop.fract.: {} %\n'.format(self.temp,
                                                                                                                   self.field,
                                                                                                                   self.dis,
                                                                                                                   self.lam, 
                                                                                                                   100.0*self.DMR))
            else:
                f.write('# Temp.: {} K, Field: {} eV,  Disorder: {} eV, Dop.fract.: {} %\n'.format(self.temp,
                                                                                               self.field,
                                                                                               self.dis, 
                                                                                               100.0*self.DMR))

            f.write('# Number of replicas for each settings configuration: {} \n'.format(self.n_r))             
            f.write('# Overview:   {} % of jobs started to run ({}).\n'.format(100.0*run/n_sys/self.n_r,run))
            if run>0:
                f.write('#             {} % of run jobs have converged ({}). \n \n'.format(round(100.0*conv/run,1),conv))                
            f.write('#\n# System size (nm)   run    non-run   conv.    non-conv.     Conv.(%)  Conv. time: Average(h)     Maximum(h)     Minimum(h)     Standarddev.(h)\n')
            for i_s in range(n_sys): f.write(analysis_data[i_s])           
        f.close()

    def plot_av_current(self,y_logmin=-15,y_logmax=15,plot_log=True,errorbar=False, only_conv=True,curve_color='b', save_to_file=True):
        # check if CurrTempSimulation already has current etc. attribute, if not check if data is already generated and if yes read it in
        if not os.path.exists(self.dest_dir+'av_current.txt'):
            self.get_av_current()   
        else:
            if np.any(self.current) == None:
                current_data = np.loadtxt(self.dest_dir+'av_current.txt',comments='#',unpack=True)
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
                    plt.errorbar(sys_size,np.log(current),yerr=std_current[1],label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                else:
                    plt.plot(sys_size,np.log(current),'-+',label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                ylabel = 'log J (A/m$^2$)'
            else:
                if errorbar:
                    plt.errorbar(sys_size,current,yerr=std_current[0],label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                else:
                    plt.plot(sys_size,current,'-+',label='Dop.fract. {} %'.format(round(self.DMR*100,2)),color=curve_color)
                ylabel = 'J (A/m$^2$)'
            if save_to_file:
                plt.xlabel('system size (nm)')
                if plot_log:
                    plt.ylim(y_logmin,y_logmax)
                plt.ylabel(ylabel)
                if self.rates == "Marcus":
                    plt.title('Temp. = {} K, Field = {} eV, disorder = {} eV, $\lambda$ = {} eV, DMR = {}.'.format(self.temp,self.field,self.dis,self.lam,self.DMR),fontsize=10)
                else:
                    plt.title('Temp. = {} K, Field = {} eV, disorder = {} eV, DMR = {}.'.format(self.temp,self.field,self.dis,self.DMR),fontsize=10)
                plt.savefig(self.dest_dir+'av_current.png')
                plt.close()

    def plot_conv_analysis(self,x_loc_scale=1.25,bar_loc=0.0,bar_width=0.25, 
                           bar_color='b',dest_dir='',save_to_file=True,figure=None,axes=None):
        n_sys      = len(self.sys_size)
        if not os.path.exists(self.dest_dir+'conv_analysis.txt'):
            self.get_av_current()   
        else:
            conv_data = np.loadtxt(self.dest_dir+'conv_analysis.txt', comments='#', unpack=True)

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
                    dmr_label = 'Dop.fract. {} %'.format(round(self.DMR*100,2))
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
            if self.rates == "Marcus":
                fig.suptitle('Temp. = {} K, Field={} V/nm,  disorder = {} eV, $\lambda$ = {} eV, DMR = {}.'.format(self.temp,self.field,self.dis,self.lam,self.DMR),size=10)
            else:
                fig.suptitle('Temp. = {} K, Field={} V/nm,  disorder = {} eV, DMR = {}.'.format(self.temp,self.field,self.dis,self.DMR),size=10)
            plt.savefig(self.dest_dir+'conv_analysis.png')
            plt.close()
        else:
            fig.suptitle('Temp. = {} K, Field = {} V/nm, disorder = {} eV.'.format(self.temp,self.field,self.dis),size=14)

class CurrTempDMRset:
    '''
    Analyse J(T)-sets for different doping molar ratios (DMR)
    and one value for field and disorder.
    The directory must contain subdirectories in the format DMR_*/temp_*/r_*
    and the J(T)-sets have to consist of the same sim. parameters.
    '''
    def __init__(self,rates="Miller",analyse_Ecoul=False,source_dir='',dest_dir='analysis'):
        '''
        Initialize attributes of DMRset: number of replicas, field, 
        disorder, doping molar ratios (DMRs), system sizes, temperatures.
        Initialize colorset for plots.
        '''
        self.source_dir = source_dir
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)        
        self.dest_dir   = dest_dir+'/'
        # get sorted list of DMR_* directories
        DMR_dir_unsorted = glob.glob(source_dir+'DMR_*')
        self.list_DMR_dirs = sorted(DMR_dir_unsorted, key = lambda x: int(x.split('_')[-1]))
        if self.list_DMR_dirs == []:
            print("SimAnalysisError: No compatible CurrTempSimulations found.")
            exit()
        else:
            print(self.list_DMR_dirs)
            self.rates = rates
            self.analyse_Ecoul = analyse_Ecoul
            for i_dmr,dir in enumerate(self.list_DMR_dirs):
                dmr_sim = CurrTempSimulation(rates = self.rates, analyse_Ecoul=self.analyse_Ecoul, source_dir=dir+'/',dest_dir='analysis/'+dir)
                if i_dmr==0:
                    self.n_r      = dmr_sim.n_r
                    self.field    = dmr_sim.field
                    self.dis      = dmr_sim.dis
                    if self.rates == "Marcus":
                        self.lam = dmr_sim.lam
                        if self.analyse_Ecoul:
                            self.Ecoul = dmr_sim.Ecoul
                    self.DMR      = [dmr_sim.DMR]
                    self.sys_size = [dmr_sim.sys_size]
                    self.temp     = dmr_sim.temp
                else:
                    if dmr_sim.n_r != self.n_r  or dmr_sim.field != self.field or dmr_sim.dis != self.dis or np.any(dmr_sim.temp != self.temp):
                        print("SimAnalysisError: default sim. parameters of CurrTemp simulations do not agree, DMR set cannot be analysed.")
                        exit()
                    else:
                        self.DMR.append(dmr_sim.DMR)
                        self.sys_size.append(dmr_sim.sys_size)

        self.act_energy = None 
        self.colorset = np.array(['salmon','orangered','darkred','violet','blueviolet','darkblue','royalblue','dodgerblue','mediumseagreen',
                                  'green','black','navajowhite','blue','brown','burlywood','cadetblue','chartreuse','chocolate', 
                                  'coral','cornflowerblue','crimson','cyan','darkblue','darkcyan','darkgoldenrod','darkgray', 
                                  'darkgreen','darkgrey','darkkhaki','darkmagenta','darkolivegreen','darkorange','darkorchid', 
                                  'darksalmon','darkseagreen','darkslateblue','darkslategray','darkturquoise','darkviolet', 
                                  'deeppink','deepskyblue','dimgray','firebrick','forestgreen','fuchsia','gainsboro','gold', 
                                  'goldenrod','greenyellow','grey','hotpink','indianred','indigo','khaki','lawngreen', 
                                  'lemonchiffon','lightblue','lightcoral','lightgreen','lightpink','lightsalmon', 
                                  'lightseagreen','lightskyblue','lightslategrey','lightsteelblue','lime','limegreen','magenta', 
                                  'maroon','mediumaquamarine','mediumblue','mediumorchid','mediumpurple','mediumslateblue', 
                                  'mediumspringgreen','mediumturquoise','mediumvioletred','midnightblue','mintcream', 
                                  'moccasin','navy','olive','olivedrab','orange','orchid','palegoldenrod','palegreen','paleturquoise', 
                                  'palevioletred','peachpuff','peru','pink','plum','powderblue','purple','rebeccapurple', 
                                  'red','rosybrown','royalblue','saddlebrown','sandybrown','seagreen','sienna','silver', 
                                  'skyblue','slateblue','slategray','slategrey','snow','springgreen','steelblue','tan','teal', 
                                  'thistle','tomato','turquoise','wheat','yellow','yellowgreen'])

        sim_data = open(self.dest_dir+'sim_param_DMR_set.txt','w')
        sim_data.write("# n_replicas\n{}\n".format(self.n_r))
        sim_data.write("# field\n{}\n".format(self.field))
        sim_data.write("# disorder\n{}\n".format(self.dis))
        if self.rates == "Marcus":
            sim_data.write("# lambda\n{}\n".format(self.lam))
            if self.analyse_Ecoul:
                sim_data.write("# Ecoul\n{}\n".format(self.Ecoul))
        sim_data.write("# doping ratios\n")
        for dmr in self.DMR: sim_data.write("{}\n".format(dmr))
        sim_data.write("# system size\n")
        for sys in self.sys_size: sim_data.write("{}\n".format(sys))
        sim_data.write("# temperatures\n")
        for t in self.temp: sim_data.write("{}\n".format(t))



    def plot_av_current(self,Tlim_low=None,Tlim_high=None,Jlim_low=None,Jlim_high=None,errorbar=False,plot_log=True, only_conv=True):
        DMR_dir = self.list_DMR_dirs
        curr_dmr_colors  = self.colorset
        for i_dmr,dir in enumerate(DMR_dir):
            dmr_sim = CurrTempSimulation(rates = self.rates, analyse_Ecoul= self.analyse_Ecoul, source_dir=dir+'/',dest_dir=self.dest_dir+dir)
            if not os.path.exists(self.dest_dir+dir+"/av_current.txt"):    
                dmr_sim.get_av_current()
            dmr_sim.plot_av_current(Tlim_low,Tlim_high,Jlim_low,Jlim_high, plot_log, errorbar, only_conv, curr_dmr_colors[i_dmr], False)
        if plot_log:
            ylabel = 'log J (A/m$^2$)'
        else:
            ylabel = 'J (A/m$^2$)'
        plt.xlabel('1/T (1000/K)')
        if self.rates == "Marcus":
            if self.analyse_Ecoul:  
                plt.title('Field = {} V/nm, disorder = {} eV, $\lambda$ = {} eV, Ecoul = {} eV.'.format(self.field,self.dis,self.lam, self.Ecoul),fontsize=14)   
            else:
                plt.title('Field = {} V/nm, disorder = {} eV, $\lambda$ = {} eV.'.format(self.field,self.dis,self.lam),fontsize=14)      
        else:
            plt.title('Field = {} V/nm, disorder = {} eV.'.format(self.field,self.dis),fontsize=14)
        plt.ylabel(ylabel)
        plt.legend()    
        plt.savefig(self.dest_dir+'current_DMR_set.png')
        plt.close()

    def plot_av_conductivity(self,errorbar=False,plot_log=True, only_conv=True,log10=False):
        DMR_dir = self.list_DMR_dirs
        curr_dmr_colors  = self.colorset
        for i_dmr,dir in enumerate(DMR_dir):
            dmr_sim = CurrTempSimulation(rates = self.rates, analyse_Ecoul= self.analyse_Ecoul, source_dir=dir+'/',dest_dir=self.dest_dir+dir)
            if not os.path.exists(self.dest_dir+dir+"/conductivity.txt"):    
                dmr_sim.get_av_conductivity()
            dmr_sim.plot_av_conductivity( plot_log, errorbar, only_conv, curr_dmr_colors[i_dmr], False)
        if plot_log:
            ylabel = 'log $\sigma$ (A/eV/m)'
        else:
            ylabel = '$\sigma$ (A/eV/m)'
        if log10:
            plt.yscale("log")
        plt.xlabel('1/T (1000/K)')
        plt.ylabel(ylabel)
        plt.legend()
        if self.rates == "Marcus":
            if self.analyse_Ecoul:  
                plt.title('Field = {} V/nm, disorder = {} eV, $\lambda$ = {} eV, Ecoul = {} eV.'.format(self.field,self.dis,self.lam, self.Ecoul),fontsize=14)   
            else:
                plt.title('Field = {} V/nm, disorder = {} eV, $\lambda$ = {} eV.'.format(self.field,self.dis,self.lam),fontsize=14)      
        else:
            plt.title('Field = {} V/nm, disorder = {} eV.'.format(dmr_sim.field,dmr_sim.dis),fontsize=14)
        plt.savefig(self.dest_dir+'conductivity_DMR_set.png')
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
            dmr_sim = CurrTempSimulation(rates = self.rates, analyse_Ecoul= self.analyse_Ecoul, source_dir=dir+'/',dest_dir=self.dest_dir+dir)
            if not os.path.exists(self.dest_dir+dir+"/conv_analysis.txt"):    
                dmr_sim.get_conv_analysis()
            if i_dmr == 0:
                n_temp = len(self.temp)
                fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(2.5*n_temp*8/10, 10))
            dmr_sim.plot_conv_analysis(x_loc_scale=x_loc_scale,bar_loc=bar_loc[i_dmr],bar_width=width,
                                    bar_color=curr_dmr_colors[i_dmr],
                                    save_to_file=False,figure=fig,axes=axes)

        plt.savefig(self.dest_dir+'conv_analysis_DMR_set.png')
        plt.close()  

    def get_act_energy(self,Tlim_low=250,Tlim_high=350):
        DMR_dir = self.list_DMR_dirs
        n_DMR    = len(DMR_dir)  
        act_energy = np.zeros(n_DMR)
        dmr        = self.DMR * 100
        if not os.path.isdir(self.dest_dir+'act_energy'):
            os.makedirs(self.dest_dir+'act_energy')       
        for i_dmr, dir in enumerate(DMR_dir):
            if os.path.exists(self.dest_dir+dir+"/act_energy/lin_fit_data.txt"):
                old_Tlim_low = (np.loadtxt(self.dest_dir+dir+"/act_energy/lin_fit_data.txt",unpack=True))[3]
                old_Tlim_high = (np.loadtxt(self.dest_dir+dir+"/act_energy/lin_fit_data.txt",unpack=True))[4]
                # check if linear fit was made of the right T range
                if old_Tlim_low != Tlim_low or old_Tlim_high != Tlim_high:
                    dmr_sim = CurrTempSimulation(rates = self.rates, analyse_Ecoul= self.analyse_Ecoul, source_dir=dir+'/',dest_dir=self.dest_dir+dir)
                    dmr_sim.get_act_energy(Tlim_low=Tlim_low,Tlim_high=Tlim_high)      
            else:
                dmr_sim = CurrTempSimulation(rates = self.rates, analyse_Ecoul= self.analyse_Ecoul, source_dir=dir+'/',dest_dir=self.dest_dir+dir)
                dmr_sim.get_act_energy(Tlim_low=Tlim_low,Tlim_high=Tlim_high)
                if not os.path.exists(self.dest_dir+dir+"/act_energy/lin_fit_data.txt"):
                    print("Activation energy could not be determined for DMR_{}.".format(i_dmr))
                    continue              
            act_energy[i_dmr] = (np.loadtxt(self.dest_dir+dir+"/act_energy/lin_fit_data.txt",unpack=True))[5]          
        if np.count_nonzero(act_energy)>0.0:
            with open(self.dest_dir+"act_energy/act_energy_DMR_set.txt",'w') as f:
                if self.rates == "Marcus":   
                    if self.analyse_Ecoul:
                        f.write('# Field = {} V/nm, disorder = {} eV, lambda = {} eV, Ecoul = {} eV\n'.format(self.field,self.dis,self.lam,self.Ecoul))
                    else:
                        f.write('# Field = {} V/nm, disorder = {} eV, lambda = {} eV\n'.format(self.field,self.dis,self.lam))
                else:    
                    f.write('# Field = {} V/nm, disorder = {} eV\n'.format(self.field,self.dis))
                f.write("# T. range for lin. fit of log J(1000/T): {} - {} K\n".format(Tlim_low,Tlim_high))    
                f.write('#\n# DMR      Activation energy(meV)\n')
                for i in range(n_DMR):
                    if act_energy[i]==0.0:
                        continue
                    else:
                        f.write('{0:<8}   {1:<6}\n'.format(dmr[i],act_energy[i]))
            f.close()
            self.act_energy = act_energy
        else:
            print("Activation energy could not be determined for this DMRset.")

    def plot_act_energy(self):
        if not os.path.exists(self.dest_dir+"act_energy/act_energy_DMR_set.txt"):
            self.get_act_energy()
        if os.path.exists(self.dest_dir+"act_energy/act_energy_DMR_set.txt"):   
            act_energy = (np.loadtxt(self.dest_dir+"act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[1]
            dmr        = (np.loadtxt(self.dest_dir+"act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[0]
            plt.plot(dmr*100,act_energy,marker="+", linestyle="None")
            plt.xlabel('DMR (%)')
            plt.xscale("log")
            plt.ylabel('$E_A$ (meV)')
            if self.rates == "Marcus":
                if self.analyse_Ecoul:  
                    plt.title('Field = {} V/nm, disorder = {} eV, $\lambda$ = {} eV, Ecoul = {} eV.'.format(self.field,self.dis,self.lam, self.Ecoul),fontsize=14)                    
                else:
                    plt.title('Field = {} V/nm, disorder = {} eV, $\lambda$ = {} eV'.format(self.field,self.dis, self.lam),fontsize=14)
            else:
                plt.title('Field = {} V/nm, disorder = {} eV'.format(self.field,self.dis),fontsize=14)
            plt.savefig(self.dest_dir+'act_energy/act_energy_DMR_set.png')
            plt.close()        
        else:
            print("Activation energy vs. DMR cannot be plotted.")

    def get_nonconv_joblist(self):
        '''
        If simulation has not converged completely write a new joblist 
        with all temperatures that have not converged in order to restart 
        the simulation with this new joblist file.
        '''
        DMR_dir = self.list_DMR_dirs 
        n_temp  = len(self.temp) 
        if not os.path.isdir(self.dest_dir+'conv_analysis.txt'):
            self.plot_conv_analysis()      
        for i_dmr, dir in enumerate(DMR_dir):
            converged_jobs = np.loadtxt(self.dest_dir+dir+"/conv_analysis.txt",unpack=True)[3]
            nonconv_jobs = []
            for i_t in range(n_temp):
                if converged_jobs[i_t] != 50:
                    n_cpus = int(np.loadtxt(self.source_dir+"joblist_DMR_"+str(i_dmr),usecols=0)[0])
                    run_kmc_name = np.loadtxt(self.source_dir+"joblist_DMR_"+str(i_dmr),usecols=1,dtype=str)[0]
                    for i_r in range(self.n_r):
                        nonconv_jobs.append("{} {} {} {} {}\n".format(n_cpus, run_kmc_name, i_dmr, i_t, i_r))
                    print("DMR_{}/temp_{} not fully converged. Moving files to non_conv".format(i_dmr,i_t)) 
                    if not os.path.exists(self.source_dir+"non_conv"):
                        os.makedirs(self.source_dir+"non_conv")   
                    try:
                        shutil.move(self.source_dir+"DMR_{}/temp_{}".format(i_dmr,i_t),self.source_dir+"non_conv/DMR_{}/temp_{}".format(i_dmr,i_t))
                    except:
                        print("Already moved")
                    
            if len(nonconv_jobs) != 0:
                print("new_joblist_DMR_{} and new .sh file are being generated.".format(i_dmr))

                new_joblist = open(self.source_dir+"new_joblist_DMR_"+str(i_dmr),'w')
                for line in nonconv_jobs: new_joblist.write(line)
                new_joblist.close()

                old_sh_name = glob.glob(self.source_dir+"submit_EA_*DMR_{}.sh".format(i_dmr))[0]
                old_sh_file = open(old_sh_name,'r')
                sh_data = old_sh_file.readlines()
                old_sh_file.close()

                sh_data[-1] = sh_data[-1].replace("joblist","new_joblist")
                new_sh_name = old_sh_name.replace("submit","new_submit")
                new_sh_file = open(new_sh_name,"w")
                for line in sh_data:
                    new_sh_file.write(line)
                new_sh_file.close()

            



class CurrFieldDMRset:
    '''
    Analyse J(E)-sets for different doping molar ratios (DMR)
    and one value for field and disorder.
    The directory must contain subdirectories in the format DMR_*/temp_*/r_*
    and the J(E)-sets have to consist of the same sim. parameters.
    '''
    # TODO: if needed add "Marcus" rates mode with lambda
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

        sim_data = open(self.dest_dir+'sim_param_DMR_set.txt','w')
        sim_data.write("# n_replicas\n{}\n".format(self.n_r))
        sim_data.write("# temperature\n{}\n".format(self.temp))
        sim_data.write("# disorder\n{}\n".format(self.dis))
        sim_data.write("# doping ratios\n")
        for dmr in self.DMR: sim_data.write("{}\n".format(dmr))
        sim_data.write("# system size\n")
        for sys in self.sys_size: sim_data.write("{}\n".format(sys))
        sim_data.write("# fields\n")
        for f in self.field: sim_data.write("{}\n".format(f))



    def plot_av_current(self,errorbar=False,plot_log=True, only_conv=True):
        DMR_dir = self.list_DMR_dirs
        curr_dmr_colors  = self.colorset
        for i_dmr,dir in enumerate(DMR_dir):
            dmr_sim = CurrFieldSimulation(source_dir=dir+'/',dest_dir='analysis/'+dir)
            if not os.path.exists(self.dest_dir+dir+"/current.txt"):    
                dmr_sim.get_av_current()
            dmr_sim.plot_av_current( plot_log, errorbar, only_conv, curr_dmr_colors[i_dmr], False)
        if plot_log:
            ylabel = 'log J (A/m$^2$)'
        else:
            ylabel = 'J (A/m$^2$)'
        plt.xlabel('E (V/nm)')
        plt.title('Simulation for temperature = {} K, disorder = {} eV.'.format(self.temp,self.dis),fontsize=14)
        plt.ylabel(ylabel)
        plt.legend()    
        plt.savefig(self.dest_dir+'current_DMR_set.png')
        plt.close()

    def plot_av_conductivity(self,errorbar=False,plot_log=True, only_conv=True,log10=False):
        # TODO: CurrFieldSimulation.get_av_conductivity() (if necessary)
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
            if not os.path.exists(self.dest_dir+dir+"/conv_analysis.txt"):    
                dmr_sim.get_conv_analysis()
            if i_dmr == 0:
                n_field = len(self.field)
                fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(2.5*n_field*8/10, 10))
            dmr_sim.plot_conv_analysis(x_loc_scale=x_loc_scale,bar_loc=bar_loc[i_dmr],bar_width=width,
                                    bar_color=curr_dmr_colors[i_dmr],
                                    save_to_file=False,figure=fig,axes=axes)

        plt.savefig(self.dest_dir+'conv_analysis_DMR_set.png')
        plt.close()  