"""
1. extract IP from QP folders of the form "10_A_shell"
"""

from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
import re
import yaml
import scipy.constants as const
import time



def main():

    #  I/O

    my_mol_num = 0


    my_pattern_short = "_A_shell"
    my_out_file_name = "Analysis/EAIP/eaip_7b043a1428ff586ba119947b48219969_IP_summary.yml"
    mol_idxs = [775, 875]


    #  prepare
    folders = return_target_folders(my_pattern_short)
    my_radii = np.sort(get_radii(folders))

    folders_sorted = []
    for i, my_radius in enumerate(my_radii):
        folders_sorted.append((str(my_radius)+my_pattern_short))

    #print(folders_sorted)  # sorted names

    # check if out_file_name exists
    updated_folders = return_folders_with_output(folders_sorted, my_out_file_name)
    print("updated folders:", updated_folders)
    my_radii_updated = np.sort(get_radii(updated_folders))
    print("updated radii:", my_radii_updated)


############ TEST
    fid = open("all_energies.txt", "w")
    n_r = len(my_radii_updated)
    n_r = 4

    ip_sd_stepwise = np.zeros([7, n_r])
    ip_fe_stepwise = np.zeros([7, n_r])
    ip_fe_final = np.zeros(n_r) # for mol_idx = 775
    ip_sd_final = np.zeros(n_r)

    ############################################### 0 ####################
    for i in range(n_r):
        fid.write("R = {} A\n".format(my_radii_updated))
        for step_num in range(7):
            my_time = time.time()
            my_ip = IP(updated_folders[i], my_radii_updated[i])
            my_ip.extract_IP()
            my_ip.extract_IP_stepwise(num_step=step_num, mol_idx = mol_idxs[my_mol_num])
            #print("Energy Charged = {} eV".format(my_ip.ip_charged))
            #print("Energy Uncharged = {} eV".format(my_ip.ip_uncharged))
            #print("Energy Charged = {} eV".format(my_ip.ip_charged))
            print("step = {}".format(step_num))
            print("IP s-d= {} eV".format(my_ip.ip_sd_stepwise))
            print("IP full_env = {} eV".format(my_ip.ip_fe_stepwise))

            ip_sd_stepwise[step_num, i] = my_ip.ip_sd_stepwise
            ip_fe_stepwise[step_num, i] = my_ip.ip_fe_stepwise

            fid.write("Charged = {} eV\tUncharged = {} eV\tIP = {} eV \n".format(my_ip.ip_charged, my_ip.ip_charged, my_ip.ip_sd_stepwise))
        ip_sd_final[i] = my_ip.raw_data[my_mol_num].single_delta
        ip_fe_final[i] = my_ip.raw_data[my_mol_num].full_env
        print("R = {} A.\tTime: {} sec".format(my_radii_updated[i], time.time() - my_time))
    ############################################### 0 ######################

    # plt.plot(range(7), ip_sd_stepwise, marker='o', label='single-delta')
    # plt.plot(6*np.ones(n_r), -ip_sd_final, marker='*', label='true sd', markersize=12)
    plt.plot(range(7), ip_fe_stepwise, marker='x', label='full env + V_C')
    plt.plot(6*np.ones(n_r), -ip_fe_final, marker='*', label='true fe')


    ############################################### 1 ####################
    my_mol_num = 1
    for i in range(n_r):
        fid.write("R = {} A\n".format(my_radii_updated))
        for step_num in range(7):
            my_time = time.time()
            my_ip = IP(updated_folders[i], my_radii_updated[i])
            my_ip.extract_IP()
            my_ip.extract_IP_stepwise(num_step=step_num, mol_idx = mol_idxs[my_mol_num])
            #print("Energy Charged = {} eV".format(my_ip.ip_charged))
            #print("Energy Uncharged = {} eV".format(my_ip.ip_uncharged))
            #print("Energy Charged = {} eV".format(my_ip.ip_charged))
            print("step = {}".format(step_num))
            print("IP s-d= {} eV".format(my_ip.ip_sd_stepwise))
            print("IP full_env = {} eV".format(my_ip.ip_fe_stepwise))

            ip_sd_stepwise[step_num, i] = my_ip.ip_sd_stepwise
            ip_fe_stepwise[step_num, i] = my_ip.ip_fe_stepwise

            fid.write("Charged = {} eV\tUncharged = {} eV\tIP = {} eV \n".format(my_ip.ip_charged, my_ip.ip_charged, my_ip.ip_sd_stepwise))
        ip_sd_final[i] = my_ip.raw_data[my_mol_num].single_delta
        ip_fe_final[i] = my_ip.raw_data[my_mol_num].full_env
        print("R = {} A.\tTime: {} sec".format(my_radii_updated[i], time.time() - my_time))
    ############################################### 1 ######################



    print("nothing1")
    fid.close()
    print(ip_sd_stepwise)

    # IP stepwise sd
    # plt.plot(range(7), ip_sd_stepwise, marker='o', label='single-delta')
    # #plt.plot(range(7), ip_fe_stepwise, marker='x', label='full env + V_C')
    # #plt.plot(6*np.ones(n_r), -ip_fe_final, marker='*', label='true fe')
    # plt.plot(6*np.ones(n_r), -ip_sd_final, marker='*', label='true sd', markersize=12)
    # plt.ylabel("IP, eV")
    # plt.xlabel("step number")
    # plt.legend()
    # plt.ylim([5,6.5])
    # plt.savefig("stepwise_{}.png".format(mol_idxs[my_mol_num]))
    # plt.close()
    #

    # IP fe
    # plt.plot(range(7), ip_sd_stepwise, marker='o', label='single-delta')
    plt.plot(range(7), ip_fe_stepwise, marker='x', label='full env + V_C')
    plt.plot(6*np.ones(n_r), -ip_fe_final, marker='*', label='true fe')
    # plt.plot(6*np.ones(n_r), -ip_sd_final, marker='*', label='true sd')
    plt.ylabel("IP, eV")
    plt.xlabel("step number")
    plt.legend()
    plt.ylim([5,6.5])
    plt.savefig("stepwise_fe_{}.png".format(mol_idxs[my_mol_num]))
    plt.close()


    # IP stepwise
    # plt.plot(range(7), ip_sd_stepwise, marker='o', label='single-delta')
    # plt.plot(range(7), ip_fe_stepwise, marker='x', label='full env + V_C')
    # plt.plot(6*np.ones(n_r), -ip_fe_final, marker='*', label='true fe')
    # plt.plot(6*np.ones(n_r), -ip_sd_final, marker='*', label='true sd')
    #
    # plt.ylabel("IP, eV")
    # plt.yscale('log')
    # plt.xlabel("step number")
    # plt.legend()
    # plt.ylim([5,6.5])
    # plt.savefig("stepwise_log_{}.png".format(mol_idxs[my_mol_num]))
    # plt.close()
    #
    np.savetxt('ip_sd_stepwise.txt', ip_sd_stepwise)
    np.savetxt('ip_fe_stepwise.txt', ip_fe_stepwise)


    exit()

########### TEST
    # read IPs from all existing folders
    IPs = []

    for i, folders in enumerate(updated_folders):
        print(i)
        my_ip = IP(updated_folders[i], my_radii_updated[i])
        my_ip.extract_IP()
        IPs.append(my_ip)


    # plot all single-delta
    all_IP = -np.array([IPs[i].mean_single_delta for i, j in enumerate(my_radii_updated)])
    # < old

    all_r = np.array([IPs[i].radius for i, j in enumerate(my_radii_updated)])
    A = -np.array([[IPs[i].mean_single_delta, IPs[i].mean_full_env, IPs[i].mean_single_delta_corrected, IPs[i].mean_full_env_corrected] for i, j in enumerate(my_radii_updated)])
    my_all_ips = AllIPs(A[:, 0], A[:, 1], A[:, 2], A[:, 3])
    print("my_all_ips.sd:", np.round(my_all_ips.sd, 2))

    # plot all vs R
    legend = ["single-delta", "full-envirement", "s-d corrected", "f-e corrected"]
    plt.plot(all_r, np.transpose([my_all_ips.sd, my_all_ips.fe, my_all_ips.sd_c, my_all_ips.fe_c]), marker='o', LineStyle='-')
    plt.legend(legend)
    plt.xlabel("$R, \AA$")
    plt.ylabel("IP, eV")
    #plt.xlim([np.min(all_r)-5, np.max(all_r)+5])
    plt.savefig("all_IP_vs_R.png")
    plt.close()


    # plot all vs 10/R
    legend = ["single-delta", "full-envirement", "s-d corrected", "f-e corrected"]
    plt.plot(10/all_r, np.transpose([my_all_ips.sd, my_all_ips.fe, my_all_ips.sd_c, my_all_ips.fe_c]), marker='o', LineStyle='-')
    plt.legend(legend)
    plt.xlabel("10/R, $10/\AA^{-1}$")
    plt.ylabel("IP, eV")
    #plt.xlim([np.min(all_r)-5, np.max(all_r)+5])
    plt.savefig("all_IP_inv.png")
    plt.close()


    # decrease the effective eps
    all_r -= 0
    #
    # plot single-delta vs 10/R together with corrections INV
    EPS = np.linspace(1, 5, 6)
    plt.plot(10.0/np.array(all_r), all_IP, LineStyle='--', marker='o', label='QP')
    for i, eps in enumerate(EPS):
        ip_corr = [correct_IP(IP_ini=all_IP[i], R=all_r[i], eps=eps) for i in range(len(my_radii_updated))]
        plt.plot(10.0/np.array(all_r), ip_corr, LineStyle=':', marker='x', label="$\epsilon$ = "+str(round(eps, 2)))
    plt.xlabel("10/R, $10/\AA^{-1}$")
    plt.ylabel("IP, eV")
    plt.legend()
    plt.savefig("mean_single_delta_INV.png")
    plt.close()

    # plot single-delta together with corrections vs R
    plt.plot(np.array(all_r), all_IP, LineStyle='--', marker='o', label='QP')
    for i, eps in enumerate(EPS):
        ip_corr = [correct_IP(IP_ini=all_IP[i], R=all_r[i], eps=eps) for i in range(len(my_radii_updated))]
        plt.plot(all_r, ip_corr, LineStyle=':', marker='x', label="$\epsilon$ = "+str(round(eps, 2)))
    plt.xlabel("$R, \AA$")
    plt.ylabel("IP, eV")
    plt.legend()
    plt.savefig("mean_single_delta.png")
    plt.close()

    # plot corrections for different R_cut
    for i, eps in enumerate(EPS):
        only_corr = [u_ana(R=all_r[i], eps=eps) for i in range(len(my_radii_updated))]
        plt.plot(10.0/np.array(all_r), only_corr, LineStyle=':', marker='x', label="$\epsilon$ = "+str(round(eps, 2)))
    plt.legend()
    plt.xlabel("10/R, Angstromes")
    plt.ylabel("Monopole correction, eV")
    plt.ylim([0, 0.4])
    plt.savefig("Monopole_correction.png")


    plt.close()


    # print correction for a typical (eps; R)
    eps = 3.5
    R = 20
    U_ana = (1 - 1. / eps) * const.value("Hartree energy") * const.value("Bohr radius") / ((R * 10 ** (-10)) * 2)
    print(U_ana*const.value("joule-electron volt relationship"))
    print(my_u_ana(R, eps))

    # fit the slope --> get the eps
    C = 7.1997176734999995  # coefficient in the formula
    coef_poly = np.polyfit(1.0/all_r, my_all_ips.sd, 1)
    print("coef_poly = ", coef_poly)
    print("Extracted dielectric permittivity:", C/(coef_poly[0]-C))

    #plot with the polynom INV
    EPS = np.linspace(1, 5, 6)
    plt.plot(10.0/np.array(all_r), all_IP, LineStyle='--', marker='o', label='QP')
    #for i, eps in enumerate(EPS):
    #    ip_corr = [correct_IP(IP_ini=all_IP[i], R=all_r[i], eps=eps) for i in range(len(my_radii_updated))]
    #    plt.plot(10.0/np.array(all_r), ip_corr, LineStyle=':', marker='x', label="$\epsilon$ = "+str(round(eps, 2)))
    plt.plot(10.0/all_r, coef_poly[1]+coef_poly[0]*1/all_r, label='poly', color = 'red')
    plt.xlabel("10/R, $10/\AA^{-1}$")
    plt.ylabel("IP, eV")
    plt.legend()
    plt.savefig("mean_single_delta_INV_compare.png")
    plt.close()


    pass

########################################### FUNCTIONS #######################################################

def return_target_folders1(pattern):  # to be deleted
    """
    return all folders of the kind "\d{2}{pattern}"
    :param: pattern (string)
    :returns: names of folders
    """

    r = re.compile(pattern)
    folder_list = os.listdir(".")
    print("Initial List: \n", folder_list, '\n')  # comment
    filtered = [folder for folder in folder_list if r.match(folder)]
    print("Filtered List: \n", filtered, '\n')  # comment
    return filtered

def return_target_folders(pattern):
    """
    return all folders of the kind "d{2}"+pattern"
    :param: pattern (string)
    :returns: names of folders
    """

    r = re.compile('\d{2}'+pattern)
    folder_list = os.listdir(".")
    #print("Initial List: \n", folder_list, '\n')  # comment
    filtered = [folder for folder in folder_list if r.match(folder)]
    #print("Filtered List: \n", filtered, '\n')  # comment
    return filtered



def get_radii(folders):
    """

    :param folders:
    :return: radii; type: int
    """
    radii = []
    for i, item in enumerate(folders):
        radius = folders[i]
        radii.append(int(radius[0:2]))
    #print("Radii:", radii)  # comment
    return(radii)

def return_folders_with_output(folders, out_file_name):
    folders_with_output = []
    for i, folder in enumerate(folders):
        path = folder+"/"+out_file_name
        #print(path)
        if os.path.isfile(path):
            folders_with_output.append(folder)

    return(folders_with_output)


class IP:
    def __init__(self, path, radius):
        self.path = path
        self.radius = radius
        fid = open(path+'/Analysis/EAIP/eaip_7b043a1428ff586ba119947b48219969_IP_summary.yml', 'r')
        self.yaml = yaml.load(fid)
        self.full_e_uncharged = 0
        self.full_e_charged = 0
        self.ip_sd_stepwise = 0
        self.ip_fe_stepwise = 0



    def extract_IP(self):
        self.mean_single_delta = self.yaml["mean_single_delta"]
        self.mean_full_env = self.yaml["mean_full_env"]
        self.mean_single_delta_corrected = self.yaml["mean_single_delta_corrected"]
        self.mean_full_env_corrected = self.yaml["mean_full_env_corrected"]
        self.raw_data = []

        for i, mol_idx in enumerate(self.yaml["raw_data"]):
            a = list(self.yaml["raw_data"].keys())[i]
            b = self.yaml["raw_data"][mol_idx]["full_env"]
            c = self.yaml["raw_data"][mol_idx]["single_delta"]
            self.raw_data.append(IPMethods(a, b, c))

    def extract_IP_stepwise(self, num_step, mol_idx):
        # charged
        fid = open(self.path+'/quantumpatch_runtime_files/Mol_{}_C_1/{}/energies.ene.yml'.format(mol_idx, num_step+9), 'r')
        self.a1 = yaml.load(fid)
        fid.close()
        self.ip_charged = self.a1[str(mol_idx)]["total"]
        # uncharged
        fid = open(self.path+'/quantumpatch_runtime_files/uncharged/{}/energies.ene.yml'.format(num_step+1), 'r')
        self.a2 = yaml.load(fid)
        fid.close()
        self.ip_uncharged = self.a2[str(mol_idx)]["total"]

        for mol_num in list(self.a1.keys()):
            self.full_e_charged += self.a1[mol_num]["total"]
            self.full_e_uncharged += self.a2[mol_num]["total"]
            diff = self.a1[mol_num]["total"] - self.a2[mol_num]["total"]
            if np.abs(diff) > 100:
                self.ip_fe_stepwise += 0  # if DFTB
            else:
                self.ip_fe_stepwise += diff
            #print('mol_num', mol_num)
            #print("e_diff charged/uncgarged", self.full_e_diff)
        self.ip_sd_stepwise = self.ip_charged - self.ip_uncharged
        #self.ip_fe_stepwise = self.full_e_charged - self.full_e_uncharged
        print("I am done")


class IPMethods:
    def __init__(self, mol_idx, full_env, single_delta):
        self.mol_idx = mol_idx
        self.full_env = full_env
        self.single_delta = single_delta


def correct_IP(IP_ini, R, eps):
    U_ana = (1 - 1. / eps) * const.value("Hartree energy") * const.value("Bohr radius") / ((R * 10 ** (-10)) * 2)
    U_ana  = U_ana/1.6E-19
    return IP_ini - U_ana


def u_ana(R, eps):
    p = (1 - 1. / eps) * const.value("Hartree energy") * const.value("Bohr radius") / ((R * 10 ** (-10)) * 2)
    p  = p*const.value("joule-electron volt relationship")
    return p


def my_u_ana(R, eps):
    a_0 = 0.529177  # A (Bohr R)
    E_H = 27.211  # eV (Hartree Energy)
    coef = E_H/2*a_0
    print("coef next to the monopole correction = ", coef)
    P_mon = coef*(1 - 1/eps)/R
    return P_mon


class AllIPs:
    def __init__(self, sd, fe, sd_c, fe_c):
        self.sd = sd
        self.fe = fe
        self.sd_c = sd_c
        self.fe_c = fe_c

if __name__ == '__main__':
    main()