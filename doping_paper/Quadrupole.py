import matplotlib.pyplot as plt
import scipy.constants
import numpy as np
import yaml
import util.xyz as xyz

EPS_0 = scipy.constants.epsilon_0
print(EPS_0)
Q_0 = scipy.constants.elementary_charge
print(Q_0)
a_0 = scipy.constants.physical_constants['Bohr radius'][0]
print(a_0)

# set paths
PATH_TO_QUADRUPOLE_MOMENTS_DICT = 'data/QM/QM.yaml'  # 'mol_name': {xx: .., yy: .., zz: ..}
PATH_TO_MOLS_XYZS = '/home/artem/soft/postprocessing/doping_paper/data/mols_xyzs'
mol_name = 'F6TCNNQ'
xyz_name = 'F6TCNNQ_zentriert_relaxiert.xyz'

# linear_sizes_of_the_molecule = {x: , y:, z:}
# end: set paths


def main():

    xyz_object = xyz.XYZ.from_file(PATH_TO_MOLS_XYZS + '/' + xyz_name)
    offsets = xyz_object.compute_box_size(offset=2.5)

    with open(PATH_TO_QUADRUPOLE_MOMENTS_DICT) as yaml_file:
        dict_with_quadrupole_moments = yaml.safe_load(yaml_file)

    qm_along_direction_in_au = dict_with_quadrupole_moments[mol_name]  # xx, not x!

    # test one distance
    distance_in_A = 3  # A
    eps_r = 2.7
    q_in_au = 50
    q_in_si = convert_from_au_to_si(q_in_au)

    qm = energy_quadrupole_monopole_interaction(distance_in_A, q_in_si, eps_r)
    print(f'Monopole-quadrupole interaction at distance {distance_in_A} [A]: {qm/Q_0} eV')
    mm = monopole_monopole_interaction(distance_in_A, eps_r)
    print(f'monopole-monopole interaction at distance {distance_in_A} [A]: {mm/Q_0} eV')
    # end: test one distance

    # distances
    directions = ('x', 'y', 'z')
    distances_along_direction = {}

    qm_along_direction_in_si = {}
    for direction in directions:
        qm_along_direction_in_si[direction] = convert_from_au_to_si(qm_along_direction_in_au[direction*2])  # x <-- xx

    for i, direction in enumerate(directions):
        distances_along_direction[direction] = np.linspace(offsets[i], 30, 100)  # in A!

    m_q_interactions = {}
    m_m_interactions = {}
    for direction in directions: 
        m_q_interactions[direction] = np.array([energy_quadrupole_monopole_interaction(distance_in_A=distance, q_zz=qm_along_direction_in_si[direction], eps_r=eps_r)
                            for distance in distances_along_direction[direction]])
        m_m_interactions[direction] = np.array([monopole_monopole_interaction(distance_in_A=distance, eps_r=eps_r)
                            for distance in distances_along_direction[direction]])

    plot_along_directions(directions, distances_along_direction, m_q_interactions, m_m_interactions)


def plot_along_directions(directions, distances_along_direction, m_q_interactions, m_m_interactions, ylim=[-1.0, 1.0]):
    plt.plot()
    for i, direction in enumerate(reversed(directions)):
        plt.plot(distances_along_direction[direction], m_q_interactions[direction]/Q_0, label='m-q', color=f'C{i}', linestyle=':')
        plt.plot(distances_along_direction[direction], m_m_interactions[direction]/Q_0, label='m-m', color=f'C{i}', linestyle='--')
        plt.plot(distances_along_direction[direction], (m_m_interactions[direction] + m_q_interactions[direction])/Q_0, label='m-m + m-q', color=f'C{i}')
    plt.xlabel('distance, $\AA$')
    plt.ylabel('V, eV')
    plt.ylim([-1.0, 1.0])
    plt.legend()
    plt.grid()
    plt.show()
    plt.close()


def convert_from_au_to_si(q_component_in_au):
    """
    Takes Q.M. in a.u., returns it in SI units
    :param q_component_in_au: zz-component of the quadrupole moment tensor
    :return: q_component_in_SI
    """
    global Q_0, a_0
    return q_component_in_au * Q_0 * a_0 * a_0


def energy_quadrupole_monopole_interaction(distance_in_A, q_zz, eps_r, charge_in_au=1):
    """
    Returns energy of the quadrupole-monopole interaction.
    :param distance_in_A: distance from the quadrupole to a given point
    :param q_zz: zz-component of the quadrupole
    :param eps_r: dielectric permittivity of the medium
    :param charge_in_au: charge of the monopole in a.u. 1 == elementary charge, abs.
    :return: energy of the monopole-quadrupole interaction at a given point in eV
    """
    global Q_0, a_0, EPS_0
    charge_in_SI = charge_in_au * Q_0
    distance_in_SI = distance_in_A*1.0E-10
    coef = charge_in_SI / (8 * np.pi * EPS_0 * eps_r)
    return coef * q_zz / distance_in_SI**3


def monopole_monopole_interaction(distance_in_A, eps_r, q1_in_au=1, q2_in_au=-1):
    """
    Returns energy of the quadrupole-monopole interaction.
    :param distance_in_A: distance from the quadrupole to a given point
    :param q_zz: zz-component of the quadrupole
    :param eps_r: dielectric permittivity of the medium
    :param charge_in_au: charge of the monopole in a.u. 1 == elementary charge, abs.
    :return: energy of the monopole-quadrupole interaction at a given point in eV
    """
    global Q_0, a_0, EPS_0
    q1_in_SI = q1_in_au * Q_0
    q2_in_SI = q2_in_au * Q_0
    distance_in_SI = distance_in_A*1.0E-10
    coef = q1_in_SI / (4 * np.pi * EPS_0 * eps_r)
    return coef * q2_in_SI / distance_in_SI

if __name__ == '__main__':
    main()

# todo. plot total coulomn interaction energy as a function of the distance.
# todo. how to deal with the dipole moement, if deal at all