# computes number density
import numpy as np
from scipy import constants as c

molar_weight = 588.74  #g/mole aNPD
density = 1.1327 # g/cm**3 aNPD

#molar_weight =  740.89# g/mole TCTA
#density = 1.1272 # g/cm^2

density = 1.7073 #C60
molar_weight = 720.66 #C60

m_molecule = molar_weight/c.N_A

number_density = density / m_molecule  # molecules / cm^3


number_density_m3 = number_density*(100)**3 # 1/m**3
number_density_nm3 = number_density_m3*1E9**(-3)

print('number_density_nm3 = ', number_density_nm3)


#SpiroOMeTAD

number_density_nm3 = 0.48 # mol/nm**3
molar_weight_g = 1225.4 # g/mole
m_molecule_g = molar_weight_g/c.N_A

density_g_nm3 = m_molecule_g*number_density_nm3

print('density SpiroOMeTAD {} [g/nm**3]'.format(density_g_nm3))
density_g_cm3 = density_g_nm3*(1E-9/1E-2)**-3
print('density SpiroOMeTAD {} [g/cm**3]'.format(density_g_cm3))
