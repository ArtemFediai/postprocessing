from scipy import constants
import numpy as np

# constant in the expression eps = c/(c - slope).
# c = e^2/(8*pi*eps_0)


e = constants.elementary_charge
eps0 = constants.epsilon_0

constant = e*e / (8 * np.pi * eps0)

print(constant)

constant = constant/constants.electron_volt/constants.angstrom

print(constant)
