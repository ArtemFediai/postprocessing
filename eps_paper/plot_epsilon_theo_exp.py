#taken from Jonas
import yaml
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.patches as mpatches

with open("exp_theo_paper.yml") as fid:
    data = yaml.load(fid, Loader=yaml.FullLoader)

# #comparison to experiment
# fig, ax = plt.subplots()
# ax.errorbar(data["exp_vals"], data["theo_vals"], xerr=data["exp_stds"], yerr=data["theo_stds"], marker="o", linestyle="none")
# plt.ylim([0,5])
# plt.xlim([0,5])
# ax.plot([0, 1], [0, 1], transform=ax.transAxes)
# plt.xlabel("Experimental relative dielectric permittivity")
# plt.ylabel("Theoretical relative dielectric permittivity")
# plt.show()

# comparison to clausius - mossotti equation
fig, ax = plt.subplots(dpi=150)
x = np.arange(len(data["theo_vals"]))
plt.xticks(x, data["molecules"])
ax.plot(x, data["exp_vals"], marker="x", markersize=8, color="black", linestyle="none", zorder=5,label="Experiment")
ax.plot(x, data["theo_vals"], marker="o", markersize=6, linestyle="none", label="Simulation") #yerr=data["theo_stds"],, fmt='o', capsize=2, capthick=5
ax.plot(x, data["claussius_mossotti_vals"], marker="^", color="darkorange", linestyle="None",zorder=5, label="CM relation")
mean_CM = np.mean(data["claussius_mossotti_vals"])
#ax.plot([0,1,2], [mean_CM]*3, marker="None", color="darkorange", linestyle="-")

#rect
rect_height = np.max(data["claussius_mossotti_vals"]) - np.min(data["claussius_mossotti_vals"])
rect = mpatches.Rectangle([0, np.min(data["claussius_mossotti_vals"])], 2, rect_height, color="lightyellow", ec="none")
ax.add_patch(rect)
plt.ylim([0,5])


plt.xlabel("Molecule")
plt.ylabel("Relative dielectric permittivity")
plt.legend()
fig.savefig("epsilon_exp_theo.png")
plt.show()