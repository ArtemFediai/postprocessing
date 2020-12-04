import seaborn as sns
sns.set_theme(style="whitegrid")
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches

#sns.set_theme(style="whitegrid")

with open("exp_theo_paper.yml") as fid:
    data = yaml.load(fid, Loader=yaml.FullLoader)

df = pd.DataFrame(data=data)

fig, ax1 = plt.subplots(figsize=(4,4))

rows4panda = []
#-> work with data
eps_kinds = ['theo_vals', 'exp_vals']
for i, molecule in enumerate(data['molecules']): #
    for j, eps_kind in enumerate(eps_kinds):
        rows4panda.append([molecule, eps_kind, data[eps_kind][i]])
#<- work with data

df1 = pd.DataFrame(data=rows4panda, columns=['molecules', 'kind_of_eps', 'eps_value'])

g = sns.lineplot(x=data['molecules'],
                 y=data['claussius_mossotti_vals'],
                 color='r',
                 marker="o",
                 ax=ax1,
                 label='Claussius-Mossotti', legend=False)

sns.barplot(
    data=df1,
    x='molecules', y='eps_value', hue='kind_of_eps',
    ci='sd', palette='dark', alpha=.6,
    ax=ax1)

plt.ylim(top=5.5)
plt.ylabel('Dielectric permittivity, $\mathrm{\epsilon}_\mathrm{r}$')
plt.xlabel(None)

labels = ["this work", "experiment", "Claussius-Mossotti"]

h, l = g.get_legend_handles_labels()
g.legend([h[1], h[2], h[0]], labels, title=None)

plt.savefig('eps_comp.svg')
plt.savefig('eps_comp.eps')
plt.savefig('eps_comp.png', dpi=600)
#plt.show()
#sns.set()