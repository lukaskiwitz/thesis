import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from plotting_rc import rc_ticks

sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

mean = 5e3
variance_range = [3 * mean, mean, mean/2, mean/10]

colors = ["cornflowerblue", "salmon", "lightgreen", "xkcd:orangeish"]

for var, variance in enumerate(variance_range):
    tmp_sigma = np.sqrt(np.log(variance ** 2 / mean ** 2 + 1))
    log_mean = np.log(mean) - 1 / 2 * tmp_sigma ** 2
    draw = np.clip(np.random.lognormal(log_mean, tmp_sigma, 5000),None, 30000)
    max = np.max(draw)
    min = np.min(draw)
    bin_width = 200
    bins = (max-min)/bin_width

    sns.histplot(data = draw, bins=int(bins), label=variance, color=colors[var], log_scale=False, alpha=0.7, shrink=.94)

handles, labels = fig.axes[0].get_legend_handles_labels()
labels = [3, 1, 0.5, 0.01]
ax.set(xscale="linear")
plt.xlim((4e2,14000))

plt.legend(handles=reversed(handles), labels = reversed(labels), title="Heterog.")
plt.xlabel("Receptors")
plt.ylabel("Frequency")
saving_string =r"/home/brunner/Documents/Current work/2022_01_21/" + "lognorm_dist"  + ".svg"
fig.savefig(saving_string, bbox_inches='tight')
plt.tight_layout()
plt.show()