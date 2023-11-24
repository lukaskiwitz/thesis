import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks

rc_ticks['figure.figsize'] = (1.7, 0.7)
sns.set_theme(context = "talk", style = "ticks", rc = rc_ticks)
fig,ax = plt.subplots()

mean = 1e3
variance_range = [10 * mean, mean, mean/10]

color = "black"
alphas = [1, 0.4, 0.2]
for var, variance in enumerate(variance_range):
    tmp_sigma = np.sqrt(np.log(variance ** 2 / mean ** 2 + 1))
    log_mean = np.log(mean) - 1 / 2 * tmp_sigma ** 2
    draw = np.clip(np.random.lognormal(log_mean, tmp_sigma, 5000 if var != 2 else 4000),None, 30000)
    max = np.max(draw)
    min = np.min(draw)
    bin_width = 48
    bins = (max-min)/bin_width

    sns.histplot(data = draw, bins=int(bins), label=variance, color=color, alpha=alphas[var], log_scale=False, shrink=.94)

handles, labels = fig.axes[0].get_legend_handles_labels()
labels = [10, 1, 0.01]
ax.set(xscale="linear")
plt.xlim((0,2000))
plt.ylim((0,800))
plt.xticks([0,1000,2000])
plt.yticks([0,400,800])

plt.xlabel("receptors")
plt.ylabel("Frequency")
saving_string = r"/home/brunner/Documents/Current work/2023_10_20/" + "lognorm_dist"  + ".pdf"
fig.savefig(saving_string, bbox_inches='tight', transparent=True)
plt.tight_layout()
plt.show()