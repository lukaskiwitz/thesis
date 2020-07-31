from thesis.main.MyPlotter import Plotter
# from parameters import path
import matplotlib.pyplot as plt
import seaborn as sns


path = '/extra/lukas/example_min/test/'
plotter = Plotter(path)

sns.set_context("paper", font_scale=2, rc={
            "lines.markersize": 4,
            "lines.linewidth": 5
        })

# fig = plt.figure()
# ax = fig.gca()


# plotter.global_time_series_plot(fig,ax,"Concentration","Concentration")
# plt.show()
# plotter.steady_state_plot(
#     fig, ax,
#     "scan_index",
#     "Concentration",
#     hue= "field_name",
#     style="scan_index")
fig = plt.figure()
ax = fig.gca()


plotter.count_plot(fig,ax,relative=True,ylim=(0,1))
# plotter.cells_time_series_plot(fig,ax,"abs_score_init","score")
plt.show()

