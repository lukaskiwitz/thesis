# from thesis.cellBehaviourUtilities.cell_solver import kineticSolver
#
# c_naive = (0.9058823529411765, 0.1607843137254902, 0.5411764705882353)
# c_treg = (0.9019607843137255, 0.6705882352941176, 0.00784313725490196)
#
# s = kineticSolver()
# x = np.linspace(0,8e3,1000)
# y = s.EC50_calculation(E_max=125e-12, E_min=0, k=860, N=1, R=x)
# plt.plot(x,y, color = c_naive)
# ylim = plt.gca().get_ylim()
# plt.vlines(100,*ylim, color = c_naive)
#
# y = s.EC50_calculation(E_max=125e-12, E_min=0, k=0.25 * 860, N=1, R=x)
# plt.plot(x,y, color = c_treg)
# plt.vlines(1000,*ylim, color = c_treg)
#
# plt.legend(["naive","treg"])
# plt.show()
#
#
#

# x = [100,200]
# X = np.meshgrid(x,x,x)
# X = np.array([np.ravel(i) for i in X])
#
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(X[0,:],X[1,:],X[2,:])
# plt.show()


# plotter = Plotter("/extra/kiwitz/Treg_competition_parallel_low_sec/200_3D/20220201_run_3_replicats/blub/")
# df = plotter.cell_df
#
# def get_steady_state(df):
#     def f(df):
#         return df.loc[(df.time_index == df.time_index.max()) | (df.time_index == 0)]
#
#     ind = pd.MultiIndex.from_frame(
#         df.loc[:,[plotter.model_index_key,plotter.scan_index_key,plotter.scan_name_key,plotter.replicat_index_key]]
#     )
#     df = df.set_index(ind,drop=True)
#     df = df.drop(columns=df.index.names)
#     t_max = df.groupby(axis=0, level= df.index.names,group_keys=False).apply(f)
#     t_max = t_max.reset_index()
#
#     return t_max
#


# df = plotter.cell_df
# df = get_steady_state(df)
# df.to_hdf(os.path.join("/extra/kiwitz/Treg_competition_parallel_low_sec/200_3D/20220201_run_3_replicats/blub/", "ss_cell_df.h5"), key="df", mode="w")
#
# df = plotter.global_df
# df = get_steady_state(df)
# df.to_hdf(os.path.join("/extra/kiwitz/Treg_competition_parallel_low_sec/200_3D/20220201_run_3_replicats/blub/", "ss_global_df.h5"), key="df", mode="w")
