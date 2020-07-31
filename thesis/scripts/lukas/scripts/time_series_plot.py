import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def my_filter(x: pd.DataFrame):

    # x = x.loc[
    #     (
    #             (x["scan_index"] == 0) |
    #             (x["scan_index"] == 4) |
    #             (x["scan_index"] == 2)
    #     )
    #     ]

    return x

def my_lineplot(**kwargs):

    x = kwargs["data"]
    y = kwargs["y_name"]

    ids = x.loc[
        (
                (x["type_name"] == "Tfh")|
                (x["type_name"] == "Th1")
        )&
        (x["t"] == 1)
    ]["id"].unique()
    # message(ids)

    x = x.loc[~x["id"].isin(ids)]
    x = my_filter(x)
    # x = x.loc[
    #     (
    #             x["type_name"] != "Tn"
    #     )
    # ]

    sns.lineplot(x="t", y=y, data=x, hue="type_name", style="scan_index", sort=False,ci=None)
    plt.xlim([1,10])
    plt.ylim([0.5,1.5])

def my_count_plot(**kwargs):
    x = kwargs["data"]
    x = my_filter(x)

    sns.lineplot(x="t",y="n",data=x,hue="type_name",style="scan_index")
    plt.xlim([1, 50])

def my_global_plot(**kwargs):

    x = kwargs["data"]
    x = my_filter(x)

    sns.lineplot(
        x="timeIndex",
        y="concentration",
        data=x,
        hue="Cytokine",
        style="scan_index",
        palette=sns.color_palette(["red","seagreen","navy"])
    )

    plt.ylim([0,0.1])



IMGPATH = "/home/kiwitz/scan_plots/"
path = "/extra/kiwitz/scan_example_large/"

# cell_df = pd.read_hdf(path+"cell_df.h5",key="df")
# global_df = pd.read_hdf(path+"global_df.h5",mode="r")
# global_df = global_df.rename(columns={"field_name":"Cytokine"})
# global_df = global_df.rename({"scanIndex":"scan_index"},axis=1)

sns.set_context("paper",font_scale=1,rc={
            "lines.markersize":0,
            "lines.linewidth":1
            }
            )



# cell_df = pd.read_hdf(path+"cell_df.h5",mode="r")

# counts = cell_df.groupby(["type_name","t","scan_index"],as_index=False).count()
# counts = counts.drop(counts.columns.drop(["t","type_name","n","scan_index"]),axis=1)
#
# # master_test_df = normalize(cell_df)
# # my_lineplot(data=master_test_df,y_name="Th1_score_norm")
# # plt.show()
# # my_lineplot(data=master_test_df,y_name="Tfh_score_norm")
# plt.show()
# my_count_plot(data=counts)
# plt.show()
# my_global_plot(data=global_df)
# plt.show()

# plt.figure()
# facetgrid = sns.FacetGrid(master_test_df, col_wrap=3, legend_out=True)
# grid = facetgrid.map_dataframe(my_lineplot,y_name="Tfh_score_norm")
# grid.add_legend()
# # plt.savefig(IMGPATH+"Tfh_score.pdf",dpi=600)
#
# plt.figure()
# facetgrid = sns.FacetGrid(master_test_df, col_wrap=3, legend_out=True)
# grid = facetgrid.map_dataframe(my_lineplot,y_name="Th1_score_norm")
# grid.add_legend()
# # plt.savefig(IMGPATH+"Th1_score.pdf",dpi=600)
# plt.figure()
#
# facetgrid = sns.FacetGrid(counts, col_wrap=3, legend_out=True)
# grid = facetgrid.map_dataframe(my_count_plot)
# grid.add_legend()
# # plt.savefig(IMGPATH+"cell_number.pdf",dpi=600)
#
# plt.figure()
# facetgrid = sns.FacetGrid(global_df, col="scan_name", col_wrap=3, legend_out=True)
# grid = facetgrid.map_dataframe(my_global_plot)
# grid.add_legend()
# # plt.savefig(IMGPATH+"concentration.pdf",dpi=600)
# plt.show()



# plt.subplot(2,2,1)
# sns.lineplot(x="t",y="n",data=g.rename(columns={"type_name":"Cell Type"}), hue="Cell Type",ci="sd",sort=False)
# plt.ylabel("Number of cells")
# plt.xlabel("Time (a. u.)")
#
# plt.legend(["Tfh","Th1","Tn"],loc=1)
# plt.subplot(2,2,2)
# sns.lineplot(
#     x="timeIndex",y="concentration",hue="Cytokine",data=global_df.loc[
#         # (global_df["field_name"] == "il2")&
#         (global_df["scanIndex"] == scan_index)&
#         (global_df["scan_name"] == subfolders[0])
#         ],palette=sns.color_palette(["red","seagreen","navy"])
# )
# plt.legend(["IL-6","IL-2",r'INFG$_{\gamma}$'],loc=1)
# plt.ylabel("avg. Concentration (nM)")
# plt.xlabel("Time (a. u.)")
#
# plt.subplot(2,2,4)
# sns.lineplot(
#     x="timeIndex",y="gradient",hue="Cytokine",data=global_df.loc[
#         # (global_df["field_name"] == "il2")&
#         (global_df["scanIndex"] == scan_index)&
#         (global_df["scan_name"] == subfolders[0])
#         ],palette=sns.color_palette(["red","seagreen","navy"])
# )
# # plt.tight_layout()
# plt.legend(["IL-6","IL-2",r'INFG$_{\gamma}$'],loc=1)
# plt.ylabel("avg. Gradient (nM/dm)")
# plt.xlabel("Time (a. u.)")
# plt.tight_layout()
# plt.savefig("/home/kiwitz/scan_default.pdf",dpi=600)
# plt.show()

