import glob

import pandas as pd


def load_csv(path):
    files = []
    for file_path in glob.glob(path + "*.csv"):
        with open(file_path) as f:
            df = pd.read_csv(f, names=["x", "y"])

            type_name = file_path.split("/")[-1].split(".")[0].split("_")[0]

            series = pd.Series([type_name] * len(df))
            series.index = df.index

            df["type_name"] = series

            series = pd.Series([0] * len(df))
            series.index = df.index

            df["time_index"] = series

            series = pd.Series([0] * len(df))
            series.index = df.index

            df["z"] = series

            files.append(df)

    cells = pd.concat(files)
    return cells
