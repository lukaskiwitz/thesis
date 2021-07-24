import os.path
from parameters import path
import numpy as np
import pandas as pd

df = pd.read_hdf(os.path.join(path,"global_df.h5"), mode="r")
error_norm_limit = 1

e = np.array(df.error)


def is_converging(error):

    error = np.abs(error)
    return len(np.where(np.gradient(error) >= 0)[0]) == 0

assert is_converging(e)
assert abs(e[-1]) < error_norm_limit
