import glob
import os
import shutil
from typing import List

PATH = "/extra/kiwitz/sensitivity/"


def check_dir(path: str, pattern: str) -> List[str]:
    os.chdir(path)
    return glob.glob("**/"+pattern, recursive=True)


def make_copy(path_list: List[str], old_prefix: str, new_prefix: str) -> None:
    for file in path_list:
        # print(new_prefix+"".join(file.split("/")[0:-1]))
        os.makedirs(new_prefix+"".join(file.split("/")[0:-1]), exist_ok=True)
        shutil.copyfile(old_prefix+file, new_prefix+file)


NEW_PREFIX = "/extra/kiwitz/sensitivity_result/"
make_copy(check_dir(PATH, "*.xml"), PATH, NEW_PREFIX)
make_copy(check_dir(PATH, "global_dataframe.h5"), PATH, NEW_PREFIX)
make_copy(check_dir(PATH, "dataframe.h5"), PATH, NEW_PREFIX)
