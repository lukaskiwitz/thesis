#!/home/kiwitz/anaconda3/envs/fenics/bin/python
import glob
import os
import shutil
import sys
from typing import List
from my_debug import message

if len(sys.argv) > 1:
    PATH = sys.argv[1]
    if not PATH[-1] == "/":
        PATH = PATH + "/"
else:
    exit(1)

def check_dir(path: str, patterns: List[str]) -> List[str]:
    os.chdir(path)
    result = []
    for pattern in patterns:
        for g in glob.glob(pattern):
            result.append(g)
    return result



def make_copy(path: str, files: List[str]) -> None:

    path_split = PATH.split("/")
    new_path = "{top}_{name}".format(top="/".join(path_split[0:-1]),name=path_split[-1]+"pp/")
    os.makedirs(new_path,exist_ok=True)
    for file in files:
        message("coping file {f} from {old} to {new}".format(f=file,old=PATH,new=new_path))
        shutil.copyfile(PATH+file, new_path+file)


r = check_dir(PATH, ["*.scan","*.xml","*.h5"])
make_copy(PATH,r)

