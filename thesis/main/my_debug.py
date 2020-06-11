import os
import time as t

import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def message(text: str):
    if rank == 0:
        log_path = os.environ.get("LOG_PATH") if os.environ.get("LOG_PATH") else "./"
        l = t.localtime(t.time())
        text = "{h}:{m}:{s}: {message}".format(h=l.tm_hour, m=l.tm_min, s=l.tm_sec, message=text)

        file = log_path + "/debug.log"
        os.makedirs(log_path, exist_ok=True)
        try:
            if not os.path.isfile(file):
                with open(file, "x") as f:
                    pass
        except Exception as e:
            pass
        try:
            with open(file, "a") as f:
                f.write(text + "\n")
            print(text)
        except Exception as e:
            pass


def total_time(time: float, pre="", post=""):
    l = t.gmtime(time)

    text = "{h}h:{m}m:{s}s".format(
        h=l.tm_hour,
        m=l.tm_min,
        s=l.tm_sec
    )
    message(pre + text + post)

def debug(text: str):
    if "DEBUG" in  globals() and globals()["DEBUG"]:
        message("DEBUG: "+text)


def alert(text: str):
    text = str(text)
    message("---------------" + text)
