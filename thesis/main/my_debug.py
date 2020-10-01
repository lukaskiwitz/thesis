import os
import time as t
import logging
import mpi4py.MPI as MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def message(text: str):
    if rank == 0:
        logging.info(text)
        print(get_cli_format(text))


def debug(text):
    if rank == 0:
        logging.debug(text)
        # print(get_cli_format(text))


def get_cli_format(text):

    l = t.localtime(t.time())
    text = "{h}:{m}:{s}: {message}".format(h=l.tm_hour, m=l.tm_min, s=l.tm_sec, message=text)
    return text

def total_time(time: float, pre="", post=""):
    l = t.gmtime(time)

    text = "{h}h:{m}m:{s}s".format(
        h=l.tm_hour,
        m=l.tm_min,
        s=l.tm_sec
    )
    message(pre + text + post)


def warning(text: str):
    if rank == 0:
        logging.warning(text)
        print(get_cli_format(text))
