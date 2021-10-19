import logging
import time as t

import mpi4py.MPI as MPI
from tqdm import tqdm

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def message(text: str):
    if rank == 0:
        logging.info(text)
        tqdm.write(get_cli_format(text))


def debug(text):
    if rank == 0:
        logging.debug(text)


def warning(text: str):
    if rank == 0:
        logging.warning(text)
        tqdm.write(get_cli_format(text))


def critical(text):
    if rank == 0:
        logging.critical(text)
        tqdm.write(get_cli_format(text))


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
