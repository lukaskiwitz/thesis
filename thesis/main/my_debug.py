import os
import time as t
import logging
import mpi4py.MPI as MPI
from time import sleep
import contextlib
import sys
from tqdm import tqdm
from tqdm.contrib import DummyTqdmFile

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

import os

# so = se = open(os.environ["LOG_PATH"]+"full.log", 'w', 0)
# # re-open stdout without buffering
# sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
# # redirect stdout and stderr to the log file opened above
# os.dup2(so.fileno(), sys.stdout.fileno())
# os.dup2(se.fileno(), sys.stderr.fileno())


# @contextlib.contextmanager
# def std_out_err_redirect_tqdm():
#     orig_out_err = sys.stdout, sys.stderr
#     try:
#         sys.stdout, sys.stderr = map(DummyTqdmFile, orig_out_err)
#         return orig_out_err[0]
#     # Relay exceptions
#     except Exception as exc:
#         raise exc
#     # Always restore sys.stdout/err if necessary
#     finally:
#         sys.stdout, sys.stderr = orig_out_err


import os
# class OutputFilter:
#
#     def __init__(self, stream):
#         self.stream = stream
#     def __getattr__(self, attr_name):
#         return getattr(self.stream, attr_name)
#
#     def clear(self):
#         os.system("clear")
#
#     def write(self, data):
#
#         self.stream.write(data)
#         self.stream.flush()
#
#     def flush(self):
#         self.stream.flush()

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


