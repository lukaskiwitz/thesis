import logging
import logging.handlers
import os
import time as t

import mpi4py.MPI as MPI
from tqdm import tqdm

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def setup_loggers(path: str, debug=False):
    if os.path.exists(os.path.join(path, "sim.log")): os.remove(os.path.join(path, "sim.log"))

    logging.getLogger().handlers = []
    logger = logging.getLogger('thesis.main')

    if debug:
        if os.path.exists(os.path.join(path, "debug.log")): os.remove(os.path.join(path, "debug.log"))
        logger.setLevel(logging.DEBUG)
        debug_handler = logging.handlers.RotatingFileHandler(os.path.join(path, "debug.log"), maxBytes=int(1e8),
                                                             backupCount=10, mode="w")
        debug_handler.setLevel(logging.DEBUG)
        debug_handler.setFormatter(
            logging.Formatter('%(asctime)2s %(name)-40s %(levelname)-5s %(message)s', datefmt='%m-%d %H:%M'))
        logger.addHandler(debug_handler)
    else:
        logger.setLevel(logging.INFO)

    info_handler = logging.handlers.RotatingFileHandler(os.path.join(path, "sim.log"), maxBytes=int(1e6),
                                                        backupCount=10, mode="w")
    info_handler.setLevel(logging.INFO)
    info_handler.setFormatter(
        logging.Formatter('%(asctime)2s %(name)-40s %(levelname)-5s %(message)s', datefmt='%H:%M'))
    logger.addHandler(info_handler)


def message(text: str, logger=None):
    if rank == 0:
        if logger is not None:
            logger.info(text)
        else:
            logging.info(text)
        try:
            tqdm.write(get_cli_format(text))
        except BlockingIOError:
            pass


def info(text, logger=None):
    if rank == 0:
        if logger is not None:
            logger.info(text)
        else:
            logging.info(text)


def debug(text, logger=None):
    if rank == 0:
        if logger is not None:
            logger.debug(text)
        else:
            logging.debug(text)


def warning(text: str, logger=None):
    if rank == 0:
        if logger is not None:
            logger.warning(text)
        else:
            logging.info(text)
        try:
            tqdm.write(get_cli_format(text))
        except BlockingIOError:
            pass


def critical(text, logger=None):
    if rank == 0:
        if logger is not None:
            logger.critical(text)
        else:
            logging.info(text)
        try:
            tqdm.write(get_cli_format(text))
        except BlockingIOError:
            pass


def get_cli_format(text):
    l = t.localtime(t.time())
    text = "{h}:{m}:{s}: {message}".format(h=l.tm_hour, m=l.tm_min, s=l.tm_sec, message=text)
    return text


def total_time(time: float, pre="", post="", logger=None):
    l = t.gmtime(time)

    text = "{h}h:{m}m:{s}s".format(
        h=l.tm_hour,
        m=l.tm_min,
        s=l.tm_sec
    )
    message(pre + text + post, logger)
