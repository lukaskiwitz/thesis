import multiprocessing.pool

"""
https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
"""


class NoDaemonProcess(multiprocessing.Process):
    """mainly for compatibility with github CI"""

    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


class Pool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess
