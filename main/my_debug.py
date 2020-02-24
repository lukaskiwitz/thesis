import time as t


def message(text: str):
    l = t.localtime(t.time())
    print("{h}:{m}:{s}: {message}".format(
        h=l.tm_hour,
        m = l.tm_min,
        s = l.tm_sec,
        message = text
    ))


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
