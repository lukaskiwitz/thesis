import os
import time as t


def message(text: str):
    log_path = os.environ.get("LOG_PATH") if os.environ.get("LOG_PATH") else "./"
    l = t.localtime(t.time())
    text = "{h}:{m}:{s}: {message}".format(h=l.tm_hour, m=l.tm_min, s=l.tm_sec, message=text)

    file = log_path + "/debug.log"
    os.makedirs(log_path, exist_ok=True)
    if not os.path.isfile(file):
        with open(file, "x") as f:
            pass

    with open(file, "a") as f:
        f.write(text + "\n")
    print(text)


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
