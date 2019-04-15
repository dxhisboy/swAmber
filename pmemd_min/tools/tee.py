import subprocess
import os
def tee(filename, *streams):
    tee = subprocess.Popen(["tee", filename], stdin=subprocess.PIPE)
    for stream in streams:
        os.dup2(tee.stdin.fileno(), stream.fileno())


