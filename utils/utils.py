from subprocess import PIPE, Popen
from os import mkdir, path
import warnings

def custom_formatwarning(msg, *a):
    # ignore everything except the message
    return str(msg) + '\n'
warnings.formatwarning = custom_formatwarning

def run(cline):
    p = Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    return out, err

def printif(string):
    if len(string) > 0:
        print string

def makedir(new_dir, force = False):
    if not path.exists(new_dir):
        mkdir(new_dir)
        msg = "New directory created at %s" % path.abspath(new_dir)
        warnings.warn(msg)
    elif force:
        rmtree(new_dir)
        mkdir(new_dir)
        msg = "Directory overwritten at %s" % path.abspath(new_dir)
        warnings.warn(msg)
    else:
        msg = "Directory already exists at %s" % path.abspath(new_dir)
        warnings.warn(msg)