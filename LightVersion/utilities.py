import os
import sys

def readConfig():
    args={}
    with open("config.ini", "r") as conf:
        configFiles=conf.readlines()
    for line in configFiles:
        if not line.startswith("#"):
            try:
                line=line.rstrip("\n").split("=")
                args[line[0]]=line[1]
            except:
                print("ERROR: invalid config file", file=sys.stderr)
                sys.exit()  
    return args