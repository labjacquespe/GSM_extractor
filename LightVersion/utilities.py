import os
import sys
from collections import OrderedDict

def ReadConfig():
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

def GetDictio(name):
    commonDictio={}
    with open(name,"r") as f:
        lines=f.readlines()
    for line in lines:
        line=line.rstrip("\n").split(",")
        commonDictio[line[0]]=line[1]
    return OrderedDict(commonDictio)