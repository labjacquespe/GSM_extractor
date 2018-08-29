import sys
import re

file=[]
with open(sys.argv[1], "r") as inf:
    file=inf.readlines()
copy=[]
for l in file:
    l=l.rstrip("\n").split(",")
    copy.append(l[0])


for l in file:
    l=l.rstrip("\n").split(",")
    name=l[0]
    regex=l[1]
    joined=[x.rstrip("\n").replace("-","") for x in copy]
    matches=[]
    matches=re.findall(regex, " - ".join(joined).lower())
    if len(matches)>2:
        print(name)
    joined=""