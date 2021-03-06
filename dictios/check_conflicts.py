import sys
import re

file=[]
conflicted=[]
good=[]
with open(sys.argv[1], "r") as inf:
    file=inf.readlines()
copy=[]
for l in file:
    l=l.rstrip("\n").split(",")
    copy.append(l[0])

joined=[x.rstrip("\n").replace("-","") for x in copy]
for l in file:
    l=l.rstrip("\n").split(",")
    name=l[0]
    regex=l[1]
    #joined=[x.rstrip("\n").replace("-","") for x in copy]
    matches=[]
    matches=re.findall(regex, " - ".join(joined).lower())
    if len(matches)>1:
        conflicted.append(name)
    else:
        good.append(",".join(l))

if conflicted:
    print("INFO: Done. Conflicts detected. Removed from output and saved in conflicted.txt", file=sys.stderr)
    with open("conflicted_targets.ignore","w") as outf:
        outf.write("\n".join(conflicted))
else:
    print("INFO: Done. No conflicts detected", file=sys.stderr)
print(*good, sep='\n')
