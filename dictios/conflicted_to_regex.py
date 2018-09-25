import re
import sys

lines=[]
with open(sys.argv[1],"r") as inf:
  lines=inf.readlines()

for name in lines:
  name=name.replace("\n","")
  if re.search(r'\d+$', name):
    line="{},((?<![a-z]){}(\D+|$))".format(name.upper(),"-?".join(re.split('(\d+)',name.replace("-","-?"))).rstrip("-?").lower())
  else:
    line="{},((?<![a-z]){}(?![a-z]))".format(name.upper(),"-?".join(re.split('(\d+)',name.replace("-","-?"))).rstrip("-?").lower())
  print(line)
