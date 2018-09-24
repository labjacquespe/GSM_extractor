import re
import sys

lines=[]
with open(sys.argv[1],"r") as inf:
  lines=inf.readlines()

for name in lines:
  name=name.replace("\n","")
  if re.search(r'\d+$', name):
    line="{},({}(\D+|$))".format(name.upper(),"-?".join(re.split('(\d+)',name.replace("-","-?"))).rstrip("-?").lower())
  else:
    line="{},({})".format(name.upper(),"-?".join(re.split('(\d+)',name.replace("-","-?"))).rstrip("-?").lower())
  print(line)
