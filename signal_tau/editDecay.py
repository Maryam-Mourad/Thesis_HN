import sys

mass = sys.argv[1]

with open("HN"+mass+"/signal_tau_mass"+mass+".decay","r") as decFile:
  content = decFile.read()
  
# replace the name of HN
content = content.replace("HN","HN"+mass)

# Overwrite the file
with open("HN"+mass+"/signal_tau_mass"+mass+".decay","w") as decFile:
  decFile.write(content)