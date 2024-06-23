import os
import sys

ID = sys.argv[1]
mass = sys.argv[2]
ctau = sys.argv[3]

parAddress="../../config/particles.dat"
already_there = False
appending_text = ID + "  HN" + mass + " --- " + mass + " 0 0 0.5 --- " + ctau + "\n"

#if not os.path.exists("HN"+mass):
#  os.mkdir("HN"+mass)
#  print("HN directory created.")



with open(parAddress,"r") as parFile:
  current_content = parFile.read()
  if ("HN"+mass in current_content):
    print("HN with this mass exists in particles.dat")
    already_there = True


if (not already_there):
  with open(parAddress,"a") as parFile:
    if(current_content[-1] == "\n"):
      parFile.write(appending_text)
    else:
      parFile.write("\n"+appending_text)
  print("Added HN to particles.dat")
    