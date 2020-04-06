from shutil import move
import pandas as pd
import os
import re
import glob

# makes final_stats file

os.mkdir("final")
output_file = "final/fubar.txt"

logs = glob.glob("families/*.log")

with open(output_file, "w") as out:
  out.write("family number_SUS\n")
  for log in logs:
    with open(log) as f:
      for line in f:
          if "## FUBAR" in line:
              if "no" not in line:
                  result = re.search("inferred(.*)sites", line)
                  out.write(log.split("/")[-1].split(".")[0] + " " + result.group(1) + "\n")

# move all positive files to new folder

fams = pd.read_csv("final/fubar.txt", r"\s+",
                   index_col=False)

families = fams["family"]

families_in_dir = os.listdir("families")

os.mkdir("families_fubar/")

for i in range(0, len(families_in_dir)):
    if families_in_dir[i].split(".")[0] in list(families):
        move("families/" + families_in_dir[i],
                 "families_fubar/" + families_in_dir[i])
