# # Remove families not in famsUnderSelection

import os
from shutil import move
import pandas as pd
import glob
import re

os.mkdir("final")
output_file = "final/absrel.txt"

logs = glob.glob("families_fubar/*.log")

with open(output_file, "w") as out:
    out.write("family num_SUS\n")
    for currentFile in input.log:
        with open(currentFile) as f:
            for line in f:
                if "Likelihood" in line:
                    if "**0**" not in line:
                        result = re.search('found(.*)branches', line)
                        out.write(currentFile.split(
                            "/")[-1].split(".")[0] + " " + result.group(1) + "\n")

# move all positive files to new folder

fams = pd.read_csv("final/absrel.txt", r"\s+",
                   index_col=False)

families = fams["family"]

families_in_dir = os.listdir("families_fubar")

os.mkdir("families_absrel/")

for i in range(0, len(families_in_dir)):
    if families_in_dir[i].split(".")[0] in list(families):
        move("families_fubar/" + families_in_dir[i],
                 "families_absrel/" + families_in_dir[i])
