# # Remove families not in famsUnderSelection

import os
from shutil import copyfile
import pandas as pd

fams = pd.read_csv("final_results/famsUnderSelection.txt", r"\s+",
                   index_col=False)

families = fams["family"]

families_in_dir = os.listdir("families")

families_in_dir[2].split(".")[0]

for i in range(0, len(families_in_dir)):
    if families_in_dir[i].split(".")[0] in list(families):
        copyfile("families/" + families_in_dir[i],
                 "families_fubar/" + families_in_dir[i])
