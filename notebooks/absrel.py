#!/usr/bin/env python
# coding: utf-8

# # Remove families not in famsUnderSelection

# In[31]:


import os
from shutil import copyfile
import pandas as pd


# In[11]:


fams = pd.read_csv("final_results/famsUnderAbsrel.txt","\s+",index_col=False)


# In[14]:


fams.head()


# In[16]:


fams["numSitesUnderSelection"].value_counts().plot()


# In[17]:


families = fams["family"]


# In[20]:


families_in_dir = os.listdir("families_fubar")


# In[29]:


families_in_dir[2].split(".")[0]


# In[30]:


len(families_in_dir)


# In[ ]:


for i in range(0,len(families_in_dir)):
     if families_in_dir[i].split(".")[0] in list(families):
         copyfile("families_fubar/"+families_in_dir[i], "families_absrel/"+families_in_dir[i])

