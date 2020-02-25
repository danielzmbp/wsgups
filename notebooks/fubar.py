#!/usr/bin/env python
# coding: utf-8

# # Remove families not in famsUnderSelection

# In[1]:


import os
from shutil import copyfile
import pandas as pd


# In[2]:


fams = pd.read_csv("final_results/famsUnderSelection.txt","\s+",index_col=False)


# In[3]:


fams.head()


# In[4]:


fams["numSitesUnderSelection"].value_counts().plot()


# In[5]:


families = fams["family"]


# In[6]:


families_in_dir = os.listdir("families")


# In[7]:


families_in_dir[2].split(".")[0]


# In[8]:


len(families_in_dir)


# In[9]:


for i in range(0,len(families_in_dir)):
     if families_in_dir[i].split(".")[0] in list(families):
         copyfile("families/"+families_in_dir[i], "families_fubar/"+families_in_dir[i])


# In[ ]:




