# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# + Collapsed="false"
# %load_ext ferretmagic

# + Collapsed="false"
conda env list

# + Collapsed="false"
# %%ferret
use c2_sst.nc
sho data
sha sst[k=1,l=1]

# + Collapsed="false"
# %%ferret
set reg/y=-15:15/x=130:300
show data

# + Collapsed="false"
# %ferret_getdata sstdict = sst 

# + Collapsed="false"
# %%ferret?

# + Collapsed="false"
sstdict.keys()

# + Collapsed="false"
sst=sstdict['data']
sst.shape

# + Collapsed="false"

