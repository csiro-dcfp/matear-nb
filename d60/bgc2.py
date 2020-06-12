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

# + [markdown] Collapsed="false"
# # Extract the BGC fields for the cafe60 for Tyler
#
# ## Southern Ocean data for SOTS
#
#

# + Collapsed="false"
import xarray as xr
#import netCDF4 as netCDF4
import numpy as np

import matplotlib.pyplot as plt
from scipy.io import netcdf
import scipy.stats as stats
import scipy.fftpack as fft
import scipy.signal as signal
import sys

#import  Cartopy
import cartopy.crs as ccrs
from cartopy import config


# + Collapsed="false"
def climatology(dsx,TIME1):
    clim = dsx.groupby(TIME1+'.month').mean(dim=TIME1)
    anom = dsx.groupby(TIME1+'.month') - clim
    season=dsx.groupby(TIME1+'.season').mean(dim=TIME1)
    return(clim,season,anom)


# + Collapsed="false"
#file='../big/data/csiro-dcfp-cafe-d60/ocean_month.zarr'
#file1='../big/data/csiro-dcfp-cafe-d60/ocean_ens_mean_at_analysis.zarr'
#file3='../big/data/csiro-dcfp-cafe-d60/atmos_isobaric_month.zarr'
file2='/OSM/CBR/OA_DCFP/data3/scratch_backup/06-Apr-2020/csiro-dcfp-cafe-d60/ocean_bgc_month.zarr'
file2='/OSM/CBR/OA_DCFP/data/model_output/CAFE/data_assimilation/d60-zarr/ocean_bgc_month.zarr'
#file1='../big/data/csiro-dcfp-cafe-d60/ocean_grid.nc'
file1='/home/mat236/area.nc'

dgrid=xr.open_dataset(file1)
dbgc = xr.open_zarr(file2,consolidated=True)


# + Collapsed="false"
#file3='../big/data/csiro-dcfp-cafe-d60/atmos_isobaric_month.zarr'
#datm = xr.open_zarr(file3,consolidated=True)
#datm.t_ref[0:10,0,45,74].plot()
def region(ds,lat1,lat2):
    dsr=ds.sel(lat=slice(lat1,lat2))
    return (dsr)


# + [markdown] Collapsed="false"
# ## extract the desired data from the zarr file

# + Collapsed="false"
# The minimum space bin I'd want to look at (where SOTS is) would be:
# Lat: -45:-50 | Lon: 140E: 145E.  
#The big box I've been looking at for regional context is:
# Lat: -35:-60 | Lon: 90E: 180E.  
dd=dbgc.sel(yt_ocean=slice(-60,-35),st_ocean=slice(0,500),xt_ocean=slice(90-360,180-360),
            time=slice('2000','2019'))

pclim,pseason,panom=climatology(dd,'time')

# + Collapsed="false"
pclim.phy.to_netcdf('/scratch1/mat236/tr/phy.nc')
pclim.zoo.to_netcdf('/scratch1/mat236/tr/zoo.nc')
pclim.no3.to_netcdf('/scratch1/mat236/tr/no3.nc')
pclim.pprod_gross.to_netcdf('/scratch1/mat236/tr/pprod_gross.nc')

# + Collapsed="false"

