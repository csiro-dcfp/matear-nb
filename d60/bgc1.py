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
# # Extract the BGC fields for the cafe60 reanalysis paper
#
# ## sections and integrate flux

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
file3='../big/data/csiro-dcfp-cafe-d60/atmos_isobaric_month.zarr'
file2='/OSM/CBR/OA_DCFP/data3/scratch_backup/06-Apr-2020/csiro-dcfp-cafe-d60/ocean_bgc_month.zarr'
file1='../big/data/csiro-dcfp-cafe-d60/ocean_grid.nc'
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
# ## Fluxes

# + Collapsed="false"
stf10_mean = dbgc.stf10.mean(axis=1)
stf10_std = dbgc.stf10.std(axis=1)
stf07_mean = dbgc.stf07.mean(axis=1)
stf07_std = dbgc.stf07.std(axis=1)

# + Collapsed="false"
# southern ocean south of 35s
sdic=stf07_mean.sel(yt_ocean=slice(-90,-35))
sadic=stf10_mean.sel(yt_ocean=slice(-90,-35))
ssum=sdic*dgrid.area_t*86400*365*12E-3
totsd=ssum.sum(axis=(1,2)).load()
ssum=sadic*dgrid.area_t*86400*365*12E-3
totsa=ssum.sum(axis=(1,2)).load()

# global ocean
sdic= stf07_mean *dgrid.area_t*86400*365*12E-3
sadic= stf10_mean *dgrid.area_t*86400*365*12E-3
totd=sdic.sum(axis=(1,2)).load()
tota=sadic.sum(axis=(1,2)).load()

# + Collapsed="false"
# Ensemble estimate of the flux from which the standard deviation can be calculated
gdic= dbgc.stf07*dgrid.area_t*86400*365*12E-3
gadic= dbgc.stf10*dgrid.area_t*86400*365*12E-3
sdic=dbgc.stf07.sel(yt_ocean=slice(-90,-35))*dgrid.area_t*86400*365*12E-3
sadic=dbgc.stf10.sel(yt_ocean=slice(-90,-35))*dgrid.area_t*86400*365*12E-3
sstd_d=sdic.sum(axis=(2,3))
sstd_a=sadic.sum(axis=(2,3))
gstd_d=gdic.sum(axis=(2,3))
gstd_a=gadic.sum(axis=(2,3))


# + Collapsed="false"
time1=np.arange(1960.0833,2019.5,.08333)

# save some output
ddout='/scratch1/mat236/work/'
totd.to_netcdf(ddout+'totd.nc')
tota.to_netcdf(ddout+'tota.nc')
gstd_d.to_netcdf(ddout+'gstd_d.nc')
gstd_a.to_netcdf(ddout+'gstd_a.nc')

totsd.to_netcdf(ddout+'totsd.nc')
totsa.to_netcdf(ddout+'totsa.nc')
sstd_d.to_netcdf(ddout+'sstd_d.nc')
sstd_a.to_netcdf(ddout+'sstd_a.nc')

# + [markdown] Collapsed="false"
# ## Flux Map

# + Collapsed="false"
#file4='/OSM/CBR/OA_DCFP/work/mat236/obs/taka_2007.nc'
file4='/OSM/CBR/OA_DCFP/work/mat236/obs/spco2_clim_1985-2015_MPI_SOM-FFN_v2016.nc'
dfco2=xr.open_dataset(file4)
oflx=dfco2.fgco2_clim.mean(axis=0).load()*12  # mol/m2/y to gC/m2/y


# + Collapsed="false"
# save some output
ddout='/scratch1/mat236/work/'
oflx.to_netcdf(ddout+'oflx.nc')

# + Collapsed="false"
# averaged flux 
aflx=dbgc.stf10.sel(time=slice('1985-01-16','2015-12-31')).mean(
    axis=(0,1)).load()*86400*365*12e-3 # mmol/m2/s to g C/m2/y

aflx.to_netcdf(ddout+'aflx.nc')

# + [markdown] Collapsed="false"
# ## Sections

# + Collapsed="false"
dir='bowen/data2/observations/GLODAPv2/mapped/GLODAPv2.2016b_MappedClimatologies/'
file1='GLODAPv2.2016b.PO4.nc' 
file2='GLODAPv2.2016b.TCO2.nc'
file3='GLODAPv2.2016b.Cant.nc'
dfco2=xr.open_dataset(dir+file2)
dpo4=xr.open_dataset(dir+file1)
daco2=xr.open_dataset(dir+file3)

Depth=dfco2.Depth
otco2=dfco2.TCO2.sel(lon=180.5)*1.025
opo4 =dpo4.PO4.sel(lon=180.5)*1.025
oaco2 =daco2.Cant.sel(lon=180.5)*1.025


# + Collapsed="false"
Depth.to_netcdf(ddout+'odepth.nc')
otco2.to_netcdf(ddout+'otco2.nc')
opo4.to_netcdf(ddout+'opo4.nc')
oaco2.to_netcdf(ddout+'oaco2.nc')

# + Collapsed="false"
tmp=dbgc.adic.mean(axis=1)  # mean of the ensemble
adic=tmp.sel(time=slice('2002-01-16','2002-12-31'),
             xt_ocean=-179.5).mean(axis=0)
adic.to_netcdf(ddout+'adic.nc')

tmp=dbgc.dic.mean(axis=1)  # mean of the ensemble
dic=tmp.sel(time=slice('2002-01-16','2002-12-31'),
            xt_ocean=-179.5).mean(axis=0)
dic.to_netcdf(ddout+'dic.nc')

tmp=dbgc.no3.mean(axis=1)  # mean of the ensemble
po4=1/16.*tmp.sel(time=slice('2002-01-16','2002-12-31'),
                  xt_ocean=-179.5).mean(axis=0)
po4.to_netcdf(ddout+'po4.nc')

# + [markdown] Collapsed="false"
# # Finished !!
### rjm 
