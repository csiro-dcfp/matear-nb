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
# # Make BGC figures for the cafe60 reanalysis paper
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
sstd_d=sdic.sum(axis=(2,3)).load()
sstd_a=sadic.sum(axis=(2,3)).load()
gstd_d=gdic.sum(axis=(2,3)).load()
gstd_a=gadic.sum(axis=(2,3)).load()


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


# + Collapsed="false"
# Plot figure
def flxplt(a,b,da,db):
    a1=(a+da).rolling(time=12, center=True).mean()
    a2=(a-da).rolling(time=12, center=True).mean()
    b1=(b+db).rolling(time=12, center=True).mean()
    b2=(b-db).rolling(time=12, center=True).mean()
    plt.plot(time1,a.rolling(time=12, center=True).mean(),label='Natural')
    plt.plot(time1,b.rolling(time=12, center=True).mean(),label='Anthropogenic')
    plt.plot(time1,a*0,color='black')
    plt.fill_between(time1,a1,a2,color='cyan',alpha=.75)
    plt.fill_between(time1,b1,b2,color='orange',alpha=0.75)
#   plt.xlim(-15,15)
    plt.title(t1)
    plt.xlabel('Year')
    plt.ylabel('Sea to Air flux (Pg C y$^{-1}$)')
    plt.legend()
    
# Plot co2 fluxes
plt.figure(figsize=(10, 10))
plt.subplot(2,2,1)
# global
t1='Global'
flxplt(totd*(-1e-15),tota*(-1e-15),gstd_d.std(axis=1)*1e-15,gstd_a.std(axis=1)*1e-15)
plt.subplot(2,2,2)
t1='Southern Ocean'
flxplt(totsd*(-1e-15),totsa*(-1e-15),sstd_d.std(axis=1)*1e-15,sstd_a.std(axis=1)*1e-15)

plt.subplot(2,2,3)
plt.xlim(1985,2020)
# global
t1='Global'
flxplt(totd*(-1e-15),tota*(-1e-15),gstd_d.std(axis=1)*1e-15,gstd_a.std(axis=1)*1e-15)
plt.subplot(2,2,4)
plt.xlim(1985,2020)
t1='Southern Ocean'
flxplt(totsd*(-1e-15),totsa*(-1e-15),sstd_d.std(axis=1)*1e-15,sstd_a.std(axis=1)*1e-15)

plt.savefig('fco2.pdf',dpi=600)

# + [markdown] Collapsed="true"
# ## Flux Map

# + Collapsed="false"
#file4='/OSM/CBR/OA_DCFP/work/mat236/obs/taka_2007.nc'
file4='/OSM/CBR/OA_DCFP/work/mat236/obs/spco2_clim_1985-2015_MPI_SOM-FFN_v2016.nc'
dfco2=xr.open_dataset(file4)
oflx=dfco2.fgco2_clim.mean(axis=0).load()*12  # mol/m2/y to gC/m2/y

# + Collapsed="false"
# averaged flux 
aflx=dbgc.stf10.sel(time=slice('1985-01-16','2015-12-31')).mean(axis=(0,1)).load()*86400*365*12e-3 # mmol/m2/s to g C/m2/y


# + Collapsed="false"
def mflxplt(a,y,x):
    plt.xlabel('Longitude East')
    plt.ylabel('Latitude')
    plt.title(t1)
    cv=np.arange(-50,60,5)
    ff=plt.contourf(x,y,a,cmap='RdBu_r',levels=cv ,extend='both')
    cbar = plt.colorbar(ff)
    cbar.set_label('Sea to Air flux (g C m$^{-2}$ y$^{-1}$)')
    
plt.figure(figsize=(12,5 ))
plt.subplot(1,2,1)
t1='Observations'
mflxplt(oflx,oflx.lat,oflx.lon)

plt.subplot(1,2,2)
t1='CAFE60'
mflxplt(aflx*(-1),aflx.yt_ocean,aflx.xt_ocean)


# + [markdown] Collapsed="false"
# ## Sections

# + Collapsed="false"
dir='/home/mat236/dcfp/d60/bowen/data2/observations/GLODAPv2/mapped/GLODAPv2.2016b_MappedClimatologies/'
file1='GLODAPv2.2016b.PO4.nc' 
file2='GLODAPv2.2016b.TCO2.nc'
file3='GLODAPv2.2016b.Cant.nc'
dfco2=xr.open_dataset(dir+file2)
dpo4=xr.open_dataset(dir+file1)
daco2=xr.open_dataset(dir+file3)

otco2=dfco2.TCO2.sel(lon=180.5).load()*1.025
opo4 =dpo4.PO4.sel(lon=180.5).load()*1.025
oaco2 =daco2.Cant.sel(lon=180.5).load()*1.025

otco2.to_netcdf('/scratch1/mat236/work/otco2.nc')
opo4.to_netcdf('/scratch1/mat236/work/opo4.nc')
oaco2.to_netcdf('/scratch1/mat236/work/oaco2.nc')

# + Collapsed="false"
tmp=dbgc.adic.mean(axis=1)  # mean of the ensemble
adic=tmp.sel(time=slice('2002-01-16','2002-12-31'),xt_ocean=-179.5).mean(axis=0).load()
adic.to_netcdf('/scratch1/mat236/work/adic.nc')

tmp=dbgc.dic.mean(axis=1)  # mean of the ensemble
dic=tmp.sel(time=slice('2002-01-16','2002-12-31'),xt_ocean=-179.5).mean(axis=0).load()
dic.to_netcdf('/scratch1/mat236/work/dic.nc')

tmp=dbgc.no3.mean(axis=1)  # mean of the ensemble
po4=1/16.*tmp.sel(time=slice('2002-01-16','2002-12-31'),xt_ocean=-179.5).mean(axis=0).load()
po4.to_netcdf('/scratch1/mat236/work/po4.nc')


# + Collapsed="false"
def splt(a,y,x):
    plt.xlabel('Latitude')
    plt.ylabel('Depth')
    plt.title(t1)
    plt.gca().patch.set_color('.25')
 #   plt.ylim(5500,0)
    ff=plt.contourf(x,y,a,cmap='RdBu_r',levels=cv ,extend='both')
    ax = plt.gca()
    ax.invert_yaxis()
    cbar = plt.colorbar(ff)
    cbar.set_label('mmol m$^{-3}$ ')

plt.figure(figsize=(11,10 ))

cv=np.arange(1800,2500,50)
plt.subplot(2,1,1)
t1='Observations'
splt(otco2,dfco2.Depth,otco2.lat)  

plt.subplot(2,1,2)
t1='CAFE60'
#tco2=adic.sel(xt_ocean=-179.5).load()
splt(oaco2,adic.st_ocean,tco2.yt_ocean)  


# + Collapsed="false"
# nutrients
plt.figure(figsize=(12,15 ))

cv=np.arange(0,4,.25)
plt.subplot(2,1,1)
t1='Observations'
splt(opo4,dpo4.Depth,opo4.lat)  

plt.subplot(2,1,2)
t1='CAFE60'
splt(po4,po4.st_ocean,po4.yt_ocean)  




# + Collapsed="false"

