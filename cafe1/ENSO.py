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

# + [markdown] Collapsed="false" slideshow={"slide_type": "notes"} toc-hr-collapsed=false
# # ENSO analysis of observation and models
#
# working ...

# + Collapsed="false" slideshow={"slide_type": "skip"}
# #%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import scipy.stats as stats
import scipy.fftpack as fft
import scipy.signal as signal
import sys
import xarray as xr
import netCDF4 as nc

import pandas as pd

#
#import  Cartopy
import cartopy.crs as ccrs
from cartopy import config

from eofs.xarray import Eof


# + Collapsed="false"
# useful procedures 

def climatology(dsx,TIME1):
    clim = dsx.groupby(TIME1+'.month').mean(dim=TIME1)
    anom = dsx.groupby(TIME1+'.month') - clim
    return(clim,anom)
    
def region(ds,lat1,lat2,lon1,lon2):
#    dsr=ds.sel(XT_OCEAN=slice(lon1,lon2),YT_OCEAN=slice(lat1,lat2))
    dsr=ds.sel(lon=slice(lon1,lon2),lat=slice(lat1,lat2))
    return (dsr)

def regresst(a,b):
    astd=a.std(dim=('time'))
    bstd=b.std(dim=('time'))
    adev=a-a.mean(dim=('time'))
    bdev=b-b.mean(dim=('time'))
    ab=adev*bdev
    abstd=ab.mean(dim='time')
    slope=abstd/(astd*astd)
    intercep=b.mean(dim=('time')) - slope*a.mean(dim=('time'))
    bhat = slope * a + intercep 
#r =abstd/astd/bstd
#slope= abstd*abstd/astd/astd
#slope per std of a = r * bstd
# intercept =
# qhat = slope*a + b.mean - slope*a.mean
#
    return astd,bstd,abstd,bhat


# + [markdown] Collapsed="false" slideshow={"slide_type": "notes"} toc-hr-collapsed=false
# ## Observations

# + Collapsed="false"
# for observations

dobs = xr.open_mfdataset('HadISST.nc', decode_times=False)
units, reference_date = dobs.time.attrs['units'].split('since')
reference_date = ' 1870-1-1 00:00:00 0:0'
print(units, reference_date)
print(dobs.sizes['time'])

dobs['time'] = pd.date_range(start=reference_date, periods=dobs.sizes['time'], freq='MS')
obs_clim,obs_anom=climatology(dobs,'time')

obs_mean=obs_clim.sst.mean(dim='month')

# + [markdown] Collapsed="false"
# ### detrend

# + Collapsed="false"
# get time and put it into a variable of float as days since initial time
p1=obs_anom.sst.coords['time']
p=(p1 - np.datetime64('1800-01-01 00:00:00')).astype(int) / (86400*1e9)

astd,bstd,abstd,bhat=regresst(p,obs_anom)
# correlation
r=abstd/(astd*bstd)
# SSTanom slope normalised 
slope=r*bstd
slope=abstd/(astd*astd)

r1=slope.sst
print(r1[90,180].values)

obs_dtrend=obs_anom.sst - bhat.sst

#plt.plot(dtrend[:,90,180])
#plt.plot(bhat.sst[90,180,:])

# + [markdown] Collapsed="false"
# ### ENSO power spectrum

# + Collapsed="false"
# nino34 
obs_nino34reg=region(obs_dtrend,-5,5,360-170,360-120)
obs_nino34 = obs_nino34reg.mean(axis=(1,2)).load()


# + Collapsed="false"
# compute power spectra
def power_spec(dt,ns):
    N=ns.size
    f,t,ss=signal.spectrogram(ns,fs=dt,scaling='density')  
    f1,s1=signal.welch(ns, fs=dt, scaling='density')
    return f1,s1,f,t,ss


obs_f1,obs_s1,f,t,ss=power_spec(1,obs_nino34)
obs_sd=ss.std(axis=1)
np1=26
#plt.plot(f[0:np1],ss[0:np1].mean(axis=1))
plt.plot(obs_f1[0:np1],obs_s1[0:np1])
sd=obs_sd[0:np1]
plt.plot(obs_f1[0:np1],obs_s1[0:np1]+sd, obs_f1[0:np1],obs_s1[0:np1]-sd)


# + [markdown] Collapsed="false"
# ### ENSO composites

# + Collapsed="false"
# Select DJF period and compute composite
def enso_composite(nino34,dtrend):
# mask other months with nan
    ds_DJF = nino34.where(nino34['time.season'] == 'DJF')
# rolling mean -> only Jan is not nan
# however, we loose Jan/ Feb in the first year and Dec in the last
    ds_DJF = ds_DJF.rolling(min_periods=3, center=True, time=3).mean()
# make annual mean DJF index
    ds_DJF = ds_DJF.groupby('time.year').mean('time')
# compute map of DJF mean
    map_DJF= dtrend.where(dtrend['time.season'] == 'DJF')
    map_DJF = map_DJF.rolling(min_periods=3, center=True, time=3).mean()
    map_DJF = map_DJF.groupby('time.year').mean('time')
# compute composites
    elnino=map_DJF.where(ds_DJF > 1).mean('year')
    lanina=map_DJF.where(ds_DJF < -1).mean('year')
    return elnino, lanina

obs_el,obs_la= enso_composite(obs_nino34,obs_dtrend)

(obs_el-obs_la).plot()


# + Collapsed="false"
# save output in Dataset
def mk_dataset(str,f1,s1,el,la,nino34):
    dtmp1= xr.DataArray(s1, dims=(str+'f1'), coords={str+'f1': f1})
    dtmp = dtmp1.to_dataset(name = str+'s1')
    dtmp[str+'el']=el
    dtmp[str+'la']=la
    dtmp[str+'nino34']=nino34
    return dtmp

obsc=mk_dataset('obs_',obs_f1,obs_s1,obs_el,obs_la,obs_nino34)
obsc.to_netcdf('obs_enso.nc')

# + [markdown] Collapsed="false" slideshow={"slide_type": "notes"} toc-hr-collapsed=false
# ## CAFE

# + Collapsed="false"
# for CAFE (use the time axis from access for getting a pandas compatible version for detrending)
daccess = xr.open_mfdataset('../cmip6/tos_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc')
dcafe = xr.open_mfdataset('c2_sst.nc')
dcafe = dcafe.rename({'XT_OCEAN':'lon','YT_OCEAN':'lat','TIME1':'time'})
cafe_clim,cafe_anom= climatology(dcafe,'time')

cafe_mean=cafe_clim.SST[:,0,:,:].mean(dim='month')
cafe_anom1=cafe_anom.SST[:,0,:,:]



# + [markdown] Collapsed="false"
# ### detrend

# + Collapsed="false"
# get time and put it into a variable of float as days since initial time
anom=cafe_anom1.assign_coords(time=daccess.time[0:1200].values)
p1=anom.coords['time']
p=(p1 - np.datetime64('1800-01-01 00:00:00')).astype(int) / (86400*1e9)

astd,bstd,abstd,bhat=regresst(p,anom)
# correlation
r=abstd/(astd*bstd)
# SSTanom slope normalised 
slope=r*bstd
slope=abstd/(astd*astd)

r1=slope
print(r1[90,180].values)

cafe_dtrend=anom - bhat

#plt.plot(cafe_dtrend[:,90,180])
#plt.plot(bhat[90,180,:])

# + [markdown] Collapsed="false"
# ### ENSO power spectrum

# + Collapsed="false"
# nino34 
cafe_nino34reg=region(cafe_dtrend,-5,5,-170,-120)
cafe_nino34 = cafe_nino34reg.mean(axis=(1,2)).load()


# + Collapsed="false"
# compute power spectra
def power_spec(dt,ns):
    N=ns.size
    f,t,ss=signal.spectrogram(ns,fs=dt,scaling='density')  
    f1,s1=signal.welch(ns, fs=dt, scaling='density')
    return f1,s1,f,t,ss


cafe_f1,cafe_s1,f,t,ss=power_spec(1,cafe_nino34)
cafe_sd=ss.std(axis=1)
np1=26
#plt.plot(f[0:np1],ss[0:np1].mean(axis=1))
plt.plot(cafe_f1[0:np1],cafe_s1[0:np1])
sd=cafe_sd[0:np1]
plt.plot(cafe_f1[0:np1],cafe_s1[0:np1]+sd, cafe_f1[0:np1],cafe_s1[0:np1]-sd)


# + [markdown] Collapsed="false"
# ### ENSO composites

# + Collapsed="false"
# Select DJF period and compute composite
def enso_composite(nino34,dtrend):
# mask other months with nan
    ds_DJF = nino34.where(nino34['time.season'] == 'DJF')
# rolling mean -> only Jan is not nan
# however, we loose Jan/ Feb in the first year and Dec in the last
    ds_DJF = ds_DJF.rolling(min_periods=3, center=True, time=3).mean()
# make annual mean DJF index
    ds_DJF = ds_DJF.groupby('time.year').mean('time')
# compute map of DJF mean
    map_DJF= dtrend.where(dtrend['time.season'] == 'DJF')
    map_DJF = map_DJF.rolling(min_periods=3, center=True, time=3).mean()
    map_DJF = map_DJF.groupby('time.year').mean('time')
# compute composites
    elnino=map_DJF.where(ds_DJF > 1).mean('year')
    lanina=map_DJF.where(ds_DJF < -1).mean('year')
    return elnino, lanina

cafe_el,cafe_la= enso_composite(cafe_nino34,cafe_dtrend)

(cafe_el-cafe_la).plot()


# + Collapsed="false"
# save output in Dataset
cafec=mk_dataset('cafe_',cafe_f1,cafe_s1,cafe_el,cafe_la,cafe_nino34)
cafec.to_netcdf('cafe_enso.nc')

# + [markdown] Collapsed="false"
# ## ACCESS-ESM

# + Collapsed="false"
# for ACCESS-ESM (use the cafe infor for lat and lon axis)
dcafe = xr.open_mfdataset('c2_sst.nc')
dcafe = dcafe.rename({'XT_OCEAN':'lon','YT_OCEAN':'lat','TIME1':'time'})

daccess = xr.open_mfdataset('../cmip6/tos_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc')
access_clim,access_anom= climatology(daccess,'time')
access_mean=access_clim.tos.mean(dim='month')


# + [markdown] Collapsed="false"
# ### detrend

# + Collapsed="false"

# get time and put it into a variable of float as days since initial time
p1=access_anom.tos.coords['time']
p=(p1 - np.datetime64('1800-01-01 00:00:00')).astype(int) / (86400*1e9)

astd,bstd,abstd,bhat=regresst(p,access_anom)
# correlation
r=abstd/(astd*bstd)
# SSTanom slope normalised 
slope=r*bstd
slope=abstd/(astd*astd)

r1=slope.tos
print(r1[90,180].values)

access_dtrend=access_anom.tos - bhat.tos

plt.plot(access_dtrend[:,90,180])
plt.plot(bhat.tos[90,180,:])

# + [markdown] Collapsed="false"
# ### ENSO power spectrum

# + Collapsed="false"
# change access grid values
dtrend_cafe_gd=access_dtrend.copy()
dtrend_cafe_gd=dtrend_cafe_gd.assign_coords(j=dcafe.lat.values,i=dcafe.lon.values)
dtrend_cafe_gd=dtrend_cafe_gd.rename({'i':'lon','j':'lat',})
# nino34 

access_nino34reg=region(dtrend_cafe_gd,-5,5,-170,-120)
access_nino34 = access_nino34reg.mean(axis=(1,2)).load()


# + Collapsed="false"
# compute power spectra
def power_spec(dt,ns):
    N=ns.size
    f,t,ss=signal.spectrogram(ns,fs=dt,scaling='density')  
    f1,s1=signal.welch(ns, fs=dt, scaling='density')
    return f1,s1,f,t,ss


access_f1,access_s1,f,t,ss=power_spec(1,access_nino34)
access_sd=ss.std(axis=1)
np1=26
#plt.plot(f[0:np1],ss[0:np1].mean(axis=1))
plt.plot(access_f1[0:np1],access_s1[0:np1])
sd=access_sd[0:np1]
plt.plot(access_f1[0:np1],access_s1[0:np1]+sd, access_f1[0:np1],access_s1[0:np1]-sd)


# + [markdown] Collapsed="false"
# ### ENSO composites

# + Collapsed="false"
# Select DJF period and compute composite
def enso_composite(nino34,dtrend):
# mask other months with nan
    ds_DJF = nino34.where(nino34['time.season'] == 'DJF')
# rolling mean -> only Jan is not nan
# however, we loose Jan/ Feb in the first year and Dec in the last
    ds_DJF = ds_DJF.rolling(min_periods=3, center=True, time=3).mean()
# make annual mean DJF index
    ds_DJF = ds_DJF.groupby('time.year').mean('time')
# compute map of DJF mean
    map_DJF= dtrend.where(dtrend['time.season'] == 'DJF')
    map_DJF = map_DJF.rolling(min_periods=3, center=True, time=3).mean()
    map_DJF = map_DJF.groupby('time.year').mean('time')
# compute composites
    elnino=map_DJF.where(ds_DJF > 1).mean('year')
    lanina=map_DJF.where(ds_DJF < -1).mean('year')
    return elnino, lanina

access_el,access_la= enso_composite(access_nino34,dtrend_cafe_gd)

(access_el-access_la).plot()


# + [markdown] Collapsed="false"
# ### save key variables
#

# + Collapsed="false"
# save output in Dataset
accessc=mk_dataset('access_',access_f1,access_s1,access_el,access_la,access_nino34)
accessc.to_netcdf('access_enso.nc')


# + [markdown] Collapsed="false"
# # Plot all fields

# + Collapsed="false"
# Generic Plot Map script using AxesGrid
def map3(z1, z2, z3, tit1,tit2,tit3):
    import cartopy
    from mpl_toolkits.axes_grid1 import AxesGrid
    from cartopy.mpl.geoaxes import GeoAxes
    proj=ccrs.Mollweide()
    proj=ccrs.Robinson(central_longitude=180)
    dproj=ccrs.PlateCarree()
    cv=np.arange(-4,5,1)
    vm=4
    
    axes_class = (GeoAxes, dict(map_projection=proj))
    axgr = AxesGrid(fig, 111, nrows_ncols=(3, 1), axes_pad=0.25, cbar_mode='none', 
               cbar_location='bottom', cbar_pad=0.2, axes_class=axes_class,
               cbar_size='5%', label_mode='')  # note the empty label_mode

#
    for i, ax in enumerate(axgr):
#    ax = plt.axes(projection=ccrs.Mollweide())
#    ax = plt.axes(projection=proj)
        ax.coastlines(resolution='110m')
        ax.gridlines()
        ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black',color='white')
        
        if (i == 0):
            lons=z1.lon
            lats=z1.lat
            p=ax.contourf(lons, lats, z1,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
            ax.set_title(tit1,loc='left')
            
        if (i == 1):
            lons=z2.lon
            lats=z2.lat
            p=ax.contourf(lons, lats, z2,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
            ax.set_title(tit2, loc= 'left')
            
        if (i == 2):
            lons=z3.lon
            lats=z3.lat
            p=ax.contourf(lons, lats, z3,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
            ax.set_title(tit3, loc= 'left')
    
    plt.colorbar(p, shrink=.8, orientation='horizontal', label='SST anomaly (C)')
    plt.title('Bias')
    
    axgr[0].cax.colorbar(p)
#    axis.label.set_text("SST Difference")

# Generic Plot Map script using subplot
def map3o(z1, z2, z3,tit1,tit2,tit3):
    import cartopy
    from mpl_toolkits.axes_grid1 import AxesGrid
    from cartopy.mpl.geoaxes import GeoAxes
    proj=ccrs.Mollweide()
    proj=ccrs.Robinson(central_longitude=180)
    dproj=ccrs.PlateCarree()
    cv=np.arange(-4,5,1)
    vm=4
    
#
    for  i in range(3):
#    ax = plt.axes(projection=ccrs.Mollweide())
#    ax = plt.axes(projection=proj)
        ax = plt.subplot(2,2,i+2,projection=proj)
        ax.coastlines(resolution='110m')
        ax.gridlines()
        ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black',color='white')
        if (i == 0):
            lons=z1.lon
            lats=z1.lat
            p=ax.contourf(lons, lats, z1,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
            plt.title(tit1)
            
        if (i == 1):
            lons=z2.lon
            lats=z2.lat
            p=ax.contourf(lons, lats, z2,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
            plt.title(tit2)
            
        if (i == 2):
            lons=z3.lon
            lats=z3.lat
            p=ax.contourf(lons, lats, z3,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
            plt.title(tit3)
            
            
#    plt.colorbar(p)
    plt.colorbar(p, shrink=.8, orientation='horizontal', label='SST difference')

    return
    

# + Collapsed="false"

fig=plt.figure(figsize=(10,16))
obs=obs_el-obs_la
cafe=cafe_el-cafe_la
access=access_el - access_la
#map3(obs,cafe,access,'Observations','CAFE','ACCESS')


# + Collapsed="false"

# power spectrum
fig=plt.figure(figsize=(16,12))


ax=plt.subplot(2,2,1)
plt.plot(obs_f1[0:np1],obs_s1[0:np1])
plt.plot(access_f1[0:np1],access_s1[0:np1])
plt.plot(cafe_f1[0:np1],cafe_s1[0:np1])

plt.xlabel('Frequency (cycles per month)')
plt.ylabel('Power')
plt.title('a)',loc='left')

map3o(obs,cafe,access,'Observations','CAFE','ACCESS')

plt.savefig('enso.pdf',dpi=600)


# + Collapsed="false"

