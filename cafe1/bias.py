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

# + [markdown] Collapsed="false" toc-hr-collapsed=false
# # Bias calculation
#

# + [markdown] Collapsed="false"
# ## setup the modules  

# + Collapsed="false"
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
    dsr=ds.sel(lon=slice(lon1,lon2),lat=slice(lat1,lat2))
    return (dsr)

def regresst(a,b):
    astd=a.std(dim=('time'))
    bstd=b.std(dim=('time'))
    adev=a-a.mean(dim=('time'))
    bdev=b-b.mean(dim=('time'))
    ab=adev*bdev
    abstd=ab.mean(dim='time')
    slope=abstd*abstd/astd/astd
    intercep=b.mean(dim=('time')) - slope*a.mean(dim=('time'))
    bhat = slope * a + intercep
#r =abstd/astd/bstd
#slope= abstd*abstd/astd/astd
#slope per std of a = r * bstd
# intercept =
# qhat = slope*a + b.mean - slope*a.mean
#
    return astd,bstd,abstd,bhat

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x,0,0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)



# + Collapsed="false"
def power1(anti,sym,len):
# do 2d fft calculation but split time into sub bits
    da=ar_reshape(anti,len)
    ds=ar_reshape(sym,len)
    Nt=da.shape[1]
    Nx=anti.lon.size
    print(Nt,Nx)
# get frequency
    xfw=fft.fftfreq(Nt,1)
    xfx=fft.fftfreq(Nx,1) *Nx  # puts it into wave number
#
    amp_anti= fft.fft2(da,axes=[1,3])
    amp_sym= fft.fft2(ds,axes=[1,3])
    xfw=fft.fftshift(xfw)
    xfx=fft.fftshift(xfx)
    amp_anti=fft.fftshift(amp_anti)
    amp_sym=fft.fftshift(amp_sym)
#
    return(xfw,xfx,amp_anti,amp_sym)

def power(anti,sym):
    Nt=anti.time.size
    Nx=anti.lon.size
# get frequency
    xfw=fft.fftfreq(Nt,1)
    xfx=fft.fftfreq(Nx,1) *Nx  # puts it into wave number
#
    amp_anti= fft.fft2(anti,axes=[0,2])
    atemp=abs(amp_anti)**2
    apower_anti= atemp.mean(axis=1)
#
    amp_sym= fft.fft2(sym,axes=[0,2])
    atemp=abs(amp_sym)**2
    apower_sym= atemp.mean(axis=1)
#
    return(xfw,xfx,apower_anti,apower_sym,amp_anti,amp_sym)

# separate signal into symmetric and anti-symmetric parts for both hemispheres
def split1(xr):
    north=xr.sel(lat=slice(-15,15))
    south=xr.sel(lat=slice(-15,15))
#south=south.reindex(lat=south.lat[::-1]) ## not needed
    south1=south.assign_coords(lat=south.lat*(-1))
    print(south1.lat.values)
    print(north.lat.values)
#
    anti = (north - south1)/2
    sym = (north+south1)/2
#
    total=anti+sym
    check=north- total
    print(check.max(axis=(0,1,2)).values)
    return(anti,sym,total)


def ar_reshape(ar,len):
    Nt=anti.time.size
    Nrep=np.int(Nt/len)
    Nlon=anti.lon.size
    Nlat=anti.lat.size
    print(Nrep,len)
    ain=np.ma.masked_invalid(ar)
    aa=np.ma.filled(ain,0)
#    aa=np.copy(ar)
    bb=aa[0:Nrep*len,:,:]
    cc=bb.reshape(Nrep,len,Nlat,Nlon)
    dd=cc*0
    from scipy.signal import blackman
    window=blackman(len)
#    window=np.ones(len)
    from itertools import product
    for l,i,j in product(range(Nrep), range(0, 14), range(0,144) ):  
        dd[l,0:len,i,j]= cc[l,0:len,i,j]*window[0:len]
#
    return(dd)



# + [markdown] Collapsed="false"
# ## Obs

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
# ## ACCESS-ESM

# + Collapsed="false"
# for ACCESS-ESM
daccess = xr.open_mfdataset('../cmip6/tos_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc')
access_clim,access_anom= climatology(daccess,'time')
access_mean=access_clim.tos.mean(dim='month')


# + Collapsed="false"


# + [markdown] Collapsed="false"
# ## CAFE

# + Collapsed="false"
# for CAFE
dcafe = xr.open_mfdataset('c2_sst.nc')
dcafe = dcafe.rename({'XT_OCEAN':'lon','YT_OCEAN':'lat','TIME1':'time'})
cafe_clim,cafe_anom= climatology(dcafe,'time')

cafe_mean=cafe_clim.SST[:,0,:,:].mean(dim='month')


# + Collapsed="false"
obs_mean.to_netcdf('sst_obs.nc')
access_mean.to_netcdf('access_sst.nc')
cafe_mean.to_netcdf('cafe_sst.nc')

# + [markdown] Collapsed="false"
# ## Final figures with all output

# + Collapsed="false"
# put all model output on the obs grid
lons=obs_mean.lon
lats=obs_mean.lat
sst_obs=obs_mean

ll=np.copy(cafe_mean.lon)
cafe_mean=cafe_mean.assign_coords(lon=np.where(ll<0,ll+360,ll))

# change access grid values
cafe_grid=access_mean.copy()
cafe_grid=cafe_grid.assign_coords(j=cafe_mean.lat.values,i=cafe_mean.lon.values)

# interpolate on to obs grid
sst_cafe=cafe_mean.interp(lon=lons, lat=lats)
sst_access=cafe_grid.interp(i=lons, j=lats)

dsst_cafe=sst_cafe - sst_obs
dsst_access=sst_access - sst_obs

# + Collapsed="false"


# + Collapsed="false"
# Generic Plot Map script
def map1(z1,lons,lats):
    import cartopy
    from mpl_toolkits.axes_grid1 import AxesGrid
    proj=ccrs.Mollweide()
    dproj=ccrs.PlateCarree()
    cv=np.arange(-4,5,1)
    vm=4
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.coastlines(resolution='110m')
    ax.gridlines()
    plt.contourf(lons, lats, z1,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
#    plt.pcolormesh(lons, lats, z1,vmin=-1*vm,vmax=vm, cmap='RdBu_r',transform=dproj)
#plt.contourf(lons, lats, sst, levels=cv, extend='both', cmap='RdBu', transform=dproj)
    plt.colorbar()
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black',color='white')
    return


# Generic Plot Map script using subplot
def map2o(z1, z2, lons,lats,tit1,tit2):
    import cartopy
    from mpl_toolkits.axes_grid1 import AxesGrid
    from cartopy.mpl.geoaxes import GeoAxes
    proj=ccrs.Mollweide()
    proj=ccrs.Robinson(central_longitude=180)
    dproj=ccrs.PlateCarree()
    cv=np.arange(-4,5,1)
    vm=4
    
    fig = plt.figure(figsize=(10, 16))
    ax = plt.subplot(2,1,1,projection=proj)
#
    for  i in range(2):
#    ax = plt.axes(projection=ccrs.Mollweide())
#    ax = plt.axes(projection=proj)
        ax = plt.subplot(2,1,i+1,projection=proj)
        ax.coastlines(resolution='110m')
        ax.gridlines()
        ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black',color='white')
        if (i == 0):
            p=ax.contourf(lons, lats, z1,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
            plt.title(tit1)
            
        if (i == 1):
            p=ax.contourf(lons, lats, z2,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
            plt.title(tit2)
            
#    plt.colorbar(p)
    plt.colorbar(p, shrink=.8, orientation='horizontal', label='SST difference')

    return
    

# Generic Plot Map script using AxesGrid
def map2(z1, z2, lons,lats,tit1,tit2):
    import cartopy
    from mpl_toolkits.axes_grid1 import AxesGrid
    from cartopy.mpl.geoaxes import GeoAxes
    proj=ccrs.Mollweide()
    proj=ccrs.Robinson(central_longitude=180)
    dproj=ccrs.PlateCarree()
    cv=np.arange(-4,5,1)
    vm=4
    
    fig = plt.figure(figsize=(10, 16))
    axes_class = (GeoAxes, dict(map_projection=proj))
    axgr = AxesGrid(fig, 111, nrows_ncols=(2, 1), axes_pad=0.25, cbar_mode='none', 
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
            p=ax.contourf(lons, lats, z1,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
            ax.set_title(tit1,loc='left')
            
        if (i == 1):
            p=ax.contourf(lons, lats, z2,levels=cv, cmap='RdBu_r',transform=dproj,extend='both')
            ax.set_title(tit2, loc= 'left')
    
    plt.colorbar(p, shrink=.8, orientation='horizontal', label='SST difference')
    plt.title('Bias')
    
    axgr[0].cax.colorbar(p)
#    axis.label.set_text("SST Difference")


    return
    


# + [markdown] Collapsed="false" toc-hr-collapsed=true
# ## extra

# + Collapsed="false"
map2(dsst_cafe,dsst_access,lons,lats,'CAFE','ACCESS')
plt.savefig('sst_bias.pdf',dpi=600)


# + Collapsed="false"
map2o(dsst_cafe,dsst_access,lons,lats,'CAFE','ACCESS')

# + Collapsed="false"

