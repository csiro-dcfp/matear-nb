#!/usr/bin/env python
# coding: utf-8

# # Make coherence plot combining jra and cafe control run analysis on mush with CAFE60 data

# In[1]:



#get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
import scipy.stats as stats
import scipy.fftpack as fft
import scipy.signal as signal
import sys
import xarray as xr
import netCDF4 as nc

#
#import  Cartopy
import cartopy.crs as ccrs
from cartopy import config

#from eofs.xarray import Eof


# In[2]:


def climatology(dsx,TIME1):
    clim = dsx.groupby(TIME1+'.month').mean(dim=TIME1)
    anom = dsx.groupby(TIME1+'.month') - clim
    season=dsx.groupby(TIME1+'.season').mean(dim=TIME1)
    return(clim,season,anom)


# In[3]:


#import h5netcdf
import pandas as pd
import xarray as xr
import xarray.ufuncs as xu
import matplotlib as mpl
import numpy as np
import scipy
import datetime
import matplotlib.pyplot as plt

#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
#import cartopy
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import dask as dask
#import seaborn as sns
import os
import re

from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab

#import proplot as plot
#import cmocean
#import cmocean.cm as cmo

import cftime as cftime


# In[4]:


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
    xlat=15
    north=xr.sel(lat=slice(xlat*-1,xlat))
    south=xr.sel(lat=slice(xlat*-1,xlat))
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
    print(ar.shape)
    Nt=ar.time.size
    Nrep=np.int(Nt/len)
    Nlon=ar.lon.size
    Nlat=ar.lat.size
    print(Nrep,len,Nlat,Nlon)
    ain=np.ma.masked_invalid(ar)
    aa=np.ma.filled(ain,0)
#    aa=np.copy(ar)
    print(aa.shape)
    bb=aa[0:Nrep*len,:,:]
    cc=bb.reshape(Nrep,len,Nlat,Nlon)
    dd=cc*0
    from scipy.signal import blackman
    window=blackman(len)
#    window=np.ones(len)
    from itertools import product
    for l,i,j in product(range(Nrep), range(0, Nlat), range(0,Nlon) ):  
        dd[l,0:len,i,j]= cc[l,0:len,i,j]*window[0:len]
#
    return(dd)

# detrend xarray 
def dtrend_t(anom):
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
    dtrend=anom - bhat
    return dtrend,bhat


# In[5]:


# coherence calculation of anti and symmetric parts of 3d fields
def coherence(twindow,dat1,dat2):
    anti,sym,total=split1(dat1)
    anti_u,sym_u,total_u=split1(dat2)
    
# spectra
    xfw,xfx,ampa,amps = power1(anti,sym,twindow)
    powera=(ampa*ampa.conj()).mean(axis=(0,2))
    powers=(amps*amps.conj()).mean(axis=(0,2))

    xfw_u,xfx_u,ampa_u,amps_u = power1(anti_u,sym_u,twindow)
    powera_u=(ampa_u*ampa_u.conj()).mean(axis=(0,2))
    powers_u=(amps_u*amps_u.conj()).mean(axis=(0,2))
# now the cross signal
    anti_cross=(ampa_u*ampa.conj()).mean(axis=(0,2))
    sym_cross= (amps_u*amps.conj()).mean(axis=(0,2))

# coherence
    coh_anti = abs(anti_cross)**2 / (powera_u*powera)
    coh_sym = abs(sym_cross)**2 / (powers_u*powers)

    print(coh_anti.max(axis= (0,1)))
    print(coh_anti.min(axis= (0,1)))

    return coh_anti,coh_sym,xfw,xfx


# In[6]:


#xarray option
#xr.set_options(display_style="html")

tmpdir = '/scratch1/mat236/work/'

dmush = xr.open_mfdataset(tmpdir+'mush_coherence.nc',parallel=True,combine='by_coords')
dcafe60 = xr.open_mfdataset(tmpdir+'1cafe60_coherence.nc',parallel=True,combine='by_coords')


# In[13]:


def anti_plot(cv,xfx,xfw,anti,tit1):
    plt.xlim(-15,15)
    plt.title(tit1+':  Anti')
    plt.xlabel('Zonal Wavenumber')
    plt.ylabel('Time Frequency (day$^{-1}$)')
    ff=plt.contourf(xfx,xfw,(abs(anti)),cmap='PuRd',levels=cv ,extend='max')
    cbar = plt.colorbar(ff)
    return
                    
def sym_plot(cv,xfx,xfw,sym):
    plt.xlim(-15,15)
    plt.xlabel('Zonal Wavenumber')
    plt.ylabel('Time Frequency (day$^{-1}$)')
    plt.title('Symmetric')
    ff=plt.contourf(xfx,xfw,(abs(sym)),cmap='PuRd',levels=cv ,extend='max')
    cbar= plt.colorbar(ff)
    return

def sym_plot1(cv,xfx,xfw,sym):
    plt.xlim(-15,15)
    plt.xlabel('Zonal Wavenumber')
    plt.ylabel('Time Frequency (day$^{-1}$)')
    plt.title('Symmetric')
    ff=plt.pcolormesh(xfx,xfw,(abs(sym)),cmap='PuRd',levels=cv ,extend='max')
    cbar= plt.colorbar(ff)
    return


# save output in Dataset
def mk_dataset(str,x,y,anti,asym):
    dtmp1= xr.DataArray(abs(anti), dims=(str+'xfw', str+'xfx'), coords={str+'xfw': x, str+'xfx': y})
    dtmp2= xr.DataArray(abs(asym), dims=(str+'xfw', str+'xfx'), coords={str+'xfw': x, str+'xfx': y})
    dtmp = dtmp1.to_dataset(name = str+'coh_anti')
    dtmp[str+'coh_sym']=dtmp2
    return dtmp


# In[8]:


# load up the data
obs_coh_anti=dmush.obs_coh_anti
obs_coh_sym =dmush.obs_coh_sym
obs_xfx=dmush.obs_xfx
obs_xfw=dmush.obs_xfw


# In[9]:


# load up the data
cafe_coh_anti=dmush.cafe_coh_anti
cafe_coh_sym =dmush.cafe_coh_sym
cafe_xfx=dmush.cafe_xfx
cafe_xfw=dmush.cafe_xfw


# In[10]:


cafe60_coh_anti=dcafe60.cafe_coh_anti
cafe60_coh_sym =dcafe60.cafe_coh_sym
cafe60_xfx=dcafe60.cafe_xfx
cafe60_xfw=dcafe60.cafe_xfw


# In[29]:


nwindow=128
cv=np.arange(0.05,.5,.05)
nt=int(nwindow/2)

plt.figure(figsize=(15, 20))

plt.subplot(3,2,1)
anti_plot(cv,obs_xfx,obs_xfw[nt:],obs_coh_anti[nt:,:],'JRA')
plt.subplot(3,2,2)
sym_plot(cv,obs_xfx,obs_xfw[nt:],obs_coh_sym[nt:,:])

plt.subplot(3,2,3)
anti_plot(cv,cafe_xfx,cafe_xfw[nt:],cafe_coh_anti[nt:,:],'CAFE')
plt.subplot(3,2,4)
sym_plot(cv,cafe_xfx,cafe_xfw[nt:],cafe_coh_sym[nt:,:])

plt.subplot(3,2,5)
anti_plot(cv,cafe60_xfx,cafe60_xfw[nt:],cafe60_coh_anti[nt:,:],'CAFE60')
plt.subplot(3,2,6)
sym_plot(cv,cafe60_xfx,cafe60_xfw[nt:],cafe60_coh_sym[nt:,:])

plt.savefig('MJO.pdf',dpi=600)


# In[ ]:


#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 
#cb = plt.colorbar(ax1, cax = cbaxes)  
#The numbers in the square brackets of add_axes refer to [left, bottom, width, height], 
#where the coordinates are just fractions that go from 0 to 1 of the plotting area.

