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

# + [markdown] Collapsed="false" toc-hr-collapsed=false Collapsed="false"
# # PDO 
#

# + Collapsed="false"

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

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x,0,0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)



# + Collapsed="false"

# compute PDO
def PDO(dat1,dat2):
    solver = Eof(dat1)
    pc1 = solver.pcs(npcs=1, pcscaling=1)
    sigs=solver.eigenvalues()
    eofs = solver.eofs()
    eof1=eofs[0]
# normalise sigs
    sigs1=sigs/sigs.sum()
# not used eof1 = solver.eofsAsCorrelation(neofs=1)
#    pcs=solver.pcs()
    
# filter first pdf
    Nt=dat1.time.size
    fp=fft.fft(pc1[:,0])
    x=fft.fftfreq(Nt,1/12.)  # cycles per year
# lowpass filter at 0.1 per year
    i=abs(x)<=.1
    fp_fil=fp*i
#plt.plot(x,abs(fp))
#plt.plot(x,abs(fp_fil))
    pfil=fft.ifft(fp_fil)
#
    pc1_fil=pc1[:,0]*0+np.real(pfil)
#plt.plot(pc1_fil)
#print(pc1)
    print(pc1_fil)
    tmp=np.imag(pfil)
    print(tmp.max)
# correlate with sst field
    
    astd,bstd,abstd,bhat=regresst(pc1_fil,dat2)
    r=abstd/(astd*bstd)
    slope=r*bstd
    
    return sigs1,eof1,slope,pc1,pc1_fil,x,fp_fil



# + Collapsed="false"
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



# + [markdown] Collapsed="false" toc-hr-collapsed=false Collapsed="false"
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

anom=obs_anom.sst.copy()
# detrend data
dtrend,trend=dtrend_t(anom)

plt.plot(dtrend[:,90,180])
plt.plot(trend[90,180,:])

obs_dtrend=dtrend.copy()

dsub=region(obs_dtrend,20,70,360-240,360-110)
sigs1,eof1,slope,pc1,pc1_fil,x,fp_fil=PDO(dsub,obs_dtrend)
#cafe_sigs,cafe_slope,cafe_pc1,cafe_pc1_fil,cafe_x,cafe_fp_fil=PDO(dsub,cafe_dtrend)
#slope.plot(levels=40)
obs_sig = sigs1.copy()
obs_eof1 = eof1.copy()
obs_slope = slope.copy()
obs_pc1 = pc1.copy()  # pc1
obs_pc1_fil = pc1_fil.copy()   # filtered pc1
obs_x=x.copy()
obs_fp_fil=fp_fil.copy()


# + [markdown] Collapsed="false"
# ### test plots

# + Collapsed="false"
plt.figure(figsize=(7,10))
plt.subplot(3,1,1)
plt.plot(pc1)
plt.subplot(3,1,2)
sigs1[0:10].plot(color='b', linewidth=2)
plt.subplot(3,1,3)
eof1.plot.contourf(levels=20)

# + Collapsed="false"



# + [markdown] Collapsed="false" toc-hr-collapsed=false Collapsed="false"
# ## ACCESS

# + Collapsed="false"
# for ACCESS-ESM (use the cafe infor for lat and lon axis)
dcafe = xr.open_mfdataset('c2_sst.nc')
dcafe = dcafe.rename({'XT_OCEAN':'lon','YT_OCEAN':'lat','TIME1':'time'})

daccess = xr.open_mfdataset('../cmip6/tos_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc')
access_clim,access_anom= climatology(daccess,'time')
access_mean=access_clim.tos.mean(dim='month')
access_anom1=access_anom.tos

# change access grid values
cafe_gd=access_anom1.copy()
cafe_gd=cafe_gd.assign_coords(j=dcafe.lat.values,i=dcafe.lon.values)
cafe_gd=cafe_gd.rename({'i':'lon','j':'lat',})

anom=cafe_gd.copy()
# detrend data
dtrend,trend=dtrend_t(anom)

plt.plot(dtrend[:,90,180])
plt.plot(trend[90,180,:])

access_dtrend=dtrend.copy()

dsub=region(access_dtrend,20,70,-240,-110)
sigs1,eof1,slope,pc1,pc1_fil,x,fp_fil=PDO(dsub,access_dtrend)
#cafe_sigs,cafe_slope,cafe_pc1,cafe_pc1_fil,cafe_x,cafe_fp_fil=PDO(dsub,cafe_dtrend)
#slope.plot(levels=40)
access_sig = sigs1.copy()
access_eof1 = eof1.copy()
access_slope = slope.copy()
access_pc1 = pc1.copy()  # pc1
access_pc1_fil = pc1_fil.copy()   # filtered pc1
access_x=x.copy()
access_fp_fil=fp_fil.copy()


# + [markdown] Collapsed="false"
# ### test plots

# + Collapsed="false"
plt.figure(figsize=(7,10))
plt.subplot(3,1,1)
plt.plot(pc1)
plt.subplot(3,1,2)
sigs1[0:10].plot(color='b', linewidth=2)
plt.subplot(3,1,3)
eof1.plot.contourf(levels=20)

# + [markdown] Collapsed="false" toc-hr-collapsed=false Collapsed="false"
# ## CAFE

# + Collapsed="false"
# for CAFE (use the time axis from access for getting a pandas compatible version for detrending)
daccess = xr.open_mfdataset('../cmip6/tos_Omon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc')
dcafe = xr.open_mfdataset('c2_sst.nc')
dcafe = dcafe.rename({'XT_OCEAN':'lon','YT_OCEAN':'lat','TIME1':'time'})
cafe_clim,cafe_anom= climatology(dcafe,'time')

cafe_mean=cafe_clim.SST[:,0,:,:].mean(dim='month')
cafe_anom1=cafe_anom.SST[:,0,:,:]

# get time and put it into a variable of float as days since initial time
anom=cafe_anom1.assign_coords(time=daccess.time[0:1200].values)
# detrend data
dtrend,trend=dtrend_t(anom)

plt.plot(dtrend[:,90,180])
plt.plot(trend[90,180,:])
cafe_dtrend=dtrend.copy()

dsub=region(cafe_dtrend,20,70,-240,-110)
sigs1,eof1,slope,pc1,pc1_fil,x,fp_fil=PDO(dsub,cafe_dtrend)
#cafe_sigs,cafe_slope,cafe_pc1,cafe_pc1_fil,cafe_x,cafe_fp_fil=PDO(dsub,cafe_dtrend)
#slope.plot(levels=40)
cafe_sig = sigs1.copy()
cafe_eof1 = eof1.copy()
cafe_slope = slope.copy()
cafe_pc1 = pc1.copy()  # pc1
cafe_pc1_fil = pc1_fil.copy()   # filtered pc1
cafe_x=x.copy()
cafe_fp_fil=fp_fil.copy()


# + [markdown] Collapsed="false"
# ### test plots

# + Collapsed="false"
plt.figure(figsize=(7,10))
plt.subplot(3,1,1)
plt.plot(pc1)
plt.subplot(3,1,2)
sigs1[0:10].plot(color='b', linewidth=2)
plt.subplot(3,1,3)
eof1.plot.contourf(levels=20)


# + [markdown] Collapsed="false"
# # Final plots

# + Collapsed="false"
# Generic Plot Map script using subplot
def map3o(z1, z2, z3,tit1,tit2,tit3):
    import cartopy
    from mpl_toolkits.axes_grid1 import AxesGrid
    from cartopy.mpl.geoaxes import GeoAxes
    proj=ccrs.Mollweide()
    proj=ccrs.Robinson(central_longitude=180)
    dproj=ccrs.PlateCarree()
    cv=np.arange(-1,1.1,0.1)
    vm=4
    
#
    for  i in range(3):
#    ax = plt.axes(projection=ccrs.Mollweide())
#    ax = plt.axes(projection=proj)
        ax = plt.subplot(3,2,i+2,projection=proj)
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
    
def fill_plot(a):
    time=a.coords['time'].values
    plt.plot(time,a, color='black')
    ax.fill_between(time,a, where=a>0, facecolor='red', interpolate=True)
    ax.fill_between(time,a, where=a<0, facecolor='blue', interpolate=True)
    return



# + Collapsed="false"

plt.figure(figsize=(16, 25))

ax=plt.subplot(3,2,1)
plt.plot(obs_x[0:20],abs(obs_fp_fil[0:20]))
plt.plot(cafe_x[0:20],abs(cafe_fp_fil[0:20]))
plt.plot(access_x[0:20],abs(access_fp_fil[0:20]))

plt.xlim(0,.11)
plt.xlabel('Frequency (cycles per month)')
plt.ylabel('Power')
plt.title('a)',loc='left')


map3o(obs_slope*-1,cafe_slope,access_slope*-1,'Observations','CAFE','ACCESS')

ax = plt.subplot(9,1,7)
fill_plot(obs_pc1_fil*-1)

ax = plt.subplot(9,1,8)
fill_plot(cafe_pc1_fil)

ax = plt.subplot(9,1,9)
fill_plot(access_pc1_fil*-1)

plt.savefig('PDO.pdf',dpi=600)


# + Collapsed="false"
plt.figure(figsize=(15, 15))
ax = plt.subplot(3,2,1)
fill_plot(obs_pc1_fil*-1)
ax = plt.subplot(9,1,7)
fill_plot(obs_pc1_fil*-1)
ax = plt.subplot(9,1,8)
fill_plot(obs_pc1_fil*-1)
ax = plt.subplot(9,1,9)
fill_plot(obs_pc1_fil*-1)





# + Collapsed="false"

