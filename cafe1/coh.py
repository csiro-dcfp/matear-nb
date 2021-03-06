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
# # MJO calculations
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



# + Collapsed="false"
# set some global parameters
# time window for calculation
nwindow=128   
xlat=15  # set in split1

# + [markdown] Collapsed="false"
# ## Obs

# + Collapsed="false"
# for observations
dobs = xr.open_mfdataset('../obs/jra/jra.55.ulwrf.toa.1958010100_2017063021.nc')
dobs = dobs.rename({'g0_lon_2':'lon','g0_lat_1':'lat','initial_time0_hours':'time'})
clim,anom= climatology(dobs,'time')
olr_obs=anom.ULWRF_GDS0_NTAT_ave3h[0:21550,:,:]
olr_obs=olr_obs.reindex(lat=olr_obs.lat[::-1])

dobs2 = xr.open_mfdataset('../obs/jra/jra.55.ugrd.850.1958010100_2016123118.nc')
dobs2 = dobs2.rename({'g0_lon_3':'lon','g0_lat_2':'lat','initial_time0_hours':'time'})
clim,anom= climatology(dobs2,'time')
u850_obs=anom.UGRD_GDS0_ISBL[:,0,:,:]
u850_obs=u850_obs.reindex(lat=u850_obs.lat[::-1])


# + Collapsed="false"
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

coh_anti,coh_sym,xfw,xfx=coherence(nwindow,olr_obs,u850_obs)
# save the coherence values in unique variables
obs_coh_anti=coh_anti.copy()
obs_coh_sym=coh_sym.copy()
obs_xfw=xfw.copy()
obs_xfx=xfx*-1


# + Collapsed="false"
cv=np.arange(0.05,.75,.05)
nt=int(nwindow/2)
plt.figure(figsize=(18, 10))
plt.subplot(2,2,1)
plt.xlim(-15,15)
plt.title('Anti')
ff=plt.contourf(xfx*-1,xfw[nt:],(abs(coh_anti[nt:,:]) ),cmap='BuPu',levels=cv,extend='max' )
cbar = plt.colorbar(ff)
plt.subplot(2,2,2)
plt.xlim(-15,15)
plt.xlabel('Zonal Wavenumber')
plt.ylabel('Time Frequency (day$^{-1}$)')
plt.title('Symmetric')
ff=plt.contourf(xfx*-1,xfw[nt:],(abs(coh_sym[nt:,:]) ),cmap='BuPu',levels=cv,extend='max' )
cbar= plt.colorbar(ff)

# + [markdown] Collapsed="false"
# ## ACCESS-ESM

# + Collapsed="false"
# for ACCESS-ESM
daccess = xr.open_mfdataset('../cmip6/rlut_day_ACCESS-ESM1-5_historical_r1i1p1f1_gn_19500101-19991231.nc')
clim1,anom1= climatology(daccess,'time')
olr_access=anom1.rlut
#anti,sym,total=split1(olr_access)

daccess2 = xr.open_mfdataset('../cmip6/ua_day*19*.nc',combine='by_coords')
lon=daccess.lon
lat=daccess.lat
clim2,anom2= climatology(daccess2,'time')
u850_access=anom2.ua[:,1,:,:]
# put on same grid as olr
dnew=u850_access.interp(lon=lon, lat=lat)
#did not work u850_access =dtmp[:,:,:].interp_like(olr_access)
#anti_u,sym_u,total_u=split1(dnew)
coh_anti,coh_sym,xfw,xfx=coherence(nwindow,olr_access,dnew)

# save the coherence values in unique variables
access_coh_anti=coh_anti.copy()
access_coh_sym=coh_sym.copy()
access_xfw=xfw*1
access_xfx=xfx*-1

# + Collapsed="false"
cv=np.arange(0.05,.75,.05)
nt=int(nwindow/2)
plt.figure(figsize=(18, 10))
plt.subplot(2,2,1)
plt.xlim(-15,15)
plt.title('Anti')
ff=plt.contourf(xfx*-1,xfw[nt:],(abs(coh_anti[nt:,:]) ),cmap='BuPu',levels=cv ,extend='max')
cbar = plt.colorbar(ff)
plt.subplot(2,2,2)
plt.xlim(-15,15)
plt.xlabel('Zonal Wavenumber')
plt.ylabel('Time Frequency (day$^{-1}$)')
plt.title('Symmetric')
ff=plt.contourf(xfx*-1,xfw[nt:],(abs(coh_sym[nt:,:]) ),cmap='BuPu',levels=cv ,extend='max')
cbar= plt.colorbar(ff)

# + Collapsed="false"



# + [markdown] Collapsed="false"
# ## CAFE

# + Collapsed="false"
# for CAFE
dcafe = xr.open_mfdataset('cafe/atmos_daily_040*.nc',combine='by_coords')
clim1,anom1= climatology(dcafe,'time')

olr_cafe=anom1.olr
u850_cafe=anom1.ucomp[:,15,:,:]

nwindow=128
coh_anti,coh_sym,xfw,xfx=coherence(nwindow,olr_cafe,u850_cafe)
# save the coherence values in unique variables
cafe_coh_anti=coh_anti.copy()
cafe_coh_sym=coh_sym.copy()
cafe_xfw=xfw*1
cafe_xfx=xfx*-1


# + Collapsed="false"
cv=np.arange(0.05,.75,.05)
nt=int(nwindow/2)
plt.figure(figsize=(18, 10))
plt.subplot(2,2,1)
plt.xlim(-15,15)
plt.title('Anti')
ff=plt.contourf(xfx*-1,xfw[nt:],(abs(coh_anti[nt:,:]) ),cmap='BuPu',levels=cv ,extend='max')
cbar = plt.colorbar(ff)
plt.subplot(2,2,2)
plt.xlim(-15,15)
plt.xlabel('Zonal Wavenumber')
plt.ylabel('Time Frequency (day$^{-1}$)')
plt.title('Symmetric')
ff=plt.contourf(xfx*-1,xfw[nt:],(abs(coh_sym[nt:,:]) ),cmap='BuPu',levels=cv ,extend='max')
cbar= plt.colorbar(ff)

# + Collapsed="false"



# + [markdown] Collapsed="false"
# ## Final figures with all output

# + Collapsed="false"
def anti_plot(cv,xfx,xfw,anti):
    plt.xlim(-15,15)
    plt.title('Anti')
    plt.xlabel('Zonal Wavenumber')
    plt.ylabel('Time Frequency (day$^{-1}$)')
    ff=plt.contourf(xfx,xfw,(abs(anti)),cmap='BuPu',levels=cv ,extend='max')
    cbar = plt.colorbar(ff)
    return
                    
def sym_plot(cv,xfx,xfw,sym):
    plt.xlim(-15,15)
    plt.xlabel('Zonal Wavenumber')
    plt.ylabel('Time Frequency (day$^{-1}$)')
    plt.title('Symmetric')
    ff=plt.contourf(xfx,xfw,(abs(sym)),cmap='BuPu',levels=cv ,extend='max')
    cbar= plt.colorbar(ff)
    return

# save output in Dataset
def mk_dataset(str,x,y,anti,asym):
    dtmp1= xr.DataArray(abs(anti), dims=(str+'xfw', str+'xfx'), coords={str+'xfw': x, str+'xfx': y})
    dtmp2= xr.DataArray(abs(asym), dims=(str+'xfw', str+'xfx'), coords={str+'xfw': x, str+'xfx': y})
    dtmp = dtmp1.to_dataset(name = str+'coh_anti')
    dtmp[str+'coh_sym']=dtmp2
    return dtmp



# + Collapsed="false"
cv=np.arange(0.05,.75,.05)
nt=int(nwindow/2)

plt.figure(figsize=(15, 20))

plt.subplot(3,2,1)
anti_plot(cv,obs_xfx,obs_xfw[nt:],obs_coh_anti[nt:,:])
plt.subplot(3,2,2)
sym_plot(cv,obs_xfx,obs_xfw[nt:],obs_coh_sym[nt:,:])


plt.subplot(3,2,3)
anti_plot(cv,cafe_xfx,cafe_xfw[nt:],cafe_coh_anti[nt:,:])
plt.subplot(3,2,4)
sym_plot(cv,cafe_xfx,cafe_xfw[nt:],cafe_coh_sym[nt:,:])

plt.subplot(3,2,5)
anti_plot(cv,access_xfx,access_xfw[nt:],access_coh_anti[nt:,:])
plt.subplot(3,2,6)
sym_plot(cv,access_xfx,access_xfw[nt:],access_coh_sym[nt:,:])

dobsc=mk_dataset('obs_',obs_xfw,obs_xfx,obs_coh_anti,obs_coh_sym)
daccessc=mk_dataset('access_',access_xfw,access_xfx,access_coh_anti,access_coh_sym)
dcafec=mk_dataset('cafe_',cafe_xfw,cafe_xfx,cafe_coh_anti,cafe_coh_sym)
dall=xr.merge([dobsc,daccessc,dcafec])
dall.to_netcdf('coherence.nc')

plt.savefig('MJO.pdf',dpi=600)



# + [markdown] Collapsed="false" toc-hr-collapsed=true Collapsed="false"
# ## extra

# + Collapsed="false"

