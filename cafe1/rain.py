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
# # Rainfall
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

def correlate_year(a,b):
    astd=a.std(dim=('year'))
    bstd=b.std(dim=('year'))
    adev=a-a.mean(dim=('year'))
    bdev=b-b.mean(dim=('year'))
    ab=adev*bdev
    abstd=ab.mean(dim='year')
    r =abstd/astd/bstd
    return r

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x,0,0)) 

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


# + [markdown] Collapsed="false"
# ## Rainfall data from awap, cafe and access

# + Collapsed="false"
# for observations
dobs=xr.open_mfdataset('../obs/pearcey/awap_rain.nc',parallel=True,combine='by_coords')
dcafe=xr.open_mfdataset('../obs/pearcey/cafe_rain.nc',parallel=True,combine='by_coords')
daccess = xr.open_mfdataset('../obs/pearcey/access_rain.nc',combine='by_coords')

obs_rain=dobs.rain_month

cafe_rain=dcafe.precip
cafe_oz = region(cafe_rain,-50,-10, 110,160) *1e-3 *1e3 *86400*365/12.

access_rain=daccess.pr
access_oz = region(access_rain,-50,-10, 110,160) *1e-3 *1e3 *86400*365/12.



# + Collapsed="false"
print(daccess)
print(dobs)
print(dcafe)

# + Collapsed="false"
# interpolate on to obs grid
lons=obs_rain.longitude
lats=obs_rain.latitude

cafe=cafe_oz.interp(lon=lons, lat=lats)
access= access_oz.interp(lon=lons, lat=lats) +obs_rain*0

# + Collapsed="false"
fig=plt.figure(figsize=(20,14))
seasons=['DJF','MAM','JJA','SON']
cv=np.arange(0,200,10)
cmap='viridis_r'

lons=cafe.lon.values
lats=cafe.lat.values

for i in range(12):
    j=int((i)/4)
    l=i-j*4
    if j == 0 :
        zar=obs_rain.copy()
    if j == 1 :
        zar=cafe.copy()
    if j == 2:
        zar=access.copy()
    ax=plt.subplot(3,4,1+i)
    print(i,l)
    p=plt.contourf(lons,lats,zar[l,:,:],levels=cv,cmap=cmap,extend='max')
    if j == 0:
        plt.title(seasons[l])
    if l== 0:
        plt.ylabel('Latitude')
    if j==2 :
        plt.xlabel('Longitude')
    
#    plt.colorbar(p)
#plt.colorbar(p, shrink=.8, orientation='horizontal', label='Monthly averaged Rainfall (mm)')


plt.savefig('rain_seaons.pdf',dpi=600)
    
    
    

# + Collapsed="false"
for i in range(4):
    print(seasons[i],obs_rain[i,:,:].mean(axis=(0,1)).values, cafe[i,:,:].mean(axis=(0,1)).values,access[i,:,:].mean(axis=(0,1)).values)

# + Collapsed="false"
# define some regions following Watterson \cite{Watterson:2013vz}


# + [markdown] Collapsed="false"
# ## Compute correlation with NINO34

# + Collapsed="false"
dcafe=xr.open_mfdataset('../obs/pearcey/cafe_anom_rain.nc',parallel=True,combine='by_coords')
daccess = xr.open_mfdataset('../obs/pearcey/access_anom_rain.nc',combine='by_coords')

#ocn_cafe1=ocn_cafe.assign_coords(lon=np.where(ll<0,ll+360,ll))
dland=xr.open_mfdataset('../cmip6/gpp_Lmon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc',
    parallel=True,combine='by_coords')
    
cafe_rain=dcafe.precip
cafe_oz = region(cafe_rain,-50,-10, 110,160) *1e-3 *1e3 *86400*365/12.

access_rain=daccess.pr
access_oz1 = region(access_rain,-50,-10, 110,160) *1e-3 *1e3 *86400*365/12.

lons=obs_rain.longitude
lats=obs_rain.latitude

#o=ocean_oz.fillna(-1)
#omask= o.where(o == -1)*-1
landmask1=dland.gpp[0,:,:]
landmask = landmask1 *0 + 1
land_oz = region(landmask,-50,-10, 110,160)
access_oz=access_oz1*land_oz

#cafe=cafe_oz.interp(lon=lons, lat=lats)
#access= access_oz.interp(lon=lons, lat=lats) +obs_rain*0

# + Collapsed="false"
# get nino3.4 index
dnino_cafe=xr.open_mfdataset('cafe_enso.nc',parallel=True,combine='by_coords')
dnino_access=xr.open_mfdataset('access_enso.nc',parallel=True,combine='by_coords')
cafe_nino34=dnino_cafe.cafe_nino34
access_nino34=dnino_access.access_nino34
# put cafe rain data on the same time axis as its nino index
cafe_oz=cafe_oz.assign_coords(time=dnino_cafe.time[0:1200].values)

# + Collapsed="false"
#astd,bstd,abstd,bhad=regresst(cafe_nino34,cafe_oz)
astd,bstd,abstd,bhat=regresst(access_nino34,access_oz)
r =(abstd/astd/bstd).copy()
ra=r*land_oz

astd1,bstd1,abstd1,bhat1=regresst(cafe_nino34,cafe_oz)
rc =(abstd1/astd1/bstd1).copy()

# mean of the correlation
print(rc.mean(axis=(0,1)).values,ra.mean(axis=(0,1)).values)

mm1 = cafe_oz.mean(axis=(1,2))
mm2 = access_oz.mean(axis=(1,2))

r1=correlate_year(cafe_nino34.groupby('time.year').mean('time'),
                 mm1.groupby('time.year').mean('time'))

r2=correlate_year(access_nino34.groupby('time.year').mean('time'),
                 mm2.groupby('time.year').mean('time'))

print(r1.values,r2.values)


# + Collapsed="false"
# Extract a given season for time series of index and maps
def season_extract(season,nino34,rain):
# mask other months with nan
    ds_SS = nino34.where(nino34['time.season'] == season)
# rolling mean -> only Jan is not nan
# however, we loose Jan/ Feb in the first year and Dec in the last
    ds_SS = ds_SS.rolling(min_periods=3, center=True, time=3).mean()
# make annual mean SS index
    ds_SS = ds_SS.groupby('time.year').mean('time')
# compute map of SS mean
    map_SS= rain.where(rain['time.season'] == season)
    map_SS = map_SS.rolling(min_periods=3, center=True, time=3).mean()
    map_SS = map_SS.groupby('time.year').mean('time')
# compute correlation 
    r = correlate_year(ds_SS,map_SS)
    m1= map_SS.mean(axis=(1,2))
    r1=correlate_year(ds_SS,m1)
    
    
    return ds_SS,map_SS,r,r1


# + Collapsed="false"
rc=cafe_oz[0:4,:,:].copy()
ra=access_oz[0:4,:,:].copy()
rc.load()
ra.load()
rm1=cafe_nino34[0:4].copy()
rm2=access_nino34[0:4].copy()
rm1.load()
rm2.load()
for i in range(4):
    tsc,mapsc,rsc,rm=season_extract(seasons[i],cafe_nino34,cafe_oz)
    rc[i,:,:]=rsc.copy()
    rm1[i]=rm.copy()
    tsa,mapsa,rsa,rm=season_extract(seasons[i],access_nino34,access_oz)
    ra[i,:,:]=rsa.copy()
    rm2[i]=rm.copy()
    

# + Collapsed="false"
from itertools import product
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# interpolate on to obs grid
lons=obs_rain.longitude
lats=obs_rain.latitude

z2=ra*land_oz
z1=rc

z1i=z1.interp(lon=lons, lat=lats)
z2i=z2.interp(lon=lons, lat=lats)

fig=plt.figure(figsize=(20,40))
seasons=['DJF','MAM','JJA','SON']
cv=np.arange(-1,1.05,.1)

cmap='RdBu'

for i,j in product (range(4), range(2)) :
    ax=plt.subplot(4,2,i*2+j+1)
    if j == 0 :
        p=plt.contourf(lons,lats,z1i[i,:,:],levels=cv,cmap=cmap)
        plt.contour(lons,lats,z1i[i,:,:],levels=(-.3,.3))
#        rc[i].plot(levels=cv,cmap=cmap)
    if j == 1 :
        p=plt.contourf(lons,lats,z2i[i,:,:],levels=cv,cmap=cmap)
        plt.contour(lons,lats,z2i[i,:,:],levels=(-.3,.3))
#        z2[i].plot(levels=cv,cmap=cmap)
    plt.title(seasons[i])
    if j == 0 : 
        plt.ylabel('Latitude')
    if i == 3 :
        plt.xlabel('Longitude')
        
axins = inset_axes(ax,
                   width="5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.05, 0., 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
fig.colorbar(p,cax=axins)
    
#    plt.colorbar(p)
#plt.colorbar(p, shrink=.8, orientation='horizontal', label='Monthly averaged Rainfall (mm)')
plt.savefig('rain_enso.pdf',dpi=600)


# + Collapsed="false"
# Australia rainfall with ENSO3.4 by season
for i in range(4) :
#    print( seasons[i], rc[i,:,:].mean(axis=(0,1)).values, z2[i,:,:].mean(axis=(0,1)).values )
    print( seasons[i], rm1[i].values, rm2[i].values )
    

# + Collapsed="false"
print( rc[i,:,:].mean(axis=(0,1)).values)

# + Collapsed="false"


# + Collapsed="false"
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

