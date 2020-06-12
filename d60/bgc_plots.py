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
# # Make the BGC figures for the cafe60 paper

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
# scratch is the logical link to the data
ddout='scratch/'
#dgrid=xr.open_dataset(file1)

# + [markdown] Collapsed="false"
# ## integrated fluxes

# + Collapsed="false"
def readn(ff):
    f1=xr.open_dataset(ddout+ff+'.nc')
    ff=f1.__xarray_dataarray_variable__.rename(ff).load()
    return(ff)

def reads1(ff,f2):
    f1=xr.open_dataset(ddout+ff+'.nc')
    ff=f1[f2]
    return(ff)

def region(ds,lat1,lat2):
    dsr=ds.sel(lat=slice(lat1,lat2))
    return (dsr)


# + Collapsed="false"
time1=np.arange(1960.0833,2019.5,.08333)

totd=readn('totd')
tota=readn('tota')
gstd_d=readn('gstd_d')
gstd_a=readn('gstd_a')

totsd=readn('totsd')
totsa=readn('totsa')
sstd_d=readn('sstd_d')
sstd_a=readn('sstd_a')

# + Collapsed="false"
oflx=reads1('oc_v1.7_daily','co2flux_ocean')
sflx=oflx.sel(lat=slice(-90,-35))
time2=np.arange(1982,2019.026,.00274)


# + Collapsed="false"
# Plot figure
def flxplt(a,b,da,db):
    a1=(a+da).rolling(time=12, center=True).mean()
    a2=(a-da).rolling(time=12, center=True).mean()
    b1=(b+db).rolling(time=12, center=True).mean()
    b2=(b-db).rolling(time=12, center=True).mean()
    plt.plot(time1,a.rolling(time=12, center=True).mean(),label=tl1,color='blue')
    plt.plot(time1,b.rolling(time=12, center=True).mean(),label=tl2,color='orange')
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
tl1='Natural'
tl2="Total"
t1='Global'
flxplt(totd*(-1e-15),tota*(-1e-15),gstd_d.std(axis=1)*1e-15,gstd_a.std(axis=1)*1e-15)
plt.subplot(2,2,2)
t1='Southern Ocean'
flxplt(totsd*(-1e-15),totsa*(-1e-15),sstd_d.std(axis=1)*1e-15,sstd_a.std(axis=1)*1e-15)

plt.subplot(2,2,3)
plt.xlim(1982,2020)
plt.ylim(-4,0)
# global
t1='Global'
tl1='Observation'
tl2="Total"
flxplt(totd*(-1e-15)*0,tota*(-1e-15)+.5,gstd_d.std(axis=1)*1e-15*0,gstd_a.std(axis=1)*1e-15)
plt.plot(time2,oflx.sum(axis=(1,2)).rolling(mtime=365, center=True).mean(),color='blue')

plt.subplot(2,2,4)
plt.xlim(1982,2020)
plt.ylim(-2.5,0)
t1='Southern Ocean'
flxplt(totsd*(-1e-15)*0,totsa*(-1e-15),sstd_d.std(axis=1)*1e-15*0,sstd_a.std(axis=1)*1e-15)

plt.plot(time2,sflx.sum(axis=(1,2)).rolling(mtime=365, center=True).mean(),color='blue')

plt.savefig('fco2.pdf',dpi=600)


# + [markdown] Collapsed="false"
# ## flux maps
# + Collapsed="false"
def reads1(ff,f2):
    f1=xr.open_dataset(ddout+ff+'.nc')
    ff=f1[f2]
    return(ff)


# + Collapsed="false"
oflx=reads1('oflx','fgco2_clim')
aflx=reads1('aflx','stf10')


# + Collapsed="false"
# Generic Plot Map script using subplot
def map1(a,y,x):
    import cartopy
    from mpl_toolkits.axes_grid1 import AxesGrid
    from cartopy.mpl.geoaxes import GeoAxes
#    proj=ccrs.Mollweide()
#    proj=ccrs.Robinson(central_longitude=180)
    dproj=ccrs.PlateCarree()
    ax.coastlines(resolution='110m')
    ax.gridlines()
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black',
                   color='grey') #'white')
    p=ax.contourf(x, y, a,levels=cv, cmap=cmap,transform=dproj,extend='both')
    
#    
    plt.colorbar(p, shrink=.8, orientation='horizontal',
                 label=cbtit)
    plt.title(t1)

    return
    


# + Collapsed="false"
plt.figure(figsize=(12,5 ))
cmap='RdBu_r'
cbtit='Sea to Air flux (g C m$^{-2}$ y$^{-1}$)'
cv=np.arange(-50,55,5)
proj=ccrs.Robinson(central_longitude=180)

ax = plt.subplot(1,2,1,projection=proj)
t1='Observations'
map1(oflx,oflx.lat,oflx.lon)

ax = plt.subplot(1,2,2,projection=proj)
t1='CAFE60'
map1(aflx*(-1),aflx.yt_ocean,aflx.xt_ocean)

plt.savefig('flx_map.pdf',dpi=600)


# + [markdown] Collapsed="false"
# ## sections

# + Collapsed="false"
def reads1(ff,f2):
    f1=xr.open_dataset(ddout+ff+'.nc')
    ff=f1[f2]
    return(ff)


# + Collapsed="false"
Depth=reads1('odepth','Depth')
otco2=reads1('otco2','TCO2')
opo4=reads1('opo4','PO4')
oaco2=reads1('oaco2','Cant')

adic=reads1('adic','adic')
dic=reads1('dic','dic')
po4=reads1('po4','no3')
aco2=adic-dic


# + Collapsed="false"
def splt(a,y,x):
    plt.xlabel('Latitude')
    plt.ylabel('Depth')
    plt.title(t1)
    plt.gca().patch.set_color('.25')
 #   plt.ylim(5500,0)
    ff=plt.contourf(x,y,a,cmap=cmap,levels=cv ,extend='both')
    ax = plt.gca()
    ax.invert_yaxis()
    cbar = plt.colorbar(ff)
    cbar.set_label('mmol m$^{-3}$ ')



# + Collapsed="false"
plt.figure(figsize=(11,10 ))

cv=np.arange(1800,2500,50)
cmap='RdBu_r'
plt.subplot(2,1,1)
t1='Observations'
splt(otco2,Depth,otco2.lat)  

plt.subplot(2,1,2)
t1='CAFE60'
#tco2=adic.sel(xt_ocean=-179.5).load()
splt(adic,adic.st_ocean,adic.yt_ocean)  

plt.savefig('adic_section.pdf',dpi=600)

# + Collapsed="false"
plt.figure(figsize=(11,10 ))

cmap='BuGn'
cv=np.arange(0,4,.25)
plt.subplot(2,1,1)
t1='Observations'
splt(opo4,Depth,opo4.lat)  

plt.subplot(2,1,2)
t1='CAFE60'
splt(po4,po4.st_ocean,po4.yt_ocean)  

plt.savefig('po4_section.pdf',dpi=600)

# + Collapsed="false"
plt.figure(figsize=(11,10 ))

cmap='PuRd'
cv=np.arange(0,60,5)
plt.subplot(2,1,1)
t1='Observations'
splt(oaco2,Depth,oaco2.lat)  

plt.subplot(2,1,2)
t1='CAFE60'
splt(aco2,aco2.st_ocean,aco2.yt_ocean)  

plt.savefig('aco2_section.pdf',dpi=600)


# + [markdown] Collapsed="false"
# ## integrate aco2 maps

# + Collapsed="false"
def reads1(ff,f2):
    f1=xr.open_dataset(ddout+ff+'.nc')
    ff=f1[f2]
    return(ff)


# + Collapsed="false"
mcant=reads1('mcinv','a')
ocant=reads1('oinv','a')

# + Collapsed="false"
plt.figure(figsize=(12,5 ))
cmap='Reds'
#cmap='PuRd'
cbtit='Anthropogenic CO$_2$ (mol m$^{-2}$)'
cv=np.arange(0,250,20)
proj=ccrs.Robinson(central_longitude=180)

ax = plt.subplot(1,2,1,projection=proj)
t1='Observations'
map1(ocant,ocant.lat,ocant.lon)

ax = plt.subplot(1,2,2,projection=proj)
t1='CAFE60'
map1(mcant,mcant.yt_ocean,mcant.xt_ocean)

plt.savefig('cant_map.pdf',dpi=600)

# + Collapsed="false"

