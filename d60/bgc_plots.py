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

# + [markdown] Collapsed="true"
# ## integrated fluxes

# + Collapsed="false"
def readn(ff):
    f1=xr.open_dataset(ddout+ff+'.nc')
    ff=f1.__xarray_dataarray_variable__.rename(ff).load()
    return(ff)


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


# + Collapsed="false" jupyter={"outputs_hidden": true}
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

# + [markdown] Collapsed="false"
# ## flux maps
# -



# +
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
# ## sections

# + Collapsed="false"
def reads1(ff,f2):
    f1=xr.open_dataset(ddout+ff+'.nc')
    ff=f1[f2]
    return(ff)


# + Collapsed="false"
Depth=np.arange(0,3300,100)
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

# + Collapsed="false"

