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
# # Calculations for WMO
# ## tsurface 
#
#

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

file2='/home/mat236/dcfp/big/data/csiro-dcfp-cafe-f5/atmos_isobaric_month.zarr'
file1='/home/mat236/area.nc'

dgrid=xr.open_dataset(file1)
ddat = xr.open_zarr(file2,consolidated=True)

# + Collapsed="false"
file3='/home/mat236/dcfp/big/data/csiro-dcfp-jra55/surface_month_cafe-grid.zarr'
#file3='/home/mat236/dcfp/big/data/csiro-dcfp-jra55/derived_surface_month_cafe-grid.zarr'

dobs = xr.open_zarr(file3,consolidated=True)


# + Collapsed="false"
# change from t_surf (skin temperature) to t_ref (2m Temperature)
ddat.t_ref[50,0,1,:,:].plot()

# + Collapsed="false"
dobs.TMP_GDS0_HTGL[50*12,:,:].plot()


# + [markdown] Collapsed="false"
# # Prediction 

# + [markdown] Collapsed="false"
# ## Maps

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
ptemp1=ddat.t_ref
ptemp2=ptemp1 - ptemp1.mean(dim=('init_date','ensemble'))

# + Collapsed="false"
plt.figure(figsize=(12,5 ))
cmap='RdBu_r'
cbtit='Ensemble Mean ($^o$C)'
cv=np.arange(-1.5,1.7,.25)
proj=ccrs.Robinson(central_longitude=180)

ax = plt.subplot(1,2,1,projection=proj)
t1='Year 1'
ptemp3= ptemp2[59,0:12,:,:,:].mean(dim=('ensemble'))
oflx=ptemp3.mean(axis=0)
ptemp3= ptemp2[59,0:12,:,:,:].mean(dim=('lead_time'))
oflx=ptemp3.mean(dim=('ensemble'))
map1(oflx,oflx.lat,oflx.lon)

ax = plt.subplot(1,2,2,projection=proj)
t1='Years 2-5'
ptemp3= ptemp2[59,12:60,:,:,:].mean(dim=('ensemble'))
oflx=ptemp3.mean(axis=0)
ptemp3= ptemp2[59,12:60,:,:,:].mean(dim=('lead_time'))
oflx=ptemp3.mean(dim=('ensemble'))
map1(oflx,oflx.lat,oflx.lon)

plt.savefig('predict_mean.pdf',dpi=600)



# + Collapsed="false"
# Ensemble spread
plt.figure(figsize=(12,5 ))
cmap='Reds'
cmap='PuRd'
cbtit='Ensemble Standard Deviation ($^o$C)'
cv=np.arange(0.25,2.25,.25)
proj=ccrs.Robinson(central_longitude=180)

ax = plt.subplot(1,2,1,projection=proj)
t1='Year 1'
ptemp3= ptemp2[59,0:12,:,:,:].mean(dim=('lead_time'))
oflx=ptemp3.std(dim=('ensemble'))
map1(oflx,oflx.lat,oflx.lon)

ax = plt.subplot(1,2,2,projection=proj)
t1='Years 2-5'
ptemp3= ptemp2[59,12:60,:,:,:].mean(dim=('lead_time'))
oflx=ptemp3.std(dim=('ensemble'))
map1(oflx,oflx.lat,oflx.lon)

plt.savefig('predict_spread.pdf',dpi=600)


# + [markdown] Collapsed="false"
# ## Global Timeseries

# + Collapsed="false"
def weighting(a):
    # weightings 
    weights = np.cos(np.deg2rad(a.lat))
    weights.name = "weights"
    
    aw=a*weights
    norm=(aw*0+1)*weights
    ntmp=norm.sum(dim=('lat','lon'))
    nn=ntmp[0,0]
    am=aw.sum(dim=("lon", "lat"))/nn
    
    return(am)

def weighting0(a):
    # weightings 
    weights = np.cos(np.deg2rad(a.lat))
    weights.name = "weights"
    
    aw=a*weights
    norm=(aw*0+1)*weights
    ntmp=norm.sum(dim=('lat','lon'))
    nn=ntmp[0]
    am=aw.sum(dim=("lon", "lat"))/nn
    
    return(am)


# Plot figure
def line1plt(a,amax,amin):
    a2=(amax).rolling(lead_time=12, center=True).mean()
    a3=(amin).rolling(lead_time=12, center=True).mean()
    plt.plot(time1,a.rolling(lead_time=12, center=True).mean(),label=tl1,color='blue')
    plt.plot(time1,a*0,color='black')
    plt.fill_between(time1,a2,a3,color='cyan',alpha=.75)
#   plt.xlim(-15,15)
    plt.title(t1)
    plt.xlabel('Year')
    plt.ylabel('Temperature ($^o$C)')
    plt.legend()
    


# + Collapsed="false"
# extract ensemble of forecasts 
ncast=59
# anomaly relative to the mean of all forecasts
tmod_all=ddat.t_ref.mean(dim=('ensemble','init_date'))
tsurf=ddat.t_ref[ncast]-tmod_all
tmean=weighting(tsurf)
time1=np.arange(2019+11/12,2025,1/12)
time1=time1[0:60]

plt.figure(figsize=(10, 10))
# global
tl1='Mean'
t1='Global Surface Temperature'
gm=tmean.mean('ensemble')
gmax=tmean.max('ensemble')
gmin=tmean.min('ensemble')

line1plt(gm,gmax,gmin)

plt.savefig('predict_global.pdf',dpi=600)

#(tmean.mean('ensemble')).plot()
#tsurf.mean(dim=('lat','lon','ensemble')).plot()



# + [markdown] Collapsed="false"
# # Model - Observation Comparison 

# + Collapsed="false"
# Generic Plot Map script using subplot
# with additional hatching
def map2(a,b,y,x):
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
    
#    p1=ax.contourf(x, y, b, 3 ,cmap='Greys_r',hatches=['',"."],transform=dproj,alpha=0.0)
    
    bf=np.copy(b)
    blon,blat=np.meshgrid(b.lon,b.lat)
    plt.scatter(blon.flatten(),blat.flatten(),bf.flatten(),c='k',marker='.',transform=dproj)
    
#    
    plt.colorbar(p, shrink=.8, orientation='horizontal',label=cbtit)
    plt.title(t1)

    return


# + [markdown] Collapsed="false"
# ## extract the desired data from the zarr file and compute anomaly
#

# + Collapsed="false"
# compute the seasonal climatology based on selected years and remove it from the whole dataset
ds=dobs.sel( initial_time0_hours=slice('1971','2000'))
pclim,pseason,panom=climatology(ds,'initial_time0_hours')

oanom = dobs.TMP_GDS0_HTGL.groupby('initial_time0_hours'+'.month') - pclim.TMP_GDS0_HTGL
oclim=pclim.TMP_GDS0_HTGL
# oanom has the observed 3d anomaly field 
# oclim has the observed 2d climatology


# + Collapsed="false"
# forecast anomaly
ptemp1=ddat.t_ref
ptemp2=ptemp1 - ptemp1.mean(dim=('init_date','ensemble'))

# + Collapsed="false"
# Match the observations with the selected forecast 
# and compute the new time axis based on the JRA file
ncast=54
ncast=50
tobs=dobs.initial_time0_hours
to = ddat.init_date[ncast]
i=tobs['initial_time0_hours.year'] == to['init_date.year']
list=np.where(i)
index=list[0][10]  # November of the matching year

osurf= oanom[index:index+60,:,:]

# ptemp2 model and osurf is the obs


# + [markdown] Collapsed="false"
# ## Global Maps

# + Collapsed="false"
print(to)
plt.figure(figsize=(12,5 ))
cmap='RdBu_r'
cbtit='Observed Mean ($^o$C)'
cv=np.arange(-1.5,1.7,.25)
proj=ccrs.Robinson(central_longitude=180)

t1='Year 1: Observations'
ax = plt.subplot(1,2,1,projection=proj)
oflx1=osurf[0:12,:,:].mean(axis=(0))
map1(oflx1,oflx1.lat,oflx1.lon)

t1='Year 1: Forecasts'
cbtit='Ensemble Mean ($^o$C)'
ax = plt.subplot(1,2,2,projection=proj)
ptemp3= ptemp2[ncast,0:12,:,:,:].mean(dim=('lead_time'))
oflx=ptemp3.mean(dim=('ensemble'))
omax=ptemp3.max(dim=('ensemble'))
omin=ptemp3.min(dim=('ensemble'))
test1= oflx1 < omax
test2 = oflx1 > omin
test3=(test1 & test2)
test=np.invert(test3) 
map2(oflx,test,oflx.lat,oflx.lon)


plt.savefig('year1.pdf',dpi=600)



# + Collapsed="false"
from scipy.stats import pearsonr
a=np.copy(oflx)
b=np.copy(oflx1)
corr =pearsonr(a.flatten(),b.flatten())
print(corr) #corr=pearsonr(a,b)

# + Collapsed="false"
plt.figure(figsize=(12,5 ))
cmap='RdBu_r'
cbtit='Observed Mean ($^o$C)'
cv=np.arange(-1.5,1.7,.25)
proj=ccrs.Robinson(central_longitude=180)

t1='Years 2-5: Observations'
ax = plt.subplot(1,2,1,projection=proj)
oflx1=osurf[12:60,:,:].mean(axis=(0))
map1(oflx1,oflx1.lat,oflx1.lon)

t1='Years 2-5: Forecasts'
cbtit='Ensemble Mean ($^o$C)'
ax = plt.subplot(1,2,2,projection=proj)
ptemp3= ptemp2[ncast,12:60,:,:,:].mean(dim=('lead_time'))
oflx=ptemp3.mean(dim=('ensemble'))
omax=ptemp3.max(dim=('ensemble'))
omin=ptemp3.min(dim=('ensemble'))
test1= oflx1 < omax
test2 = oflx1 > omin
test3=(test1 & test2)
test=np.invert(test3)   #.load()
map2(oflx,test,oflx.lat,oflx.lon)

plt.savefig('year5.pdf',dpi=600)



# + Collapsed="false"
from scipy.stats import pearsonr
a=np.copy(oflx)
b=np.copy(oflx1)
corr =pearsonr(a.flatten(),b.flatten())
print(corr) #corr=pearsonr(a,b)


# + [markdown] Collapsed="false"
# ## Global Mean

# + Collapsed="false"
# Plot figure
def line2plt(a,amax,amin,b):
    a2=(amax).rolling(lead_time=12, center=True).mean()
    a3=(amin).rolling(lead_time=12, center=True).mean()
    plt.plot(time1,a.rolling(lead_time=12, center=True).mean(),label=tl1,color='blue')
    plt.plot(time1,b.rolling(initial_time0_hours=12, center=True).mean(),label=tl2,color='red')
    plt.plot(time1,a*0,color='black')
    plt.fill_between(time1,a2,a3,color='cyan',alpha=.75)
#   plt.xlim(-15,15)
    plt.title(t1)
    plt.xlabel('Year')
    plt.ylabel('Temperature ($^o$C)')
    plt.legend()
    


# + Collapsed="false"
# extract ensemble of forecasts 
print(ncast)
ptemp3= ptemp2[ncast,:,:,:,:]
tmean=weighting(ptemp3)
omean=weighting0(osurf).load()
yy=to['init_date.year']
time1=np.arange(yy+11/12,2025,1/12)
time1=time1[0:60]

plt.figure(figsize=(10, 10))
# global
tl1='Mean'
tl2='Observed'
t1='Global Surface Temperature'
gm=tmean.mean('ensemble')
gmax=tmean.max('ensemble')
gmin=tmean.min('ensemble')

line2plt(gm,gmax,gmin,omean)

plt.savefig('timeseries.pdf',dpi=600)

#(tmean.mean('ensemble')).plot()
#tsurf.mean(dim=('lat','lon','ensemble')).plot()



# + [markdown] Collapsed="true"
# # Original

# + Collapsed="false"
# unweighted
oanom.mean(dim=('lat','lon')).plot()
panom.TMP_GDS0_HTGL.mean(axis=(1,2)).plot()  #anomaly over the period of the climatology
# weighted
weights = np.cos(np.deg2rad(oanom.lat))
weights.name = "weights"

oair_weighted = oanom*(weights)

norm = (oanom*0+1)*weights
ntmp=norm.mean(('lat','lon'))
nn=ntmp[0]

omean = oair_weighted.mean(("lon", "lat"))/nn
omean.load()

omean.plot()

# + [markdown] Collapsed="true"
# # Junk

# + Collapsed="false"
# extract ensemble of forecasts 
# and compute the new time axis based on the JRA file
ncast=54
ncast=50
tobs=dobs.initial_time0_hours
to = ddat.init_date[ncast]
i=tobs['initial_time0_hours.year'] == to['init_date.year']
list=np.where(i)
index=list[0][10]  # November of the matching year

ttnew= dobs['initial_time0_hours'][index:index+60]

tmod_all=ddat.t_ref.mean(dim=('ensemble','init_date'))
tsurf=ddat.t_ref[ncast]-tmod_all


# + Collapsed="false"
#plt.plot(ttnew,ddat.t_surf[ncast,:,1,:,:].mean(dim=('lon','lat')))
#plt.plot(tobs[index:index+60],dobs.TMP_GDS0_HTGL[index:index+60,:,:].mean(axis=(1,2)))

tmp = oanom[index:index+60,:,:].load()
omeanr=omean[index:index+60].copy()
tmp1=tmp.expand_dims({'ensemble':10})
tmp3=tmp1.copy()

tmp2=tsurf[:,:,:,:]
aa=np.copy(tmp2)
A=np.einsum('jikl->ijkl',aa)
tmp3 = tmp3*0 + A

manom = (tmp3.groupby('initial_time0_hours'+'.month') - oclim).load()
manom=tmp3.load()

manom.mean(dim=('lat','lon','ensemble')).plot()
omeanr.plot()

# + Collapsed="false"
mair_weighted = manom*(weights)

norm = (manom*0+1)*weights
ntmp=norm.sum(dim=('lat','lon'))
nn=ntmp[0,0]

mmean = (mair_weighted.sum(dim=("lon", "lat"))/nn).load()

(mmean.mean('ensemble')).plot()
omeanr.plot()
manom.mean(dim=('lat','lon','ensemble')).plot()

# + Collapsed="false"
nr=12
omeanr.rolling(initial_time0_hours=nr, center=True).mean().plot(color='orange')
mmean.mean(dim=('ensemble')).rolling(initial_time0_hours=nr, center=True).mean().plot(color='black',linestyle='solid')
#manom.mean(dim=('lat','lon','ensemble')).rolling(initial_time0_hours=nr, center=True).mean().plot()
mmean.max(axis=0).rolling(initial_time0_hours=nr, center=True).mean().plot(color='cyan')
mmean.min(axis=0).rolling(initial_time0_hours=nr, center=True).mean().plot(color='cyan')



# + [markdown] Collapsed="true"
# # extra

# + Collapsed="false"

norm.mean(('lat','lon')).plot()

# + Collapsed="false"

