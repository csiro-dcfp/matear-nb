# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Calculate climate indices on Zarr files
#
#
#

# +
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


# -

def climatology(dsx,TIME1):
    clim = dsx.groupby(TIME1+'.month').mean(dim=TIME1)
    anom = dsx.groupby(TIME1+'.month') - clim
    season=dsx.groupby(TIME1+'.season').mean(dim=TIME1)
    return(clim,season,anom)


# +
import sys

print("Starting")
try:
    get_ipython
    syear='2015'
    ci='nino4'
except:
    print( 'Number of arguments:', len(sys.argv), 'arguments.')
    print( 'Argument List:', str(sys.argv) )
    syear=(sys.argv[1])
    ci=(sys.argv[2])
    print(syear)
    quit

print(syear,ci)



# +

bdir='/g/data/v14/vxk563/CAFE/forecasts/f6/WIP/'
dyear='c5-d60-pX-f6-'+syear+'1101/ZARR/'
file1='ocean_month.zarr.zip'
file2='ocean_bgc_month.zarr.zip'


docn = xr.open_zarr(bdir+dyear+file1,consolidated=True)
#dbgc = xr.open_zarr(dir_fcst+file2,consolidated=True)
# -

# ## Observations

# +
fobs='/g/data/v14/observations/sst/sst.mnmean.v4.nc'
fobs='/g/data/v14/rxm599/obs/sst.mnmean.nc'
fobs='/g/data/v14/rxm599/obs/sst_mask.nc'
dobs=xr.open_dataset(fobs)
d1=dobs.rename(lat='yt_ocean',lon='xt_ocean')
d2=d1.assign_coords(xt_ocean=d1.xt_ocean-360)

# compute the climatology and anomaly
oclim,oseason,oanom=climatology(d2,'time')

# -


# ## F5

f5='/g/data/v14/vxk563/CAFE/forecasts/f5/WIP/ZARR/'
df5=xr.open_zarr(f5+file1)


# ## Climate Index

# use weighting to compute spatial mean
def v_weight(var,weights):
    weights.name = "weights"
    v_weighted = var.weighted(weights)
    w_mean = v_weighted.mean(("xt_ocean", "yt_ocean"))
    return(w_mean)


# +
# convert index calculation to a procedure

def CI():
    tmp = docn.sel(yt_ocean=slice(y2,y1),xt_ocean=slice(x1,x2))
    ftmp= df5.sel(yt_ocean=slice(y2,y1),xt_ocean=slice(x1,x2))
    otmp=d2.sel(yt_ocean=slice(y1,y2),xt_ocean=slice(x1,x2))
# my new mask file fixed the lat direction
    otmp=d2.sel(yt_ocean=slice(y2,y1),xt_ocean=slice(x1,x2))
# observations raw, clim, anomaly
    otmp1=oclim.sel(yt_ocean=slice(y2,y1),xt_ocean=slice(x1,x2))
    otmp2=oanom.sel(yt_ocean=slice(y2,y1),xt_ocean=slice(x1,x2))

    oweights = np.cos(np.deg2rad(d2.sst.yt_ocean))
    ovar=otmp.sst
    onino3r=v_weight(ovar,oweights)
    ovar=otmp1.sst
    onino3c=v_weight(ovar,oweights)
    ovar=otmp2.sst
    onino3=v_weight(ovar,oweights)
        
# f6            
    weights=tmp.area_t
    var=tmp.sst
    nino3r=v_weight(var,weights)
# f5
    var=ftmp.sst
    nino3fr=v_weight(var,weights)
# corrected with observations
    nino3=nino3r.groupby('time.month') - onino3c
# corrected with F5
    corr= nino3fr.mean(('init_date','ensemble'))
    nino3r=nino3r.rename(time='lead_time')
    nino3r=nino3r.assign_coords(lead_time=corr.lead_time)
    acorr= nino3r-corr
    return(nino3,acorr,onino3)


# +
# set reg/y=0:-10/x=90W:80W  ! nino1.2

if ci == 'tcoral' :
# #!T coral sea  coral= sst[k=1,x=142:166@ave,y=-18:-10@ave]
    y1=-10; y2=-18; x1=-216; x2=-194
elif ci == 'nino34' :    
# nino3.4  #set reg/y=-5:5/x=170W:120W  
    y1=5; y2=-5; x1=-170; x2= -120
elif ci == 'nino4' :
# set reg/y=5N:5S/x=160E:150W !nino4
    y1=5; y2=-5; x1=-200; x2= -150
elif ci == 'nino3' :
# nino3 calculation for forecasts and observations
    y1=5; y2=-5; x1=-150; x2=-90
else :
    y1=0; y2=-1; x1=-200; x2= -150


# +
nino3,acorr,onino3=CI()
nino3.load()  # f6 corrected with obs
acorr.load()  # f6 corrected with f5
onino3.load() # obs

acorr.to_netcdf(ci+'_'+syear+'.nc')

# +
# find observation time that matches forecasts
to=docn.time[0]
fyear=to['time.year']
fmonth=to['time.month']
i=dobs.time['time.year']== fyear
list=np.where(i)
if list[0].size < 10: 
    index=441
else:
    index=list[0][10] 


plt.plot(onino3.time[index:],onino3[index:])


# +
# Plot figure
def line2plt(a,amax,amin,b):
    a2=(amax).rolling(time=1, center=True).mean()
    a3=(amin).rolling(time=1, center=True).mean()
    plt.plot(time1,a.rolling(time=1, center=True).mean(),label=tl1,color='blue')
    plt.plot(time2,b.rolling(time=1, center=True).mean(),label=tl2,color='red')
    plt.plot(time1,a*0,color='black')
    plt.fill_between(time1,a2,a3,color='cyan',alpha=.75)
#   plt.xlim(-15,15)
    plt.title(t1)
    plt.xlabel('Year')
    plt.ylabel('Temperature ($^o$C)')
    plt.legend()
    
# Plot figure
def line3plt(a,amax,amin,b,bmax,bmin):
    a2=(amax).rolling(lead_time=1, center=True).mean()
    a3=(amin).rolling(lead_time=1, center=True).mean()
    b2=(bmax).rolling(lead_time=1, center=True).mean()
    b3=(bmin).rolling(lead_time=1, center=True).mean()
    plt.plot(time1,a.rolling(lead_time=1, center=True).mean(),label=tl1,color='blue')
    plt.plot(time1,b.rolling(lead_time=1, center=True).mean(),label=tl2,color='red')
    plt.plot(time1,a*0,color='black')
    plt.fill_between(time1,a2,a3,color='cyan',alpha=.95)
    plt.fill_between(time1,b2,b3,color='orange',alpha=.95)
    
#   plt.xlim(-15,15)
    plt.title(t1)
    plt.xlabel('Lead Time (months)')
    plt.ylabel('Index ($^o$C)')
    plt.legend()


# +
time1=acorr.lead_time
time2=np.arange(onino3[index:].size)
plt.figure(figsize=(10, 10))
# global
tl1='Mean'
tl2='F5 mean'
t1='F5 Corrected: '+ci+' '+str(np.int(fyear))+'/'+str(np.int(fmonth))
gm=acorr.mean('ensemble')
gmax=acorr.max('ensemble')
gmin=acorr.min('ensemble')

b=acorr[:,0:10].mean('ensemble')
bmax=acorr[:,0:10].max('ensemble')
bmin=acorr[:,0:10].min('ensemble')


line3plt(gm,gmax,gmin,b,bmax,bmin)
plt.plot(onino3[index:],color='black',linewidth=4,label='Obs')

plt.savefig(ci+'_'+syear+'.pdf',dpi=600)
# +
# january pdf
a=acorr[3:120:12,:]

plt.figure(figsize=(10, 10))
plt.hist(np.ravel(a),bins=30,label='f6')
plt.hist(np.ravel(a[:,0:10]),color='orange',label='f5')

plt.hist(onino3[0:464:12],color='green',label='Obs') 
plt.title(ci+'-'+syear+': Jan')
plt.xlabel('Index ')
plt.ylabel('Number')

plt.text(-3,80,'Samples='+str(a.size))
plt.legend()
plt.savefig(ci+'_'+syear+'_hist.pdf',dpi=600)
# -




# +
time1=acorr.lead_time
time2=np.arange(onino3[index:].size)
plt.figure(figsize=(10, 10))
# global
tl1='Mean'
tl2='Observed'
t1='Corrected with Observations'
gm=nino3.mean('ensemble')
gmax=nino3.max('ensemble')
gmin=nino3.min('ensemble')

omean=onino3[index:]

#line2plt(gm,gmax,gmin,omean)
# -



