# Plot 3 maps of the ofam output

def map3(z1,z2,x1,y1,cfig,title,title1,title2,cval1,cval2,cval3,pal1,pal2):
#	 matplotlib inline
#	 
	 import numpy as np
	 import matplotlib.pyplot as plt
	 from scipy.io import netcdf
	 import scipy.stats as stats
	 import sys
	 import netCDF4 as nc
	 import xarray as xr
#
	 import cartopy.crs as ccrs
	 from cartopy import config
	 from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

	 plt.figure(figsize=(14, 12))


	 ax = plt.subplot(2, 1, 1,projection=ccrs.PlateCarree())
	 p1=plt.contourf(x1,y1,z1,cval1,extend='both',cmap=pal1,
			   transform=ccrs.PlateCarree())

	 ax.set_xticks([120, 130, 140, 150, 160, 170, 180], crs=ccrs.PlateCarree())
	 ax.set_yticks([-15, -10, -5, 0, 5, 10, 15], crs=ccrs.PlateCarree())
	 ax.coastlines()

	 lon_formatter = LongitudeFormatter(zero_direction_label=True)
	 lat_formatter = LatitudeFormatter()
	 ax.xaxis.set_major_formatter(lon_formatter)
	 ax.yaxis.set_major_formatter(lat_formatter)
#	 cbar = plt.colorbar(p1, pad=0.25, extend='both')
#	 cbar.ax.set_ylabel('Eady Growth Rate (days) ', fontsize=14) 
	 cbar=plt.colorbar(p1,label=title1 )

	 xx=[128,128,138,138,128]
	 yy=[5,10,10,5,5]
	 plt.plot(xx,yy,'k')
	 plt.title(title)
	 plt.xlabel('Longitude')
	 plt.ylabel('Latitude') 

# second panel
	 ax2 = plt.subplot(2, 1, 2,projection=ccrs.PlateCarree())
	 p1=plt.contourf(x1,y1,z2,cval2,extend='both',cmap=pal2,
			   transform=ccrs.PlateCarree())
	 ax2.set_xticks([120, 130, 140, 150, 160, 170, 180], crs=ccrs.PlateCarree())
	 ax2.set_yticks([-15, -10, -5, 0, 5, 10, 15], crs=ccrs.PlateCarree())

	 ax2.coastlines()
	 plt.contour(x1,y1,z2,[0.],
			   transform=ccrs.PlateCarree())

	 lon_formatter = LongitudeFormatter(zero_direction_label=True)
	 lat_formatter = LatitudeFormatter()
	 ax2.xaxis.set_major_formatter(lon_formatter)
	 ax2.yaxis.set_major_formatter(lat_formatter)
	 cbar=plt.colorbar(p1,label=title2)

	 xx=[128,128,138,138,128]
	 yy=[5,10,10,5,5]
	 plt.plot(xx,yy,'k')


	 plt.title(title)
	 plt.xlabel('Longitude')
	 plt.ylabel('Latitude') 

	 plt.savefig(cfig+'a.pdf',dpi=600)


	 plt.figure(figsize=(9, 12))

	 clev = cval3
	 rdsst=z2[200:251,79:181]
	 ax = plt.subplot(2, 1, 1,projection=ccrs.PlateCarree())
	 p1=plt.contourf(x1[79:181],y1[200:251],rdsst,clev,
			   extend='both',cmap=pal2,
			   transform=ccrs.PlateCarree())
	 ax.set_xticks([128, 130, 132, 134, 136, 138], crs=ccrs.PlateCarree())
	 ax.set_yticks([5, 6, 7, 8, 9, 10], crs=ccrs.PlateCarree())
	 ax.coastlines()
	 plt.contour(x1[79:181],y1[200:251],rdsst,[0.],
			   transform=ccrs.PlateCarree())

	 lon_formatter = LongitudeFormatter(zero_direction_label=True)
	 lat_formatter = LatitudeFormatter()
	 ax.xaxis.set_major_formatter(lon_formatter)
	 ax.yaxis.set_major_formatter(lat_formatter)
	 cbar=plt.colorbar(p1,label=title2,orientation='vertical',shrink=.5)
	 plt.title(title)
	 plt.xlabel('Longitude')
	 plt.ylabel('Latitude') 
	 plt.savefig(cfig+'b.pdf',dpi=600)


