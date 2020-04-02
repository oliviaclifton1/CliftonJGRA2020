# oeclifton
# this script calculates the 2010s to 2090s change in a seasonally averaged variable
# and plots map for northern hemisphere; plots land only values 

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import MFDataset, Dataset
from mpl_toolkits.basemap import shiftgrid
from mpl_toolkits.basemap import Basemap
import seaborn as sns

# define some variables for plotting 
fontsize = 7
# inspired by http://nipunbatra.github.io/2014/08/latexify/
params = {
    'axes.labelsize': fontsize, # fontsize for x and y labels (was 10)
    'axes.titlesize': fontsize,
    'font.size': fontsize, # was 10
    'legend.fontsize': fontsize, # was 10
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
}
matplotlib.rcParams.update(params)

# define colorbar
mycmap = sns.diverging_palette(240,10,as_cmap=True)

# define season
seas = 'DJF'
month1 = 0
month2 = 1
month3 = 11

# define time period
time = ['2010','2090']
extension = ['_newest_final','_newest_final']
nyears = [10,10]
time10years = ['201','209']

# choose variable 
var = 'LAI'
section = 'land'
factor = 1.
letter = ''
var4label = r'$\Delta$ LAI'
units4label = r'(m$^2$ m$^{-2}$)'
vmin = -5
vmax = 5

# define output file name
filename = 'northhemi_trend_xactive_'+seas+'_'+var+'_rcp85.pdf' #name output file·

# load data 
xactive_seas = np.empty([90,144,2])
for t in range(0, 2):
    # load xactive ozone
    temp_files = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+time[t]+'_am3dd'+ \
        extension[t]+'/pp/'+section+'/ts/monthly/1yr/'+section+'.'+time10years[t]+'*.'+var+'.nc')
    if var == 'sphum' or var == 'O3':
        temp1 = np.squeeze(np.array(temp_files.variables[var][:])[:,47,:,:])
    else:
        temp1 = np.squeeze(np.array(temp_files.variables[var][:])[:,:,:]) #get the i-th tracer from netcdf file·
    # now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
    temp_o3 = np.squeeze(temp1.reshape(-1,nyears[t],12,90,144))
    lon = temp_files.variables['lon'][:] #grab model lat and lons
    lat = temp_files.variables['lat'][:] #need last part to actually get the data
    temp_files.close()

    # average across season
    xactive_seas[:,:,t] = np.squeeze(np.mean(temp_o3[:,[month1,month2,month3],:,:],axis=(0,1)))

# calc change from 2010s to 2090s 
trend_xactive = (np.squeeze(xactive_seas[:,:,1]-xactive_seas[:,:,0]))*factor

# need land mask for static file·
staticfile = Dataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd_newest_final/pp/tracer_level/tracer_level.static.nc')
land = np.array(staticfile.variables['land_mask'][:])[:,:]
staticfile.close()

# remove oceanic values 
trend_xactive[land<0.5] = np.nan

# plot northern hemisphere map of seasonal changes from 2010s to 2090s
fig=plt.figure()
ax = plt.subplot(5,1,1)
m = Basemap(projection='merc',resolution = 'l',area_thresh=10000,llcrnrlon=-180,urcrnrlon=180,llcrnrlat=22,urcrnrlat=65)
title = '2090s-2010s'
temp,lons_out = shiftgrid(178.75,trend_xactive,lon,start=False)
x,y = m(*np.meshgrid(lons_out,lat))
pc = m.pcolormesh(x,y,temp,vmin=vmin,vmax=vmax,cmap=mycmap)
plt.title(seas+' '+title,pad=1.5,fontweight='bold')
m.drawcoastlines(linewidth=0.25)
m.drawmapboundary(linewidth=0.5)
m.drawmeridians(np.arange(-170.,170.,50.),labels=[0,0,0,1],linewidth=0.05,fontsize=fontsize)
m.drawparallels(np.arange(0,90.,20),labels=[1,0,0,0],linewidth=0.05,fontsize=fontsize)
xpt,ypt = m(-165,30)
ax.text(xpt,ypt, letter)
plt.subplots_adjust(wspace=0.1,hspace=0)
# add colorbar.
cax = fig.add_axes([0.25,0.68,0.5, 0.015])
cb = fig.colorbar(pc,cax=cax,pad=0.1,  orientation='horizontal')
cb.ax.set_xlabel(var4label+' '+units4label,labelpad=1)
cb.ax.tick_params(direction='out',pad=1.4,width=0.4)
cb.outline.set_linewidth(0.4)
fig.savefig(filename) #save figure to pdf·
