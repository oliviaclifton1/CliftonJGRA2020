# oeclifton
# this script plots map of correlation coefficients for daily values between 1) surface ozone & ozone deposition velocity
# and 2) effective stomatal conductance & effective wet cuticular conductance
# also creates netcdf files for maps of correlation coefficients 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,shiftgrid
from matplotlib import ticker
from netCDF4 import MFDataset, Dataset

# set some variables related to plotting
fontsize = 7
params = {
    'axes.labelsize': fontsize, # fontsize for x and y labels (was 10)
    'axes.titlesize': fontsize,
    'font.size': fontsize, # was 10
    'legend.fontsize': fontsize, # was 10
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
}
matplotlib.rcParams.update(params)

# set some variables related to colorbar
cmap = plt.cm.get_cmap("PuBuGn")
cmap.set_under('red')
met4plotlabel = 'correlation coefficient'

# select season 
seas = 'JJA'

if seas == 'JJA':
    days_of_year = np.arange(1,366,1)
    ind = np.logical_and(days_of_year > 151, days_of_year < 244)
    ndays = 31+31+30
elif seas == 'DJF':
    days_of_year = np.arange(1,366,1)
    ind = np.logical_or(days_of_year > 334, days_of_year < 60)
    ndays = 31+31+28
elif seas == 'July':
    days_of_year = np.arange(1,366,1)
    ind = np.logical_and(days_of_year > 181, days_of_year < 213)
    ndays = 31

# set some variables related to simulation and plotting 
year = '2010'
extension = year+'_am3dd_newest_final'
nyears = 10
title = [(seas+' '+year+'s'+ r' -r(v$_d$,surface ozone)'),(seas+' '+year+'s'+ ' r(stomatal uptake,wet-cuticular uptake)')]
nrows = 4
ncols = 1

#define min and max for plot
datamin = 0
datamax = 1

# need land mask (don't want to show correlations over the ocean)
staticfile = Dataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+extension+'/pp/tracer_level/tracer_level.static.nc')
land = np.array(staticfile.variables['land_mask'][:])[:,:]
staticfile.close()

# load variables for correlation between vd and surface ozone 
met = 'O3_tot_con'
metclass = 'drydep'
f = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+extension+'/pp/'+metclass+'/ts/daily/1yr/'+metclass+'.201*.'+met+'.nc')
temp = np.array(f.variables[met][:])[:,:,:] #get the i-th tracer from netcdf file 
mmvaluesf = np.squeeze(temp.reshape(-1,nyears,365,90,144)) # adds an extra dimension (not sure why)
lon = f.variables['lon'][:] #grab model lat and lons
lat = f.variables['lat'][:] #need last part to actually get the data
f.close()
mmvaluesf = mmvaluesf[:,ind,:,:]

met = 'O3'
metclass = 'atmos_level'
b = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+extension+'/pp/'+metclass+'/ts/daily/1yr/'+metclass+'.201*.'+met+'.nc')
temp = np.squeeze(np.array(b.variables[met][:])[:,47,:,:]) #get the i-th tracer from netcdf file·
mmvaluesb = np.squeeze(temp.reshape(-1,nyears,365,90,144)) # adds an extra dimension (not sure why)
b.close()
mmvaluesb = mmvaluesb[:,ind,:,:]

# get rid of oceanic vals 
mmvaluesf[:,:,land<0.5]=np.nan
mmvaluesb[:,:,land<0.5]=np.nan

# do correlation 
mmvalues = np.empty([90,144])
for i in range(0,144):
    for j in range(0,90):
        temp = np.corrcoef(np.ravel(np.squeeze(mmvaluesb[:,:,j,i])),np.ravel(np.squeeze(mmvaluesf[:,:,j,i])))
        # multiply correlation by -1 b/c only show colorbar for postive values 
        mmvalues[j,i] = np.multiply(-1.,temp[0,1])

# plotting 
fig, axes = plt.subplots(nrows,ncols) #open the figure·
filename = ("map_northhemi_correls_paper_"+seas+"_"+year+"s.pdf") #name output file·
m = Basemap(projection='merc',resolution = 'c',area_thresh=10000, \
    llcrnrlon=-180,urcrnrlon=180,llcrnrlat=22,urcrnrlat=65,ax=axes[0])
m.drawcoastlines(linewidth=0.25)
m.drawmeridians(np.arange(-170.,170.,50.),labels=[0,0,0,1],linewidth=0.05,fontsize=fontsize)
m.drawparallels(np.arange(0,90.,20),labels=[1,0,0,0],linewidth=0.05,fontsize=fontsize)
mmvalues_shift,lons_out = shiftgrid(178.75,mmvalues,lon,start=False)
x, y =  m(*np.meshgrid(lons_out,lat)) # compute map proj coordinates.
pc = m.pcolormesh(x,y,mmvalues_shift,vmin=datamin, vmax=datamax,cmap=cmap,ax=axes[0])
# Draw a border around the whole
m.drawmapboundary(color='k', linewidth=0.5)
axes[0].set_title(title[0],pad=1.5,fontweight='bold')
xpt,ypt = m(-165,30)
axes[0].text(xpt,ypt, '(g)')

# save correlation info to netcdf 
# need to rename lat and lon
lat1 = lat
lon1 = lon
# multiply correlation by -1 
mmvalues = np.multiply(-1.,mmvalues)
# write a netcdf file with correlation
root_grp = Dataset('/home/oeclifton/python/AM3DD/correl_sfc_o3_o3ddv.jja.nh.2010s.nc', 'w', format='NETCDF3_CLASSIC')
root_grp.description = 'Correlation (r) between surface ozone and ozone deposition velocity on daily basis during June-August 2010s '
# dimensions
root_grp.createDimension('lat', 90)
root_grp.createDimension('lon', 144)
# variables
lon = root_grp.createVariable('lon', 'f4', ('lon',))
lat = root_grp.createVariable('lat', 'f4', ('lat',))
field = root_grp.createVariable('r', 'f4', ('lat', 'lon',))
field.units = 'unitless'
lat.cartesian_axis = "Y"
lon.units =  "degrees_E"
lon.long_name = "longitude"
lon.cartesian_axis = "X"
lon[:] = lon1
lat[:] = lat1
field[:] = mmvalues
root_grp.close()

# load variables for correlation between stomatal uptake and cuticular uptake 
met = 'O3_econ_stom'
metclass = 'drydep'
f = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+extension+'/pp/'+metclass+'/ts/daily/1yr/'+metclass+'.201*.'+met+'.nc')
temp = np.array(f.variables[met][:])[:,:,:] #get the i-th tracer from netcdf file·
mmvaluesf = np.squeeze(temp.reshape(-1,nyears,365,90,144)) # adds an extra dimension (not sure why)
f.close()
mmvalues_st = mmvaluesf[:,ind,:,:]

met = 'O3_econ_cu_wet'
metclass = 'drydep'
f = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+extension+'/pp/'+metclass+'/ts/daily/1yr/'+metclass+'.201*.'+met+'.nc')
temp = np.array(f.variables[met][:])[:,:,:] #get the i-th tracer from netcdf file·
mmvaluesf = np.squeeze(temp.reshape(-1,nyears,365,90,144)) # adds an extra dimension (not sure why)
f.close()
mmvalues_cu = mmvaluesf[:,ind,:,:]

# get rid of oceanic vals·
mmvalues_st[:,:,land<0.5]=np.nan
mmvalues_cu[:,:,land<0.5]=np.nan

# calculate correlation coefficient 
mmvalues = np.empty([90,144])
for i in range(0,144):
    for j in range(0,90):
        temp = np.corrcoef(np.ravel(np.squeeze(mmvalues_st[:,:,j,i])),np.ravel(np.squeeze(mmvalues_cu[:,:,j,i])))
        mmvalues[j,i] = temp[0,1]

# plotting·
m = Basemap(projection='merc',resolution = 'c',area_thresh=10000,llcrnrlon=-180,urcrnrlon=180,llcrnrlat=22,urcrnrlat=65,ax=axes[1])
m.drawcoastlines(linewidth=0.25)
m.drawmeridians(np.arange(-170.,170.,50.),labels=[0,0,0,1],linewidth=0.05,fontsize=fontsize)
m.drawparallels(np.arange(0,90.,20),labels=[1,0,0,0],linewidth=0.05,fontsize=fontsize)
mmvalues_shift,lons_out = shiftgrid(178.75,mmvalues,lon1,start=False)
x, y =  m(*np.meshgrid(lons_out,lat1)) # compute map proj coordinates.
pc = m.pcolormesh(x,y,mmvalues_shift,vmin=datamin, vmax=datamax,cmap=cmap,ax=axes[1])
# Draw a border around the whole
m.drawmapboundary(color='k', linewidth=0.5)
axes[1].set_title(title[1],pad=1.5,fontweight='bold')
xpt,ypt = m(-165,30)
axes[1].text(xpt,ypt, '(h)')

# adjust figure a bit
plt.subplots_adjust(wspace=0.1,hspace=0.3)
axes[2].set_visible(False)
axes[3].set_visible(False)

# add colorbar
cax = fig.add_axes([0.25,0.465,0.5, 0.015])
cbar = fig.colorbar(pc,cax=cax,pad=0.1,orientation='horizontal')
cbar.ax.tick_params(labelsize=fontsize,pad=1.2,width=0.25)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.set_xlabel(met4plotlabel,labelpad=1)
cbar.outline.set_linewidth(0.25)
fig.savefig(filename) #save figure to pdf

# write a netcdf file with correlation
root_grp = Dataset('/home/oeclifton/python/AM3DD/correl_O3_econ_stom_O3_econ_cu_wet.jja.nh.2010s.nc', 'w', format='NETCDF3_CLASSIC')
root_grp.description = 'Correlation (r) between surface ozone (O3_econ_stom) and wet cuticular uptake (O3_econ_cu_wet) on daily basis during June-August 2010s '
# dimensions
root_grp.createDimension('lat', 90)
root_grp.createDimension('lon', 144)
# variables
lon = root_grp.createVariable('lon', 'f4', ('lon',))
lat = root_grp.createVariable('lat', 'f4', ('lat',))
field = root_grp.createVariable('r', 'f4', ('lat', 'lon',))
field.units = 'unitless'
lat.cartesian_axis = "Y"
lon.units =  "degrees_E"
lon.long_name = "longitude"
lon.cartesian_axis = "X"
lon[:] = lon1
lat[:] = lat1
field[:] = mmvalues
root_grp.close()
