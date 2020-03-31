# oeclifton, 3/3/2020
# plot ozone vertical profile difference between two sims
# land only; one region

import numpy as np
from scipy.io import netcdf_file
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import MFDataset, Dataset
from mpl_toolkits.basemap import shiftgrid
from mpl_toolkits.basemap import Basemap
import seaborn as sns


def geo_idx(dd, dd_array):
   """
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
    """
   geo_idx = (np.abs(dd_array - dd)).argmin()
   return geo_idx


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


# define season
seas = 'DJF'
month1 = 0
month2 = 1
month3 = 11

# define time period
time = '2010'

# define output file name
filename = 'northhemi_profile_diff_xactive_static_ozone_'+seas+'_rcp85_'+time+'.pdf' #name output file·

# load xactive ozone
temp_files = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+time+'_am3dd_newest_final/pp/tracer_level/ts/monthly/1yr/tracer_level.201*.O3.nc')
temp1 = np.squeeze(np.array(temp_files.variables['O3'][:])[::,:,:]) #get the i-th tracer from netcdf file·
# now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
temp2 = np.squeeze(temp1.reshape(-1,10,12,48,90,144))
temp = np.squeeze(temp2.mean(axis=0)) # average across years
lon = temp_files.variables['lon'][:] #grab model lat and lons
lat = temp_files.variables['lat'][:] #need last part to actually get the data
pressure = temp_files.variables['pfull'][:]
temp_files.close()
# average across season
xactive_seas = np.squeeze(np.mean(temp[[month1,month2,month3],:,:,:],axis=0))


# load static ozone
temp_files = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+time+'_am3dd_static_o3dd_newest_final/pp/tracer_level/ts/monthly/1yr/tracer_level.*.O3.nc')
temp1 = np.squeeze(np.array(temp_files.variables['O3'][:])[:,:,:,:]) #get the i-th tracer from netcdf file·
# now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
temp2 = np.squeeze(temp1.reshape(-1,10,12,48,90,144))
temp = np.squeeze(temp2.mean(axis=0)) # average across years
temp_files.close()
# average across season
static_seas = np.squeeze(np.mean(temp[[month1,month2,month3],:,:,:],axis=0))

# difference
diff = xactive_seas-static_seas
# convert to desired units
diff = diff*1e9

# need land mask (don't want to show correlations over the ocean)
staticfile = Dataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd_new/pp/tracer_level/tracer_level.static.nc')
land = np.array(staticfile.variables['land_mask'][:])[:,:]
staticfile.close()

# mask out values with less than 1/2 grid cell land
diff[:,land<0.5] = np.nan

# average over northern hemisphere (mid lats + boreal only)
reg_minlat = np.array([40])
reg_maxlat = np.array([65])
reg_minlon = np.array([0])
reg_maxlon = np.array([360])
reg_minlon[reg_minlon<0.] = reg_minlon[reg_minlon<0.]+360.
reg_maxlon[reg_maxlon<0.] = reg_maxlon[reg_maxlon<0.]+360.
reg_minlat_idx = geo_idx(reg_minlat,lat)
reg_maxlat_idx = geo_idx(reg_maxlat,lat)
reg_minlon_idx = geo_idx(reg_minlon,lon)
reg_maxlon_idx = geo_idx(reg_maxlon,lon)
data = np.nanmean(diff[:,reg_minlat_idx:reg_maxlat_idx,reg_minlon_idx:reg_maxlon_idx],axis=(1,2))

# plot ozone vertical profile
fig=plt.figure()
plt.plot(data,pressure,color='black')
plt.ylabel('pressure (hPa')
plt.xlabel('ozone (ppb)')
plt.title(seas+' xactive-static')
plt.xlim(-15,0)
plt.ylim(400,1000)
# invert y-axis
plt.gca().invert_yaxis()
fig.savefig(filename) #save figure to pdf·
