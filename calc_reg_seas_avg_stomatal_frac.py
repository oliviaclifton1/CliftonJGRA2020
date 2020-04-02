# oeclifton
# this script calculates the regional average stomatal fraction of 
# ozone dry deposition for a given season
# only includes grid cells where leaf area index (LAI)> 2 m2/m2

import numpy as np
from scipy.io import netcdf_file
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import MFDataset, Dataset

def geo_idx(dd, dd_array):
   """
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
    """
   geo_idx = (np.abs(dd_array - dd)).argmin()
   return geo_idx

def create_reg_mean_txtfile(mmvalues,filename,lat,lon):
    """
     create text file with regional averages
    """

    plotregions = ['northern mid-lats']
    reg_minlat = np.array([40.])
    reg_maxlat = np.array([55.])
    reg_minlon = np.array([0.])
    reg_maxlon = np.array([365])
    nregions = 1

    outF = open("./txt/"+filename,"w")
    for r in range(0,nregions):
       # get indices for lat/lon bounds for a given region
       reg_minlat_idx = geo_idx(reg_minlat[r],lat)
       reg_maxlat_idx = geo_idx(reg_maxlat[r],lat)
       reg_minlon_idx = geo_idx(reg_minlon[r],lon)
       reg_maxlon_idx = geo_idx(reg_maxlon[r],lon)
       data = np.nanmean(mmvalues[reg_minlat_idx:reg_maxlat_idx,reg_minlon_idx:reg_maxlon_idx])
       outF.write(plotregions[r]+" ")
       outF.write(str(round(data,3)))
       outF.write("\n")
    outF.close()

# simulation/season/variable info 
sim = "dynamic"
extension = '_newest_final'
section = 'drydep'
nyears = 10
year = '2010'
seas = 'JJA'
month1 = 5
month2 = 6
month3 = 7

# load ozone deposition velocity
slp_files = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+year+'_am3dd'+extension+'/pp/'+section+ '/ts/monthly/1yr/'+section+'.201*.O3_tot_con.nc')
slp1 = np.array(slp_files.variables['O3_tot_con'][:])[:,:,:] #get the i-th tracer from netcdf file
slp1 = np.squeeze(slp1[:,:,:])
# now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
slp2 = np.squeeze(slp1.reshape(-1,nyears,12,90,144))
slp = np.squeeze(slp2.mean(axis=0)) # average across years
lon = slp_files.variables['lon'][:] #grab model lat and lons
lat = slp_files.variables['lat'][:] #need last part to actually get the data
slp_files.close()

# load effective stomatal conductance 
slp_files = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+year+'_am3dd'+extension+'/pp/'+section+ '/ts/monthly/1yr/'+section+'.201*.O3_econ_stom.nc')
slp1 = np.array(slp_files.variables['O3_econ_stom'][:])[:,:,:] #get the i-th tracer from netcdf file·
slp1 = np.squeeze(slp1[:,:,:])
# now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
slp2 = np.squeeze(slp1.reshape(-1,nyears,12,90,144))
slp2 = np.squeeze(slp2.mean(axis=0)) # average across years
slp_files.close()

# calculate stomatal fraction
slp = np.divide(slp2,slp)
# average across months 
slp = np.mean(slp[[month1,month2,month3],:,:],axis=0)

# load LAI 
slp_files = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+year+'_am3dd'+extension+'/pp/land/ts/monthly/1yr/land.201*.LAI.nc')
slp1 = np.array(slp_files.variables['LAI'][:])[:,:,:] #get the i-th tracer from netcdf file·
slp1 = np.squeeze(slp1[:,:,:])
# now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
slp2 = np.squeeze(slp1.reshape(-1,nyears,12,90,144))
slp2 = np.squeeze(slp2.mean(axis=0)) # average across years
slp_files.close()
# average across months
slp2 = np.mean(slp2[[month1,month2,month3],:,:],axis=0)

# only keep stomatal fractions with less than 2 m2/m2 LAI 
slp[slp2<2] = np.nan

# write textfile
create_reg_mean_txtfile(slp,'reg_mean_stomatal_frac_lai_ge_2_at_2010s_jja.txt',lat,lon)

