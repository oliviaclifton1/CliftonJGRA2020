# oeclifton
# plot northern hemisphere maps of effective conductances for a given season
# left panel of plots are 2010s values
# right panel of plots are 2090s-2010s values

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import MFDataset, Dataset
from mpl_toolkits.basemap import shiftgrid
from mpl_toolkits.basemap import Basemap
import seaborn as sns

def create_reg_mean_txtfile(mmvalues,filename,lat,lon):
    """
     create text file with regional averages
    """

    plotregions = ['NE US', 'SE US', 'IMW US', \
               'C. Europe', 'S. Asia', 'E. Asia','midlatplusboreal']
    reg_minlat = np.array([36., 30., 36., 45.,25., 25.,40])
    reg_maxlat = np.array([46., 38., 46., 54.,40., 40.,65])
    reg_minlon = np.array([-80., -95.,-110., 0., 65,105.,0])
    reg_maxlon = np.array([-70., -80.,-100., 25.,90, 120.,360])
    nregions = 7

    reg_minlon[reg_minlon<0.] = reg_minlon[reg_minlon<0.]+360.
    reg_maxlon[reg_maxlon<0.] = reg_maxlon[reg_maxlon<0.]+360.

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

def geo_idx(dd, dd_array):
   """
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
    """
   geo_idx = (np.abs(dd_array - dd)).argmin()
   return geo_idx

# define some variables related to plotting
fontsize = 6
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

# defines some variables related to simulations 
sim = 'dynamic'
extension = ['_newest_final','_newest_final']
nyears =[10,10]

# pick season 
seas = 'DJF'
month1 = 0
month2 = 1
month3 = 11

if seas == 'DJF':
    variables = ['O3_econ_cu_dry','O3_econ_bk','O3_econ_cu_frz','O3_econ_g_dry']
    vars4label = ['dry and wet cuticular uptake','stem uptake','snow uptake to cuticles and ground','dry ground uptake']
    nvars = 4
    datamax = 0.6
elif seas == 'JJA' or seas == 'MAM':
    variables = ['O3_econ_cu_dry','O3_econ_cu_wet','O3_econ_stom','O3_econ_g_dry']
    vars4label = ['dry cuticular uptake','wet cuticular uptake','stomatal uptake','dry ground uptake']
    nvars = 4
    datamax = 0.6

# define output file name
filename = "northhemi_egx_"+seas+"_rcp85_"+sim+extension[0]+"_2010s_trend.pdf" #name output file·

# read in data 
mmvaluesb = np.empty([nvars,12,90,144])
mmvaluesf = np.empty([nvars,12,90,144])
for n in range(0,nvars):
    b = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd'+extension[0]+'/pp/drydep/ts/monthly/1yr/drydep.201*12.'+variables[n]+'.nc')
    # load all years
    mmvaluesb1 = np.array(b.variables[variables[n]][:])[:,:,:] #get the i-th tracer from netcdf file
    # now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
    mmvaluesb1 = np.squeeze(mmvaluesb1.reshape(-1,nyears[0],12,90,144))
    mmvaluesb[n,:,:,:] = np.squeeze(mmvaluesb1.mean(axis=0)) # average across years
    lon = b.variables['lon'][:] #grab model lat and lons
    lat = b.variables['lat'][:] #need last part to actually get the data
    b.close()
    f = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_2090_am3dd'+extension[1]+'/pp/drydep/ts/monthly/1yr/drydep.209*12.'+variables[n]+'.nc')
    # load all years
    mmvaluesf1 = np.array(f.variables[variables[n]][:])[:,:,:] #get the i-th tracer from netcdf file
    # now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
    mmvaluesf1 = np.squeeze(mmvaluesf1.reshape(-1,nyears[1],12,90,144))
    mmvaluesf[n,:,:,:] = np.squeeze(mmvaluesf1.mean(axis=0)) # average across years
    f.close()

# also include frozen ground uptake & wet cuticular uptake for DJF
if seas == 'DJF':
    b = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd'+extension[0]+'/pp/drydep/ts/monthly/1yr/drydep.201*12.O3_econ_g_frz.nc')
    mmvaluesb1 = np.array(b.variables['O3_econ_g_frz'][:])[:,:,:] #get the i-th tracer from netcdf file
    mmvaluesb1 = np.squeeze(mmvaluesb1.reshape(-1,nyears[0],12,90,144))
    mmvaluesb[2,:,:,:] = mmvaluesb[2,:,:,:]+np.squeeze(mmvaluesb1.mean(axis=0)) # average across years
    b.close()
    f = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_2090_am3dd'+extension[1]+'/pp/drydep/ts/monthly/1yr/drydep.209*12.O3_econ_g_frz.nc')
    mmvaluesf1 = np.array(f.variables['O3_econ_g_frz'][:])[:,:,:] #get the i-th tracer from netcdf file
    # now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
    mmvaluesf1 = np.squeeze(mmvaluesf1.reshape(-1,nyears[1],12,90,144))
    # add frozen ground uptake to frozen cuticular uptake 
    mmvaluesf[2,:,:,:] = mmvaluesf[2,:,:,:]+np.squeeze(mmvaluesf1.mean(axis=0)) # average across years
    f.close()
    b = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd'+extension[0]+'/pp/drydep/ts/monthly/1yr/drydep.201*12.O3_econ_cu_wet.nc')
    mmvaluesb1 = np.array(b.variables['O3_econ_cu_wet'][:])[:,:,:] #get the i-th tracer from netcdf file
    mmvaluesb1 = np.squeeze(mmvaluesb1.reshape(-1,nyears[0],12,90,144))
    mmvaluesb[0,:,:,:] = mmvaluesb[0,:,:,:]+np.squeeze(mmvaluesb1.mean(axis=0)) # average across years
    b.close()
    f = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_2090_am3dd'+extension[1]+'/pp/drydep/ts/monthly/1yr/drydep.209*12.O3_econ_cu_wet.nc')
    mmvaluesf1 = np.array(f.variables['O3_econ_cu_wet'][:])[:,:,:] #get the i-th tracer from netcdf file
    # now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
    mmvaluesf1 = np.squeeze(mmvaluesf1.reshape(-1,nyears[1],12,90,144))
    # add wet cuticular uptake to dry cuticular uptake 
    mmvaluesf[0,:,:,:] = mmvaluesf[0,:,:,:]+np.squeeze(mmvaluesf1.mean(axis=0)) # average across years
    f.close()

# need land mask (don't want to show values over the ocean)
staticfile = Dataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd'+extension[0]+'/pp/tracer_level/tracer_level.static.nc')
land = np.array(staticfile.variables['land_mask'][:])[:,:]
staticfile.close()

# average across season
mmvaluesb_seas = np.squeeze(np.mean(mmvaluesb[:,[month1,month2,month3],:,:],axis=1))
mmvaluesf_seas = np.squeeze(np.mean(mmvaluesf[:,[month1,month2,month3],:,:],axis=1))

# calculate difference between time periods
mmvaluesf_seas = mmvaluesf_seas-mmvaluesb_seas

# convert to desired units
mmvaluesb_seas = mmvaluesb_seas*100.
mmvaluesf_seas = mmvaluesf_seas*100.

# mask out values with less than 1/2 grid cell land
mmvaluesb_seas[:,land<0.5] = np.nan
mmvaluesf_seas[:,land<0.5] = np.nan

# print regional average snow uptake for 2010s
if seas == 'DJF':
   mmvalues = np.squeeze(mmvaluesb_seas[2,:,:])
   create_reg_mean_txtfile(mmvalues,'snow_uptake_2010s_DJF.txt',lat,lon)

#shift  to -180 to 180 for longitudes
mmvaluesb_seas_shift,lons_out = shiftgrid(178.75,mmvaluesb_seas,lon,start=False)
mmvaluesf_seas_shift,lons_out = shiftgrid(178.75,mmvaluesf_seas,lon,start=False)

# plot
fig=plt.figure()
# plot 2010s 
cy = [0,0.55,0.3,1]
cx = [0.1,0.1,0.1,0.1]
subplot_number = [1,3,5,7]
letter = ['(a)','(c)','(e)','(g)']
mycmap = "cubehelix_r"
for n in range(0, nvars):
   ax = plt.subplot(nvars,2,subplot_number[n]) 
   m = Basemap(projection='merc',resolution = 'c',area_thresh=10000,llcrnrlon=-180,urcrnrlon=180,llcrnrlat=22,urcrnrlat=65)
   x,y = m(*np.meshgrid(lons_out,lat))
   pc = m.pcolormesh(x,y,mmvaluesb_seas_shift[n,:,:],vmin=0,vmax=datamax,cmap=mycmap)
   plt.title(seas+' 2010s '+vars4label[n],pad=1.5,fontweight='bold')
   cbar = plt.colorbar(shrink=0.8,orientation='horizontal')
   cbar.outline.set_linewidth(0.5)
   cbar.set_label('effective conductance (cm s$^{-1}$)')
   cbar.ax.tick_params(width=0.5)
   if n != nvars-1:
       cbar.remove()
   m.drawcoastlines(linewidth=0.25)
   m.drawmapboundary(linewidth=0.5)
   m.drawmeridians(np.arange(-170.,170.,50.),labels=[0,0,0,1],linewidth=0.05,fontsize=5)
   m.drawparallels(np.arange(0,90.,20),labels=[1,0,0,0],linewidth=0.05,fontsize=5)
   xpt,ypt = m(-165,30)
   ax.text(xpt,ypt, letter[n])
   if n < 2:
       ax.set_anchor((cx[n],cy[n]))

# plot 2090s-2010s
cx = [0.6,0.6,0.6,0.6,0.6]
subplot_number = [2,4,6,8]
letter = ['(b)','(d)','(f)','(h)']
mycmap = sns.diverging_palette(240,10,as_cmap=True)
for n in range(0, nvars):
   ax = plt.subplot(nvars,2,subplot_number[n])
   m = Basemap(projection='merc',resolution = 'c',area_thresh=10000,llcrnrlon=-180,urcrnrlon=180,llcrnrlat=22,urcrnrlat=65)
   x,y = m(*np.meshgrid(lons_out,lat))
   pc = m.pcolormesh(x,y,mmvaluesf_seas_shift[n,:,:],vmin=-0.3,vmax=0.3,cmap=mycmap)
   plt.title(seas+' 2090s-2010s '+ vars4label[n],pad=1.5,fontweight='bold')
   cbar = plt.colorbar(shrink=0.8,orientation='horizontal')
   cbar.outline.set_linewidth(0.5)
   cbar.set_label('change in effective conductance (cm s$^{-1}$)')
   cbar.ax.tick_params(width=0.5)
   if n != nvars-1:
       cbar.remove()
   m.drawcoastlines(linewidth=0.25)
   m.drawmapboundary(linewidth=0.5)
   m.drawmeridians(np.arange(-170.,170.,50.),labels=[0,0,0,1],linewidth=0.05,fontsize=5)
   m.drawparallels(np.arange(0,90.,20),labels=[1,0,0,0],linewidth=0.05,fontsize=5)
   xpt,ypt = m(-165,30)
   ax.text(xpt,ypt, letter[n])
   if n < 2:
       ax.set_anchor((cx[n],cy[n]))

plt.subplots_adjust(wspace=0.1,hspace=0)
fig.savefig(filename) #save figure to pdf·

