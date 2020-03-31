# oeclifton
# plots the fractional variance explained by each depositional pathway in daily variations in
# ozone deposition velocity across the northern hemisphere for a given season

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,shiftgrid
from matplotlib import ticker
from netCDF4 import MFDataset, Dataset
from matplotlib.colors import ListedColormap
import seaborn as sns

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

    plotregions = ['NE US', 'SE US', 'IMW US', \
               'C. Europe', 'C. Asia', 'E. Asia' ]
    reg_minlat = np.array([36., 30., 36., 45.,25., 25.])
    reg_maxlat = np.array([46., 38., 46., 54.,40., 40.])
    reg_minlon = np.array([-80., -95.,-110., 0., 65,105.])
    reg_maxlon = np.array([-70., -80.,-100., 25.,90, 120.])
    latlon_reg = ['[36-46ºN,70-80ºW]','[30-38ºN,80-95ºW]', \
    '[36-46ºN,100-110ºW]', \
    '[45-54ºN,0-25ºE]', \
    '[25-40ºN,65-90ºE]',\
    '[25-40ºN,105-120ºE]']
    nregions = 6

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

# define some variables related to plotting
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

# pick colorbar
mycmap = ListedColormap(sns.color_palette("Blues",n_colors=20))

#select dep pathways
variable = ['O3_econ_cu_dry','O3_econ_cu_wet','O3_econ_stom','O3_econ_g_dry']
variable4plotting = ['dry-cuticular uptake','wet-cuticular uptake','stomatal uptake','ground uptake']
npathways = 4

# define time period & simulation info 
year = '2010'
nyears = 10
extension = '_newest_final'
sim = "dynamic"

#define month(s) to average over & season name
seas = 'JJA'
startday = 152-1
endday = 244-1
ndays = endday-startday+1

#define min and max for plot
vmin = 0
vmax = 1

# land mask (don't want to plot oceanic values)
staticfile = Dataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd'+extension+'/pp/tracer_level/tracer_level.static.nc')
land = np.array(staticfile.variables['land_mask'][:])[:,:]
staticfile.close()

# loop over depositional pathways in loading data
mmvaluesf_new = np.empty([npathways,nyears,ndays,90,144])
for s in range(0,npathways): # loop through sims
    f = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+year+'_am3dd'+extension+'/pp/drydep/ts/daily/1yr/drydep.201*.'+variable[s]+'.nc')
    lon1 = f.variables['lon'][:] #grab model lat and lons
    lat1 = f.variables['lat'][:] #need last part to actually get the data
    # load all years
    temp = np.array(f.variables[variable[s]][:])[:,:,:] #get the i-th tracer from netcdf file·
    mmvaluesf = np.squeeze(temp.reshape(-1,nyears,365,90,144)) # adds an extra dimension (not sure why)
    # remove oceanic values 
    mmvaluesf[:,:,land<0.5] = np.nan
    f.close()
    # select only relevant season
    mmvaluesf_new[s,:,:,:,:] = mmvaluesf[:,startday:endday+1:1,:,:]

# add wet ground to dry ground
f = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+year+'_am3dd'+extension+'/pp/drydep/ts/daily/1yr/drydep.201*.O3_econ_g_wet.nc')
temp = np.array(f.variables['O3_econ_g_wet'][:])[:,:,:]
mmvaluesf = np.squeeze(temp.reshape(-1,nyears,365,90,144))
mmvaluesf[:,:,land<0.5] = np.nan
f.close()
mmvaluesf_new[npathways-1,:,:,:,:] = mmvaluesf_new[npathways-1,:,:,:,:]+mmvaluesf[:,startday:endday+1:1,:,:]

# calculate covariances and variances
var_eachpathway = np.empty([90,144,npathways])
cov01 = np.empty([90,144])
cov12 = np.empty([90,144])
cov23 = np.empty([90,144])
cov03 = np.empty([90,144])
cov02 = np.empty([90,144])
cov13 = np.empty([90,144])
for i in range(0,144):
    for j in range(0,90):
    	for s in range(0,npathways): # loop through sims
    	    temp = np.var(np.ravel(np.squeeze(mmvaluesf_new[s,:,:,j,i])))
    	    var_eachpathway[j,i,s] = temp
    	temp = np.cov(np.ravel(np.squeeze(mmvaluesf_new[0,:,:,j,i])),np.ravel(np.squeeze(mmvaluesf_new[1,:,:,j,i])))
    	cov01[j,i] = temp[0,1]
    	temp = np.cov(np.ravel(np.squeeze(mmvaluesf_new[1,:,:,j,i])),np.ravel(np.squeeze(mmvaluesf_new[2,:,:,j,i])))
    	cov12[j,i] = temp[0,1]
    	temp = np.cov(np.ravel(np.squeeze(mmvaluesf_new[2,:,:,j,i])),np.ravel(np.squeeze(mmvaluesf_new[3,:,:,j,i])))
    	cov23[j,i] = temp[0,1]
    	temp = np.cov(np.ravel(np.squeeze(mmvaluesf_new[0,:,:,j,i])),np.ravel(np.squeeze(mmvaluesf_new[3,:,:,j,i])))
    	cov03[j,i] = temp[0,1]
    	temp = np.cov(np.ravel(np.squeeze(mmvaluesf_new[0,:,:,j,i])),np.ravel(np.squeeze(mmvaluesf_new[2,:,:,j,i])))
    	cov02[j,i] = temp[0,1]
    	temp = np.cov(np.ravel(np.squeeze(mmvaluesf_new[1,:,:,j,i])),np.ravel(np.squeeze(mmvaluesf_new[3,:,:,j,i])))
    	cov13[j,i] = temp[0,1]
vartotal = np.sum(var_eachpathway,axis=2) + 2*cov01 + 2*cov12 + 2*cov23 + 2*cov03 + 2*cov02 + 2*cov13

# plot maps of fractional variance
filename = ("map_northhemi_frac_var_o3ddv_"+seas+"_"+sim+"_"+year+"s.pdf") #name output file
fig=plt.figure()
letter = ['(i)','(j)','(k)','(l)']
for s in range(0,npathways):
    # plot map for each pathway 
    ax = plt.subplot(npathways,1,s+1)
    m = Basemap(projection='merc',resolution = 'c',area_thresh=10000,llcrnrlon=-180,urcrnrlon=180,llcrnrlat=22,urcrnrlat=65)
    mmvalues = np.divide(np.squeeze(var_eachpathway[:,:,s]),vartotal)
    create_reg_mean_txtfile(mmvalues,'reg_mean_frac_var_o3ddv_'+variable[s]+'_at_2010s_jja.txt',lat1,lon1)
    temp_shift,lons_out = shiftgrid(178.75,mmvalues,lon1,start=False)
    x,y = m(*np.meshgrid(lons_out,lat1))
    pc = m.pcolormesh(x,y,temp_shift,vmin=vmin,vmax=vmax,cmap=mycmap)
    plt.title(seas+" "+year+"s "+variable4plotting[s],pad=1.5,fontweight='bold')
    m.drawcoastlines(linewidth=0.25)
    m.drawmapboundary(linewidth=0.5)
    m.drawmeridians(np.arange(-170.,170.,50.),labels=[0,0,0,1],linewidth=0.05,fontsize=fontsize)
    m.drawparallels(np.arange(0,90.,20),labels=[1,0,0,0],linewidth=0.05,fontsize=fontsize)
    xpt,ypt = m(-165,30)
    ax.text(xpt,ypt, letter[s])
    # write netcdf of fractional variance for each pathway 
    root_grp = Dataset('/home/oeclifton/python/AM3DD/frac_var_'+variable[s]+'.jja.nh.2010s.nc', 'w', format='NETCDF3_CLASSIC')
    root_grp.description = 'Fraction of variance that '+variable[s]+'contributes to ozone deposition velocity during June-August 2010s '
    # dimensions
    root_grp.createDimension('lat', 90)
    root_grp.createDimension('lon', 144)
    # variables
    lon = root_grp.createVariable('lon', 'f4', ('lon',))
    lat = root_grp.createVariable('lat', 'f4', ('lat',))
    field = root_grp.createVariable('frac_var', 'f4', ('lat', 'lon',))
    field.units = 'unitless'
    lat.cartesian_axis = "Y"
    lon.units =  "degrees_E"
    lon.long_name = "longitude"
    lon.cartesian_axis = "X"
    lon[:] = lon1
    lat[:] = lat1
    field[:] = mmvalues
    root_grp.close()

plt.subplots_adjust(wspace=0.1,hspace=0.3)
# add colorbar.
cax = fig.add_axes([0.25,0.065,0.5, 0.015])
cb = fig.colorbar(pc,cax=cax,pad=0.1,  orientation='horizontal')
cb.ax.set_xlabel(r'variance explained in v$_d$',labelpad=1)
cb.ax.tick_params(pad=1.2,width=0.25)
cb.outline.set_linewidth(0.25)
fig.savefig(filename) #save figure to pdf·
