# oeclifton
# evaluate model surface ozone with TOAR surface ozone dataset
# look at winter and summer northern hemisphere
# plot bias (simulated-observed) and change in bias with new model simulation

import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import ticker
from netCDF4 import MFDataset

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

# define this function for sampling model at a given location 
def geo_idx(dd, dd_array):
   """
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
   """
   geo_idx = (np.abs(dd_array - dd)).argmin()
   return geo_idx

# choose baseline simulation to compare to·
sim = 'static'
extension = '_static_o3dd_newest_final'
nyears=10

#define three months to average over & season name
seas = ["DJF","JJA"]
nseas = 2 

# define some variables related to plotting 
fig, axes = plt.subplots(3,2) #open the figure·
filename = ("map_northhemi_seas_ozone_bias_diff_xactive_"+sim+".pdf") #name output file·
datamin = 0 # for absolute bias
datamax = 25 # for absolute bias
datamin2 = -15 # for change in bias 
datamax2 = 15 # for change in bias 

# define observational dataset
obsfilename = '/home/oeclifton/TOAR_sfc_ozone_seasonal_global_2008-2015_aggregated.csv'

# plot bias
for i in range(0,nseas): # loop over seasons 
	if seas[i] == 'DJF':
	    # define variables related to reading obs dataset
	    ind_data_capture = 25
	    ind_smean = 26
	    ind_count = 35
	    month1 = 0
	    month2 = 1
	    month3 = 11
	    # define plotting variables 
	    letter = ['(a)','(c)','(e)']
	elif seas[i] == 'JJA':
	    # define variables related to reading obs dataset
	    ind_data_capture = 47
	    ind_smean = 48
	    ind_count = 47
	    month1 = 5
	    month2 = 6
	    month3 = 7
	    # define plotting variables
	    letter = ['(b)','(d)','(f)']
	# load TOAR observational dataset
	# this file has some hashtags in it, but genfromtxt doesn't like this
	# setting comments (the character used to indicate the start of a comment) to 'OLIVIA' fixes things
	siteinfo = np.genfromtxt(obsfilename, skip_header=84, dtype=str,delimiter =';', autostrip=True,comments='OLIVIA')
	nsites,ncols = np.shape(siteinfo)
	slon = []
	slat = []
	data_capture = []
	years_count = []
	site_details = []
	site_details2 = []
	smean = []
	for s in siteinfo:
	    slon.append(s[8])
	    slat.append(s[7])
	    data_capture.append(s[ind_data_capture])
	    smean.append(s[ind_smean])
	    years_count.append(s[ind_count])
	    site_details.append(s[5])
	    site_details2.append(s[6])
	slon = np.array(slon)
	slon = slon.astype(np.float)
	slat = np.array(slat)
	slat = slat.astype(np.float)
	data_capture = np.array(data_capture)
	data_capture = data_capture.astype(np.float)
	years_count = np.array(years_count)
	years_count = data_capture.astype(np.float)
	smean = np.array(smean)
	smean = smean.astype(np.float)
	site_details2 = np.array(site_details2)
	site_details = np.array(site_details)
	# remove sites with certain conditions 
	smean[np.logical_or(np.logical_or(site_details2 == 'urban',site_details2 == 'suburban'), \
	np.logical_or(site_details == 'traffic',site_details == 'industry'))]=np.nan
	smean[np.logical_or(data_capture<0.5,years_count<0.5)]=np.nan
	
	# load model data
	xactive = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd_newest_final/pp/tracer_level/ts/monthly/1yr/tracer_level.201*12.O3.nc')
	static = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd'+extension+'/pp/tracer_level/ts/monthly/1yr/tracer_level.*12.O3.nc')
	# load all years for xactive
	mmvaluesx = np.array(xactive.variables['O3'][:])[:,47,:,:] #get the i-th tracer from netcdf file·
	# now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
	mmvaluesx = np.squeeze(mmvaluesx.reshape(-1,10,12,90,144))
	mmvaluesx = np.squeeze(mmvaluesx.mean(axis=0)) # average across years
	mmvaluesx = np.divide((mmvaluesx[month1,:,:]+mmvaluesx[month2,:,:]+mmvaluesx[month3,:,:]),3.0)
	mmvaluesx = np.squeeze(mmvaluesx)*1e9 #get rid of singleton dimension
	xactive.close()
	# load all years for static 
	mmvaluess = np.array(static.variables['O3'][:])[:,47,:,:] #get the i-th tracer from netcdf file·
	# now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
	mmvaluess = np.squeeze(mmvaluess.reshape(-1,nyears,12,90,144))
	mmvaluess = np.squeeze(mmvaluess.mean(axis=0)) # average across years
	mmvaluess = np.divide((mmvaluess[month1,:,:]+mmvaluess[month2,:,:]+mmvaluess[month3,:,:]),3.0)
	mmvaluess = np.squeeze(mmvaluess)*1e9 #get rid of singleton dimension
	lon = static.variables['lon'][:] #grab model lat and lons
	lat = static.variables['lat'][:] #need last part to actually get the data
	static.close()
	
	# sample model at observational sites
	mmvaluesxc = np.empty([nsites])
	mmvaluessc = np.empty([nsites])
	slon2 = slon
	for s in range(0,nsites):
	   if slon[s]<0:
	      temp = slon[s]+360
	   else:
	      temp = slon[s]
	   s_mod_lon_idx = geo_idx(temp,lon)
	   s_mod_lat_idx = geo_idx(slat[s],lat)
	   mmvaluesxc[s] = mmvaluesx[s_mod_lat_idx,s_mod_lon_idx]
	   mmvaluessc[s] = mmvaluess[s_mod_lat_idx,s_mod_lon_idx]
	
	# calculate model bias (model-obs)
	mmxbias = mmvaluesxc-smean
	mmsbias = mmvaluessc-smean
	# calculate change in the bias 
	mmdiffbias = mmxbias-mmsbias
	# remove any time when mmsbias is negative
	mmdiffbias[mmsbias<0.]=np.nan
	# print number of times this happens 
	mmsbias2=mmsbias
	mmsbias2=mmsbias2[mmsbias<0.]
	print(len(mmsbias2))
	cmap = plt.cm.get_cmap("YlOrBr")
	cmap.set_under('lightblue')
	for si in range(0,2): # simulations
	    m = Basemap(resolution='c',area_thresh=10000,llcrnrlon=-180,urcrnrlon=180,llcrnrlat=0,urcrnrlat=90,ax=axes[si,i])
	    m.drawcoastlines(linewidth=0.25)
	    m.drawmeridians(np.arange(-170.,170.,50.),labels=[0,0,0,1],linewidth=0.05,fontsize=5)
	    m.drawparallels(np.arange(0,90.,20),labels=[1,0,0,0],linewidth=0.05,fontsize=5)
	    x,y = m(slon, slat)
	    if si == 0:
	        bias = mmxbias
	        title = (seas[i]+' 2010s dynamic')
	    elif si == 1:
	        bias = mmsbias
	        title = (seas[i]+" 2010s "+sim)
	    pc = m.scatter(x, y, s=0.25,c=bias, cmap=cmap,vmin=datamin,vmax=datamax,ax=axes[si,i])
	    m.fillcontinents(color='lightgrey',lake_color='white',zorder=0)
	    axes[si,i].set_title(title,pad=1.5,fontweight='bold')
	    xpt,ypt = m(-165,75)
	    axes[si,i].text(xpt,ypt, letter[si])
	# plot change in bias
	cmap = plt.cm.get_cmap("RdBu_r")
	m = Basemap(resolution='c',area_thresh=10000,llcrnrlon=-180,urcrnrlon=180,llcrnrlat=0,urcrnrlat=90,ax=axes[2,i])
	m.drawcoastlines(linewidth=0.25)
	m.drawmeridians(np.arange(-170.,170.,50.),labels=[0,0,0,1],linewidth=0.05,fontsize=5)
	m.drawparallels(np.arange(0,90.,20),labels=[1,0,0,0],linewidth=0.05,fontsize=5)
	x,y = m(slon, slat)
	bias = mmdiffbias
	pc2 = m.scatter(x, y, s=0.25,c=bias, cmap=cmap,vmin=datamin2,vmax=datamax2,ax=axes[2,i])
	m.fillcontinents(color='lightgrey',lake_color='white',zorder=0)
	axes[2,i].set_title(seas[i]+' 2010s dynamic bias - static bias',pad=1.5,fontweight='bold')
	xpt,ypt = m(-165,75)
	axes[2,i].text(xpt,ypt, letter[2])

plt.subplots_adjust(wspace=0.1,hspace=0.)

# add colorbar
cbaxes = fig.add_axes([0.15, 0.075,0.3,0.01])
cbar = fig.colorbar(pc, cax=cbaxes,orientation='horizontal')
cbar.ax.tick_params(labelsize=fontsize,width=0.5)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.set_label('surface ozone bias (ppb)')
cbar.outline.set_linewidth(0.5)

# add colorbar.
cbaxes = fig.add_axes([0.55, 0.075,0.3,0.01])
cbar = fig.colorbar(pc2,cax=cbaxes,orientation='horizontal')
cbar.ax.tick_params(labelsize=fontsize,width=0.5)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.set_label('bias difference (ppb)')
cbar.outline.set_linewidth(0.5)

# save file 
fig.savefig(filename) #save figure to pdf·
