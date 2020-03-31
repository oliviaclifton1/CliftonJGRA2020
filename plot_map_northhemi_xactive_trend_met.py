# oeclifton
# plot northern hemisphere seasonal change from 2010s to 2090s from xactive/dynamic simulations
# in leaf area index (LAI), precipitation, vapor pressure deficit (VPD), and soil wetness

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import MFDataset, Dataset
from mpl_toolkits.basemap import shiftgrid
from mpl_toolkits.basemap import Basemap
import seaborn as sns

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

# define color map 
mycmap = sns.diverging_palette(240,10,as_cmap=True)

# need land mask for static file·
staticfile = Dataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd_newest_final/pp/tracer_level/tracer_level.static.nc')
land = np.array(staticfile.variables['land_mask'][:])[:,:]
staticfile.close()

# define season
seas = 'JJA'
month1 = 5
month2 = 6
month3 = 7

# define time period
time = ['2010','2090']
extension = ['_newest_final','_newest_final']
nyears = [10,10]
time4file = ['201','209']

# define output file name
filename = 'northhemi_trend_met_biophys_xactive'+extension[0]+'_2010s_'+seas+'_rcp85.pdf'

# define some variables related to variables to plot  
var = ['LAI', 'precip','sphum','theta']
section = ['land','atmos_level','tracer_level','land']
factor = [1,1e05,1,1]
letter = ['(a)','(b)','(c)','(d)']
var4label = ['LAI','precip', 'VPD','soil wetness']
units4label = [r'(m$^{2}$ m$^{-2}$)',r'(10$^{-4}$ kg m$^{-2}$ s$^{-1}$)', '(kPa)',r'(m$^{3}$ m$^{-3}$)']
vmin = [-5,-1.5,-1.5,-0.3]
vmax = [5,1.5,1.5,0.3]

# read data and plot 
fig=plt.figure()
for v in range(0,4): 
	xactive_seas = np.empty([90,144,2])
	for t in range(0, 2):
	    print('/data5/oeclifton/c48L48_am3p12_rcp85_'+time[t]+'_am3dd'+extension[t]+'/pp/'+section[v]+'/ts/monthly/1yr/'+section[v]+'.'+time4file[t]+'*.'+var[v]+'.nc')
	    temp_files = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+ \
	       time[t]+'_am3dd'+extension[t]+'/pp/'+section[v]+'/ts/monthly/1yr/'+section[v]+'.'+time4file[t]+'*.'+var[v]+'.nc')
	    if var[v] == 'sphum':
	        temp1 = np.squeeze(np.array(temp_files.variables[var[v]][:])[:,47,:,:])
	    else:
	        temp1 = np.squeeze(np.array(temp_files.variables[var[v]][:])[:,:,:]) 
	    # now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
	    temp_o3 = np.squeeze(temp1.reshape(-1,nyears[t],12,90,144))
	    lon = temp_files.variables['lon'][:] #grab model lat and lons
	    lat = temp_files.variables['lat'][:] #need last part to actually get the data
	    temp_files.close()
	
	    if var4label[v] == 'VPD':
	        f1 = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+time[t]+\
	             '_am3dd'+extension[t]+'/pp/'+section[v]+'/ts/monthly/1yr/'+section[v]+'.'+time4file[t]+'*.temp.nc')
	        temp = np.array(f1.variables['temp'][:])[:,:,:,:]
	        temp = np.squeeze(temp[:,47,:,:])
	        mmvaluesf1 = np.squeeze(temp.reshape(-1,nyears[t],12,90,144))
	        mmvaluesf1 = mmvaluesf1 - 273.15 # convert to degrees C
	        f1.close()
	        f2 = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+time[t]+\
	             '_am3dd'+extension[t]+'/pp/'+section[v]+'/ts/monthly/1yr/'+section[v]+'.'+time4file[t]+'*.ps.nc')
	        temp = np.array(f2.variables['ps'][:])[:,:,:]
	        mmvaluesf2 = np.squeeze(temp.reshape(-1,nyears[t],12,90,144))
	        mmvaluesf2 = np.divide(mmvaluesf2,100) # convert to hPa
	        f2.close()
	        # need to calculate saturation vapor pressure (es)
	        es = np.multiply(6.112,np.exp(np.divide(np.multiply(17.67,mmvaluesf1),243.5+mmvaluesf1))) # in HPa
	        ws = 0.622*np.divide(es,mmvaluesf2)
	        w = np.divide(temp_o3,1-temp_o3)
	        e = np.multiply(np.divide(w, w+0.622),mmvaluesf2)
	        temp_o3 = es-e # vpd in hPa
	        temp_o3 = temp_o3*0.1 # in kPa
	    # average across season
	    xactive_seas[:,:,t] = np.squeeze(np.mean(temp_o3[:,[month1,month2,month3],:,:],axis=(0,1)))
	
	# calc trend
	trend_xactive = (np.squeeze(xactive_seas[:,:,1]-xactive_seas[:,:,0]))*factor[v]
	# remove oceanic values 
	trend_xactive[land<0.5] = np.nan
	# plot 
	ax = plt.subplot(5,1,v+1)
	m = Basemap(projection='merc',resolution = 'c',area_thresh=10000, \
	    llcrnrlon=-180,urcrnrlon=180,llcrnrlat=22,urcrnrlat=65)
	temp,lons_out = shiftgrid(178.75,trend_xactive,lon,start=False)
	x,y = m(*np.meshgrid(lons_out,lat))
	pc = m.pcolormesh(x,y,temp,vmin=vmin[v],vmax=vmax[v],cmap=mycmap)
	if v==0: 
	    plt.title(seas+' 2090s-2010s',pad=1.5,fontweight='bold')
	m.drawcoastlines(linewidth=0.25)
	m.drawmapboundary(linewidth=0.5)
	m.drawmeridians(np.arange(-170.,170.,50.),labels=[0,0,0,1],linewidth=0.05,fontsize=fontsize)
	m.drawparallels(np.arange(0,90.,20),labels=[1,0,0,0],linewidth=0.05,fontsize=fontsize)
	xpt,ypt = m(-165,30)
	ax.text(xpt,ypt, letter[v])
	# add colorbar.
	cb = fig.colorbar(pc,pad=0.04,fraction=0.046)
	cb.ax.set_ylabel(r' $\Delta$ '+var4label[v]+'\n '+units4label[v],labelpad=1)
	cb.ax.tick_params(direction='out',pad=1.4,width=0.4)
	cb.outline.set_linewidth(0.4)
fig.savefig(filename) #save figure to pdf·
