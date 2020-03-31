# oeclifton
# plot pdfs of regional mean daily surface ozone during a given season
# for two simulations
# also include pdfs for rainy days for southeastern US 

from scipy.stats.kde import gaussian_kde
from numpy import linspace
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import MFDataset, Dataset

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
green = "#00AB66"

def geo_idx(dd, dd_array):
   """
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
    """
   geo_idx = (np.abs(dd_array - dd)).argmin()
   return geo_idx


def create_txtfile_with_daily_vals(mmvalues,filename):
    """
     create text file with daily regional mean values
     one per region 
    """
    outF = open("./txt/"+filename,"w")
    mmvalues.tolist()
    for item in mmvalues:
        outF.write(str(round(item,3)))
        outF.write("\n")
    outF.close()

# define season

#seas = "DJF"
#days_of_year = np.arange(1,366,1)
#ind = np.logical_or(days_of_year > 334, days_of_year < 60)
#ndays = 31+31+28

seas = "JJA"
days_of_year = np.arange(1,366,1)
ind = np.logical_and(days_of_year > 151, days_of_year < 244)
ndays = 31+31+30

# define some variables related to the simulation and time periods 
fileextension = ['dynamic','static']
extension = ['_newest_final','_static_o3dd_newest_final']
nsims = 2
time = '2010'
time4label = '2010s'
nyears = [10,10]

# define some variables related to plotting
factor = 1e9
ylim = [0,0.13]
yticks = [0,0.04,0.08,0.12,]
xlim = [25, 90]
hours4ticks = [30,50,70,90]
hour_names = ['30','50','70','90']

filename = 'pdf_reg_'+seas+'_dm_ozone_xactive'+extension[0]+'_static_'+time4label+'_regmean.pdf'

# need land mask (don't want to show correlations over the ocean)
staticfile = Dataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd_new/pp/tracer_level/tracer_level.static.nc')
land = np.array(staticfile.variables['land_mask'][:])[:,:]
staticfile.close()

# define regions
plotregions = ['NE US', 'SE US', 'IMW US', \
               'C. Europe', 'C. Asia', 'E. Asia' ]
fileregions = ['neus','seus','imwus','ceurope','casia','easia']
letter = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)']
reg_minlat = np.array([36., 30., 36., 45.,25., 25.])
reg_maxlat = np.array([46., 38., 46., 54.,40., 40.])
reg_minlon = np.array([-80., -95.,-110., 0., 65,105.])
reg_maxlon = np.array([-70., -80.,-100., 25.,90, 120.])
latlon_reg = ['[36-46º N,70-80º W]','[30-38º N,80-95º W]', \
'[36-46º N,100-110º W]', \
'[45-54º N,0-25º E]', \
'[25-40º N,65-90º E]',\
'[25-40º N,105-120º E]']
nregions = 6

reg_minlon[reg_minlon<0.] = reg_minlon[reg_minlon<0.]+360.
reg_maxlon[reg_maxlon<0.] = reg_maxlon[reg_maxlon<0.]+360.

# load ozone from xactive/dynamic simulations 
xactive = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+time+'_am3dd'+extension[0]+ \
                        '/pp/atmos_level/ts/daily/1yr/atmos_level.201*.O3.nc')
# load all years
temp = np.array(xactive.variables['O3'][:])[:,:,:,:] #get the i-th tracer from netcdf file
temp = np.squeeze(temp[:,47,:,:])
temp = np.squeeze(temp.reshape(-1,nyears[0],365,90,144)) # adds an extra dimension (not sure why)
dmvaluesx = temp[:,ind,:,:] # go from startday to endday with 1 step size
xactive.close()

# load precip from xactive/dynamic simulations 
b = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd_newest_final/pp/atmos_level/ts/daily/1yr/atmos_level.201*.precip.nc')
# load all years
temp1 = np.array(b.variables['precip'][:])[:,:,:] #get the i-th tracer from netcdf file·
temp1 = np.squeeze(temp1.reshape(-1,10,365,90,144)) # adds an extra dimension (not sure why)
precip_xactive = temp1[:,ind,:,:] # go from startday to endday with 1 step size
b.close()

# load ozone from static simulations
static = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+time+ \
                       '_am3dd'+extension[1]+'/pp/atmos_level/ts/daily/1yr/atmos_level.*.O3.nc')
# load all years
temp = np.array(static.variables['O3'][:])[:,:,:,:] #get the i-th tracer from netcdf file
temp = np.squeeze(temp[:,47,:,:])
temp = np.squeeze(temp.reshape(-1,nyears[1],365,90,144))
dmvaluess = temp[:,ind,:,:]
lon = static.variables['lon'][:] #grab model lat and lons
lat = static.variables['lat'][:] #need last part to actually get the data
static.close()
    
# load precip from static simulations 
b = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd_static_o3dd_newest_final/pp/atmos_level/ts/daily/1yr/atmos_level.*.precip.nc')
# load all years
temp1 = np.array(b.variables['precip'][:])[:,:,:] #get the i-th tracer from netcdf file·
temp1 = np.squeeze(temp1.reshape(-1,10,365,90,144)) # adds an extra dimension (not sure why)
precip_static = temp1[:,ind,:,:] # go from startday to endday with 1 step size
b.close()

# mask out values with less than 1/2 grid cell land
dmvaluess[:,:,land<0.5] = np.nan
dmvaluesx[:,:,land<0.5] = np.nan
precip_xactive[:,:,land<0.5] = np.nan
precip_static[:,:,land<0.5] = np.nan 

# convert precip kg/m2/s to mm/day
# kg/m2 = 1 mm
# given that 1 kg of rain water spread over 1 square meter of surface is 1 mm in thickness
precip_xactive = np.multiply( precip_xactive, 86400.)
precip_static = np.multiply( precip_static, 86400.)

# a wet day is defined is precip higher than threshold in mm/day
threshold = [1.,6.,1.,1.,1.,1.]

# calculate & plot pdfs for each region
fig = plt.figure() #open the figure
plot_number = 1
linestyle = '-'
color_wet = ['olivedrab','grey']
for r in range(0,nregions):
   print(plotregions[r])
   # get indices for lat/lon bounds for a given region
   reg_minlat_idx = geo_idx(reg_minlat[r],lat)
   reg_maxlat_idx = geo_idx(reg_maxlat[r],lat)
   reg_minlon_idx = geo_idx(reg_minlon[r],lon)
   reg_maxlon_idx = geo_idx(reg_maxlon[r],lon)
   ax = plt.subplot(4,3,r+1)
   for s in range(0,nsims):
           if s == 1:
               color = 'black'
               precip1_xactive = np.ravel(np.squeeze(np.nanmean(precip_xactive[:,:,reg_minlat_idx:reg_maxlat_idx,reg_minlon_idx:reg_maxlon_idx],axis=(2,3))))
               temp = np.nanmean(dmvaluesx[:,:,reg_minlat_idx:reg_maxlat_idx,reg_minlon_idx:reg_maxlon_idx],axis=(2,3))
               data1 = np.ravel(np.squeeze(temp))
               mask = ~np.isnan(data1)
               data1 = data1[mask]
               precip1_xactive = precip1_xactive[mask]
               if fileregions[r] == 'seus':
                  indices=np.array(precip1_xactive>threshold[r])
                  data_wet = np.extract(indices,data1)
           elif s == 0:
               color = green
               precip1_static = np.ravel(np.squeeze(np.nanmean(precip_static[:,:,reg_minlat_idx:reg_maxlat_idx,reg_minlon_idx:reg_maxlon_idx],axis=(2,3))))
               temp = np.nanmean(dmvaluess[:,:,reg_minlat_idx:reg_maxlat_idx,reg_minlon_idx:reg_maxlon_idx],axis=(2,3))
               data1 = np.ravel(np.squeeze(temp))
               mask = ~np.isnan(data1)
               data1 = data1[mask]
               precip1_static = precip1_static[mask]
               if fileregions[r] == 'seus':
                  indices=np.array(precip1_static>=threshold[r])
                  data_wet = np.extract(indices,data1)
           data1 = np.multiply(data1,factor)
           # this create the kernel, given an array it will estimate the probability over that values
           kde = gaussian_kde( data1 )
           # these are the values over wich your kernel will be evaluated
           dist_space = linspace( min(data1), max(data1), 100 )
           # plot the results
           plt.plot( dist_space, kde(dist_space),color=color,linestyle=linestyle,lw=2)
           if fileregions[r] == 'seus':
              data_wet = np.multiply(data_wet,factor)
              kde_wet = gaussian_kde ( data_wet )
              dist_space_wet = linspace( min(data_wet), max(data_wet), 100 )
              plt.plot( dist_space_wet, kde_wet(dist_space_wet),color=color_wet[s],linestyle=linestyle,lw=1)
              create_txtfile_with_daily_vals(data_wet,'reg_mean_daily_'+seas+'_sfc_o3_wet_days_ppb_'+time4label+'_'+fileextension[s]+'_'+fileregions[r]+'.txt')
           create_txtfile_with_daily_vals(data1,'reg_mean_daily_'+seas+'_sfc_o3_ppb_'+time4label+'_'+fileextension[s]+'_'+fileregions[r]+'.txt')
   ax.set_ylim(ylim)
   ax.get_yaxis().set_ticks(yticks)
   ax.set_xlim(xlim)
   ax.get_xaxis().set_ticks(hours4ticks)
   if r == nregions-1 or r == nregions-2 or r == nregions-3:
      ax.set_xticklabels(hour_names)
      ax.set_xlabel('surface ozone (ppb)')
   else:
      ax.set_xticklabels([])
   if plot_number == 1 or plot_number == 4 :
      ax.set_yticklabels(yticks)
      ax.set_ylabel(r'ppb$^{-1}$')
   else:
      ax.set_yticklabels([])
   ax.text(30,0.102, letter[r])
   plt.title(plotregions[r]+' '+latlon_reg[r],fontsize=fontsize,fontweight='bold')
   plot_number = plot_number + 1

plt.figtext((0.18)+0.6, 0.84, ('dynamic '+time4label), color='black',fontsize=fontsize)
plt.figtext((0.185)+0.6, 0.812, ('static '+ time4label), color=green,fontsize=fontsize)
plt.figtext((0.19)+0.6, 0.784, ('dynamic wet'), color=color_wet[1],fontsize=fontsize)
plt.figtext((0.195)+0.6, 0.756, ('static wet'), color=color_wet[0],fontsize=fontsize)

plt.subplots_adjust(wspace=0.1,hspace=0.4)
fig.savefig(filename) #save figure to pdf
