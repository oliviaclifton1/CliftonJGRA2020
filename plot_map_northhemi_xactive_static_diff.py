# oeclifton
# plot northern hemisphere seasonal averages for differences in ozone deposition velocity (vd) and surface ozone
# between xactive/dynamic simulations and static simulations, as well as change in vd and surface ozone from 2010s and 2090s

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import MFDataset, Dataset
from mpl_toolkits.basemap import shiftgrid
from mpl_toolkits.basemap import Basemap
import seaborn as sns
from matplotlib.patches import Polygon

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
               'C. Europe', 'S. Asia', 'E. Asia','boreal','n mid lats']
    reg_minlat = np.array([36., 30., 36., 45.,25., 25.,55,40])
    reg_maxlat = np.array([46., 38., 46., 54.,40., 40.,65,55])
    reg_minlon = np.array([-80., -95.,-110., 0., 65,105.,0,0,])
    reg_maxlon = np.array([-70., -80.,-100., 25.,90, 120.,360,360])
    nregions = 8

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
# define color map
mycmap = sns.diverging_palette(240,10,as_cmap=True)

# define locations where to draw boxes denoting regions on first map
reg_minlat = np.array([36., 30., 36., 45.,25., 25.])
reg_maxlat = np.array([46., 38., 46., 54.,40., 40.])
reg_minlon = np.array([-80., -95.,-110., 0., 65, 105.]) 
reg_maxlon = np.array([-70., -80.,-100., 25.,90, 120.])
nregions=6

# define season
seas = ['DJF','JJA']
month1 = [11,5]
month2 = [0,6]
month3 = [1,7]

# define time period and simulation info 
time = ['2010','2090']
extension = ['_newest_final','_newest_final']
extension_static = ['_static_o3dd_newest_final','_static_o3dd_newest_final']
nyears_static =[10,10]
nyears = [10,10]
time4file = ['201','209']

# define max and min for change in ozone plots 
ozone_min = -15
ozone_max = 15

# define output file name
filename = 'northhemi_diff_trend_xactive_static'+extension[0]+'_ozone_vd_rcp85.pdf' #name output file·


# load data
xactive_seas = np.empty([90,144,2,2])
vd_xactive_seas = np.empty([90,144,2,2])
static_seas = np.empty([90,144,2,2])
for t in range(0, 2):
    # load xactive ozone
    temp_files = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+time[t]+'_am3dd'+\
        extension[t]+'/pp/tracer_level/ts/monthly/1yr/tracer_level.'+time4file[t]+'*.O3.nc')
    temp1 = np.squeeze(np.array(temp_files.variables['O3'][:])[:,47,:,:]) #get the i-th tracer from netcdf file·
    # now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
    temp2 = np.squeeze(temp1.reshape(-1,nyears[t],12,90,144))
    temp_o3 = np.squeeze(temp2.mean(axis=0)) # average across years
    lon = temp_files.variables['lon'][:] #grab model lat and lons
    lat = temp_files.variables['lat'][:] #need last part to actually get the data
    temp_files.close()

    # load xactive vd
    temp_files = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+time[t]+'_am3dd'+\
        extension[t]+'/pp/drydep/ts/monthly/1yr/drydep.'+time4file[t]+'*.O3_tot_con.nc')
    temp1 = np.squeeze(np.array(temp_files.variables['O3_tot_con'][:])[:,:,:]) #get the i-th tracer from netcdf file·
    # now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
    temp2 = np.squeeze(temp1.reshape(-1,nyears[t],12,90,144))
    temp_vd = np.squeeze(temp2.mean(axis=0)) # average across years
    temp_files.close()

    # load static ozone
    temp_files = MFDataset('/data5/oeclifton/c48L48_am3p12_rcp85_'+time[t]+'_am3dd'+\
        extension_static[t]+'/pp/tracer_level/ts/monthly/1yr/tracer_level.*.O3.nc')
    temp1 = np.squeeze(np.array(temp_files.variables['O3'][:])[:,47,:,:]) #get the i-th tracer from netcdf file·
    # now shape is nyears*12,  90, 144, reshape to nyears,12,90,144
    temp2 = np.squeeze(temp1.reshape(-1,nyears_static[t],12,90,144))
    temp_o3_static = np.squeeze(temp2.mean(axis=0)) # average across years
    temp_files.close()

    # save data to arrays and average across seasons 
    for s in range(0,2):
        # average across season
        xactive_seas[:,:,t,s] = np.squeeze(np.mean(temp_o3[[month1[s],month2[s],month3[s]],:,:],axis=0))
        vd_xactive_seas[:,:,t,s] = np.squeeze(np.mean(temp_vd[[month1[s],month2[s],month3[s]],:,:],axis=0))
        static_seas[:,:,t,s] = np.squeeze(np.mean(temp_o3_static[[month1[s],month2[s],month3[s]],:,:],axis=0))

# load static vd 
# need to load that archived by GFDL model because grid is the same as xactive vd 
vd_static_seas = np.empty([90,144,2])
temp_file = Dataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd_static_o3dd_newest_final/pp/atmos_diurnal/ts/monthly/1yr/atmos_diurnal.201101-201112.O3_dvel.nc')
temp1 = np.array(temp_file.variables['O3_dvel'])[:,:,:,:]
temp_file.close()
for s in range(0,2):
    # average across season
    vd_static_seas[:,:,s] = np.squeeze(np.mean(temp1[[month1[s],month2[s],month3[s]],:,:,:],axis=(0,1)))

# calculate difference at 2010s between xactive and static 
diff_o3_at_2010s = np.squeeze(xactive_seas[:,:,0,:]-static_seas[:,:,0,:])*1e9
diff_vd_at_2010s = (np.squeeze(vd_xactive_seas[:,:,0,:])-vd_static_seas[:,:,:])*100

# calculate 2090s-2010s changes in xactive and static
trend_vd = (np.squeeze(vd_xactive_seas[:,:,1,:]-vd_xactive_seas[:,:,0,:]))*100
trend_o3_xactive = (np.squeeze(xactive_seas[:,:,1,:]-xactive_seas[:,:,0,:]))*1e9
trend_o3_static = (np.squeeze(static_seas[:,:,1,:]-static_seas[:,:,0,:]))*1e9

# need land mask for static file·
staticfile = Dataset('/data5/oeclifton/c48L48_am3p12_rcp85_2010_am3dd_new/pp/tracer_level/tracer_level.static.nc')
land = np.array(staticfile.variables['land_mask'][:])[:,:]
staticfile.close()

# remove oceanic values for vd
diff_vd_at_2010s[land<0.5,:] = np.nan
trend_vd[land<0.5,:] = np.nan

# archive regional mean wintertime vd for static (O3_dvel) and xactive (O3_tot_con)
mmvalues = np.squeeze(vd_static_seas[:,:,0])*100
mmvalues[land<0.5] = np.nan
create_reg_mean_txtfile(mmvalues,'O3_dvel_DJF.txt',lat,lon)

mmvalues = np.squeeze(vd_xactive_seas[:,:,0,0])*100
mmvalues[land<0.5] = np.nan
create_reg_mean_txtfile(mmvalues,'O3_tot_con_DJF.txt',lat,lon)

# plot 
plt.figure(figsize=(7.2,4.45), dpi=300)
fig=plt.figure()
for i in range(0,5):
    for s in range(0,2):
        temp1 = 0.
        temp = 0.
        if s == 0:
             subplot_number = [1,3,5,7,9]
             letter = ['(a)','(c)','(e)','(g)','(i)']
        elif s == 1:
             subplot_number = [2,4,6,8,10]
             letter = ['(b)','(d)','(f)','(h)','(j)']
        ax = plt.subplot(5,2,subplot_number[i])
        m = Basemap(projection='merc',resolution = 'c',area_thresh=10000,llcrnrlon=-180,urcrnrlon=180,llcrnrlat=22,urcrnrlat=65)
        if i == 0:
             temp = np.squeeze(diff_o3_at_2010s[:,:,s])
             temp1 = temp
             temp1[land<0.5] = np.nan
             create_reg_mean_txtfile(temp1,'reg_mean_diff_o3_at_2010s'+seas[s]+'.txt',lat,lon)
             title = '2010s dynamic-static surface ozone'
             vmin = ozone_min
             vmax = ozone_max
        elif i == 1:
             temp = np.squeeze(diff_vd_at_2010s[:,:,s])
             create_reg_mean_txtfile(temp,'reg_mean_diff_vd_at_2010s'+seas[s]+'.txt',lat,lon)
             title = r'2010s dynamic-static v$_d$'
             vmin = -0.4
             vmax = 0.4
        elif i == 3:
             temp = np.squeeze(trend_o3_xactive[:,:,s])
             temp1 = temp
             temp1[land<0.5] = np.nan
             create_reg_mean_txtfile(temp1,'reg_mean_xactive_trend_o3'+seas[s]+'.txt',lat,lon)
             title = '2090s-2010s dynamic surface ozone'
             vmin = ozone_min
             vmax = ozone_max
        elif i == 2:
             temp = np.squeeze(trend_o3_static[:,:,s])
             temp1 = temp
             temp1[land<0.5] = np.nan
             create_reg_mean_txtfile(temp1,'reg_mean_static_trend_o3'+seas[s]+'.txt',lat,lon)
             title = '2090s-2010s static surface ozone'
             vmin = ozone_min
             vmax = ozone_max
        elif i == 4:
             temp = np.squeeze(trend_vd[:,:,s])
             create_reg_mean_txtfile(temp,'reg_mean_xactive_trend_vd'+seas[s]+'.txt',lat,lon)
             title = r'2090s-2010s dynamic v$_d$'
             vmin = -0.4
             vmax = 0.4
        temp_shift,lons_out = shiftgrid(178.75,temp,lon,start=False)
        x,y = m(*np.meshgrid(lons_out,lat))
        pc = m.pcolormesh(x,y,temp_shift,vmin=vmin,vmax=vmax,cmap=mycmap)
        plt.title(seas[s]+' '+title,pad=1.5,fontweight='bold')
        m.drawcoastlines(linewidth=0.25)
        m.drawmapboundary(linewidth=0.5)
        m.drawmeridians(np.arange(-170.,170.,50.),labels=[0,0,0,1],linewidth=0.05,fontsize=5)
        m.drawparallels(np.arange(0,90.,20),labels=[1,0,0,0],linewidth=0.05,fontsize=5)
        xpt,ypt = m(-165,30)
        ax.text(xpt,ypt, letter[i])
        # add regional boundaries onto first plot 
        if i == 0 and s == 0:
           for n in range(0,nregions): # loop through regions
                x1,y1 = m(reg_minlon[n],reg_maxlat[n])
                x2,y2 = m(reg_maxlon[n],reg_maxlat[n])
                x3,y3 = m(reg_maxlon[n],reg_minlat[n])
                x4,y4 = m(reg_minlon[n],reg_minlat[n])
                poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],facecolor='none',edgecolor='black',linewidth=1)
                ax.add_patch(poly)
plt.subplots_adjust(wspace=0.1,hspace=0)
# add colorbar.
cax = fig.add_axes([0.25,0.05,0.5, 0.015])
cb = fig.colorbar(pc,cax=cax,pad=0.1,  orientation='horizontal')
cb.ax.set_xlabel(r'change in v$_d$ (cm s$^{-1}$)',labelpad=1)
cb.ax.tick_params(direction='in',pad=1.2,width=0.25)
ax2 = plt.twiny(ax=cax)
ax2.set_xlim(ozone_min, ozone_max)
ax2.set_xlabel('change in surface ozone (ppb)')
ax2.tick_params(direction='in',pad=1.2,width=0.25)
cb.outline.set_linewidth(0.25)
fig.savefig(filename,dpi=300) #save figure to pdf·

