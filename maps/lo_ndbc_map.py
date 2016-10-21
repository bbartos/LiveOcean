# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 18:50:02 2016

For creation of a map of NDBC buoy stations used in the LiveOcean model.

@author: Bradley
"""
import os
import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages
from mpl_toolkits.basemap import Basemap

plt.close('all')

# full station list with locations (station, lat, lon)
sn_loc = np.array([[46088,48.334,-123.165],[46087,48.494,-124.728],[46041,47.353,-124.731],
          [46029,46.159,-124.514],[46089,45.893,-125.819],[46005,45.95,-131],
          [46050,44.656,-124.526],[46015,42.764,-124.832],[46002,42.614,-130.49]])

# set time limits to plot
t0 = datetime(2013,1,2)
t1 = datetime(2016,8,31)
date_string0 = t0.strftime('%Y.%m.%d')
date_string1 = t1.strftime('%Y.%m.%d')

# Set directories from LiveOcean
gridname = 'cascadia1'
tag = 'base'
ex_name = 'lobio1'
list_type = 'low_pass' # backfill, low_pass
Ldir = Lfun.Lstart(gridname, tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name
indir = (Ldir['LOo'] + 'moor/' + Ldir['gtagex'] + '_' + 
        list_type + '_' + date_string0 + '_' + 
        date_string1 + '/')

# Set directories and station list from process_ndbc.py
# Choose time filter
tf = 'w' # 'm', 'w', or 'd'

# Find home directory
which_home = os.environ.get("HOME")
if which_home == '/Users/PM5': # Mac
    dirname = which_home + '/Documents/LiveOcean_data/ndbc/'
    dirname2 = which_home + '/Documents/tools_data/obs_data/ndbc/'
    savname = which_home + '/Documents/LiveOcean_output/ndbc/maps/'
elif which_home == '/home/parker': # Fjord
    dirname = which_home + '/LiveOcean_data/ndbc/'
    dirname2 = which_home + '/tools_data/obs_data/ndbc/'
    savname = which_home + '/LiveOcean_output/ndbc/maps/'
elif which_home == '/home/bbartos': # Bradley's Fjord
    dirname = which_home + '/LiveOcean_data/ndbc/'
    dirname2 = which_home + '/tools_data/obs_data/ndbc/'
    savname = which_home + '/LiveOcean_output/ndbc/maps/'
elif which_home == None: # Windows version
    which_home = os.path.expanduser("~")
    dirname = which_home.replace('\\','/') + '/Documents/Research Work/Parker/LiveOcean_data/ndbc/'
    dirname2 = which_home.replace('\\','/') + '/Documents/Research Work/Parker/tools_data/obs_data/ndbc/'
    savname = which_home.replace('\\','/') + '/Documents/Research Work/Parker/LiveOcean_output/ndbc/maps/'
else:
    print('Trouble filling out environment variables')
Lfun.make_dir(savname)

# open ndbc data dictionary
fn = open((dirname + tf + '_ndbc_lo_data.txt'),'rb')
ndbc_lo_df = pickle.load(fn)

# list of stations to plot
sn_ltp = ['46088','46087','46089','46041','46050','46029'] # choose stations here

sn_df = pd.DataFrame(index = sn_ltp, columns = ['x','y','u_summ','v_summ','T_summ','u_wint','v_wint','T_wint','u_lo_summ','v_lo_summ','u_lo_wint','v_lo_wint'])
num_df = pd.DataFrame(index = sn_ltp, columns = ['x','y','u_summ','v_summ','T_summ','u_wint','v_wint','T_wint','u_lo_summ','v_lo_summ','u_lo_wint','v_lo_wint'])
c = 0
for j in range(len(sn_loc)):
    sn = str(int(sn_loc[j,0]))
    if sn in sn_ltp:
        ndbc_df = ndbc_lo_df[sn]
        # choose weeks to include in each season
        # full year
        summ = list(np.arange(14,40))
        wint = list(np.arange(1,14)) + list(np.arange(40,54))
        # half of year
#        summ = list(np.arange(23,40))
#        wint = list(np.arange(1,9)) + list(np.arange(45,54))
        
        summ_df = pd.DataFrame(index = ndbc_df.index, columns = ndbc_df.columns)
        wint_df = pd.DataFrame(index = ndbc_df.index, columns = ndbc_df.columns)
        
        # separate data seasonally
        for t in range(len(ndbc_df)):
            if ndbc_df.index.week[t] in summ:
                summ_df.ix[t] = ndbc_df.ix[t]
            elif ndbc_df.index.week[t] in wint:
                wint_df.ix[t] = ndbc_df.ix[t]
            else:
                pass
        
        # fill station averages
        sn_df['x'].ix[sn] = sn_loc[j,2]
        sn_df['y'].ix[sn] = sn_loc[j,1]
        sn_df['u_summ'].ix[sn] = np.nanmean(list(summ_df['Uwind']))
        sn_df['v_summ'].ix[sn] = np.nanmean(list(summ_df['Vwind']))
        sn_df['T_summ'].ix[sn] = np.nanmean(list(summ_df['WTMP']))
        sn_df['u_wint'].ix[sn] = np.nanmean(list(wint_df['Uwind']))
        sn_df['v_wint'].ix[sn] = np.nanmean(list(wint_df['Vwind']))
        sn_df['T_wint'].ix[sn] = np.nanmean(list(wint_df['WTMP']))
        sn_df['u_lo_summ'].ix[sn] = np.nanmean(list(summ_df['lobio1_Uwind']))
        sn_df['v_lo_summ'].ix[sn] = np.nanmean(list(summ_df['lobio1_Vwind']))
        sn_df['u_lo_wint'].ix[sn] = np.nanmean(list(wint_df['lobio1_Uwind']))
        sn_df['v_lo_wint'].ix[sn] = np.nanmean(list(wint_df['lobio1_Vwind']))
        
        # fill number of data points at each station
        num_df['x'].ix[sn] = sn_loc[j,2]
        num_df['y'].ix[sn] = sn_loc[j,1]
        num_df['u_summ'].ix[sn] = sum(summ_df['Uwind'].notnull())
        num_df['v_summ'].ix[sn] = sum(summ_df['Vwind'].notnull())
        num_df['T_summ'].ix[sn] = sum(summ_df['WTMP'].notnull())
        num_df['u_wint'].ix[sn] = sum(wint_df['Uwind'].notnull())
        num_df['v_wint'].ix[sn] = sum(wint_df['Vwind'].notnull())
        num_df['T_wint'].ix[sn] = sum(wint_df['WTMP'].notnull())
        num_df['u_lo_summ'].ix[sn] = sum(summ_df['lobio1_Uwind'].notnull())
        num_df['v_lo_summ'].ix[sn] = sum(summ_df['lobio1_Vwind'].notnull())
        num_df['u_lo_wint'].ix[sn] = sum(wint_df['lobio1_Uwind'].notnull())
        num_df['v_lo_wint'].ix[sn] = sum(wint_df['lobio1_Vwind'].notnull())
        c += 1

# map edges
lat_min = int(float(min(sn_df['y']))) - 1
lon_min = int(float(min(sn_df['x']))) - 2
lat_max = int(float(max(sn_df['y']))) + 2
lon_max = int(float(max(sn_df['x']))) + 1

# figure setup
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=[20,20])
axes[0].set_title('Summer (Apr-Sep)', fontsize='x-large')
axes[1].set_title('Winter (Oct-Mar)', fontsize='x-large')

m = Basemap(ax=axes[0], llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max, projection='cass', lon_0=(lon_max+lon_min)/2, lat_0=(lat_max+lat_min)/2, resolution='i')
m2 = Basemap(ax=axes[1], llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max, projection='cass', lon_0=(lon_max+lon_min)/2, lat_0=(lat_max+lat_min)/2, resolution='i')

# mask with map image
#m.etopo()

# draw coastlines, country boundaries, fill land
m.drawcoastlines(linewidth=1, color='darkgrey', zorder=1)
m2.drawcoastlines(linewidth=1, color='darkgrey', zorder=1)
m.drawstates(linewidth=1, linestyle='--', color='darkgrey', zorder=1)
m2.drawstates(linewidth=1, linestyle='--', color='darkgrey', zorder=1)
m.drawcountries(linewidth=1, linestyle='--', color='darkgrey', zorder=1)
m2.drawcountries(linewidth=1, linestyle='--', color='darkgrey', zorder=1)
m.fillcontinents(color='lightgrey', zorder=0)
m2.fillcontinents(color='lightgrey', zorder=0)

# draw the edge of the map and fill the ocean
m.drawmapboundary()
m2.drawmapboundary()

# draw lat/lon grid lines 
m.drawmeridians(np.arange(lon_min,lon_max,2), labels=[0,0,0,1], fontsize='x-large')
m2.drawmeridians(np.arange(lon_min,lon_max,2), labels=[0,0,0,1], fontsize='x-large')
m.drawparallels(np.arange(lat_min,lat_max,2), labels=[1,0,0,0], fontsize='x-large')
m2.drawparallels(np.arange(lat_min,lat_max,2), labels=[1,0,0,0], fontsize='x-large')

# variables
x = list(sn_df['x'])
y = list(sn_df['y'])
u_summ = list(sn_df['u_summ'])
v_summ = list(sn_df['v_summ'])
T_summ = list(sn_df['T_summ'])
u_wint = list(sn_df['u_wint'])
v_wint = list(sn_df['v_wint'])
T_wint = list(sn_df['T_wint'])
u_lo_summ = list(sn_df['u_lo_summ'])
v_lo_summ = list(sn_df['v_lo_summ'])
u_lo_wint = list(sn_df['u_lo_wint'])
v_lo_wint = list(sn_df['v_lo_wint'])

# station scatterplots
loc1 = m.scatter(x, y, s=75, latlon=True, c='none')
loc2 = m2.scatter(x, y, s=75, latlon=True, c='none')

# wind vectors
a_scale = m.quiver(x, y, u_summ, v_summ, scale=18, units='width', latlon=True, color='blue')
arrow1 = m2.quiver(x, y, u_wint, v_wint, scale=18, units='width', latlon=True, color='blue')
m.quiver(x, y, u_lo_summ, v_lo_summ, scale=18, units='width', latlon=True, color='red')
arrow2 = m2.quiver(x, y, u_lo_wint, v_lo_wint, scale=18, units='width', latlon=True, color='red')

# vector scale
plt.quiverkey(a_scale, 0.25, 0.05, 2, '2 ms$^{-1}$', coordinates='axes', labelpos='E', color='k')
plt.quiverkey(arrow1, 0.25, 0.05, 2, '2 ms$^{-1}$', coordinates='axes', labelpos='E', color='k')

# vector legend
plt.quiverkey(arrow1, 0.75, 0.09, 1, 'NDBC Buoys', coordinates='axes', labelpos='E')
plt.quiverkey(arrow2, 0.75, 0.05, 1, 'LiveOcean', coordinates='axes', labelpos='E')

# label stations
labels = sn_df.index
for j in range(len(sn_df)):
    x_coord, y_coord = m(x[j]-0.15, y[j]-0.1)
    xn_coord, yn_coord = m(x[j]-0.15, y[j]-0.2)
    axes[0].text(x_coord, y_coord, int(labels[j]), ha='right', fontsize=12, fontweight='bold')
    axes[0].text(xn_coord, yn_coord, ('N = '+str(num_df.ix[j]['u_summ'])), ha='right', fontsize=12, fontweight='bold')
    axes[1].text(x_coord, y_coord, int(labels[j]), ha='right', fontsize=12, fontweight='bold')
    axes[1].text(xn_coord, yn_coord, ('N = '+str(num_df.ix[j]['u_wint'])), ha='right', fontsize=12, fontweight='bold')

plt.show()

# save basemap
plt.savefig(savname + 'Wind_Vectors.png', bboxinches='tight')