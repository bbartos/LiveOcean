# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 18:59:53 2016

@author: Bradley

Plotting both ndbc buoy data and moorings of the LiveOcean model for comparison.
Use ptools/ndbc/get_ndbc.py for retrieval of ndbc data.
Use ptools/ndbc/process_ndbc.py for processing ndbc data.
On Fjord, use LiveOcean/moor/moor_ndbc.py for retrieval of LiveOcean data.
On personal computer, use LiveOcean/get_LiveOcean.py for retrieval of LO data 
    from online, then use LiveOcean/moor/moor_ndbc.py.
Use LiveOcean/stats_LiveOcean_ndbc.py to process LO data and retrieve stats.
"""

import os
import sys
alp = os.path.abspath('./alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from warnings import filterwarnings
filterwarnings('ignore') # skip some warning messages

# set time limits to plot
t0 = datetime(2013,1,2)
t1 = datetime(2016,8,8)
date_string0 = t0.strftime('%Y.%m.%d')
date_string1 = t1.strftime('%Y.%m.%d')
t02 = datetime(2013,1,2)
t12 = datetime(2016,8,31)
date_string02 = t02.strftime('%Y.%m.%d')
date_string12 = t12.strftime('%Y.%m.%d')

# Set directories from LiveOcean
gridname = 'cascadia1'
tag = 'base'
ex_name = 'lo1'
ex_name2 = 'lobio1'
list_type = 'low_pass' # backfill, low_pass
Ldir = Lfun.Lstart(gridname, tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name
Ldir['gtagex2'] = Ldir['gtag'] + '_' + ex_name2
indir = (Ldir['LOo'] + 'moor/' + Ldir['gtagex'] + '_' + 
        list_type + '_' + date_string0 + '_' + date_string1 + '/')
indir2 = (Ldir['LOo'] + 'moor/' + Ldir['gtagex2'] + '_' + 
        list_type + '_' + date_string02 + '_' + date_string12 + '/')

# Set directories and station list from process_ndbc.py
# Choose time filter
tf = 'w' # 'm', 'w', or 'd'

# Choose directories
dirname = Ldir['data'] + '/ndbc/'
dirname2 = Ldir['parent'] + '/tools_data/obs_data/ndbc/'
savname = Ldir['LOo'] + '/ndbc/lo1_lobio1/'
Lfun.make_dir(savname)

# open ndbc data dictionary
fn = open((dirname + tf + '_ndbc_lo_data.txt'),'rb')
ndbc_lo_df = pickle.load(fn)

# open climatology dictionary
fn_clim = open((dirname2 + 'ndbc_clim_df.txt'),'rb')
ndbc_clim_df = pickle.load(fn_clim)
fn_clim_std = open((dirname2 + 'ndbc_clim_std_df.txt'),'rb')
ndbc_clim_std_df = pickle.load(fn_clim_std)

# open stats dictionary
fn_skill = open((dirname + tf + '_skill_dictionary.txt'), 'rb')
sn_skill = pickle.load(fn_skill)

# ndbc stations for standard LiveOcean
sn_list = ['46088','46087','46041','46029','46089','46050']
sn_list_1 = sn_list[0:2]
sn_list_2 = sn_list[2:]

# For variable stations
#sn_list = list(ndbc_df.keys())
#sn_list_1 = sn_list[0:int(0.5*len(sn_list)+1)]
#sn_list_2 = sn_list[int(0.5*len(sn_list)+1):]

# Construct figures with subplots
plt.close('all')
NC = 3
# Time Series Figures
NR_1 = len(sn_list_1)
fig_1, axes_1 = plt.subplots(nrows=NR_1, ncols=NC, figsize=(30,18))
plt.suptitle('Salish Sea', fontsize='x-large')
plt.subplots_adjust(left=0.05, right=0.95, top=0.925)

NR_2 = int(len(sn_list_2)/2)
fig_2, axes_2 = plt.subplots(nrows=NR_2, ncols=NC, figsize=(30,18))
plt.suptitle('Northern Stations', fontsize='x-large')
plt.subplots_adjust(left=0.05, right=0.95, top=0.925)

NR_3 = len(sn_list_2) - NR_2
fig_3, axes_3 = plt.subplots(nrows=NR_3, ncols=NC, figsize=(30,18))
plt.suptitle('Southern Stations', fontsize='x-large')
plt.subplots_adjust(left=0.05, right=0.95, top=0.925)

cc_time = 0

# Property_Property Figures
fig_4,axis_4 = plt.subplots(nrows=2,ncols=3,figsize=(30,18))
plt.suptitle('Water Temp (\u00B0C)', fontsize='x-large')
plt.subplots_adjust(left=0.05, right=0.95, top=0.925)

fig_5,axis_5 = plt.subplots(nrows=2,ncols=3,figsize=(30,18))
plt.suptitle('U Wind (ms$^{-1}$)', fontsize='x-large')
plt.subplots_adjust(left=0.05, right=0.95, top=0.925)

fig_6,axis_6 = plt.subplots(nrows=2,ncols=3,figsize=(30,18))
plt.suptitle('V Wind (ms$^{-1}$)', fontsize='x-large')
plt.subplots_adjust(left=0.05, right=0.95, top=0.925)

cc_prop = 0

# Set minimums and maximums
mintemp = 5
maxtemp = 20
minwind = -15
maxwind = 15

# Begin station loop
for sn in sn_list:
    if sn in sn_list_1:
        NR = NR_1
        fig = fig_1
        axes = axes_1
    elif sn in sn_list_2[0:NR_2]:
        NR = NR_2
        fig = fig_2
        axes = axes_2
    else:
        NR = NR_3
        fig = fig_3
        axes = axes_3

#%% load Data
    
# LiveOcean
    try:
        fn = open((indir + tf + 'filter_' + Ldir['gtagex'] + '_' + sn + '_' + 
            list_type + '_' + date_string0 + '_' + 
            date_string1 + '.nc'), 'rb')
        V = pickle.load(fn)
    except FileNotFoundError:
        print('LiveOcean data for Station ' + sn + ' and time ' + date_string0 + 
            '_' + date_string1 + ' not found.')
        continue
    try:
        fn2 = open((indir2 + tf + 'filter_' + Ldir['gtagex2'] + '_' + sn + '_' + 
            list_type + '_' + date_string02 + '_' + 
            date_string12 + '.nc'), 'rb')
        V2 = pickle.load(fn2)
    except FileNotFoundError:
        print('LiveOceanBio data for Station ' + sn + ' and time ' + date_string0 + 
            '_' + date_string1 + ' not found.')
        continue
    
    # Unit dictionary
    V_units = dict()
    V_units['Uwind'] = 'ms$^{-1}$'
    V_units['Vwind'] = 'ms$^{-1}$'
    V_units['temp'] = u'\N{DEGREE SIGN}C'
    
    # Dates for plotting
    mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])
    mdt = mdates.num2date(mdays) # list of datetimes of data
    mdays2 = Lfun.modtime_to_mdate_vec(V2['ocean_time'])
    mdt2 = mdates.num2date(mdays2)
    zero_line = np.zeros(len(mdt))
    
# ndbc data
    try:
        DFF = ndbc_lo_df[sn]
    except FileNotFoundError:
        print('ndbc data for Station ' + sn + ' and time ' + date_string0 + 
            '_' + date_string1 + ' not found.')
        continue

# Climatology data
    clim_yr = ndbc_clim_df[sn] # dataframe of ordinal weeks
    clim_DFF = pd.DataFrame(columns=clim_yr.columns)
    for yr in pd.date_range('2013',t1,freq='AS'):
        clim_yr.index = pd.date_range(start=yr,periods=len(clim_yr),freq='W')
        clim_DFF = clim_DFF.append(clim_yr)
    # standard deviation
    clim_std_DFF = ndbc_clim_std_df[sn] # dataframe with ordinal week index
    clim_std_yr = ndbc_clim_std_df[sn] # dataframe for later years
    for yr in range(t0.year,t1.year+1):
        if yr == t0.year:
            clim_std_DFF.index = pd.date_range(start=yr,periods=len(clim_std_DFF),freq='W')
        else:
            clim_std_yr.index = pd.date_range(start=yr,periods=len(clim_std_yr),freq='W')
            clim_std_DFF = clim_std_DFF.append(clim_std_yr)            

# Skill tests
    skill_df = sn_skill[sn]

# Time Series Figures
    # Set subplot index
    ir_time = cc_time
    cc_time += 1
    if sn in sn_list_1 and cc_time == len(sn_list_1):
        cc_time = 0
    elif sn == sn_list_2[NR_2-1] and cc_time == NR_2:
        cc_time = 0
    for ic_time in range(NC):
        ax = axes[ir_time,ic_time]

    # WTMP column
        if ic_time==0: 
            # LiveOcean
            vn2 = 'temp'
            ax.plot(mdt, V[vn2], 'g', linewidth=2)
            ax.plot(mdt2, V2[vn2], 'darkblue', linewidth=2)
            # ndbc
            vn = 'WTMP'
            DFF[vn].plot(ax=ax, xlim=(t0,t1), ylim=(mintemp,maxtemp), color='r', linewidth=2)
            ax.text(.05, .9, 'NDBC Station ' + str(sn), transform=ax.transAxes)
            # climatology
            clim_DFF[vn].plot(ax=ax, xlim=(t0,t1), ylim=(mintemp,maxtemp), color='k', linewidth=1)
            clim_min = clim_DFF[vn].as_matrix()-clim_std_DFF[vn].as_matrix()
            clim_max = clim_DFF[vn].as_matrix()+clim_std_DFF[vn].as_matrix()
            ax.fill_between(clim_DFF.index, list(clim_min), list(clim_max), color='grey', alpha=0.3)
            ax.axhline(linestyle='--', color='k', linewidth=1)
            ax.axvline(x=datetime(2014,1,1), color='k', linewidth=0.75)
            ax.axvline(x=datetime(2015,1,1), color='k', linewidth=0.75)
            ax.axvline(x=datetime(2016,1,1), color='k', linewidth=0.75)
            if ir_time==0:
                ax.set_xlabel('')
                ax.set_xticklabels([])
                ax.set_title(vn2+' ['+V_units[vn2]+']')
            elif ir_time==(NR-1):
                pass
            else:
                ax.set_xlabel('')
                ax.set_xticklabels([])

    # Uwind column
        if ic_time==1:
            vn = 'Uwind'
            # LiveOcean
            ax.plot(mdt, V[vn], 'g', linewidth=2)
            ax.plot(mdt2, V2[vn], 'darkblue', linewidth=2)
            # year and zero lines
            ax.axhline(linestyle='--', color='k', linewidth=1)
            ax.axvline(x=datetime(2014,1,1), color='k', linewidth=0.75)
            ax.axvline(x=datetime(2015,1,1), color='k', linewidth=0.75)
            ax.axvline(x=datetime(2016,1,1), color='k', linewidth=0.75)
            # ndbc
            DFF[vn].plot(ax=ax, xlim=(t0,t1), ylim=(minwind,maxwind), color='r', linewidth=2)
            # climatology
            clim_DFF[vn].plot(ax=ax, xlim=(t0,t1), ylim=(minwind,maxwind), color='k', linewidth=1)
            clim_min = clim_DFF[vn].as_matrix()-clim_std_DFF[vn].as_matrix()
            clim_max = clim_DFF[vn].as_matrix()+clim_std_DFF[vn].as_matrix()
            ax.fill_between(clim_DFF.index, list(clim_min), list(clim_max), color='grey', alpha=0.3)
            if ir_time==0:
                ax.set_xlabel('')
                ax.set_xticklabels([])
                ax.set_title(vn+' ['+V_units[vn]+']')
            elif ir_time==(NR-1):
                pass
            else:
                ax.set_xlabel('')
                ax.set_xticklabels([])

    # Vwind column
        if ic_time==2:
            vn = 'Vwind'
            # LiveOcean
            ax.plot(mdt, V[vn], 'g', label='LiveOcean', linewidth=2)
            ax.plot(mdt2, V2[vn], 'darkblue', label='LiveOcean Bio', linewidth=2)
            # year and zero lines
            ax.axhline(linestyle='--', color='k', linewidth=1)
            ax.axvline(x=datetime(2014,1,1), color='k', linewidth=0.75)
            ax.axvline(x=datetime(2015,1,1), color='k', linewidth=0.75)
            ax.axvline(x=datetime(2016,1,1), color='k', linewidth=0.75)
            # ndbc
            DFF[vn].plot(ax=ax, xlim=(t0,t1), ylim=(minwind,maxwind), color='r', linewidth=2, label='ndbc buoys')
            # climatology
            clim_DFF[vn].plot(ax=ax, xlim=(t0,t1), ylim=(minwind,maxwind), color='k', linewidth=1, label='Climatology')
            clim_min = clim_DFF[vn].as_matrix()-clim_std_DFF[vn].as_matrix()
            clim_max = clim_DFF[vn].as_matrix()+clim_std_DFF[vn].as_matrix()
            ax.fill_between(clim_DFF.index, list(clim_min), list(clim_max), color='grey', alpha=0.3, label='1 Std of Clim')
            if ir_time==0:
                ax.set_xlabel('')
                ax.set_xticklabels([])
                ax.set_title(vn+' ['+V_units[vn]+']')
            elif ir_time==(NR-1):
                ax.legend(bbox_to_anchor=(0.825,1.25), ncol=2)
            else:
                ax.set_xlabel('')
                ax.set_xticklabels([])

    # Stats
        ax.text(.575, .15, 'MB lo1 = %.2f' %skill_df['lo1_'+vn].ix['MB'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['MB'],transform=ax.transAxes)
        ax.text(.575, .1, 'SDE lo1 = %.2f' %skill_df['lo1_'+vn].ix['SDE'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['SDE'],transform=ax.transAxes)
        ax.text(.575, .05, 'CC lo1 = %.2f' %skill_df['lo1_'+vn].ix['CC'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['CC'], transform=ax.transAxes)

#Property-Property Figures
    # set subplot number
    ir_prop = int(cc_prop/3)
    ic_prop = int(cc_prop-3*ir_prop)
    
    # Re-pull data
    DFF = ndbc_lo_df[sn]
    skill_df = sn_skill[sn]
    
    # regression line data
    x_temp = np.arange(mintemp,maxtemp+1)
    y_lo1_temp = skill_df['lo1_WTMP'].ix['slope'] * x_temp + skill_df['lo1_WTMP'].ix['yint']
    y_lobio1_temp = skill_df['lobio1_WTMP'].ix['slope'] * x_temp + skill_df['lobio1_WTMP'].ix['yint']
    x_wind = np.arange(minwind,maxwind+1)
    y_lo1_Uwind = skill_df['lo1_Uwind'].ix['slope'] * x_wind + skill_df['lo1_Uwind'].ix['yint']
    y_lobio1_Uwind = skill_df['lobio1_Uwind'].ix['slope'] * x_wind + skill_df['lobio1_Uwind'].ix['yint']
    y_lo1_Vwind = skill_df['lo1_Vwind'].ix['slope'] * x_wind + skill_df['lo1_Vwind'].ix['yint']
    y_lobio1_Vwind = skill_df['lobio1_Vwind'].ix['slope'] * x_wind + skill_df['lobio1_Vwind'].ix['yint']
    
    # WTMP
    vn = 'WTMP'
    ax=axis_4[ir_prop,ic_prop]
    ax.scatter(DFF[vn], DFF['lo1_'+vn], color='g', s=15, label='LiveOcean')
    ax.scatter(DFF[vn], DFF['lobio1_'+vn], color='darkblue', s=5, label='LiveOcean Bio')
    ax.plot(x_temp, y_lo1_temp, 'g') # model best fit line
    ax.plot(x_temp, y_lobio1_temp, 'darkblue')
    ax.plot([mintemp,maxtemp], [mintemp,maxtemp], '--k') # line y=x for comparison
    ax.set_xlim([mintemp,maxtemp])
    ax.set_ylim([mintemp,maxtemp])
    ax.set_xlabel('NDBC Data')
    ax.set_ylabel('Model Data')
    ax.set_title('NDBC Station ' + sn)
    ax.text(.05, .9, 'R$^2$ lo1 = %.2f' %skill_df['lo1_'+vn].ix['R2'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['R2'],transform=ax.transAxes)
    ax.text(.05, .85, 'Slope Offest lo1 = %.2f' %abs(skill_df['lo1_'+vn].ix['slope']-1) + ' lobio1 = %.2f' %abs(skill_df['lobio1_'+vn].ix['slope']-1),transform=ax.transAxes)
    ax.text(.5, .1, 'RMSE lo1 = %.2f' %skill_df['lo1_'+vn].ix['RMSE'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['RMSE'],transform=ax.transAxes)
    ax.text(.5, .05, 'SS lo1 = %.2f' %skill_df['lo1_'+vn].ix['SS'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['SS'],transform=ax.transAxes)
    if ic_prop != 0:
        ax.set_ylabel('')
    if ir_prop == 0:
        ax.set_xlabel('')
    if ir_prop == 1 and ic_prop == 2:
        ax.legend(bbox_to_anchor=(1.05,1.175))

    # Uwind
    vn = 'Uwind'
    ax=axis_5[ir_prop,ic_prop]
    ax.scatter(DFF[vn], DFF['lo1_'+vn], color='g', s=15, label='LiveOcean')
    ax.scatter(DFF[vn], DFF['lobio1_'+vn], color='darkblue', s=5, label='LiveOcean Bio')
    ax.plot(x_wind, y_lo1_Uwind, 'g') # model best fit line
    ax.plot(x_wind, y_lobio1_Uwind, 'darkblue')
    ax.plot([minwind,maxwind], [minwind,maxwind], '--k') # line y=x for comparison
    ax.set_xlim([minwind,maxwind])
    ax.set_ylim([minwind,maxwind])
    ax.set_xlabel('NDBC Data')
    ax.set_ylabel('Model Data')
    ax.set_title('NDBC Station ' + sn)
    ax.text(.05, .9, 'R$^2$ lo1 = %.2f' %skill_df['lo1_'+vn].ix['R2'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['R2'],transform=ax.transAxes)
    ax.text(.05, .85, 'Slope Offest lo1 = %.2f' %abs(skill_df['lo1_'+vn].ix['slope']-1) + ' lobio1 = %.2f' %abs(skill_df['lobio1_'+vn].ix['slope']-1),transform=ax.transAxes)
    ax.text(.5, .1, 'RMSE lo1 = %.2f' %skill_df['lo1_'+vn].ix['RMSE'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['RMSE'],transform=ax.transAxes)
    ax.text(.5, .05, 'SS lo1 = %.2f' %skill_df['lo1_'+vn].ix['SS'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['SS'],transform=ax.transAxes)
    if ic_prop != 0:
        ax.set_ylabel('')
    if ir_prop == 0:
        ax.set_xlabel('')
    if ir_prop == 1 and ic_prop == 2:
        ax.legend(bbox_to_anchor=(1.05,1.175))

    # Vwind
    vn = 'Vwind'
    ax=axis_6[ir_prop,ic_prop]
    ax.scatter(DFF[vn], DFF['lo1_'+vn], color='g', s=15, label='LiveOcean')
    ax.scatter(DFF[vn], DFF['lobio1_'+vn], color='darkblue', s=5, label='LiveOcean Bio')
    ax.plot(x_wind, y_lo1_Vwind, 'g') # model best fit line
    ax.plot(x_wind, y_lobio1_Vwind, 'darkblue')
    ax.plot([minwind,maxwind], [minwind,maxwind], '--k') # line y=x for comparison
    ax.set_xlim([minwind,maxwind])
    ax.set_ylim([minwind,maxwind])
    ax.set_xlabel('NDBC Data')
    ax.set_ylabel('Model Data')
    ax.set_title('NDBC Station ' + sn)
    ax.text(.05, .9, 'R$^2$ lo1 = %.2f' %skill_df['lo1_'+vn].ix['R2'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['R2'],transform=ax.transAxes)
    ax.text(.05, .85, 'Slope Offest lo1 = %.2f' %abs(skill_df['lo1_'+vn].ix['slope']-1) + ' lobio1 = %.2f' %abs(skill_df['lobio1_'+vn].ix['slope']-1),transform=ax.transAxes)
    ax.text(.5, .1, 'RMSE lo1 = %.2f' %skill_df['lo1_'+vn].ix['RMSE'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['RMSE'],transform=ax.transAxes)
    ax.text(.5, .05, 'SS lo1 = %.2f' %skill_df['lo1_'+vn].ix['SS'] + ' lobio1 = %.2f' %skill_df['lobio1_'+vn].ix['SS'],transform=ax.transAxes)
    if ic_prop != 0:
        ax.set_ylabel('')
    if ir_prop == 0:
        ax.set_xlabel('')
    if ir_prop == 1 and ic_prop == 2:
        ax.legend(bbox_to_anchor=(1.05,1.175))
    cc_prop += 1



# Save figures
fig_1.savefig(savname + 'Time_Series_1_lo1_lobio1.png', bbox_inches='tight')
fig_2.savefig(savname + 'Time_Series_2_lo1_lobio1.png', bbox_inches='tight')
fig_3.savefig(savname + 'Time_Series_3_lo1_lobio1.png', bbox_inches='tight')
fig_4.savefig(savname + 'Property_1_lo1_lobio1.png', bbox_inches='tight')
fig_5.savefig(savname + 'Property_2_lo1_lobio1.png', bbox_inches='tight')
fig_6.savefig(savname + 'Property_3_lo1_lobio1.png', bbox_inches='tight')

