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
    savname = which_home + '/Documents/LiveOcean_output/ndbc/lobio1/'
elif which_home == '/home/parker': # Fjord
    dirname = '/data1/parker/LiveOcean_data/ndbc/'
    dirname2 = '/data1/parker/tools_data/obs_data/ndbc/'
    savname = '/data1/parker/LiveOcean_output/ndbc/lobio1/'
elif which_home == '/home/bbartos': # Bradley's Fjord
    dirname = which_home + '/LiveOcean_data/ndbc/'
    dirname2 = which_home + '/tools_data/obs_data/ndbc/'
    savname = which_home + '/LiveOcean_output/ndbc/lobio1/'
elif which_home == None: # Windows version
    which_home = os.path.expanduser("~")
    dirname = which_home.replace('\\','/') + '/Documents/Research Work/Parker/LiveOcean_data/ndbc/'
    dirname2 = which_home.replace('\\','/') + '/Documents/Research Work/Parker/tools_data/obs_data/ndbc/'
    savname = which_home.replace('\\','/') + '/Documents/Research Work/Parker/LiveOcean_output/ndbc/lobio1/'
else:
    print('Trouble filling out environment variables')
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

# Property-Property Figures
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
    
    # Unit dictionary
    V_units = dict()
    V_units['Uwind'] = 'ms$^{-1}$'
    V_units['Vwind'] = 'ms$^{-1}$'
    V_units['temp'] = u'\N{DEGREE SIGN}C'
    
    # Dates for plotting
    mdays = Lfun.modtime_to_mdate_vec(V['ocean_time'])
    mdt = mdates.num2date(mdays) # list of datetimes of data
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
    ir = cc_time
    cc_time += 1
    if sn in sn_list_1 and cc_time == len(sn_list_1):
        cc_time = 0
    elif sn == sn_list_2[NR_2-1] and cc_time == NR_2:
        cc_time = 0
    for ic in range(NC):
        ax = axes[ir,ic]

    # WTMP column
        if ic==0: 
            # LiveOcean
            vn2 = 'temp'
            ax.plot(mdt, V[vn2], 'darkblue', linewidth=2)
            # ndbc
            vn = 'WTMP'
            DFF[vn].plot(ax=ax, xlim=(t0,t1), ylim=(mintemp,maxtemp), color='r', linewidth=2)
            ax.text(.05, .9, 'NDBC Station ' + str(sn), transform=ax.transAxes)
            # climatology
            clim_DFF[vn].plot(ax=ax, xlim=(t0,t1), ylim=(mintemp,maxtemp), color='k', linewidth=1)
            clim_min = clim_DFF[vn].as_matrix()-clim_std_DFF[vn].as_matrix()
            clim_max = clim_DFF[vn].as_matrix()+clim_std_DFF[vn].as_matrix()
            ax.fill_between(clim_DFF.index, list(clim_min), list(clim_max), color='grey', alpha=0.3)
            # year and zero lines
            ax.axhline(linestyle='--', color='k', linewidth=1)
            ax.axvline(x=datetime(2014,1,1), color='k', linewidth=0.75)
            ax.axvline(x=datetime(2015,1,1), color='k', linewidth=0.75)
            ax.axvline(x=datetime(2016,1,1), color='k', linewidth=0.75)
            if ir==0:
                ax.set_xlabel('')
                ax.set_xticklabels([])
                ax.set_title(vn2+' ['+V_units[vn2]+']')
            elif ir==(NR-1):
                pass
            else:
                ax.set_xlabel('')
                ax.set_xticklabels([])

    # Uwind column
        if ic==1:
            vn = 'Uwind'
            # LiveOcean
            ax.plot(mdt, V[vn], 'darkblue', linewidth=2)
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
            if ir==0:
                ax.set_xlabel('')
                ax.set_xticklabels([])
                ax.set_title(vn+' ['+V_units[vn]+']')
            elif ir==(NR-1):
                pass
            else:
                ax.set_xlabel('')
                ax.set_xticklabels([])

    # Vwind column
        if ic==2:
            vn = 'Vwind'
            # LiveOcean
            ax.plot(mdt, V[vn], 'darkblue', label='LiveOcean', linewidth=2)
            # year and zero lines
            ax.axhline(linestyle='--', color='k', linewidth=1)
            ax.axvline(x=datetime(2014,1,1), color='k', linewidth=0.75)
            ax.axvline(x=datetime(2015,1,1), color='k', linewidth=0.75)
            ax.axvline(x=datetime(2016,1,1), color='k', linewidth=0.75)
            # ndbc
            DFF[vn].plot(ax=ax, xlim=(t0,t1), ylim=(minwind,maxwind), color='r', label='ndbc buoys', linewidth=2)
            # climatology
            clim_DFF[vn].plot(ax=ax, xlim=(t0,t1), ylim=(minwind,maxwind), color='k', linewidth=1, label='Climatology')
            clim_min = clim_DFF[vn].as_matrix()-clim_std_DFF[vn].as_matrix()
            clim_max = clim_DFF[vn].as_matrix()+clim_std_DFF[vn].as_matrix()
            ax.fill_between(clim_DFF.index, list(clim_min), list(clim_max), color='grey', alpha=0.3, label='1 Std of Clim')
            if ir==0:
                ax.set_xlabel('')
                ax.set_xticklabels([])
                ax.set_title(vn+' ['+V_units[vn]+']')
            elif ir==(NR-1):
                ax.legend(bbox_to_anchor=(0.825,1.15), ncol=2)
            else:
                ax.set_xlabel('')
                ax.set_xticklabels([])

# Property-Property Figures
    # set subplot number
    ir_prop = int(cc_prop/3)
    ic_prop = int(cc_prop-3*ir_prop)
    
    # Re-pull data
    DFF = ndbc_lo_df[sn]
    skill_df = sn_skill[sn]
    
    # regression line data
    x_temp = np.arange(mintemp,maxtemp+1)
    y_temp = skill_df[ex_name+'_WTMP'].ix['slope'] * x_temp + skill_df[ex_name+'_WTMP'].ix['yint']
    
    x_wind = np.arange(minwind,maxwind+1)
    y_Uwind = skill_df[ex_name+'_Uwind'].ix['slope'] * x_wind + skill_df[ex_name+'_Uwind'].ix['yint']
    y_Vwind = skill_df[ex_name+'_Vwind'].ix['slope'] * x_wind + skill_df[ex_name+'_Vwind'].ix['yint']
    
    # WTMP
    vn = 'WTMP'
    ax=axis_4[ir_prop,ic_prop]
    ax.scatter(DFF[vn], DFF[ex_name+'_'+vn], s=15, color='darkblue', label='LiveOcean')
    ax.plot(x_temp, y_temp, 'darkblue') # model best fit line
    ax.plot([mintemp,maxtemp], [mintemp,maxtemp], '--k') # line y=x for comparison
    ax.set_xlim([mintemp,maxtemp])
    ax.set_ylim([mintemp,maxtemp])
    ax.set_xlabel('NDBC Data')
    ax.set_ylabel('Model Data')
    ax.set_title('Station ' + sn)
    # statistics
    ax.text(.05, .9, 'MB = %.2f' %skill_df[ex_name+'_'+vn].ix['MB'], transform=ax.transAxes)
    ax.text(.05, .85, 'SDE = %.2f' %skill_df[ex_name+'_'+vn].ix['SDE'], transform=ax.transAxes)
    ax.text(.05, .8, 'CC = %.2f' %skill_df[ex_name+'_'+vn].ix['CC'], transform=ax.transAxes)
    ax.text(.95, .2, 'R$^2$ = %.2f' %skill_df[ex_name+'_'+vn].ix['R2'], ha='right', transform=ax.transAxes)
    ax.text(.95, .15, 'Slope Offest = %.2f' %abs(skill_df[ex_name+'_'+vn].ix['slope']-1), ha='right', transform=ax.transAxes)
    ax.text(.95, .1, 'RMSE = %.2f' %skill_df[ex_name+'_'+vn].ix['RMSE'], ha='right', transform=ax.transAxes)
    ax.text(.95, .05, 'SS = %.2f' %skill_df[ex_name+'_'+vn].ix['SS'], ha='right', transform=ax.transAxes)
    if ic_prop != 0:
        ax.set_ylabel('')
    if ir_prop == 0:
        ax.set_xlabel('')

    # Uwind
    vn = 'Uwind'
    ax=axis_5[ir_prop,ic_prop]
    ax.scatter(DFF[vn], DFF[ex_name+'_'+vn], s=15, color='darkblue', label='LiveOcean')
    ax.plot(x_wind, y_Uwind, 'darkblue') # model best fit line
    ax.plot([minwind,maxwind], [minwind,maxwind], '--k') # line y=x for comparison
    ax.set_xlim([minwind,maxwind])
    ax.set_ylim([minwind,maxwind])
    ax.set_xlabel('NDBC Data')
    ax.set_ylabel('Model Data')
    ax.set_title('Station ' + sn)
    # statistics
    ax.text(.05, .9, 'MB = %.2f' %skill_df[ex_name+'_'+vn].ix['MB'], transform=ax.transAxes)
    ax.text(.05, .85, 'SDE = %.2f' %skill_df[ex_name+'_'+vn].ix['SDE'], transform=ax.transAxes)
    ax.text(.05, .8, 'CC = %.2f' %skill_df[ex_name+'_'+vn].ix['CC'], transform=ax.transAxes)
    ax.text(.95, .2, 'R$^2$ = %.2f' %skill_df[ex_name+'_'+vn].ix['R2'], ha='right', transform=ax.transAxes)
    ax.text(.95, .15, 'Slope Offest = %.2f' %abs(skill_df[ex_name+'_'+vn].ix['slope']-1), ha='right', transform=ax.transAxes)
    ax.text(.95, .1, 'RMSE = %.2f' %skill_df[ex_name+'_'+vn].ix['RMSE'], ha='right', transform=ax.transAxes)
    ax.text(.95, .05, 'SS = %.2f' %skill_df[ex_name+'_'+vn].ix['SS'], ha='right', transform=ax.transAxes)
    if ic_prop != 0:
        ax.set_ylabel('')
    if ir_prop == 0:
        ax.set_xlabel('')


    # Vwind
    vn = 'Vwind'
    ax=axis_6[ir_prop,ic_prop]
    ax.scatter(DFF[vn], DFF[ex_name+'_'+vn], s=15, color='darkblue', label='LiveOcean')
    ax.plot(x_wind, y_Vwind, 'darkblue') # model best fit line
    ax.plot([minwind,maxwind], [minwind,maxwind], '--k') # line y=x for comparison
    ax.set_xlim([minwind,maxwind])
    ax.set_ylim([minwind,maxwind])
    ax.set_xlabel('NDBC Data')
    ax.set_ylabel('LiveOcean Data')
    ax.set_title('Station ' + sn)
    # statistics
    ax.text(.05, .9, 'MB = %.2f' %skill_df[ex_name+'_'+vn].ix['MB'], transform=ax.transAxes)
    ax.text(.05, .85, 'SDE = %.2f' %skill_df[ex_name+'_'+vn].ix['SDE'], transform=ax.transAxes)
    ax.text(.05, .8, 'CC = %.2f' %skill_df[ex_name+'_'+vn].ix['CC'], transform=ax.transAxes)
    ax.text(.95, .2, 'R$^2$ = %.2f' %skill_df[ex_name+'_'+vn].ix['R2'], ha='right', transform=ax.transAxes)
    ax.text(.95, .15, 'Slope Offest = %.2f' %abs(skill_df[ex_name+'_'+vn].ix['slope']-1), ha='right', transform=ax.transAxes)
    ax.text(.95, .1, 'RMSE = %.2f' %skill_df[ex_name+'_'+vn].ix['RMSE'], ha='right', transform=ax.transAxes)
    ax.text(.95, .05, 'SS = %.2f' %skill_df[ex_name+'_'+vn].ix['SS'], ha='right', transform=ax.transAxes)
    if ic_prop != 0:
        ax.set_ylabel('')
    if ir_prop == 0:
        ax.set_xlabel('')
    cc_prop += 1



# Save figures
fig_1.savefig(savname + 'Time_Series_1.png', bbox_inches='tight')
fig_2.savefig(savname + 'Time_Series_2.png', bbox_inches='tight')
fig_3.savefig(savname + 'Time_Series_3.png', bbox_inches='tight')
fig_4.savefig(savname + 'Property_1.png', bbox_inches='tight')
fig_5.savefig(savname + 'Property_2.png', bbox_inches='tight')
fig_6.savefig(savname + 'Property_3.png', bbox_inches='tight')

