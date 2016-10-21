# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 08:14:12 2016

@author: Bradley

Performing a skills test on the LiveOcean and ndbc comparison.
Use ptools/ndbc/get_ndbc.py for retrieval of ndbc data.
Use ptools/ndbc/process_ndbc.py for processing ndbc data.
On Fjord, use LiveOcean/moor/moor_ndbc.py for retrieval of LiveOcean data.
On personal computer, us LiveOcean/get_LiveOcean.py for retrieval of LO data 
    from online, then use LiveOcean/moor/moor_ndbc.py.
"""

import os
import sys
alp = os.path.abspath('./alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun
import Lfun
import pickle
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime
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
        list_type + '_' + date_string0 + '_' + 
        date_string1 + '/')
indir2 = (Ldir['LOo'] + 'moor/' + Ldir['gtagex2'] + '_' + 
        list_type + '_' + date_string02 + '_' + 
        date_string12 + '/')

# Set directories and station list from process_ndbc.py
# Choose time filter
tf = 'w' # 'm', 'w', or 'd'

# Find home directory
which_home = os.environ.get("HOME")
if which_home == '/Users/PM5': # Mac
    dirname = which_home + '/Documents/tools_data/obs_data/ndbc/'
    savname = which_home + '/Documents/LiveOcean_output/ndbc/'
elif which_home == '/home/parker': # Fjord
    dirname = '/data1/parker/tools_data/obs_data/ndbc/'
    savname = '/data1/parker/LiveOcean_output/ndbc/'
elif which_home == '/home/bbartos': # Bradley's Fjord
    dirname = which_home + '/tools_data/obs_data/ndbc/'
    savname = which_home + '/LiveOcean_output/ndbc/'
elif which_home == None: # Windows version
    which_home = os.path.expanduser("~")
    dirname = which_home.replace('\\','/') + '/Documents/Research Work/Parker/tools_data/obs_data/ndbc/'
    savname = which_home.replace('\\','/') + '/Documents/Research Work/Parker/LiveOcean_data/ndbc/'
else:
    print('Trouble filling out environment variables')
Lfun.make_dir(savname)

# open ndbc data dictionary
fn = open(os.path.join(dirname,('ndbc_df_' + tf + '.txt')),'rb')
ndbc_df = pickle.load(fn)

# ndbc stations for standard LiveOcean
sn_list = ['46088','46087','46041','46029','46089','46050']
sn_list_1 = sn_list[0:2]
sn_list_2 = sn_list[2:]

# For variable stations
#sn_list = list(ndbc_df.keys())
#sn_list_1 = sn_list[0:int(0.5*len(sn_list)+1)]
#sn_list_2 = sn_list[int(0.5*len(sn_list)+1):]

# Begin station loop
tf_ndbc_lo = dict()
tf_ndbc_lo2 = dict()
sn_skill = dict()
for sn in sn_list:

# load Data

# LiveOcean
    try:
        fn = (indir + Ldir['gtagex'] + '_' + sn + '_' + 
            list_type + '_' + date_string0 + '_' + 
            date_string1 + '.nc')
        ds = nc.Dataset(fn)
        fn2 = (indir2 + Ldir['gtagex2'] + '_' + sn + '_' + 
            list_type + '_' + date_string02 + '_' + 
            date_string12 + '.nc')
        ds2 = nc.Dataset(fn2)
    except:
        print('LiveOcean data for Station ' + sn + ' and time ' + date_string0 + 
            '_' + date_string1 + ' not found.')
        pass

    # Create dictionary of variables
    V = dict()
    V2 = dict()
    V_units = dict()
    V_units['Uwind'] = 'ms$^{-1}$'
    V_units['Vwind'] = 'ms$^{-1}$'
    V_units['temp'] = u'\N{DEGREE SIGN}C'
    list_to_plot = ['Uwind','Vwind','temp']  # adjust to switch variables
    ltp = list_to_plot.copy()
    ltp.append('ocean_time')
    for vv in ltp:
        if ds.variables[vv].ndim == 2:
            V[vv] = ds.variables[vv][-1,:]  # for variables with z dimension
            V2[vv] = ds2.variables[vv][-1,:]
        else:
            V[vv] = ds.variables[vv][:]
            V2[vv] = ds2.variables[vv][:]
        
        # Filter data
        if vv == 'ocean_time':
            continue
        else:
            if list_type == 'low_pass':
                if tf == 'm':
                    V[vv] = zfun.filt_hanning(V[vv], n=30)
                    V2[vv] = zfun.filt_hanning(V2[vv], n=30)
                elif tf == 'w':
                    V[vv] = zfun.filt_hanning(V[vv], n=7)
                    V2[vv] = zfun.filt_hanning(V2[vv], n=7)
                else:  # since data is already daily
                    pass 
            if list_type == 'backfill':
                if tf == 'm':
                    V[vv] = zfun.filt_hanning(V[vv], n=720)
                    V2[vv] = zfun.filt_hanning(V2[vv], n=720)
                elif tf == 'w':
                    V[vv] = zfun.filt_hanning(V[vv], n=168)
                    V2[vv] = zfun.filt_hanning(V2[vv], n=168)
                elif tf == 'd': 
                    V[vv] = zfun.filt_hanning(V[vv], n=24)
                    V2[vv] = zfun.filt_hanning(V2[vv], n=24)


   # Saving filtered data
    savfn = open(os.path.join(indir + tf + 'filter_' + Ldir['gtagex'] + '_' + 
            sn + '_' + list_type + '_' + date_string0 + '_' + 
            date_string1 + '.nc'), 'wb')
    pickle.dump(V,savfn)
    savfn2 = open(os.path.join(indir2 + tf + 'filter_' + Ldir['gtagex2'] + '_' + 
            sn + '_' + list_type + '_' + date_string02 + '_' + 
            date_string12 + '.nc'), 'wb')
    pickle.dump(V2,savfn2)

# ndbc data
    try:
        DFF = ndbc_df[sn]
    except:
        print('ndbc data for Station ' + sn + ' and time ' + date_string0 + 
            '_' + date_string1 + ' not found.')
        pass

# Add LiveOcean data to ndbc dataset
    DFF_comp = DFF.reindex(index=pd.date_range(t0,t1))
    for vn in V:
        vn_ser = pd.Series(V[vn],index=pd.date_range(t0,t1))
        vn_ser2 = pd.Series(V2[vn],index=pd.date_range(t02,t12))
        if vn=='ocean_time':
            continue
        if vn=='temp': 
            vn_ndbc = 'WTMP'
        else: 
            vn_ndbc = vn
        DFF_comp[ex_name + '_' + vn_ndbc] = np.nan
        DFF_comp[ex_name2 + '_' + vn_ndbc] = np.nan
        try:
            for x in range(len(vn_ser)):
                if np.isnan(DFF_comp[vn_ndbc][vn_ser.index[x]]) == False:
                    DFF_comp['lo1_'+vn_ndbc][vn_ser.index[x]] = vn_ser[x]
                    DFF_comp['lobio1_'+vn_ndbc][vn_ser.index[x]] = vn_ser2[x]
        except KeyError:
            print('End of ' + vn + ' for Station ' + sn)
            pass
    tf_ndbc_lo[sn] = DFF_comp
    
# Skill Tests
    vn_ndbc = ['WTMP', 'Uwind', 'Vwind']
    vn_lo1 = ['lo1_WTMP', 'lo1_Uwind', 'lo1_Vwind']
    vn_lobio1 = ['lobio1_WTMP', 'lobio1_Uwind', 'lobio1_Vwind']
    means = dict()
    std = dict()
    anomaly = dict()
    sqrmean = dict()
    persist = dict()
    mse_r = dict()
    # mean, standard deviation, and anomaly
    for vn in (vn_ndbc+vn_lo1+vn_lobio1):
        means[vn] = np.nanmean(DFF_comp[vn])
        std[vn] = np.nanstd(DFF_comp[vn])
        anomaly[vn] = DFF_comp[vn] - means[vn]
    # persistance anomaly for reference forecast
        anom = np.ones(len(DFF_comp[vn]))*np.nan
        c = 0
        for t in range(len(DFF_comp[vn])):
            if c == 0 and np.isnan(DFF_comp[vn][t]) == True:
                continue
            elif c == 0 and np.isnan(DFF_comp[vn][t]) == False:
                x0 = DFF_comp[vn][t]
                c = 1
            elif c == 1:
                anom[t] = (x0-means[vn]) * (DFF_comp[vn][t]-means[vn])
    # persistance correlation and reference Mean Square Error
        persist[vn] = np.nanmean(anom)/std[vn]**2
        mse_r[vn] = (1-persist[vn]**2) * std[vn]**2
 
    mb = dict()
    sde = dict()
    cc = dict()
    londbc = dict()
    slope = dict()
    yint = dict()
    r2 = dict()
    mse = dict()
    sqrtmse = dict()
    ss = dict()
    skill_df = pd.DataFrame(index=['mean','std','MB','SDE','CC','slope','yint',
                'R2','MSE','RMSE','SS'], columns=[vn_ndbc+vn_lo1+vn_lobio1])
    for vn in vn_ndbc:
    # Mean Bias
        mb['lo1_'+vn] = abs(means['lo1_'+vn] - means[vn])
        mb['lobio1_'+vn] = abs(means['lobio1_'+vn] - means[vn])
    # Standard Deviation Error
        sde['lo1_'+vn] = abs(std['lo1_'+vn] - std[vn])
        sde['lobio1_'+vn] = abs(std['lobio1_'+vn] - std[vn])
    # Cross-Correlation
        cc['lo1_'+vn] = (std['lo1_'+vn]**-1 * std[vn]**-1
                    * np.nanmean(anomaly['lo1_'+vn]*anomaly[vn]))
        cc['lobio1_'+vn] = (std['lobio1_'+vn]**-1 * std[vn]**-1
                    * np.nanmean(anomaly['lobio1_'+vn]*anomaly[vn]))
    # Linear Regression
        sqrmean[vn] = np.nanmean(DFF_comp[vn]**2)
        londbc['lo1_'+vn] = np.nanmean(DFF_comp[vn] * DFF_comp['lo1_'+vn])
        londbc['lobio1_'+vn] = np.nanmean(DFF_comp[vn] * DFF_comp['lobio1_'+vn])
        slope['lo1_'+vn] = (means[vn] * means['lo1_'+vn] - londbc['lo1_'+vn]) / (means[vn]**2 - sqrmean[vn])
        slope['lobio1_'+vn] = (means[vn] * means['lobio1_'+vn] - londbc['lobio1_'+vn]) / (means[vn]**2 - sqrmean[vn])
        yint['lo1_'+vn] = means['lo1_'+vn] - slope['lo1_'+vn] * means[vn]
        yint['lobio1_'+vn] = means['lobio1_'+vn] - slope['lobio1_'+vn] * means[vn]
    # Coefficient of Determination (R^2)
        r2['lo1_'+vn] = 1 - (np.nansum((DFF_comp['lo1_'+vn] - slope['lo1_'+vn] * DFF_comp[vn] - yint['lo1_'+vn])**2) / np.nansum((DFF_comp['lo1_'+vn] - means['lo1_'+vn])**2))
        r2['lobio1_'+vn] = 1 - (np.nansum((DFF_comp['lobio1_'+vn] - slope['lobio1_'+vn] * DFF_comp[vn] - yint['lobio1_'+vn])**2) / np.nansum((DFF_comp['lobio1_'+vn] - means['lobio1_'+vn])**2))
    # Mean Square Error (from Oke, 2001)
        mse['lo1_'+vn] = (mb['lo1_'+vn]**2 + sde['lo1_'+vn]**2 + 
            2*std['lo1_'+vn]*std[vn]*(1-cc['lo1_'+vn]))
        mse['lobio1_'+vn] = (mb['lobio1_'+vn]**2 + sde['lobio1_'+vn]**2 + 
            2*std['lobio1_'+vn]*std[vn] * (1-cc['lobio1_'+vn]))
        sqrtmse['lo1_'+vn] = np.sqrt(mse['lo1_'+vn])
        sqrtmse['lobio1_'+vn] = np.sqrt(mse['lobio1_'+vn])
    # Skill Score (from Oke, 2001, and Murphy, 1992)
        ss['lo1_'+vn] = 1 - (mse['lo1_'+vn]/mse_r[vn])
        ss['lobio1_'+vn] = 1 - (mse['lobio1_'+vn]/mse_r[vn])

    skill_df.ix['mean'] = means
    skill_df.ix['std'] = std
    skill_df.ix['MB'] = mb
    skill_df.ix['SDE'] = sde
    skill_df.ix['CC'] = cc
    skill_df.ix['slope'] = slope
    skill_df.ix['yint'] = yint
    skill_df.ix['R2'] = r2
    skill_df.ix['MSE'] = mse
    skill_df.ix['RMSE'] = sqrtmse
    skill_df.ix['SS'] = ss
    sn_skill[sn] = skill_df

# save skill dictionary
savfn = open((savname + tf + '_skill_dictionary.txt'), 'wb')
pickle.dump(sn_skill, savfn)

savfn2 = open((savname + tf + '_ndbc_lo_data.txt'), 'wb')
pickle.dump(tf_ndbc_lo, savfn2)

savfn.close()
savfn2.close()