"""
Code for particle tracking of rockfish larvae.
Derivative of tracker_b.py

With rockfish experiment, larvae are released at one location 
over two months and following the lunar cycle. 
Then, experiment continues to run for four months.
Will run for years listed in yr_list, default is 2006.
If changing the years, update the lun_dt_ind as well.

Designed to read csv files. One needs to contain the columns:
"Experiment", "model", "species", "lat", "lon", and "depth (m)",.
The other must contain: "days", and "# of particles".

The experiment will be chosen by either an argument >tracker_rf.py -ex 1_1
or an input prompt after script has started to run
"""

#%% setup
import argparse
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timedelta
import time
import pickle

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()
# Note that we override parts of the Ldir logic,
# using Ldir['gtagex'] to identify a run.

from importlib import reload
import zfun 
reload(zfun)
import trackfun
reload(trackfun)

import warnings
warnings.simplefilter('error')

# some run specifications
ic_name = 'rockfish'
dir_tag = 'forward' # 'forward' or 'reverse'
method = 'rk4' # 'rk2' or 'rk4'
surface = False # Boolean, True for trap to surface
turb = True # Boolean, True to include vertical turbulence
windage = 0 # a small number >= 0
ndiv = 1 # number of divisions to make between saves for the integration
        # e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
        # for the integration (but still only report fields hourly)
bound = 'reflect' # either 'stop' for halting at coast or 'reflect to bounce off
dep_max = -100 # maximum depth, set to None if no max is required

# create possible inputs
parser = argparse.ArgumentParser()
parser.add_argument('-e', '--experiment', nargs='?', type=str)
parser.add_argument('-t', '--testing', action='store_true') # Use to limit the number of years, particles, and timesteps, default is False
args = parser.parse_args()
testing = args.testing

# retrieve experimental data
exdf = pd.read_csv(Ldir['data'] + 'tracker/rockfish_latlon.csv', index_col = 0)
pardf = pd.read_csv(Ldir['data'] + 'tracker/rockfish_partruition.csv', index_col=0)

# choose experiment if not already provided as an argument
if args.experiment != None:
    exrow = args.experiment
else:
    print('\n%s\n' % '** Choose Experiment **')
    for ind in exdf.index:
        print(ind+'  '+str(exdf.ix[ind][0])+'  '+str(exdf.ix[ind][1]))
    exrow = str(input('-- Input Experiment -- '))

# set number of particles
if testing:
    NP0 = 100000
else:
    NP0 = exdf.ix[exrow]['# of particles']

# set particle initial locations
plon00 = np.array([exdf.ix[exrow]['lon']])
plat00 = np.array([exdf.ix[exrow]['lat']])
dep_orig = exdf.ix[exrow]['depth (m)']
# depth array has full number of larvae
dep00 = np.linspace(dep_orig, dep_orig, NP0)

# set model and species
gtagex = exdf.ix[exrow]['model']
species = exdf.ix[exrow]['species']

# set number of days
if testing:
    days_to_track = 1
else:
    days_to_track = 180 # standard run will be 6 months

# make the full IC vectors, which will have equal length
NSP = len(dep00)
NXYP = len(plon00)
plon0 = plon00.reshape(NXYP,1) * np.ones((NXYP,NSP))
plat0 = plat00.reshape(NXYP,1) * np.ones((NXYP,NSP))
dep0 = dep00.reshape(1,NSP) * np.ones((NXYP,NSP))
plon0 = plon0.flatten()
plat0 = plat0.flatten()
dep0 = dep0.flatten()

# years to track
yr_list = [2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009]
# start dates for each year coincide with new moon
if species == 'canary':
    lunar_dt_dict = dict(zip(yr_list, pd.to_datetime(['2002/01/13', 
                  '2003/01/02', '2004/01/21', '2005/01/10', '2006/01/29',
                  '2007/01/18', '2008/01/08', '2009/01/25'])))
elif species == 'yelloweye':
    lunar_dt_dict = dict(zip(yr_list, pd.to_datetime(['2002/05/12', 
                  '2003/05/01', '2004/05/18', '2005/05/08', '2006/05/26', 
                  '2007/05/16', '2008/05/05', '2009/05/24'])))

# choose year for testing or full list
if testing:
    yr_list = [2006,]
else:
    pass

# MoSSea model only has one year
if gtagex == 'MoSSea':
    yr_list = [2006,]

# save some things in Ldir
Ldir['gtagex'] = gtagex
Ldir['ic_name'] = ic_name
Ldir['dir_tag'] = dir_tag
Ldir['method'] = method
Ldir['surface'] = surface
Ldir['turb'] = turb
Ldir['windage'] = windage
Ldir['ndiv'] = ndiv
Ldir['experiment'] = exrow
Ldir['days_to_track'] = days_to_track

# track in year loop
for yr in yr_list:
    print('Working on ' + str(yr))

    idt = lunar_dt_dict[yr]
    fn_list = trackfun.get_fn_list(idt, Ldir, yr=yr)

    # create release index array
    if testing:
        mask_ini = np.zeros(NP0)
    else:
        mask_ini = np.zeros(int(pardf.iloc[0]['# of particles']),dtype='int8')
        for rel in pardf.index[1:]:
            pnum = int(pardf.ix[rel]['# of particles'])
            mask_ini = np.concatenate((mask_ini, np.ones(pnum,dtype='int8')*rel))
    # create mask of staggered release
    fn_mask = trackfun.get_fn_mask(fn_list, mi=mask_ini)
    
    # pull start dataset
    ds0 = nc.Dataset(fn_list[0])

    # make vectors to feed to interpolant maker
    G, S = zfun.get_basic_info(fn_list[0], getT=False)
    R = dict()
    R['rlonr'] = G['lon_rho'][0,:].squeeze()
    R['rlatr'] = G['lat_rho'][:,0].squeeze()
    R['rlonu'] = G['lon_u'][0,:].squeeze()
    R['rlatu'] = G['lat_u'][:,0].squeeze()
    R['rlonv'] = G['lon_v'][0,:].squeeze()
    R['rlatv'] = G['lat_v'][:,0].squeeze()
    R['rcsr'] = S['Cs_r'][:]
    R['rcsw'] = S['Cs_w'][:]

    # when initial position is on the boundary, need to move one 
    # grid cell towards ocean so velocities are not zero
    pcs_temp = np.array([0])
    plon0, plat0 = trackfun.change_position(ds0, plon0, plat0, pcs_temp, R)

    # create initial pcs from depths
    ZH0 = trackfun.get_V(['zeta', 'h'], ds0, plon0, plat0, pcs_temp, R)
    Tot_Dep = ZH0[0,0] + ZH0[0,1]
    pcs0 = -dep0/Tot_Dep
    
#%% DO THE TRACKING
    tt0 = time.time()
    print(' - Starting on ' + time.asctime())

    P, G, S = trackfun.get_tracks(fn_list, plon0, plat0, pcs0, dir_tag,
                                  method, surface, turb, ndiv, windage, 
                                  bound=bound, dep_max=dep_max, fn_mask=fn_mask)

    print(' - Took %0.1f sec for %d days'
          % (time.time() - tt0, Ldir['days_to_track']))
            
#%% save the results

    # define the output directory
    outdir = (Ldir['LOo'] + 'tracks/' + Ldir['gtagex'] + '_' +
        Ldir['ic_name'] + '_' + Ldir['method'] + '_' + 'ndiv' + 
        str(Ldir['ndiv']) + '_' + Ldir['dir_tag'] + '_' + 'surface' + 
        str(Ldir['surface']) + '_' + 'turb' + str(Ldir['turb']) + '_' + 
        'windage' + str(Ldir['windage']) + '_boundary' + bound + '_max_depth' + str(abs(dep_max)) + '/')
    Lfun.make_dir(outdir)
    #Lfun.make_dir(outdir, clean=True) # use to wipe output directory

    outname = (Ldir['ic_name'] + '_' + 'Experiment_' + Ldir['experiment']
        + '_' + str(yr) + '.p')

    pickle.dump((P, G, S, Ldir), open(outdir + outname, 'wb'))
    print('Results saved to:\n' + outdir + outname)
    print(50*'*')
