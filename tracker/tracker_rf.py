"""
Code for particle tracking of rockfish larvae.
Derivative of tracker_b.py

With rockfish experiment, larvae are released at one location 
over two months and following the lunar cycle. 
Then, experiment continues to run for four months.
Will run for years listed in yr_list, default is 2002-2009.
If changing the years, update the lun_dt_ind as well.

Designed to read csv files. One needs to contain the columns:
"Experiment", "model", "species", "lat", "lon", and "depth (m)",.
The other must contain: "days", and "# of particles".
"""

#%% setup
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

# some run specifications
gtagex = 'cascadia1_base_lobio1' # 'cascadia1_base_lobio1', 'D2005_his','C2009'
ic_name = 'rockfish'
dir_tag = 'forward' # 'forward' or 'reverse'
method = 'rk4' # 'rk2' or 'rk4'
surface = False # Boolean, True for trap to surface
turb = True # Boolean, True to include vertical turbulence
windage = 0 # a small number >= 0
mort = 0.25 # percentage of daily mortality
ndiv = 1 # number of divisions to make between saves for the integration
        # e.g. if ndiv = 3 and we have hourly saves, we use a 20 minute step
        # for the integration (but still only report fields hourly)
days_to_track = 180 # standard run will be 6 months

# retrieve experimental data
exdf = pd.read_csv(Ldir['data'] + 'tracker/rockfish_latlon.csv', index_col = 0)
pardf = pd.read_csv(Ldir['data'] + 'tracker/rockfish_partruition.csv', index_col=0)

# choose experiment
print('\n%s\n' % '** Choose Experiment **')
for ind in exdf.index:
    print(ind+'  '+str(exdf.ix[ind][0])+'  '+str(exdf.ix[ind][1]))
exrow = str(input('-- Input Experiment -- '))

# set particle initial locations, model, and species
plon0 = np.array(exdf.ix[exrow]['lon'])
plat0 = np.array(exdf.ix[exrow]['lat'])
dep_orig = np.array(exdf.ix[exrow]['depth (m)'])
# depth array has full number of larvae
dep0 = np.linspace(dep_orig, dep_orig, 102750)
gtagex = exdf.ix[exrow]['model']
species = exdf.ix[exrow]['species']

# years to track
yr_list = [2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009]
# start dates for each new year, coincide with new moon
if species == 'canary':
    lunar_dt_dict = dict(zip(yr_list, pd.to_datetime(['2002/01/13', 
                  '2003/01/02', '2004/01/21', '2005/01/10', '2006/01/29',
                  '2007/01/18', '2008/01/08', '2009/01/25'])))
elif species == 'yelloweye':
    lunar_dt_dict = dict(zip(yr_list, pd.to_datetime(['2002/05/12', 
                  '2003/05/01', '2004/05/18', '2005/05/08', '2006/05/26', 
                  '2007/05/16', '2008/05/05', '2009/05/24'])))
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
Ldir['mort'] = mort
Ldir['ndiv'] = ndiv
Ldir['experiment'] = exrow
Ldir['days_to_track'] = days_to_track

# track in year loop
for yr in yr_list:
    idt = lunar_dt_dict[yr]
    fn_list = trackfun.get_fn_list(idt, Ldir, yr=yr)

    # create mask for staggered starts
    mask_start = np.zeros(pardf.iloc[0]['# of particles'])
    for rel in pardf.index:
        pnum = pardf.ix[rel]['# of particles']
        np.concatenate((mask_start, np.ones(pnum)*rel))

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

    # create initial pcs from depths
    ds0 = nc.Dataset(fn_list[0])
    pcs_temp = np.zeros(len(plon0))
    ZH0 = trackfun.get_V(['zeta', 'h'], ds0, plon0, plat0, pcs_temp, R)
    Tot_Dep = ZH0[0] + ZH0[1]
    pcs0 = -dep0/Tot_Dep
    
#%% DO THE TRACKING
    tt0 = time.time()

    P, G, S = trackfun.get_tracks(fn_list, plon0, plat0, pcs0, dir_tag,
                                  method, surface, turb, ndiv, windage, fn_mask=fn_mask)

    print(' - Took %0.1f sec for %d days'
          % (time.time() - tt0, Ldir['days_to_track']))

    # remove values for mortality
    NP = len(plon0)
    # create list of all possible days (assuming less than a year)
    dt0 = P['ot'][0]
    days = dt0 + np.arange(400)*86400
    # go through each timestep
    for i in range(len(P['ot'])):
        # each day, create a new mask
        if P['ot'][i] in days:
            non_mort = np.isfinite(P['lon'][i,:])
            # create random mask of all non-nan particles
            mort_mask = np.random.choice(np.arange(NP)[non_mort], int(sum(non_mort)*mort))
            # set all future values at the mask to nan
            for var in P:
                P[var][i:,mort_mask] = np.nan
            
#%% save the results

    # define the output directory
    outdir = (Ldir['LOo'] + 'tracks/' + Ldir['gtagex'] + '_' +
        Ldir['ic_name'] + '_' + Ldir['method'] + '_' + 'ndiv' + 
        str(Ldir['ndiv']) + '_' + Ldir['dir_tag'] + '_' + 'surface' + 
        str(Ldir['surface']) + '_' + 'turb' + str(Ldir['turb']) + '_' + 
        'windage' + str(Ldir['windage']) + '_mort' + str(Ldir['mort']))
    Lfun.make_dir(outdir)
    #Lfun.make_dir(outdir, clean=True) # use to wipe output directory

    outname = (Ldir['ic_name'] + '_' + 'Experiment:_' + Ldir['experiment']
        + '_' + yr + '.p')

    pickle.dump((P, G, S, Ldir), open(outdir + outname, 'wb'))
    print('Results saved to:\n' + outdir + outname)
    print(50*'*')
