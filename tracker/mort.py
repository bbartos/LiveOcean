"""
mort.py takes the output dictionary from tracking experiments and applies
a random, daily mortality to the particles.

The folder and filename to apply mortality to will be chosen in the run.

Output: Dictionary with the same shape as input, but with NaNs after 
particles have stopped.
"""

import pickle
import numpy as np

import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
Ldir = Lfun.Lstart()

# percentage of particles to kill each day
mort = 25

indir = Ldir['LOo'] + 'tracks/'

# choose the run directory
print('\n%s\n' % '** Choose file to apply mortality **')
d_list_raw = os.listdir(indir)
d_list = []
for d in d_list_raw:
    d_list.append(d)
Ndt = len(d_list)
for ndt in range(Ndt):
    print(str(ndt) + ': ' + d_list[ndt])
my_ndt = int(input('-- Input number -- '))
dirname = d_list[my_ndt] + '/'

# create the list of run files
m_list_raw = os.listdir(indir + dirname)
m_list = []
for m in m_list_raw:
    if m[-2:] == '.p' and 'mort' not in m:
        m_list.append(m)
Npt = len(m_list)
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_ndt = int(input('-- Input number (99 for all) -- '))
if my_ndt == 99:
    pass
else:
    m_list = [m_list[my_ndt],]

# loop through each file
for m_f in m_list:
    fn = open(indir + dirname + m_f, 'rb')
    
    # load dictionaries from file
    P, G, S, Ldir = pickle.load(fn)
        
    # remove values for mortality
    NP = P['lon'].shape[1]
    # create list of all possible days (assuming less than a year)
    dt0 = P['ot'][0]
    days = dt0 + np.arange(400)*86400
    # go through each timestep
    for i in range(len(P['ot'])):
        # each day, create a new mask
        if P['ot'][i] in days:
            non_mort = np.isfinite(P['u'][i,:])
            # create random mask of all non-nan particles
            mort_mask = np.random.choice(np.arange(NP)[non_mort], int(sum(non_mort)*(mort/100)))
            # set all future values at the mask to nan
            for var in P:
                if var == 'ot':
                    pass
                elif var in ['lon', 'lat']:
                    P[var][i:,mort_mask] = P[var][i, mort_mask]
                else:
                    P[var][i:,mort_mask] = np.nan
            
    #%% save the results
    
    outname = (m_f[:-2] + '_mort' + str(mort) + '.p')
    
    pickle.dump((P, G, S, Ldir), open(indir + dirname + outname, 'wb'))
    print('Results saved to:\n' + indir + dirname + outname)
    print(50*'*')