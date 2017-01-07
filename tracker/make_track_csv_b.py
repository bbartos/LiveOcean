# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:43:43 2016

@author: PM5
"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import pandas as pd
import pickle  # python 3

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tracks/'

# choose the run directory
print('\n%s\n' % '** Choose mooring file to plot **')
d_list_raw = os.listdir(indir)
d_list = []
for d in d_list_raw:
    if d[-4:] == 'days':
        d_list.append(d)
Ndt = len(d_list)
d_dict = dict(zip(range(Ndt), d_list))
for ndt in range(Ndt):
    print(str(ndt) + ': ' + d_list[ndt])
my_ndt = int(input('-- Input number --'))
dirname = d_dict[my_ndt] + '/'

# create the list of run files
m_list_raw = os.listdir(indir + dirname)
m_list = []
for m in m_list_raw:
    if m[-2:] == '.p':
        m_list.append(m)
Npt = len(m_list)
m_dict = dict(zip(range(Npt), m_list))
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_npt = int(input('-- Input number (99 for all) --'))

inname_list = []
if my_npt == 99:
    for ii in m_dict.keys():
        inname_list.append(m_dict[ii])   
else:
    inname_list.append(m_dict[my_npt])
    
for inname in inname_list:
      
    P, G, S, PLdir = pickle.load( open( indir + dirname + inname, 'rb' ) )
    
    # make sure the output directory exists
    outdir = Ldir['LOo'] + 'tracks/saved_csv/'
    Lfun.make_dir(outdir)
    outdir2 = outdir + dirname
    Lfun.make_dir(outdir2)
    outname = inname[:-2]
    odir = outdir2 + outname + '/'
    Lfun.make_dir(odir)
    
    for vn in P.keys():
        a = pd.DataFrame(P[vn], index=P['ot'])
        a.index.name = 'seconds_since_1_1_1970'
        a.to_csv(odir + vn + '.csv')