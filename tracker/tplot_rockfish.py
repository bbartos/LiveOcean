# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:03:12 2016

@author: PM5
"""

"""
Plot results of tracker.
"""

# setup
import os
import sys
alp = os.path.abspath('../alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import matfun
import pickle
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

plp = os.path.abspath('../plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

Ldir = Lfun.Lstart()
indir = Ldir['LOo'] + 'tracks/'
fn_coast = Ldir['data'] + 'coast/pnw_coast_combined.mat'

# choose the run directory
print('\n%s\n' % '** Choose mooring file to plot **')
d_list_raw = os.listdir(indir)
d_list = []
for d in d_list_raw:
#    if d[-4:] == 'days':
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
    if m[-2:] == '.p':
        m_list.append(m)
Npt = len(m_list)
for npt in range(Npt):
    print(str(npt) + ': ' + m_list[npt])
my_ndt = int(input('-- Input number (99 for all) -- '))
if my_ndt == 99:
    pass
else:
    m_list = [m_list[my_ndt],]

# output directory
outdir = indir + 'plots/' + dirname
Lfun.make_dir(outdir)
#Lfun.make_dir(outdir, clean=True) # use this to clear previous plots

# create plots for each run in the run directory
if my_ndt == 99:
    plt.ioff() # use this to supress plot output
    
for inname in m_list:
    
    P, G, S, PLdir = pickle.load( open( indir + dirname + inname, 'rb' ) )
    
    # set number of times and points
    NT, NP = P['lon'].shape
    # set map range
    aa = [np.min(P['lon'])-.1, np.max(P['lon'])+.1, np.min(P['lat'])-.1, 
          np.max(P['lat'])+.1]
    # set depth contour levels
    depth_levs = [100, 200, 500, 1000, 2000, 3000]
    
    # get coastline
    cmat = matfun.loadmat(fn_coast)
    
    # PLOTTING

    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(1,2,1)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(inname)
    
    # Depth Contours
#    ax.contour(G['lon_rho'], G['lat_rho'], G['h'], depth_levs, colors='g')
    # Coastline
    ax.plot(cmat['lon'], cmat['lat'], '-k', linewidth=.5)
    # Set axis limits
    ax.axis(aa)
    # Configure axis scales
    pfun.dar(ax)
    # Make lat/lon grid
    ax.grid()
    
    # tracks
    ax.plot(P['lon'][:], P['lat'][:], '-r', linewidth=0.5, alpha=0.5)
    
    # ending points
    ax.plot(P['lon'][-1,:],P['lat'][-1,:],'ob', markersize=4, label='End')
    
    # starting points
    ax.plot(P['lon'][0,:], P['lat'][0,:], '^y', markersize=10, label='Start')

    ax.legend()
    
    # contour legend
#    ax.text(.85, .25, 'Depth Contours', horizontalalignment='center', transform=ax.transAxes, color='g')
#    dd = 1
#    for d in depth_levs:
#        ax.text(.85, .25 - .03*dd, str(d), horizontalalignment='center', transform=ax.transAxes, color='g')
#        dd += 1
    
    
    # TIME SERIES
    tdays = (P['ot'] - P['ot'][0])/86400.
    # Depth
    ax = fig.add_subplot(2,2,2)
    ax.plot(tdays, P['z'],'-', alpha=0.25)
    ax.set_ylabel('Z (m)')
    
    # Depth Contours
    # range of levels from 0 to the multiple of 25 below the lowest level
#    d_levels = np.arange(int((np.min(P['z'])-25)/25)*25, 1, 25)[::-1]
#    d_prop = np.zeros((len(d_levels),len(tdays)))
#    # loop through each day and depth level
#    for tind in np.arange(1,len(tdays)):
#        for dind in np.arange(1,len(d_levels)):
#            # previous and current depth levels
#            ddt = d_levels[dind-1]
#            ddb = d_levels[dind]
#            # boolean of particles between depths
#            td_arr = (P['z'][tind,:]>ddb)&(P['z'][tind,:]<=ddt)
#            # proportion of all alive particles between depths
#            d_prop[dind-1, tind] = sum(td_arr)/sum(np.isfinite(P['z'][tind,:]))
#    ax = fig.add_subplot(3,2,4)
#    c1 = ax.pcolor(tdays, d_levels, d_prop, cmap='OrRd')
#    cb1 = plt.colorbar(c1, fraction=0.02, pad=0.01)
#    cb1.set_label('Proportion of Alive Larvae')
#    ax.set_ylabel('Z (m)')
#    ax.grid()

    # Distance from Start
    lat_dis = np.ones(P['z'].shape)
    lon_dis = np.ones(P['z'].shape)
    for tind in range(len(tdays)):
    # change in lat/lon and total distance=sqrt(lat^2+lon^2)
        lat_dis[tind,:] = (P['lat'][tind,:] - P['lat'][0,:]) * 111
        lon_dis[tind,:] = ((P['lon'][tind,:] - P['lon'][0,:]) * 111*
                        np.cos(np.nanmean(P['lat'][tind,:])*np.pi/180))
    dis = np.sqrt(lat_dis**2 + lon_dis**2)
    # remove dead larvae
    dis_alive = np.where(np.isfinite(P['z']), dis, np.nan)
    ax = fig.add_subplot(2,2,4)
    ax.plot(tdays, dis_alive, '-', alpha=0.25)
    ax.set_ylabel('Distance From Start (km)')
    ax.set_xlabel('Days')
    ax.grid()


    # save figures
    outfn = outdir + inname[:-2] + '.png'
    plt.savefig(outfn)
    
    if my_ndt != 99:
        plt.show()

if my_ndt == 99:
    plt.close('all')
    plt.ion()
    