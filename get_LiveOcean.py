# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 18:09:30 2016

@author: Bradley

Retrieve Live Ocean data from online.
"""

import os
import sys
alp = os.path.abspath('./alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import shutil
import pandas as pd
import urllib.request as U

# Choose directory
Ldir = Lfun.Lstart()
dirname = Ldir['roms'] + '/output/cascadia1_base_lo1/'
Lfun.make_dir(dirname, clean=False)

# Date range
date_list = pd.date_range('2016-7-01','2016-7-2')
for dt in date_list:

# Retrieve data
    try:
        date = dt.strftime('%Y%m%d')
        Lfun.make_dir(os.path.join(dirname, 'f' + dt.strftime('%Y.%m.%d')))
        for hh in range(2,26):
            hhhh = ('0000' + str(hh))[-4:]
            print('Attempting to get ' + date + ' ' + str(hhhh))
            url_str = ('https://pm2.blob.core.windows.net/f' + date + '/ocean_his_' + hhhh + '.nc')
            fn = dirname + 'f' + dt.strftime('%Y.%m.%d') + '/ocean_his_' + hhhh + '.nc'
            with U.urlopen(url_str) as response, open(fn, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
    except:
        print('Could not get ' + str(date))
        continue
