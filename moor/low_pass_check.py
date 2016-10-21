import sys

def make_date_list(dt0,dt1,Ldir): # a helpful function
    del_dt = timedelta(1)
    date_list = []
    dt = dt0
    while dt <= dt1:
        date_list.append(dt.strftime('%Y.%m.%d'))
        dt = dt + del_dt
    return date_list
    
dt0 = datetime.strptime(Ldir['date_string0'],'%Y.%m.%d') # first day
dt1 = datetime.strptime(Ldir['date_string1'],'%Y.%m.%d') # last day
date_list = make_date_list(dt0,dt1,Ldir)
for dd in date_list:
    try:
        sys.stdout.flush()
        indir = Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + dd + '/'
        fn = indir + 'low_passed.nc'
        lp = open(fn)
        lp.close()
    except:
        print('No low_passed file for ' + dd)