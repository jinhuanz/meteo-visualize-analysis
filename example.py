import os
import gc
import sys
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset
from utils.field_meteo import Field2D
import psutil

# log memory use
process = psutil.Process(os.getpid())
def print_memory_usage():
    nr_mbytes = process.memory_info()[0] / 1048576.0
    sys.stdout.write("{}\n".format(nr_mbytes))
    sys.stdout.flush()
    return nr_mbytes


def aqi_p_hourly(element):

    print("visualizing %s ..." % (element))
    # step 1: download aqi forecast file
    fid_aqi = 'aqi_daily.nc'
    # step 2: read
    fid_latlon = "lonlat.nc"
    f_aqi = Dataset(fid_aqi, 'r')
    f_latlon = Dataset(fid_latlon, 'r')
    lats = f_latlon.variables['lat'][:]
    lons = f_latlon.variables['lon'][:]
    aqi_data = f_aqi.variables[element][:]

    # step 3: visualize
    dt0 = datetime.now() 
    dt = dt0

    # six colors [c1, c2, c3, c4,c5, c6] according to GB
    # [(15, 196, 142), (255, 214, 44), (255, 163, 0), (240, 105, 30), (221, 54, 29), (145, 13, 49)]
    colors = [(15, 196, 142), (15, 196, 142), (255, 214, 44), (255, 163, 0),
              (240, 105, 30), (221, 54, 29), (145, 13, 49), (145, 13, 49)]
    if element == 'aqi':
        #levs_base: [0, 50, 100, 150, 200, 300, 500] according to GB
        # 0-50: c1; 50-100: c2; 100-150: c3; 150-200: c4; 200-300: c5; 300-500: c6;
        # 0-25: c1->c1; 25-75: c1->c2 linear; 75-125: c2->c3 linear; 125-175: c3->c4 linear;
        # 175-225: c4->c5 linear; 225-325: c5-c6 linear; 325-500 c6->c6;
        # using mpl.colors.LinearSegmentedColormap
        # 500- beyond, color 6 using extend="max"
        levs = [0, 25, 75, 125, 175, 225, 325, 500]
    elif element == 'pm2_5':
        #levs_base = [0, 35, 75, 115, 150, 250, 500] according to GB
        levs = [0, 20, 55, 95, 135, 200, 275, 500]
    elif element == 'pm10':
        #levs_base = [0, 50, 150, 250, 350, 420, 600] according to GB
        levs = [0, 25, 100, 200, 300, 385, 450, 600]
    elif element == 'o3':
        #levs_base = [0, 160, 200, 300, 400, 800, 1200] according to GB
        levs = [0, 80, 180, 250, 350, 600, 1000, 1200]
    elif element == 'co':
        #levs_base = [0, 2, 4, 14, 24, 36, 60] according to GB
        levs = [0, 1, 3, 9, 19, 30, 42, 60]
    elif element == 'so2':
        #levs_base = [0, 50, 150, 475, 800, 1600, 2620] according to GB
        levs = [0, 25, 100, 300, 625, 1200, 2000, 2620]
    elif element == 'no2':
        #levs_base = [0, 40, 80, 180, 280, 565, 940] according to GB
        levs = [0, 20, 60, 130, 230, 410, 700, 940]
    levels = np.linspace(0, levs[-1], 250)
    clip_set = {'shpfile': './shapfiles/country1', 'region': 'China', 'region_code': 3}

    for i in range(1):
        print(i)
        aqi = Field2D(element, aqi_data[i, :, :], lons, lats, 'curvilinear')
        png = "%s_forecast/%s/%s.png" % (element, dt.strftime("%Y%m%d"), dt.strftime("%Y%m%d%H"))
        os.system("mkdir -p %s_forecast/%s" % (element, dt.strftime("%Y%m%d")))
        aqi.filled_field(contourf_set={'levels': levels, 'extend': 'both'},
                         clip_set=clip_set,
                         drawmap=False,
                         #nobar=True,
                         #nobounding=True,
                         #nomargin=True,
                         #nobackground=True,
                         custom_cmap_set={"colors": colors, "cmap_color_num": 256,
                                          "position":[x/levs[-1] for x in levs]},
                         figsize=(9.6, 7.2),
                         savefig_set={'fname': png})
        mem_s = print_memory_usage()
        del aqi, png
        gc.collect()
        mem_collect = print_memory_usage() - mem_s
        print("memory collected after each aqi prediction time step is %s" % (mem_collect))
        dt = dt + timedelta(hours=1)




#aqi_p_hourly('aqi')
#aqi_p_hourly('pm2_5')
#aqi_p_hourly('o3')
#aqi_p_hourly('pm10')
aqi_p_hourly('co')
aqi_p_hourly('so2')
aqi_p_hourly('no2')
