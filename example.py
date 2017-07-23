from meteo-utils.field_meteo import Field2D
from netCDF4 import Dataset

# read temperature data
f_nc = "Z_NWGD_C_BABJ_20170719070513_P_RFFC_SCMOC-TMP_201707190800_24003.nc"
fh = Dataset(f_nc, mode='r')
var_name = "TMP_2maboveground"
lon_name = "longitude"
lat_name = "latitude"
vars = fh.variables[var_name][:] - 273.15
lons = fh.variables[lon_name][:]
lats = fh.variables[lat_name][:]

# instantiation -- create an object: tmp2m
tmp2m = Field2D('tmp2m',vars[0, :, :], lons, lats)

# actions of tmp2m

# smooth tmp2m using Gaussian with default sigma((sigma_x=10.0, sigma_y=10.0))
tmp2m_smooth_1 = tmp2m.smoothing_gaussian()
# change sigma
tmp2m_smooth_2 = tmp2m.smoothing_gaussian(sigma_x=30.0, sigma_y=30.0)
# take values every other grid(default)
tmp2m_rare_1, lats_rare_1, lons_rare_1 = tmp2m.rarefaction()
# take values every three grid,
tmp2m_rare_2, lats_rare_2, lons_rare_2 = tmp2m.rarefaction(interval_x=3, interval_y=3)
# visuallize - 1
tmp2m.filled_field()
# visuallize - 2
tmp2m.filled_field(drawmap=True,
                   contourf_set={'cmap':'jet', 'extend':"both", 'alpha':1},
                   savefig_set={'fname':'./test.jpg', 'other':{'dpi':100}},figsize=(16, 12))
# visuallize - 3
tmp2m.filled_field(drawmap=True,
                   contourf_set={'cmap':'jet', 'extend':"both", 'alpha':0.6},
                   savefig_set={'fname':'./test.jpg'})
# visuallize - 4
colors = [(117, 251, 253), (110, 227, 159), (125, 189, 68), (200, 230, 82), (254, 255, 145),
          (254, 255, 85), (248, 218, 73), (243, 179, 62), (238, 124, 48), (233, 51, 36),
          (187, 39, 26), (117, 21, 45), (189, 72, 110), (234, 51, 247), (239, 135, 248)]
clevs = [x / 10.0 for x in range(0, 500, 5)]
clip_set={'shpfile':'./shapfiles/country1', 'region':'China', 'region_code':3}
tmp2m.filled_field(contourf_set={'cmap':'jet', 'extend':"both", 'levels':clevs},
                   clip_set=clip_set,
                   nobar=False,
                   drawmap=False,
                   nobounding=True,
                   custom_cmap_set = {"colors":colors, "cmap_color_num":256},
                   savefig_set={'fname':"./test.png", 'other':{'dpi':100}}, figsize=(16, 12))






