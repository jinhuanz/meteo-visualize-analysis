#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2017 Jinhuan Zhu <jinhuanzhu@seniverse.com>. All rights reserved.

class Field2D(object):
    """
    Deal with meteorology 2-dimensional field with coordinate(rectilinear grid) of one variable,
    1) data smoothing;
    2) data rarefaction;
    3) field ploting;
    4) get the nearest grid of giving cities.

    1) Gaussian smoothing is provided now. (TODO, other filter/smoothing methods).
    2) Now, data rarefaction uses method of extractiong data every n point, eg: data[::2,::2].
       (TODO, method ...)
    3) Filled  data above map is providing using Matplotlib-Basemap,
       (TODO, add other filled visuallization method such as imshow, hexbin and contour ...)
       For more please see http://matplotlib.org/basemap/users/examples.html, and
       http://basemaptutorial.readthedocs.io/en/latest/plotting_data.html#
    4) The nearest grids are found by using cKDTree (function nearest_cKDTree)
       (TODO: mytree.query(cities,2) will give the first two points which is the closet to one city)
    """

    def __init__(self, variable, data, lons, lats):
        """
        set up meteorology 2-dimensional field with rectilinear grid.

        Args:
            variable: meteorology variable name.
            data: meteorology 2-dimensional data.
            lons: longitude, monotonical 1-dimensional data.
            lats: latitude, monotonical 1-dimensional data.
        """
        self.variable = variable
        self.data = data
        self.lons = lons
        self.lats = lats

    def nearest_ckdtree(self, cities):
        """
        we have m (n, 2) cities and we need to find time series of temperature from
        reanalysis/station obsevations whose coordinates are grids (n, 2)

        Args:
            cities: m cities of (lat, lon in each row)

        Returns:
            lists whose length is same as cities, each element of (indexes, dist) is the index of
            nearest point and the distance between the city and the point.
        """
        import scipy.spatial as spatial

        # points: n points of (lat, lon in each row)
        points = []
        for ix in self.lons:
            for iy in self.lats:
                points.append([iy, ix])
        mytree = spatial.cKDTree(points)
        dist, indexes = mytree.query(cities)

        return dist, indexes

    def rarefaction(self, interval_x=2, interval_y=2):
        """
        Simple data rarefaction by extracting data every 3 points

        Args:
            interval_x: the drop interval of longitude
            interval_y: the drop interval of latitude

        Returns:
            rarefaction variable with its coord(longitude, latitude)
        """

        return self.data[::interval_y, ::interval_x], self.lats[::interval_y], \
               self.lons[::interval_x]

    def smoothing_gaussian(self, sigma_x=10.0, sigma_y=10.0):
        """
        Gaussian smoothing

        Args:
            sigma_x/sigma_y: standard deviation of Gaussian kernel

        Returns:
            var_filter: the smoothed 2-D data
        """

        import scipy.ndimage
        import scipy as sp

        sigma = [sigma_y, sigma_x]
        var_filter = sp.ndimage.filters.gaussian_filter(self.data, sigma, mode='constant')
        return var_filter

    def filled_field(self, **kwargs):
        """
        Filled meteorology 2-dimensional field with rectilinear grid.

        If kwargs are not passing, the default:
            projection is web Mactor(epsg=3857),
            extent is (llcrnrlon=73, llcrnrlat=17, urcrnrlon=136, urcrnrlat=54),
            plotting method is contourf with cmap='jet'
            with mapboundary and coastlines.

        The following personal design is supported:
        1） map setting using parameter mapset(a dict to set basemap(parameters in Basemap),
            eg: mapset = {'llcrnrlon':73,'llcrnrlat':17,'urcrnrlon':136,'urcrnrlat':54,'epsg':3857}
            For more please see http://matplotlib.org/basemap/users/mapsetup.html,
        2) set size of figure using figsize(a tuple of inches of length and width, eg: (6.4, 4.8)
        3) plot setting such as (levels+cmap or levels+colors), extend and so on
           contourf_set:  eg: contourf_set={cmap='jet', extend="both", alpha=0.6};
                (For more http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.contourf)
           pcolormesh_set: eg: pcolormesh_set = {cmap='jet', vmin=-30, vmax=50}
                (For more http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.pcolormesh)
        4) figure attributes for figures used by Gaode.
           to add figure to Gaode map, the margin, background color and bounding box should be off,
           the map projection should be webMactor,
           and the (llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat) should be given to Gaode map.
           4.1) nomargin=True means no margin
           4.2) nobacground=True means off background color
           4.3) nobounding=True means off bounding box
           4.4) nobar=True means off colorbar
        5) Using clip_set(a dictionary) to clip figure accoding to shapfile (function shp2clip)
           eg: clip_set={'shpfile':'./shapfiles/country1', 'region':'China', 'region_code':3}
           eg: clip_set={'shpfile':'dijishi_2004',
                         'region':b'\xb1\xb1\xbe\xa9\xca\xd0',
                         'region_code':-1}
        6)  save figure by setting savefig_set, a dict
            to save figure to disk using savefig_set['fname']
            eg: savefig_set={'fname':'/home/jinhuanz/works/test.pdf'}
            eg: savefig_set= {'fname':'../test.png', 'other':{'facecolor':'w', 'edgecolor':'w'}}
            emf, eps, pdf, png, ps, raw, rgba, svg, svgz, jpg(should install pillow) are suppoted.
            For more, see http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.savefig
        7)  create own colormap using custom_cmap_set (use function make_cmap)
             eg: custom_cmap_set = {"colors":[(117, 251, 253), (110, 227, 159)]}
             eg: custom_cmap_set = {"colors":[(117, 251, 253), (110, 227, 159)],
                                    "cmap_color_num":1024}
             if custom_cmap_set sets, cmap in contourf_set or pcolormesh_set is replaced.
        8) draw boundary coastlines, countries on map
           if drawmap=True or drawmap is not setted, draw mapboundary, coastline, and countries.
        """

        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap

        if "figsize" in kwargs:
            figsize = kwargs["figsize"]
            if isinstance(figsize, tuple):
                fig = plt.figure(figsize=figsize)
            else:
                print("Warning: the figsize must be a tuple as (length, width), "
                      "use figure size default (automatically depending on the figure plotted).")
                fig = plt.figure()
        else:
            fig = plt.figure()
        ax = fig.add_subplot(111)
        if "nomargin" in kwargs:
            nomargin = kwargs["nomargin"]
            if nomargin:
                plt.subplots_adjust(left=0, bottom=0,
                                             right=1, top=1, wspace=0, hspace=0)
        if "mapset" not in kwargs:
            mapset = {'llcrnrlon': 73, 'llcrnrlat': 17, 'urcrnrlon': 136, 'urcrnrlat': 54,
                      'epsg': 3857}
        m = Basemap(**mapset)
        longtide, latitude = np.meshgrid(self.lons, self.lats)
        xi, yi = m(longtide, latitude)
        if "drawmap" in kwargs:
            drawmap = kwargs["drawmap"]
            if drawmap:
                m.drawmapboundary()
                m.drawcoastlines()
                m.drawcountries()
        else:
            m.drawmapboundary()
            m.drawcoastlines()
            m.drawcountries()
        my_cmap = None
        if "custom_cmap_set" in kwargs:
            my_cmap = self.make_cmap(**kwargs["custom_cmap_set"])

        if "contourf_set" in kwargs:
            if isinstance(kwargs["contourf_set"], dict):
                if my_cmap != None:
                    kwargs["contourf_set"]['cmap'] = my_cmap
                cs = m.contourf(xi, yi, self.data, **kwargs["contourf_set"])
                method = "contourf"
            else:
                print("Warning: the contourf_set must be a dictionary, "
                      "using default plotting method.")
                cs = m.contourf(xi, yi, self.data, cmap='jet', extend="both")
                method = "contourf"
        elif "pcolormesh_set" in kwargs:
            if isinstance(kwargs["pcolormesh_set"], dict):
                if my_cmap != None:
                    kwargs["pcolormesh_set"]['cmap'] = my_cmap
                cs = m.contourf(xi, yi, self.data, **kwargs["pcolormesh_set"])
                method = "pcolormesh"
            else:
                print("Warning: the pcolormesh_set must be a dictionary, "
                      "using default plotting method.")
                cs = m.contourf(xi, yi, self.data, cmap='jet', extend="both")
                method = "contourf"
        else:
            cs = m.contourf(xi, yi, self.data, cmap='jet', extend="both")
            method = "contourf"
        if "clip_set" in kwargs:
            self.shp2clip(cs, ax, m, method, **kwargs["clip_set"])
        if "nobackground" in kwargs:
            nobackground = kwargs["nobackground"]
            if nobackground:
                for item in [fig, ax]:
                    item.patch.set_visible(False)
        if "nobounding" in kwargs:
            nobounding = kwargs["nobounding"]
            if nobounding:
                ax.set_frame_on(False)
        if "nobar" in kwargs:
            nobar = kwargs["nobar"]
            if nobar:
                pass
            else:
                m.colorbar(cs)
        else:
            m.colorbar(cs)
        if "savefig_set" in kwargs:
            if isinstance(kwargs["savefig_set"], dict):
                if "other" in kwargs["savefig_set"]:
                    plt.savefig(kwargs["savefig_set"]["fname"], **kwargs["savefig_set"]["other"])
                else:
                    plt.savefig(kwargs["savefig_set"]["fname"])
                plt.close(fig)
            else:
                print("Warning: the savefig_set must be a dictionary, using default saving setting")
                figure = "./%s_filed.png" % (self.variable)
                plt.savefig(figure)
                plt.close(fig)
        else:
            plt.show()
            plt.close(fig)

    def shp2clip(self, originfig, ax, m, method, shpfile=None, region=None, region_code=None):
        """
        This module enables you to maskout the unneccessary data outside the interest region
        on a matplotlib-plotted output instance.

        Args:
            'originfig': the matplotlib instance
            'ax':  the Axes instance
            'shapefile':  the shape file used for generating a basemap A
            'region': the name of a region of on the basemap A,outside the region the data is
                      to be maskout
            'region_code': code in shape file which region name located
            'm': the Basemap instance

        Returns:
            'clip': clipped matplotlib instance

        Tips:
             "if region_code and region depends on shapfiles (import shapefile to check)
                 *2004.shp ==> shape_rec.record[-1]==b'\xb1\xb1\xbe\xa9\xca\xd0'，(gbk not unicode)
                 eg: b'\xb1\xb1\xbe\xa9\xca\xd0' (beijing)
        """

        import shapefile
        from matplotlib.path import Path
        from matplotlib.patches import PathPatch

        sf = shapefile.Reader(shpfile)
        for shape_rec in sf.shapeRecords():
            if shape_rec.record[region_code] == region:
                vertices = []
                codes = []
                pts = shape_rec.shape.points
                prt = list(shape_rec.shape.parts) + [len(pts)]
                for i in range(len(prt) - 1):
                    for j in range(prt[i], prt[i + 1]):
                        vertices.append(m(pts[j][0], pts[j][1]))
                    codes += [Path.MOVETO]
                    codes += [Path.LINETO] * (prt[i + 1] - prt[i] - 2)
                    codes += [Path.CLOSEPOLY]
                clip = Path(vertices, codes)
                clip = PathPatch(clip, transform=ax.transData)
        if method == 'contourf':
            for contour in originfig.collections:
                contour.set_clip_path(clip)
        if method == 'pcolormesh':
            originfig.set_clip_path(clip)

    def make_cmap(self, colors=None, position=None, bit=True, cmap_color_num=256):
        """
        Takes a list of tuples which contain RGB values making a cmap with equally spaced colors
        using matplotlib.colors.LinearSegmentedColormap.

        More explanation of cdict:
            Example: suppose you want red to increase from 0 to 1 over the bottom
            half, green to do the same over the middle half, and blue over the top
            half.  Then you would use:

            cdict = {'red':   ((0.0,  0.0, 0.0),
                               (0.5,  1.0, 1.0),
                               (1.0,  1.0, 1.0)),

                     'green': ((0.0,  0.0, 0.0),
                               (0.25, 0.0, 0.0),
                               (0.75, 1.0, 1.0),
                               (1.0,  1.0, 1.0)),

                     'blue':  ((0.0,  0.0, 0.0),
                               (0.5,  0.0, 0.0),
                               (1.0,  1.0, 1.0))}

            If, as in this example, there are no discontinuities in the r, g, and b
            components, then it is quite simple: the second and third element of
            each tuple, above, is the same--call it "y".  The first element ("x")
            defines interpolation intervals over the full range of 0 to 1, and it
            must span that whole range.  In other words, the values of x divide the
            0-to-1 range into a set of segments, and y gives the end-point color
            values for each segment.

            Now consider the green. cdict['green'] is saying that for
            0 <= x <= 0.25, y is zero; no green.
            0.25 < x <= 0.75, y varies linearly from 0 to 1.
            x > 0.75, y remains at 1, full green.

            If there are discontinuities, then it is a little more complicated.
            Label the 3 elements in each row in the cdict entry for a given color as
            (x, y0, y1).  Then for values of x between x[i] and x[i+1] the color
            value is interpolated between y1[i] and y0[i+1].

            Going back to the cookbook example, look at cdict['red']; because y0 !=
            y1, it is saying that for x from 0 to 0.5, red increases from 0 to 1,
            but then it jumps down, so that for x from 0.5 to 1, red increases from
            0.7 to 1.  Green ramps from 0 to 1 as x goes from 0 to 0.5, then jumps
            back to 0, and ramps back to 1 as x goes from 0.5 to 1.

            row i:   x  y0  y1
                            /
                           /
            row i+1: x  y0  y1

            Above is an attempt to show that for x in the range x[i] to x[i+1], the
            interpolation is between y1[i] and y0[i+1].  So, y0[0] and y1[-1] are
            never used.

        Args:
            colors, bit: The RGB values may either be in 8-bit [0 to 255]
                         or arithmetic [0 to 1] ((in which bit must be set to False).
            position: Position contains values from 0 to 1 to dictate the location of each color.
            cmap_color_num: bin numbers between two colors

        Returns:
            cmap: a cmap with 256 equally spaced colors.
            One can change cmap_color_num to modify the color number
            since if one have 1024 levels then need 1024 colors to plot contourf.
        """

        import sys
        import numpy as np
        import matplotlib as mpl

        bit_rgb = np.linspace(0, 1, 256)
        if position is None:
            position = np.linspace(0, 1, len(colors))
        else:
            if len(position) != len(colors):
                sys.exit("position length must be the same as colors")
            elif position[0] != 0 or position[-1] != 1:
                sys.exit("position must start with 0 and end with 1")

        if bit:
            colors_arithmetic = []
            for color in colors:
                colors_arithmetic.append((bit_rgb[color[0]], bit_rgb[color[1]], bit_rgb[color[2]]))
            colors = colors_arithmetic
        cdict = {'red': [], 'green': [], 'blue': []}
        for pos, color in zip(position, colors):
            cdict['red'].append((pos, color[0], color[0]))
            cdict['green'].append((pos, color[1], color[1]))
            cdict['blue'].append((pos, color[2], color[2]))
        cmap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, cmap_color_num)
        return cmap
