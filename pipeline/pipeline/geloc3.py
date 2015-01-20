from mpl_toolkits.basemap import Basemap, addcyclic
import numpy as np
import matplotlib.pyplot as plt
import geoloc2
import shapely.ops as ops

fig=plt.figure()
ax=fig.add_axes([0.1,0.1,0.8,0.8])
# setup me≡jedi=0, rcator map projection.≡ (llcrnrlon = None, *llcrnrlat = None*, urcrnrlon = None, urcrnrlat = None, llcrnrx = None, llcrnry = None, urcrnrx = None, urcrnry = None, width = None, height = None, projection = 'cyl', resolution = 'c', area_thresh = None, rsphere = 6370997.0, ellps = None, lat_ts = None, lat_1 = None, lat_2 = None, lat_0 = None, lon_0 = None, lon_1 = None, lon_2 = None, o_lon_p = None, o_lat_p = None, k_0 = None, no_rot = False, suppress_ticks = True, satellite_height = 35786000, boundinglat = None, fix_aspect = True, anchor = 'C', celestial = False, round = False, epsg = None, ax = None) ≡jedi≡
m = Basemap(llcrnrlon=-50.,llcrnrlat=-50.,urcrnrlon=340.,urcrnrlat=65.,\
                        rsphere=(6378137.00,6356752.3142),\
                        resolution='h',projection='merc',\
                        lat_0=40.,lon_0=-20.,lat_ts=20.)
# nylat, nylon are lat/lon of New York
nylat = 40.78; nylon = -73.98
# lonlat, lonlon are lat/lon of London.
lonlat = 51.53; lonlon = 0.08
# draw great circle route between NY and London
m.drawgreatcircle(nylon,nylat,lonlon,lonlat,linewidth=2,color='b')
m.drawcoastlines()
m.fillcontinents()
# draw parallels
m.drawparallels(np.arange(10,90,20),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
ax.set_title('Great Circle from New York to London')

c = geoloc2.load_countries(geoloc2.s)
scotland = c.subset(countries=['US'], regions=["Asia"])
for s in scotland:

    
    s.patch = ops.transform(m, s.patch)
scotland.plot(ax, False, fc='red')

plt.show()

