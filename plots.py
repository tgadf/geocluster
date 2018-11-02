
# coding: utf-8

from timeUtils import clock, elapsed
from matplotlib import pyplot as plt
try:
    import pygeohash as geohash
except:
    import geohash

class drawClusters():
    def __init__(self, soloclusters = None, geoCounts = None, bitlen = None, clusters = None):
        self.clusters  = clusters
        self.geoCounts = geoCounts
        self.bitlen    = bitlen
        self.soloclusters = soloclusters
        self.routes=[('n', 'w'), ('n', 'e'), ('s', 'e'), ('s', 'w'), ('n', 'w')]
        self.rangeLat  = None
        self.rangeLng  = None
        
        self.showCoM   = True
        self.showSoloClusters = False
        
    def setGeoCounts(self, geoCounts):
        self.geoCounts = geoCounts
        
    def setGeoClusters(self, clusters):
        self.clusters = clusters
        
    def setSoloGeoClusters(self, soloclusters):
        self.soloclusters = soloclusters
        self.showSoloClusters = True
        
    def setBitlen(self, bitlen):
        self.bitlen = bitlen
        
    def setCenterRadius(self, latCenter, lngCenter, radius):
        self.rangeLat = (latCenter-radius,latCenter+radius)
        self.rangeLng = (lngCenter-radius,lngCenter+radius)
        
    def setLatRange(self, minval,maxval):
        self.rangeLat = (minval,maxval)
    
    def setLongRange(self, minval, maxval):
        self.rangeLng = (minval, maxval)
        
    def getCenter(self, gh):
        from pygeohash import decode_exactly
        lat,long = decode_exactly(gh)[:2]
        return (lat, long)
        
    def getBB(self, gh, istuple=False):
        from pygeohash import bbox
        if istuple is False:
            bb = {'x': [], 'y': []}
        else:
            bb = []
        boundbox = bbox(gh)
        for route in self.routes:
            if istuple is False:
                bb['x'].append(boundbox[route[1]])
                bb['y'].append(boundbox[route[0]])
            else:
                bb.append(tuple([boundbox[route[0]], boundbox[route[1]]]))
        return bb

    def drawBB(self, ax1=None):
        import matplotlib.pyplot as plt
        if ax1 is None:
            fig, ax1 = plt.subplots()
            ax1.set_facecolor('white')
            
        if self.bitlen is None:
            raise ValueError("bitlen is not set for drawing")
            
        if isinstance(self.bitlen, list):
            colors=['red', 'black', 'blue', 'green', 'magenta']
            for i,blen in enumerate(self.bitlen):
                try:
                    for geo,geoCnt in self.geoCounts[blen].most_common():
                        bb = self.getBB(geo)
                        ax1.plot(bb['x'], bb['y'], linestyle=':', color=colors[i], lw=0.5)  # solid        
                except:
                    pass
        else:
            try:
                for geo,geoCnt in self.geoCounts[self.bitlen].most_common():
                    bb = self.getBB(geo)
                    ax1.plot(bb['x'], bb['y'], linestyle=':', color='black', lw=0.5)  # solid        
            except:
                pass
            
        ax1.set_xlabel("Longitude")
        ax1.set_ylabel("Latitude")

        return ax1


    def plotClusters(self, ax1=None, long=None, lat=None):
        from seaborn import color_palette
        from haversine import haversine
        if ax1 is None:
            fig, ax1 = plt.subplots()
            ax1.set_facecolor('white')
        
        if self.bitlen is not None and self.geoCounts is not None and self.showCoM is True:
            self.drawBB(ax1)
    
        cp            = color_palette('deep', len(self.clusters))
        lineStyles    = ['-' for i,_ in enumerate(self.clusters)]
        lineWidths    = [2 for i,_ in enumerate(self.clusters)]

        latRange = [None, None]
        lngRange = [None, None]

        from matplotlib.patches import Circle, Wedge, Polygon
        from matplotlib.collections import PatchCollection
        patches=[]
        pcols=[]


        if self.showCoM is True and self.clusters is not None:
            lats = []
            lngs = []
            wgts = []
            radi = []
            for cluster,clusterData in self.clusters.items():
                lats.append(clusterData['CoM'][0])
                lngs.append(clusterData['CoM'][1])
                wgts.append(clusterData['Total'])
                radi.append(clusterData['Radius']/111319.9)
                pcols.append(tuple(cp[len(patches)]))
                circle = Circle((lngs[-1], lats[-1]), radi[-1])
                patches.append(circle)
            wgts = [max([x,1.5]) for x in wgts]
            wgts = [min([25,x]) for x in wgts]
            latRange[0] = min(lats)
            latRange[1] = max(lats)
            lngRange[0] = min(lngs)
            lngRange[1] = max(lngs)
            
            
            p = PatchCollection(patches, alpha=0.5)
            from numpy import array, linspace
            p.set_array(linspace(0,1,len(pcols)))
            ax1.add_collection(p)
            ax1.scatter(lngs, lats, s=wgts, marker='x', linewidth=2, c='black', alpha=0.5)
            

            if self.showSoloClusters is True and self.soloclusters is not None:
                color = 'black'
                lineStyle = ':'
                lineWidth = 1
                lats = []
                lngs = []
                wgts = []
                pcols= []
                patches=[]
                for cluster,clusterData in self.soloclusters.items():
                    lats.append(clusterData['CoM'][0])
                    lngs.append(clusterData['CoM'][1])
                    wgts.append(clusterData['Total'])
                    pcols.append(1)
                    circle = Circle((lngs[-1], lats[-1]), 20/111319.9)
                    patches.append(circle)
                    
                latRange[0] = min([min(lats),latRange[0]])
                latRange[1] = max([max(lats),latRange[1]])
                lngRange[0] = min([min(lngs),lngRange[0]])
                lngRange[1] = max([max(lngs),lngRange[1]])
                
                p = PatchCollection(patches, alpha=0.15)
                from numpy import array, linspace
                p.set_array(linspace(0,1,len(pcols)))
                ax1.add_collection(p)
                ax1.scatter(lngs, lats, s=0.25, marker='.', linewidth=lineWidth, c=color, alpha=0.5)
                
        elif showCoM is False and self.clusters is not None:
            color = 'red'
            lineStyle = lineStyles[0]
            lineWidth = lineWidths[0]
            for cluster,clusterData in self.clusters.items():
                bb = self.getBB(geo)
                ax1.plot(bb['x'], bb['y'], linestyle=lineStyle, color=color, lw=lineWidth)
                if latRange[0] is None:
                    latRange[0] = bb['y']
                else:
                    latRange[0] = min([latRange[0], bb['y']])
                if latRange[1] is None:
                    latRange[1] = bb['y']
                else:
                    latRange[1] = max([latRange[1], bb['y']])
                if lngtRange[0] is None:
                    lngRange[0] = bb['x']
                else:
                    lngRange[0] = min([lngRange[0], bb['x']])
                if lngRange[1] is None:
                    lngRange[1] = bb['x']
                else:
                    lngRange[1] = max([lngRange[1], bb['x']])
                
        xmin = None
        xmax = None
        ymin = None
        ymax = None
        if self.rangeLng is not None:
            xmin = max([lngRange[0], self.rangeLng[0]])
            xmax = min([lngRange[1], self.rangeLng[1]])
            xdif = xmax - xmin
            xmin = xmin - 0.05*xdif
            xmax = xmax + 0.05*xdif
            ax1.set_xlim(xmin,xmax)
        if self.rangeLat is not None:
            ymin = max([latRange[0], self.rangeLat[0]])
            ymax = min([latRange[1], self.rangeLat[1]])
            ydif = ymax - ymin
            ymin = ymin - 0.05*ydif
            ymax = ymax + 0.05*ydif
            ax1.set_ylim(ymin,ymax)
            
        if all([xmin,xmax,ymin,ymax]):
            distLat = 1000*haversine((xmin,ymin), (xmin,ymax))
            distLng = 1000*haversine((xmin,ymin), (xmax,ymin))
            
            longMeters = [0, distLng]
            latMeters  = [0, distLat]
            
            ax2 = ax1.twinx()
            ax2.plot(longMeters, latMeters, color='b', lw=0)
            ax2.set_xlabel("Relative Longitude [m]")
            ax3 = ax1.twiny()
            ax3.plot(longMeters, latMeters, color='b', lw=0)
            ax3.set_ylabel("Relative Latitude [m]")

#        if long is not None and lat is not None:
#            ax1.scatter(long, lat, s=1, linewidth=0, c='black', alpha=0.25)

from timeUtils import clock
start = clock("Compiled Plotting Tool")