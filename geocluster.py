# coding: utf-8

from timeUtils import clock, elapsed
from collections import OrderedDict
from haversine import haversine
from pandas import DataFrame, Series
from pandasUtils import isSeries, isDataFrame
import geohash
    
    

##############################################################################################################################
# Geo Clusters Class
##############################################################################################################################
class geoClusters():
    def __init__(self, key, points = None, geoCnts = None, distMax = 150, debug=False):
        self.device = key
        self.points = None
        self.bitlen = None
        
        self.geoDataFrame  = None
        self.geoCntsSeries = None
        
        self.clusters    = None
        self.clusterCoMs = None
        
        self.clusterPrefix = "cl"
        self.distMax       = distMax
        
        self.summary = {}
        
        self.convertPoints(points, debug)    
        self.findGeos(debug)
        self.findGeoCounts(debug)
            
        
    #########################################################################################################
    # Getter/Setter Functions
    #########################################################################################################
    def setGeoDataFrame(self, geoDataFrame):
        self.geoDataFrame = geoDataFrame
    
    def getGeoDataFrame(self):
        return self.geoDataFrame
    
    def setGeoCntsSeries(self, geoCntsSeries):
        self.geoCntsSeries = geoCntsSeries
    
    def getGeoCntsSeries(self):
        return self.geoCntsSeries

    def setMaxDistance(self, dist):
        self.distMax = dist

    def getMaxDistance(self):
        return self.distMax

    def getClusters(self):
        return self.clusters

    def getClusterByIndex(self, idx, debug=False):
        name    = "{0}{1}".format(self.clusterPrefix, idx)
        cluster = self.clusters.get(name)
        if cluster is None:
            if debug:
                print("Could not find cluster {0} in list of [{1}] available clusters".format(name, list(self.clusters.keys())))
        return cluster

    def getClusterNameByIndex(self, idx, debug=False):
        name    = "{0}{1}".format(self.clusterPrefix, idx)
        cluster = self.clusters.get(name)
        if cluster is None:
            name = None
            if debug:
                print("Could not find cluster {0} in list of [{1}] available clusters".format(name, list(self.clusters.keys())))
        return name

    def getCluster(self, name, debug=False):
        cluster = self.clusters.get(name)
        if cluster is None:
            if debug:
                print("Could not find cluster {0} in list of [{1}] available clusters".format(name, list(self.clusters.keys())))
        return cluster
    
    def getNClusters(self):
        return self.summary.get('Clusters')
    
    def getNCells(self):
        return self.summary.get('Cells')
    
    def getNCounts(self):
        return self.summary.get('Counts')

    def getClusterRadius(self, name, debug=False):
        cluster = self.clusters.get(name)
        if cluster is None:
            if debug:
                print("Could not find cluster {0} in list of [{1}] available clusters".format(name, list(self.clusters.keys())))
            return None
        return cluster.getRadius()

    def getClusterQuantiles(self, name, debug=False):
        cluster = self.clusters.get(name)
        if cluster is None:
            if debug:
                print("Could not find cluster {0} in list of [{1}] available clusters".format(name, list(self.clusters.keys())))
            return None
        return cluster.getQuantiles()

    def getClusterCoM(self, name, debug=False):
        cluster = self.clusters.get(name)
        if cluster is None:
            if debug:
                print("Could not find cluster {0} in list of [{1}] available clusters".format(name, list(self.clusters.keys())))
            return None
        return cluster.getCoM()

    def getClusterCells(self, name, debug=False):
        cluster = self.clusters.get(name)
        if cluster is None:
            if debug:
                print("Could not find cluster {0} in list of [{1}] available clusters".format(name, list(self.clusters.keys())))
            return None
        return cluster.getNCells()

    def getClusterCellNames(self, name, debug=False):
        cluster = self.clusters.get(name)
        if cluster is None:
            if debug:
                print("Could not find cluster {0} in list of [{1}] available clusters".format(name, list(self.clusters.keys())))
            return None
        return cluster.getCellNames()
            
        
    #########################################################################################################
    # Collect Cluster Center of Masses
    #########################################################################################################
    def getClusterCoMs(self, debug=False):
        if debug:
            start, cmt = clock("Collecting Cluster CoMs")
            
        if self.clusterCoMs is None:
            coms = DataFrame([x.getCoM() for name,x in self.getClusters().items()])
            self.clusterCoMs = coms
        else:
            coms = self.clusterCoMs
            
        if debug:
            elapsed(start, cmt)
        
        return coms
        
        
    #########################################################################################################
    # Convert Points to Correct Format
    #########################################################################################################
    def convertPoints(self, points, debug=False):
        if points is None:
            raise ValueError("Points is None and cannot convert!")
            
        if debug:
            start, cmt = clock("Converting {0} Points To Correct Format".format(len(points)))
            
        from numpy import ndarray, stack
        from pandas.core.series import Series
        from pandas.core.frame import DataFrame
        if isinstance(points, ndarray):
            self.points = points
        elif isinstance(points, DataFrame):
            self.points = points.values
        elif isinstance(points, (list,set,tuple)):
            if len(points) == 2:
                x = points[0]
                y = points[1]
            else:
                raise ValueError("Not sure how to parse data points of type {0} for device {1}".format(type(points), device))

            if isinstance(x, Series):
                x = x.values
            if isinstance(y, Series):
                y = y.values
            if not all([isinstance(x, ndarray), isinstance(y, ndarray)]):
                raise ValueError("Data is not a numpy array!")

            self.points = stack((x,y), axis=-1)
        self.convertToMeters()
        
        if debug:
            print("Data has correct format with a {0} shape.".format(self.points.shape))
        
        if debug:
            elapsed(start, cmt)
    
    
    #########################################################################################################
    # Print Functions
    #########################################################################################################
    def showGeoValues(self):
        print(self.geoDataFrame)
        
    def showClusters(self):
        print("Cl\tSeed\t\tTotal\tRadius\tLatitude,Longitude\tNum Subclusters")
        print("--\t----\t\t-----\t------\t------------------\t---------------")
        for clusterName, cluster in self.clusters.items():
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(clusterName, cluster.getSeed(), cluster.getCounts(),
                                                        cluster.getRadius(), cluster.getCoM(), cluster.getNCells()))
                
    
    #########################################################################################################
    # Geohash and Geohash Counts Finders
    #########################################################################################################
    def findGeos(self, debug=False):      
        def f(x):
            prec=8
            try:
                lat = x[0]
                long = x[1]
                retval = geohash.encode(lat, long, precision=prec)
            except:
                retval = 'xxxxxxxx'
            return retval 
        
        if debug:
            start, cmt = clock("Finding Geohash (BitLen=8) Values from {0} Points".format(len(self.points)))

        from pandas.core.frame import DataFrame
        from numpy import vectorize
        geoDataFrame = DataFrame(self.points, columns=['lat', 'long'])
        geoDataFrame['geo'] = geoDataFrame[['lat', 'long']].apply(f, axis=1).values
        self.setGeoDataFrame(geoDataFrame)
        
        if debug:
            elapsed(start, cmt)
        

    def findGeoCounts(self, debug=False):
        if not isDataFrame(self.geoDataFrame):
            raise ValueError("Must call findGeos with points before trying to count geohashs!")
                
        if debug:
            start, cmt = clock("Finding Geohash (BitLen=8) Frequency Values from Geohash DataFrame")
            
        from collections import Counter
        from pandas.core.series import Series
        
        s = self.getGeoDataFrame()['geo']
        geoCntsSeries = Series(Counter(dict(s.value_counts())))
        geoCntsSeries.sort_values(ascending=False, inplace=True)
        self.setGeoCntsSeries(geoCntsSeries)
        
        if debug:
            elapsed(start, cmt)
            
        
    def convertToMeters(self):
        from numpy import stack
        try:
            longMeters = 111319.9*(self.points[:,1] - min(self.points[:,1]))
            latMeters  = 111319.9*(self.points[:,0] - min(self.points[:,0]))
            self.meters = stack((latMeters, longMeters), axis=-1)

            maxDist = max([max(longMeters), max(latMeters)])
            if maxDist > 5000:
                self.bitlen=5
            elif maxDist > 1000:
                self.bitlen=6
            elif maxDist > 100:
                self.bitlen=7
            else:
                self.bitlen=8
        except:
            self.bitlen = 5
            self.meters = None
            
            
    def getDist(self, gcode1, gcode2, units='m'):
        from haversine import haversine
        if all((isinstance(x, str) for x in [gcode1, gcode2])):
            try:
                pnt1 = geohash.decode_exactly(gcode1)[:2]
                pnt2 = geohash.decode_exactly(gcode2)[:2]
                dist = haversine(pnt1, pnt2)
            except:
                dist = None
        elif all((isinstance(x, tuple) for x in [gcode1, gcode2])):
            try:
                dist = haversine(gcode1, gcode2)
            except:
                dist = None
        else:
            raise ValueError("Did not understand types {0} and {1} for getDist() inputs in geoClustering.py".format(type(gcode1), type(gcode2)))
            
        if units == 'm':
            try:
                dist = 1000*dist
            except:
                dist = None
        return dist
        
    
    
    #########################################################################################################
    #
    # Find Geo Clusters
    #
    #########################################################################################################
    def findClusters(self, seedMin=2, addMin=2, debug=False):
        from collections import OrderedDict
        if debug:
            start, cmt = clock("Finding Clusters with at least {0} counts".format(seedMin))
            
        clusters    = OrderedDict()     
        totalCounts = 0
        totalGeos   = 0
        verydebug   = False
        
        geoCounts = self.getGeoCntsSeries()
        if not isSeries(geoCounts):
            raise ValueError("Cannot cluster because geoCounts is not a Seriers!")
        if geoCounts is None:
            raise ValueError("Cannot cluster because there are no geoCounts in findCluster!")
        

        ## Loop over geo counts (geo, geoCnt)
        if verydebug:
            print("There are {0} remaining cells".format(geoCounts.count()))
        clCount = -1
        while len(clusters) - clCount > 0 and geoCounts.count() > 0:
            clCount = len(clusters)
            
            ## Take top cell as seed
            idx     = geoCounts.index
            seed    = idx[0]
            seedCnt = geoCounts[seed]
            
            ## Check for None
            if seed is None:
                continue

            ## Apply cluster cuts
            if seedCnt < seedMin:
                break

            ## Set cluster seed
            cluster      = geoCluster(seed=seed, cnts=seedCnt, clnum=len(clusters), clusterPrefix=self.clusterPrefix, debug=debug)
            totalGeos   += 1
            totalCounts += seedCnt


            ## Loop over geos
            for geo, geoCnt in geoCounts.iteritems():
                if geo == seed:
                    continue
                dist  = round(self.getDist(seed,geo),1)
                #print("  Check: {0}  ,  Dist: {1} < {2}".format(geoChk,dist,self.distMax))
                if dist <= self.distMax:
                    cluster.addCell(geo, geoCnt)
                    totalGeos   += 1
                    totalCounts += geoCnt

            cells     = list(cluster.getCells().keys())
            geoCounts = geoCounts.drop(labels=cells)
            cluster.findCoM()
            clusters[cluster.getName()] = cluster
            
            if verydebug:
                print("Created cluster with seed {0} and {1} cells".format(geo, len(cluster.getCells())))
                print("There are {0} remaining cells".format(geoCounts.count()))

        self.clusters = clusters
        
        self.summary['Clusters']  = len(clusters)
        self.summary['Cells']     = totalGeos
        self.summary['Counts']    = totalCounts
                
        if debug:
            elapsed(start, cmt)
            

            
    #########################################################################################################
    #
    # Cluster Matching
    #
    #########################################################################################################
    def getNearestCluster(self, lat, long, debug=False):
        if debug:
            start, cmt = clock("Computing Nearest Cluster for ({0}, {1})".format(lat, long))
            
        testSubClusters = False
        from scipy import spatial
        
        ## Get Cluster Center of Masses
        Acl  = self.getClusterCoMs()
        try:
            Acl = Acl.values
        except:
            raise ValueError("Could not extract values from cluster CoMs DataFrame!")
        
        if isinstance(lat,(int,float)) and isinstance(long,(int,float)):
            pt   = [lat,long]
            try:
                distance,index = spatial.KDTree(Acl).query(pt)
                distance *= 1000 # convert to meters
            except:
                if debug:
                    print("There was an error with the cluster finding in KDTree for {0}".format(pt))
                return None, -1, -1, (lat, long)
                
            if distance < self.distMax:
                if debug:
                    print("  Found nearest cluster {0} with distance {1}".format(index, round(distance,1)))
                clusterName = self.getClusterNameByIndex(index, debug)
                if clusterName is None:
                    raise ValueError("Returned cluster is NULL for index {0}!!".format(index))
                if debug:
                    elapsed(start, cmt)
                return clusterName, index, round(1000*distance,1), Acl[index]
            else:
                if debug:
                    print("  Nearest cluster {0} is too far away: distance {1} > {2}".format(index, round(distance,1), self.distMax))
                return None, -1, -1, (lat, long)
        else:
            if debug is True:
                print("Latitude {0} and Longitude {1} are not set for device {2}!".format(type(lat), type(long), self.device))
            return None, -1, -1, (lat, long)
        
        
    def getNearestClusters(self, gpsData, debug=False):
        if debug:
            start, cmt = clock("Computing Nearest Clusters for {0}".format(gpsData.values))
        
        try:
            latlong = gpsData.values
            cluster, index, distance, position = self.getNearestCluster(latlong[0], latlong[1])
        except:
            raise ValueError("Could not get cluster for {0} for device {1}".format(latlong, self.device))
            
        if debug:
            elapsed(start, cmt)
            
        return [cluster, index, distance, position]
            
        
        
##############################################################################################################################
# Geo Cluster Class
##############################################################################################################################
class geoCluster():
    def __init__(self, seed, cnts, clnum, clusterPrefix, debug=False):
        self.seed     = seed
        self.seedCnts = cnts
        self.clnum    = clnum
        self.clusterPrefix = clusterPrefix
        self.name     = "{0}{1}".format(clusterPrefix, clnum)
        
        ## Initialize list of cells
        self.cells    = OrderedDict()

        ## Add the seed cell to the list of seeds
        self.addCell(seed, cnts)
        
        self.clDataFrame = None
        self.quantiles   = None
        self.com         = None
        self.radius      = None
        self.counts      = None
        self.geos        = None
        self.quantiles   = None
        
        if debug:
            print("  --> Creating cluster {0} with seed {1} and {2} counts".format(self.name, self.seed, self.seedCnts))
        
    def getName(self):
        return self.name
        
    def getCounts(self):
        return self.counts
    
    def setDataFrame(self, df):
        self.clDataFrame = df
        
    def getDataFrame(self):
        return self.clDataFrame
    
    def setCoM(self, com):
        self.com = com
    
    def getCoM(self):
        return self.com
    
    def setSeed(self, seed):
        self.seed = seed
    
    def getSeed(self):
        return self.seed
    
    def setRadius(self, radius):
        self.radius = radius
    
    def getRadius(self):
        return self.radius
    
    def setQuantiles(self, quantiles):
        self.quantiles = quantiles
    
    def getQuantiles(self):
        return self.quantiles
            
    def addCell(self, geo, cnts, debug=False):
        if self.cells.get(geo) is not None:
            raise ValueError("Trying to add geo {0} to cluster {1}, but it's already there".format(geo, self.name))
        
        if debug:
            print("\tAdding cell {0}/{1} with counts {2} to cluster {3}".format(geo, len(self.cells), cnts, self.name))
        
        self.cells[geo] = cnts
        
    def getCells(self):
        return self.cells
    
    def getNCells(self):
        return len(self.cells)
    
    def getCellNames(self):
        return list(self.cells.keys())
    
    
    #########################################################################################################
    # Show Functions
    #########################################################################################################
    def show(self):
        print("Information for Cluster {0}".format(self.name))
        print("  Cells:     {0}".format(self.geos))
        print("  Counts:    {0}".format(self.counts))
        print("  CoM:       {0}".format(self.com))
        print("  Radius:    {0}".format(self.radius))
        print("  Quantiles: {0}".format(self.quantiles))
        print(self.clDataFrame)
        print("")
        
        
        
    #########################################################################################################
    # Helper
    #########################################################################################################            
    def getDist(self, gcode1, gcode2, units='m'):
        from haversine import haversine
        if all((isinstance(x, str) for x in [gcode1, gcode2])):
            try:
                pnt1 = geohash.decode_exactly(gcode1)[:2]
                pnt2 = geohash.decode_exactly(gcode2)[:2]
                dist = haversine(pnt1, pnt2)
            except:
                dist = None
        elif all((isinstance(x, tuple) for x in [gcode1, gcode2])):
            try:
                dist = haversine(gcode1, gcode2)
            except:
                dist = None
        else:
            raise ValueError("Did not understand types {0} and {1} for getDist() inputs in geoClustering.py".format(type(gcode1), type(gcode2)))
            
        if units == 'm':
            try:
                dist = 1000*dist
            except:
                dist = None
        return dist
        
        
    
    #########################################################################################################
    # Compute Features
    #########################################################################################################
    def findCoM(self, debug=False):
        if debug:
            print("\tComputing Center of Mass and Radius for {0} cells".format(len(self.cells)))

        lats  = []
        lngs  = []
        wgts  = []
        dists = []
        geos  = []
        self.counts = 0
        self.geos   = 0
        for geo,cnts in self.cells.items():
            lat, long = geohash.decode_exactly(geo)[:2]
            wgts.append(cnts)
            lats.append(lat)
            lngs.append(long)
            geos.append(geo)
            self.counts += cnts
            self.geos   += 1
        
        ## Compute Center of Mass
        swgts = sum(wgts)
        latC = round(sum(wgt*lats[i] for i,wgt in enumerate(wgts))/swgts,5)
        lngC = round(sum(wgt*lngs[i] for i,wgt in enumerate(wgts))/swgts,5)
        com = (latC,lngC)
    
        ## Computer Distances from CoM
        dists = []
        for geo,cnts in self.cells.items():
            geoPnt = geohash.decode_exactly(geo)[:2]
            dist   = self.getDist(com, geoPnt, units='m')
            dists.append(dist)
            
        clDataFrame = DataFrame(list(zip(geos, lats, lngs, dists)), columns=['geo', 'lat', 'long', 'distance'])
        dists       = Series(dists)
        radius      = round(dists.mean(),0)
        quantiles   = [round(x) for x in dists.quantile(q=[0, 0.25, 0.5, 0.75, 1])]
        
        self.setRadius(radius)
        self.setCoM(com)
        self.setQuantiles(quantiles)
        self.setDataFrame(clDataFrame)
        
        if debug:
            print("\tCenter of Mass = {0}  ;  Radius = {1}".format(self.getCoM(), self.getRadius()))