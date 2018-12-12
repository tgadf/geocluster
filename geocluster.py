# coding: utf-8

from timeUtils import clock, elapsed
from collections import OrderedDict, Counter
from haversine import haversine
from pandas import DataFrame, Series
from pandasUtils import isSeries, isDataFrame
from geoUtils import getDist, applyGeo8
import geohash
from convertData import convertData
from numpy import ndarray, stack, vectorize

    
    

##############################################################################################################################
# Geo Clusters Class
##############################################################################################################################
class geoClusters():
    def __init__(self, key, points = None, geoCnts = None, distMax = 150, mergeFactor = 2.0, debug=False):
        self.device = key
        self.points = None
        self.bitlen = None
        
        self.cellsDataFrame  = None
        self.cells = None
        
        self.mergeFactor = mergeFactor
        
        self.protoCells = None
        self.perimCells = None
        
        self.protoClusters    = None
        self.seedlessClusters = None
        self.mergedClusters   = None
        self.protoclusterCoMs = None
        
        self.clusters    = None
        self.clusterCoMs = None
        
        self.clusterPrefix = "cl"
        self.distMax       = distMax
        
        self.summary = {}

        if points is not None:
            self.setInputs(points=points, debug=debug)
        else:
            if debug:
                print("Creating empty geo cluster object!")

                
    def setInputs(self, points, debug=False):        
        ## Convert Data
        cd = convertData()
        cd.setData(points)
        self.points = cd.getArray()
        #self.createCells(debug=debug)
            
        
    #########################################################################################################
    # Getter/Setter Functions
    #########################################################################################################
    def setMaxDistance(self, dist):
        self.distMax = dist

    def getMaxDistance(self):
        return self.distMax

    ######################################
    # Cells
    ######################################
    def setCellsDataFrame(self, cellsDataFrame, debug=False):
        if debug:
            print("  Setting Cells DataFrame (Raw Data) with {0} Inputs".format(cellsDataFrame.shape[0]))
        self.cellsDataFrame = cellsDataFrame

    def getCellsDataFrame(self):
        return self.cellsDataFrame
    
    def setCells(self, cells, debug=False):
        if debug:
            print("  Setting {0} Cells Data".format(cells.shape[0]))
        self.cells = cells

    def getCells(self):
        return self.cells
    
    def setPerimCells(self, cells, debug=False):
        if debug:
            print("  Setting {0} PerimCells Data".format(cells.shape[0]))
        self.perimCells = cells

    def getPerimCells(self):
        return self.perimCells
    
    def setSeedlessCells(self, cells, debug=False):
        if debug:
            print("  Setting {0} SeedlessCells Data".format(cells.shape[0]))
        self.seedlessCells = cells

    def getSeedlessCells(self):
        return self.seedlessCells
    
    def setProtoCells(self, cells, debug=False):
        if debug:
            print("  Setting {0} ProtoCells Data".format(cells.shape[0]))
        self.protoCells = cells

    def getProtoCells(self):
        return self.protoCells

    
    ######################################
    # Proto Clusters
    ######################################
    def setProtoClusters(self, protoClusters, debug=False):
        if debug:
            print("  Setting {0} ProtoCluster Data".format(len(protoClusters)))
        self.protoClusters = protoClusters

    def getProtoClusters(self):
        return self.protoClusters

    
    ######################################
    # Seedless Clusters
    ######################################
    def setSeedlessClusters(self, seedlessClusters, debug=False):
        if debug:
            print("  Setting {0} SeedlessCluster Data".format(len(seedlessClusters)))
        self.seedlessClusters = seedlessClusters

    def getSeedlessClusters(self):
        return self.seedlessClusters

    
    ######################################
    # Merged Clusters
    ######################################
    def setMergedClusters(self, mergedClusters, debug=False):
        if debug:
            print("  Setting {0} MergedCluster Data".format(len(mergedClusters)))
        self.mergedClusters = mergedClusters

    def getMergedClusters(self):
        return self.mergedClusters
    
    
    ######################################
    # Clusters
    ######################################
    def setClusters(self, clusters, debug=False):
        if debug:
            print("Setting Clusters Data")
        self.clusters = clusters

    def getClusters(self):
        return self.clusters

    def getNClusters(self):
        return len(self.clusters)

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
                
    
    #########################################################################################################################################################
    #
    # Create Cluster Cells
    #
    ######################################################################################################################################################### 
    def createCells(self, debug=False):      
        if debug:
            start, cmt = clock("Finding Geohash (BitLen=8) Values from {0} Points".format(len(self.points)))
        
        cellDataFrame = DataFrame(self.points, columns=['lat', 'long'])
        cellDataFrame['geo'] = applyGeo8(cellDataFrame)
        self.setCellsDataFrame(cellDataFrame, debug=debug)
        cells = Series(Counter(dict(cellDataFrame['geo'].value_counts())))
        cells.sort_values(ascending=False, inplace=True)
        self.setCells(cells, debug=debug)
        
        if debug:
            elapsed(start, cmt)
        
    
    
    #########################################################################################################################################################
    #
    # Create Proto Clusters
    #
    #########################################################################################################################################################
    def createProtoClusters(self, seedMin=10, debug=False):
        if debug:
            start, cmt = clock("Creating ProtoClusters with at least {0} counts".format(seedMin))
        
        cells = self.getCells().copy()
        protoCells = {}
        verydebug=False
        
        protos  = {}
        nProtos = -1
        while len(protos) != nProtos:

            ## Take top cell as seed            
            seed    = cells.index[0]
            cnt     = cells[seed]
            if cnt < seedMin:
                break
            
            protoname = "p{0}".format(len(protos))
            protos[protoname] = {"Seed": [seed,cnt], "Counts": cnt, "CoM": None, "Radius": 0, "Cells": []}
            protoCells[seed] = protoname
            for geo in geohash.neighbors(seed):
                try:
                    geocnts = cells[geo]
                except:
                    geocnts = None
                
                if geocnts is not None:
                    dist = getDist(seed, geo, units='m')
                    protos[protoname]["Radius"] = max(dist, protos[protoname]["Radius"])
                    protos[protoname]["Cells"].append([geo, geocnts])
                    protos[protoname]["Counts"] += geocnts
                    protoCells[geo] = protoname
                    
                for gh2 in geohash.neighbors(geo):
                    if protoCells.get(gh2) is not None:
                        continue
                    try:
                        geocnts2 = cells[gh2]
                    except:
                        geocnts2 = None

                    if geocnts2 is not None:
                        dist = getDist(seed, gh2, units='m')
                        protos[protoname]["Radius"] = max(dist, protos[protoname]["Radius"])
                        protos[protoname]["Cells"].append([gh2, geocnts2])
                        protos[protoname]["Counts"] += geocnts2
                        protoCells[gh2] = protoname
                    
                
                
            pcells = protos[protoname]["Cells"]+[protos[protoname]["Seed"]]
            coms = [geohash.decode_exactly(x[0])[:2] for x in pcells]
            com  = tuple([Series(x[0] for x in coms).mean(), Series(x[1] for x in coms).mean()])
            rad   = max([getDist(com,x) for x in coms])
            protos[protoname]["CoM"]    = com
            protos[protoname]["Radius"] = rad
            
            if verydebug:
                print("  Found ProtoCluster at {0}\tRadius {1}\t{2} Counts\t{3} Cells".format([round(x,4) for x in com], 
                                                                                                      round(protos[protoname]["Radius"],1),
                                                                                                      protos[protoname]["Counts"], 
                                                                                                      len(protos[protoname]["Cells"])))
            
            dropCells = [protos[protoname]["Seed"][0]] + [x[0] for x in protos[protoname]["Cells"]]
            cells = cells.drop(labels=dropCells)

                        
        self.setProtoCells(protoCells)
        self.setPerimCells(cells)
        self.setProtoClusters(protos)
            
        if debug:
            print("Found {0} proto clusters".format(len(protos)))
            print("  There are {0} remaining perimeter cells".format(self.perimCells.shape[0]))
        
    
    
    #########################################################################################################################################################
    #
    # Create Proto Clusters
    #
    #########################################################################################################################################################
    def createSeedlessClusters(self, seedMin=2, debug=False):
        if debug:
            start, cmt = clock("Creating SeedlessClusters with at least {0} counts".format(seedMin))
        
        cells = self.getPerimCells()
        seedlessCells = {}
        verydebug=False
        
        seedless  = {}
        nSeedless = -1
        
        while len(seedless) != nSeedless:

            ## Take top cell as seed
            geo = cells.index[0]
            cnt = cells[geo]
            if cnt < seedMin:
                break
        
            name = "s{0}".format(len(seedless))
            seedless[name] = {"Counts": cnt, "CoM": None, "Radius": 0, "Cells": []}
            seedless[name]["Cells"].append([geo,cnt])
            seedlessCells[geo] = name
            for gh in geohash.neighbors(geo):
                if seedlessCells.get(gh) is not None:
                    continue
                try:
                    geocnts = cells[gh]
                except:
                    geocnts = None
                
                if geocnts is not None:
                    dist = getDist(gh, geo, units='m')
                    seedless[name]["Radius"] = max(dist, seedless[name]["Radius"])
                    seedless[name]["Cells"].append([gh, geocnts])
                    seedless[name]["Counts"] += geocnts
                    seedlessCells[gh] = name
                    
                for gh2 in geohash.neighbors(gh):
                    if seedlessCells.get(gh2) is not None:
                        continue
                    try:
                        geocnts2 = cells[gh2]
                    except:
                        geocnts2 = None

                    if geocnts2 is not None:
                        dist = getDist(gh2, geo, units='m')
                        seedless[name]["Radius"] = max(dist, seedless[name]["Radius"])
                        seedless[name]["Cells"].append([gh2, geocnts2])
                        seedless[name]["Counts"] += geocnts2
                        seedlessCells[gh2] = name

            scells = seedless[name]["Cells"]
            coms = [geohash.decode_exactly(x[0])[:2] for x in scells]
            com  = tuple([Series(x[0] for x in coms).mean(), Series(x[1] for x in coms).mean()])
            rad   = max([getDist(com,x) for x in coms])
            seedless[name]["CoM"] = com
            seedless[name]["Radius"] = rad
                
            if verydebug:
                print("  Found SeedlessCluster at {0}\tRadius {1}\t{2} Counts\t{3} Cells".format([round(x,4) for x in com], 
                                                                                                      round(seedless[name]["Radius"],1),
                                                                                                      seedless[name]["Counts"], 
                                                                                                      len(seedless[name]["Cells"])))
            
            dropCells = [x[0] for x in seedless[name]["Cells"]]
            cells = cells.drop(labels=dropCells)
            
        self.setSeedlessCells(seedlessCells)
        self.setPerimCells(cells)
        self.setSeedlessClusters(seedless)
        
        if debug:
            print("Found {0} seedless clusters".format(len(seedless)))
            print("  There are {0} remaining perimeter cells".format(self.perimCells.shape[0]))
    
    
    #########################################################################################################################################################
    #
    # Split/Merge Clusters
    #
    #########################################################################################################################################################
    def mergeClusters(self, debug=False):
        if debug:
            start, cmt = clock("Merge ProtoClusters and SeedlessClusters")
            
        from copy import deepcopy
        protoClusters = deepcopy(self.getProtoClusters())
        seedlessClusters = deepcopy(self.getSeedlessClusters())
        mergedClusters = {}
        dropClusters = {}
        verydebug=False

        for pcl in list(protoClusters.keys()):
            if dropClusters.get(pcl) is not None:
                continue
            protoCluster = protoClusters[pcl]
            pCoM = tuple(protoCluster["CoM"])
            pRad = protoCluster["Radius"]
            dropClusters[pcl] = 1

            protoMerges = []
            for pcl2,protoCluster2 in protoClusters.items():
                if dropClusters.get(pcl2) is not None:
                    continue
                pCoM2 = tuple(protoCluster2["CoM"])
                pRad2 = protoCluster2["Radius"]
                dist  = getDist(pCoM,pCoM2)
                if dist < self.mergeFactor*max([pRad,pRad2]):
                    #print(dist,'\t',pCoM,pRad,'\t',pCoM2,pRad2)
                    protoMerges.append(pcl2)
                    dropClusters[pcl2] = 1
                    if verydebug:
                        print("  Merging {0} into {1}".format(pcl2,pcl))

            seedlessMerges = []
            for scl,seedlessCluster in seedlessClusters.items():
                if dropClusters.get(scl) is not None:
                    continue
                sCoM = tuple(seedlessCluster["CoM"])
                sRad = seedlessCluster["Radius"]
                dist  = getDist(pCoM,sCoM)
                if dist < self.mergeFactor*max([pRad,sRad]):
                    #print(dist,'\t',pCoM,pRad,'\t',sCoM,sRad)
                    seedlessMerges.append(scl)
                    dropClusters[scl] = 1
                    if verydebug:
                        print("  Merging {0} into {1}".format(scl,pcl))

            ## Merge everything
            cells = protoCluster["Cells"] + [protoCluster["Seed"]]
            for pcl2 in protoMerges:
                cells += protoClusters[pcl2]["Cells"] + [protoClusters[pcl2]["Seed"]]
            for scl in seedlessMerges:
                cells += seedlessClusters[scl]["Cells"]
                
            coms  = [geohash.decode_exactly(x[0])[:2] for x in cells]
            com   = tuple([Series(x[0] for x in coms).mean(), Series(x[1] for x in coms).mean()])
            cnt   = sum([x[1] for x in cells])
            rad   = max([getDist(com,x) for x in coms])
            mname = "cl{0}".format(len(mergedClusters))
            mergedClusters[mname] = {"Counts": cnt, "CoM": com, "Radius": rad, "Cells": cells}

            if debug:
                print("  Found MergedClusters {0} at {1}\tRadius {2}\t{3} Counts\t{4} Cells".format(mname,[round(x,4) for x in com], 
                                                                                                      round(mergedClusters[mname]["Radius"],1),
                                                                                                      mergedClusters[mname]["Counts"], 
                                                                                                      len(mergedClusters[mname]["Cells"])))
                
        for name in dropClusters.keys():
            if name.startswith("s"):
                del seedlessClusters[name]
            if name.startswith("p"):
                del protoClusters[name]
                        

        dropClusters = {}
        for scl in list(seedlessClusters.keys()):
            if dropClusters.get(scl) is not None:
                continue
            seedlessCluster = seedlessClusters[scl]
            sCoM = tuple(seedlessCluster["CoM"])
            sRad = seedlessCluster["Radius"]
            dropClusters[scl] = 1

            seedlessMerges = []
            for scl2,seedlessCluster2 in seedlessClusters.items():
                if dropClusters.get(scl2) is not None:
                    continue
                sCoM2 = tuple(seedlessCluster2["CoM"])
                sRad2 = seedlessCluster2["Radius"]
                dist  = getDist(sCoM,sCoM2)
                if dist < self.mergeFactor*max([sRad,sRad2]):
                    #print(dist,'\t',sCoM,sRad,'\t',sCoM2,sRad2)
                    seedlessMerges.append(scl2)
                    dropClusters[scl2] = 1
                    if verydebug:
                        print("  Merging {0} into {1}".format(scl2,scl))

            ## Merge everything
            cells = seedlessCluster["Cells"]
            for scl2 in seedlessMerges:
                cells += seedlessClusters[scl2]["Cells"]
                
            coms  = [geohash.decode_exactly(x[0])[:2] for x in cells]
            com   = tuple([Series(x[0] for x in coms).mean(), Series(x[1] for x in coms).mean()])
            cnt   = sum([x[1] for x in cells])
            rad   = max([getDist(com,x) for x in coms])
            mname = "cl{0}".format(len(mergedClusters))
            mergedClusters[mname] = {"Counts": cnt, "CoM": com, "Radius": rad, "Cells": cells}

            if debug:
                print("  Found MergedClusters {0} at {1}\tRadius {2}\t{3} Counts\t{4} Cells".format(mname,[round(x,4) for x in com], 
                                                                                                      round(mergedClusters[mname]["Radius"],1),
                                                                                                      mergedClusters[mname]["Counts"], 
                                                                                                      len(mergedClusters[mname]["Cells"])))
        for name in dropClusters.keys():
            if name.startswith("s"):
                del seedlessClusters[name]
            if name.startswith("p"):
                del protoClusters[name]
            
            
        self.setMergedClusters(mergedClusters)
        
        
        ## Set the final clusters        
        clusters = {}
        for cl,mergedCluster in mergedClusters.items():
            mcl   = geoCluster(clnum=len(clusters), clusterPrefix="cl")
            cells = mergedCluster["Cells"]
            cells = dict(zip([t[0] for t in cells], [t[1] for t in cells]))
            mcl.setCells(cells)
            mcl.analyze()
            name = mcl.getName()
            clusters[name] = mcl
            
            
        self.setClusters(clusters)
  
        if debug:
            print("Found {0} merged clusters".format(len(mergedClusters)))
            print("  There are {0} remaining proto clusters".format(len(protoClusters)))
            print("  There are {0} remaining seedless clusters".format(len(seedlessClusters)))
            print("Created {0} final clusters".format(len(clusters)))
                
        if debug:
            elapsed(start, cmt)
    
    
    #########################################################################################################################################################
    #
    # Find Geo Clusters
    #
    #########################################################################################################################################################
    def findClusters(self, seedMin=2, addMin=2, debug=False):
        raise ValueError("Don't call this anymore!")
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
            print("There are {0} remaining cells".format(geoCounts.count(axis=0)))
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

        self.setClusters(clusters)
        
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
    def __init__(self, clnum, clusterPrefix, debug=False):
        self.clnum    = clnum
        self.clusterPrefix = clusterPrefix
        self.name     = "{0}{1}".format(clusterPrefix, clnum)
        
        ## Initialize list of cells
        self.cells    = OrderedDict()
        
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
        
    def setCells(self, cells):
        self.cells = cells
    
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
    # Compute Features
    #########################################################################################################
    def analyze(self, debug=False):
            
        if self.cells is None:
            raise ValueError("Cells are None in geoCluster()!")

        if len(self.cells) == 0:
            raise ValueError("Cells are Empty in geoCluster()!")
            
        if debug:
            print("\tComputing Center of Mass, Radius, and Quantiles for {0} cells".format(len(self.cells)))

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
            dist   = getDist(com, geoPnt, units='m')
            dists.append(dist)
            
        clDataFrame = DataFrame(list(zip(geos, lats, lngs, dists)), columns=['geo', 'lat', 'long', 'distance'])
        dists       = Series(dists)
        quantiles   = [round(x) for x in dists.quantile(q=[0, 0.25, 0.5, 0.75, 1])]
        radius      = round(dists.max(),0)
        
        self.setRadius(radius)
        self.setCoM(com)
        self.setQuantiles(quantiles)
        self.setDataFrame(clDataFrame)
        
        if debug:
            print("\tCenter of Mass = {0}  ;  Radius = {1}".format(self.getCoM(), self.getRadius()))