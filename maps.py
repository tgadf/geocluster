from folium import PolyLine, CircleMarker, Circle, Marker, Icon, FeatureGroup, Map, LayerControl
import geohash


class foliumMap():
    def __init__(self, df=None, pc=None, gc=None):
        self.pc = None
        if pc is not None:
            self.setProtoClusters(pc)
        self.gc = None
        if gc is not None:
            self.setGeoClusters(gc)
        self.df = None
        if df is not None:
            self.setTripsDataFrame(df)
        
        self.m  = None
        self.colorPoints = False
        
        self.init_zoom = 10
        
        self.colors = ['red', 'blue', 'gray', 'darkred', 'lightred', 'orange', 'beige', 'green', 'darkgreen', 'lightgreen',
                       'darkblue', 'lightblue', 'purple', 'darkpurple', 'pink', 'cadetblue', 'lightgray']
        
        self.colors = ['red', 'blue', 'green', 'orange', 'purple', 'darkred', 'darkgreen', 'darkblue'] #, 'darkpurple']
        
        
    ########################################################################################
    # Setters
    ########################################################################################
    def setTripsDataFrame(self, df, showStart=True, showEnd=True):
        self.df = df.copy()
        if "geo0" not in df.columns:
            self.df["geo0"] = 0
        if "geo1" not in df.columns:
            self.df["geo1"] = 0
        
        if showStart and showEnd:
            self.df         = df[["lat0", "long0", "geo0"]]
            self.df.columns = ["lat", "long", "label"]
            pnts            = df[["lat1", "long1", "geo1"]]
            pnts.columns    = ["lat", "long", "label"]
            self.df         = self.df.append(pnts)
        elif showStart:
            self.df         = df[["lat0", "long0", "geo0"]]
            self.df.columns = ["lat", "long", "label"]
        elif showEnd:
            self.df         = df[["lat1", "long1", "geo1"]]
            self.df.columns = ["lat", "long", "label"]
        
        self.ranking = list(self.df['label'].value_counts().index)
            
        
    def setProtoClusters(self, pc):
        self.pc = pc
        
    def setGeoClusters(self, gc):
        self.gc = gc
        
        
        
    ########################################################################################
    # Getters
    ########################################################################################        
    def getMap(self):
        return self.m
    
    def saveMap(self, savename="folium.html", dirname="/Users/tgadf/Downloads"):
        from os.path import join
        filename = join(dirname, savename)
        print("Saving folium map to {0}".format(filename), end="... ")
        self.m.save(filename)
        print("Done.")
        
        
        
    ########################################################################################
    ########################################################################################
    # Create Map
    ########################################################################################
    ########################################################################################
    def createMapFromTripsDataFrame(self, location=None, zoom=None, debug=False):
        if self.df is None:
            print("There is no trips DataFrame object!")
            return
        
        if location is not None:
            try:
                lat0 = location[0]
                long0 = location[1]
            except:
                raise ValueError("Location of map is not parseable")
        else:
            try:
                lat0           = self.df['lat'].mean()
                long0          = self.df['long'].mean()
            except:
                raise ValueError("Could not get center of geo clusters and create map!")
            
        if zoom is None:
            zoom = self.init_zoom
        self.m = Map(location=[lat0, long0], zoom_start=zoom)
        
        
    def createMapFromProtoClusters(self, location=None, zoom=None, debug=False):
        if self.pc is None:
            print("There is no ProtoClusters object!")
            return
        
        if location is not None:
            try:
                lat0 = location[0]
                long0 = location[1]
            except:
                raise ValueError("Location of map is not parseable")
        else:
            try:
                lats  = Series(x["CoM"][0] for pcl,x in self.pc.items())
                lngs  = Series(x["CoM"][1] for pcl,x in self.pc.items())
                lat0  = lats.mean()
                long0 = lngs.mean()
            except:
                raise ValueError("Could not get center of geo clusters and create map!")
            
        if zoom is None:
            zoom = self.init_zoom
        self.m = Map(location=[lat0, long0], zoom_start=zoom)
        
        
    def createMapFromGeoClusters(self, location=None, zoom=None, debug=False):
        if self.gc is None:
            print("There is no GeoClusters object!")
            return
        
        if location is not None:
            try:
                lat0 = location[0]
                long0 = location[1]
            except:
                raise ValueError("Location of map is not parseable")
        else:
            try:
                maxW = 0
                com  = None
                for cl,cluster in self.gc.getClusters().items():
                    if cluster.getCounts() > maxW:
                        com  = cluster.getCoM()
                        maxW = cluster.getCounts()
                lat0  = com[0]
                long0 = com[1]
            except:
                raise ValueError("Could not get center of geo clusters and create map!")
            
        if zoom is None:
            zoom = self.init_zoom
        self.m = Map(location=[lat0, long0], zoom_start=zoom)
        
        
    def createMap(self, location=None, zoom=None, debug=False):
        if self.gc is not None:
            self.createMapFromGeoClusters(location=location, zoom=zoom, debug=debug)
        elif self.pc is not None:
            self.createMapFromProtoClusters(location=location, zoom=zoom, debug=debug)
        elif self.df is not None:
            self.createMapFromTripsDataFrame(location=location, zoom=zoom, debug=debug)
        else:
            raise ValueError("Cannot create map because there is no data object!")
            
        
        
    ########################################################################################
    ########################################################################################
    # Points/Clusters
    ########################################################################################
    ########################################################################################      
    def addPointsFromTripsDataFrame(self, debug=False):
        if self.m is None:
            print("Folium Map is None!")
            return
        
        if self.df is None:
            print("DataFrame is None!")
            return
        
        cols = ['darkblue', 'lightblue', 'pink', 'lightgray']
        
        ## Number of clusters (if known)
        from hashlib import sha256
        from numpy import log1p
        ncls  = len(set(self.df["label"].unique()))
        ncols = min([ncls, len(self.colors)])
        if ncols == 1:
            cols = ['black']
            freq = {"0": 1}            
        else:
            cols  = self.colors[:ncols]
            freq  = dict(self.df["label"].value_counts())
            self.colorPoints = True            
            

        ranking = self.gc.getClusterRanking(returnSeries=False)
        
        rad = 5
        weight = 1
        for ir,row in self.df.iterrows():            
            cl   = str(row["label"])
            try:
                rank = self.ranking.index(cl)
            except:
                rank = 0
            com  = row[["lat", "long"]].values
            hval = int(sha256(cl.encode('utf-8')).hexdigest(), 16) % ncols
            col  = cols[rank % ncols]
            try:
                rad  = log1p(freq[cl])
            except:
                rad  = 1
            Circle(com, color=col, radius=rad, fill=True, fill_color=col, weight=weight, opacity=0).add_to(self.m)
    
    def addPointsFromGeoClusterCells(self, debug=False):
        if self.m is None:
            print("Folium Map is None!")
            return
        
        if self.gc is None:
            print("GeoClusters is None!")
            return

        cells = self.gc.getCells()
        for geo,cnt in cells.iteritems():
            com    = geohash.decode_exactly(geo)[:2]
            wgt    = int(cnt)
            popup  = str(wgt)
            Circle(com, color='black', radius=5, fill=True, fill_color='black', weight=wgt, opacity=0, popup=popup).add_to(self.m)
    
    
    def addPointsFromProtoClusters(self, debug=False):
        if self.m is None:
            print("Folium Map is None!")
            return
        
        if self.gc is None:
            print("GeoClusters is None!")
            return
        
        cols = ['darkblue', 'lightblue', 'pink', 'lightgray']
        
        from pandas import Series

        for pcl,protoCluster in self.gc.getProtoClusters().items():
            com    = protoCluster["CoM"]
            rad    = max([protoCluster["Radius"], 5])
            counts = protoCluster["Counts"]
            weight = int(counts)
            name   = pcl
            popup  = "{0} : N = {1}".format(name, counts)

            Circle(com, color=cols[0], radius=rad, fill=True, fill_color=cols[0], weight=weight, opacity=0).add_to(self.m)
    
    
    def addPointsFromSeedlessClusters(self, debug=False):
        if self.m is None:
            print("Folium Map is None!")
            return
        
        if self.gc is None:
            print("GeoClusters is None!")
            return
        
        cols = ['darkgreen']
        
        from pandas import Series

        for scl,seedlessCluster in self.gc.getSeedlessClusters().items():
            com    = seedlessCluster["CoM"]
            rad    = max([seedlessCluster["Radius"], 5])
            counts = seedlessCluster["Counts"]
            weight = int(counts)
            name   = scl
            popup  = "{0} : N = {1}".format(name, counts)

            Circle(com, color=cols[0], radius=rad, fill=True, fill_color=cols[0], weight=weight, opacity=0).add_to(self.m)
    
    
    def addPointsFromMergedClusters(self, debug=False):
        if self.m is None:
            print("Folium Map is None!")
            return
        
        if self.gc is None:
            print("GeoClusters is None!")
            return
        
        cols = ['darkred']
        
        from pandas import Series

        for cl,cluster in self.gc.getMergedClusters().items():
            com    = cluster["CoM"]
            rad    = max([cluster["Radius"], 5])
            counts = cluster["Counts"]
            weight = int(counts)
            name   = cl
            popup  = "{0} : N = {1}".format(name, counts)

            Circle(com, color=cols[0], radius=rad, fill=True, fill_color=cols[0], weight=weight, opacity=0).add_to(self.m)
    
    
    def addPointsFromGeoClusters(self, debug=False):
        if self.m is None:
            print("Folium Map is None!")
            return
        
        if self.gc is None:
            print("GeoClusters is None!")
            return
        
        cols = ['darkblue', 'lightblue', 'pink', 'lightgray']
        
        ## Number of clusters (if known)
        from hashlib import sha256
        from numpy import log1p
        ncls  = len(self.gc.getClusters())
        ncols = min([ncls, len(self.colors)])
        print(ncols,ncls)
        
        from pandas import Series
        feature_group_1 = FeatureGroup(name="Top 3 (#1=Red)")
        feature_group_2 = FeatureGroup(name="Top 10")
        feature_group_3 = FeatureGroup(name="Top 25")
        feature_group_4 = FeatureGroup(name="All Else")

        weights = Series([cluster.getCounts() for cl,cluster in self.gc.getClusters().items()])
        ranking = self.gc.getClusterRanking(returnSeries=False)
        
        top3    = min(weights.nlargest(3))
        top10   = min(weights.nlargest(10))
        top25   = min(weights.nlargest(25))

        for rank, cl in enumerate(ranking):
            cluster = self.gc.getCluster(name=cl)
            com    = cluster.getCoM()
            rad    = max([cluster.getRadius(), 10])
            counts = cluster.getCounts()
            weight = float(counts)
            quant  = cluster.getQuantiles()
            name   = cl
            cells  = ",".join(cluster.getCells())
            col    = self.colors[rank % ncols]
            
            #hval = int(sha256(cl.encode('utf-8')).hexdigest(), 16) % ncols
            #col  = self.colors[hval]
                
            if counts >= top3:
                fg = feature_group_1
            elif counts >= top10:
                fg = feature_group_2
            elif counts >= top25:
                fg = feature_group_3
            else:
                fg = feature_group_4            
            
            popup  = "{0} : N = {1} : R = {2} : {3} : Rank = {4}".format(name, weight, rad, com, rank+1)
            #popup  = ""
            


            Marker(com, icon=Icon(color=col, icon_color="white", icon="car", angle=0, prefix='fa'), popup=popup).add_to(fg)
            Circle(com, color=col, radius=rad, fill=False, fill_color='white', weight=log1p(weight), opacity=1).add_to(fg)

        feature_group_1.add_to(self.m)
        feature_group_2.add_to(self.m)
        feature_group_3.add_to(self.m)
        feature_group_4.add_to(self.m)
        LayerControl().add_to(self.m)    
        
           
    def addPoints(self, debug=False):
        if self.m is None:
            raise ValueError("Cannot add points to an empty map")
        if self.pc is not None:
            self.addPointsFromProtoClusters(debug=debug)
        elif self.gc is not None:
            self.addPointsFromGeoClusters(debug=debug)
        elif self.df is not None:
            self.addPointsFromTripsDataFrame(debug=debug)
        else:
            print("Cannot add points because there is map and object!")