from numpy import pi, cos, sin, ceil, vstack, array, repeat
from numpy.random import shuffle, uniform, normal, choice
from pandas import DataFrame
from haversine import haversine
from geoUtils import convertMetersToLat, convertMetersToLong
from pandasUtils import castDateTime, castFloat64
from collections import OrderedDict
from random import choices
from pandas import DataFrame


def genCenters(n, latRange, lngRange, retList=False, debug=False):
    if debug:
        print("Generating {0} Latitude centers between ({1}, {2})".format(n, latRange[0], latRange[1]))
    latvals = uniform(low=latRange[0], high=latRange[1], size=n)
    
    if debug:
        print("Generating {0} Longitude centers between ({1}, {2})".format(n, lngRange[0], lngRange[1]))
    lngvals = uniform(low=lngRange[0], high=lngRange[1], size=n)
        
    if retList is True:
        retval = list(zip(latvals, lngvals))
    else:
        arrs = [array(latvals), array(lngvals)]
        retval = array(arrs).transpose()
    
    return retval


def genCluster(npoints, center, dist="gauss", maxrad=100, retList=False, debug=False):
    ## Generate Latitude
    if dist == "gauss":
        mulat, siglat = center[0], convertMetersToLat(maxrad) # mean and standard deviation
        if debug:
            print("Generating Latitude {0} samples of Gaussian({1}, {2})".format(npoints, mulat, round(siglat,5)))
        lat = normal(mulat, siglat, npoints)
    if dist == "uniform":
        mulat, siglat = center[0], convertMetersToLat(maxrad) # mean and standard deviation
        if debug:
            print("Generating Latitude {0} samples of Uniform({1}, {2})".format(npoints, mulat, round(siglat,5)))
        lat = uniform(mulat, -siglat/2, siglat/2, npoints)
    
    
    ## Generate Longitude
    if dist == "gauss":
        mulng, siglng = center[1], convertMetersToLong(maxrad, center[0]) # mean and standard deviation
        if debug:
            print("Generating Longitude {0} samples of Gaussian({1}, {2})".format(npoints, mulng, round(siglng,5)))
        lng = normal(mulng, siglng, npoints)
    if dist == "uniform":
        mulng, siglng = center[1], convertMetersToLong(maxrad, center[0]) # mean and standard deviation
        if debug:
            print("Generating Longitude {0} samples of Gaussian({1}, {2})".format(npoints, mulng, round(siglng,5)))
        lng = uniform(mulng, -siglng/2, siglat/2, npoints)
    
    arrs = [array(lat), array(lng)]
    arr2d = array(arrs)
    arr2d = arr2d.transpose()
    
    return arr2d

def genClusters(n, ppc, latRange, lngRange, dist="gauss", maxrad=100, mix=True, debug=False):
    if debug:
        print("Generating {0} Centers".format(n))
    centers = genCenters(n, latRange, lngRange, retList=True, debug=debug)
    
    clusters = []
    for center in centers:
        cluster = genCluster(npoints=ppc, center=center, dist=dist, maxrad=maxrad, retList=False, debug=debug)
        clusters.append(cluster)
    
    retval = vstack(clusters)
    if mix is True:
        shuffle(retval)
        
    return retval


def genImportantClusters(gc):
    clusters = gc.getClusters()
    cls   = list(clusters.keys())

    impcls  = list(set(choice(array(cls), 10)))[:6]
    homecl  = impcls[0]
    impcls  = impcls[1:]
    cls.remove(homecl)
    for cl in impcls:
        cls.remove(cl)
    return {"Home": homecl, "Imps": impcls, "Rest": cls}
    
    
def getStartEndLocation(trip, clusters, returnLoc=True):
    # Start
    cl  = trip[0]
    if returnLoc is True:
        cluster = clusters[cl]
        rad = uniform(0, convertMetersToLat(cluster.getQuantiles()[-1]))
        phi = uniform(0, 2*pi)
        com = cluster.getCoM()
        lat = com[0] + rad*cos(phi)
        lng = com[1] + rad*sin(phi)
        start = (lat,lng)
    else:
        start = cl

    # End
    cl  = trip[1]
    if returnLoc is True:
        cluster = clusters[cl]
        rad = uniform(0, convertMetersToLat(cluster.getQuantiles()[-1]))
        phi = uniform(0, 2*pi)
        com = cluster.getCoM()
        lat = com[0] + rad*cos(phi)
        lng = com[1] + rad*sin(phi)
        end = (lat,lng)
    else:
        end = cl    
        
    return [start, end, trip]
        
    

def genTripsBetweenClusters(n, gc, returnLoc=True, returnDF=False):
    clusters = gc.getClusters()
    cls   = list(clusters.keys())

    retval = genImportantClusters(gc)
    trips = []

    cls = OrderedDict()
    cls[retval['Home']] = 0.05
    pval = 0.45/len(retval['Imps'])
    for cl in retval['Imps']:
        cls[cl] = pval

    pval = (1.0 - sum(cls.values()))/len(retval['Rest'])
    for cl in retval['Rest']:
        cls[cl] = pval

    p = list(cls.values())
    a = list(cls.keys())

    dst = choice(a=a, p=p, size=1000)
    src = repeat(retval['Home'], 1000)

    trips += list(zip(src, dst))


    for cl in retval['Imps']:
        cls = OrderedDict()
        cls[retval['Home']] = 0.20
        cls[cl] = 0.01
        for cl2 in retval['Imps']:
            if cl == cl2:
                continue
            cls[cl2] = 0.1

        pval = (1.0 - sum(cls.values()))/len(retval['Rest'])
        for cl in retval['Rest']:
            cls[cl] = pval

        p = list(cls.values())
        a = list(cls.keys())

        dst = choice(a=a, p=p, size=1000)
        src = repeat(cl, 1000)

        trips += list(zip(src, dst))


    for cl in retval['Rest']:
        cls = OrderedDict()
        cls[retval['Home']] = 0.05
        cls[cl] = 0.01
        for cl2 in retval['Imps']:
            cls[cl2] = 0.1
        pval = (1.0 - sum(cls.values()))/(len(retval['Rest'])-1)
        for cl2 in retval['Rest']:
            if cl == cl2:
                continue
            cls[cl2] = pval

        p = list(cls.values())
        a = list(cls.keys())

        dst = choice(a=a, p=p, size=1000)
        src = repeat(cl, 1000)

        trips += list(zip(src, dst))    
    
    
    ## Randomize trips
    shuffle(trips)
    
    
    ## Select 'n' trips
    selectedTrips = choices(trips, k=n)
    print("Selected {0} randomized trips".format(len(selectedTrips)))
    
    genTrips = []
    for trip in selectedTrips:
        genTrips.append(getStartEndLocation(trip, clusters, returnLoc=returnLoc))
    print("Found Start/End for the {0} randomized trips".format(len(genTrips)))
        
    print(genTrips[0])
    if returnDF is True:
        genTrips = convertTripsToDataFrame(array(genTrips))
        
    return genTrips


def convertTripsToDataFrame(trips):
    print("Converting {0} trips to a DataFrame".format(trips.shape))
    getDist = lambda x: haversine((x[0], x[1]), (x[2], x[3]))
    df = DataFrame(trips.reshape(trips.shape[0], 6))    
    df.columns = ["lat0", "long0", "lat1", "long1", "cl0", "cl1"]
    df["lat0"] = castFloat64(df["lat0"])
    df["lat1"] = castFloat64(df["lat1"])
    df["long0"] = castFloat64(df["long0"])
    df["long1"] = castFloat64(df["long1"])
    df['total_miles'] = list(map(getDist, df[["lat0", "long0", "lat1", "long1"]].values))
    df['duration'] = 60.0*df['total_miles']/3
    
    import datetime as dt
    start = dt.datetime.now()
        
    df['start'] = start
    getEnd = lambda x: x[0] + dt.timedelta(0, float(x[1]))
    df['end'] = list(map(getEnd, df[['start', 'duration']].values))

    df['start'] = castDateTime(df['start'])
    df['end'] = castDateTime(df['end'])    
    return df