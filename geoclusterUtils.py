from numpy import pi, cos, sin, ceil, vstack, array
from numpy.random import shuffle, uniform, normal, choice
from pandas import DataFrame
from haversine import haversine

    


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

def genTripsBetweenClusters(n, gc, returnLoc=True, returnDF=False):
    clusters = gc.getClusters()
    cls   = list(clusters.keys())
    trips = list(zip(choice(array(cls), size=n), choice(array(cls), size=n)))
    
    genTrips = []
    for trip in trips:
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
        
        genTrips.append([start, end])
        
    if returnDF is True:
        genTrips = convertTripsToDataFrame(array(genTrips))
        
    return genTrips


def convertTripsToDataFrame(trips):
    getDist = lambda x: haversine((x[0], x[1]), (x[2], x[3]))
    df = DataFrame(trips.reshape(1000, 4))
    df.columns = ["lat0", "long0", "lat1", "long1"]
    df['total_miles'] = list(map(getDist, df[["lat0", "long0", "lat1", "long1"]].values))
    df['duration'] = 60.0*df['total_miles']/3
    return df