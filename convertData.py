from pandas import DataFrame, Series
from pandasUtils import isSeries, isDataFrame
from numpyUtils import isArray
from numpy import ndarray, asarray

##############################################################################################################################
# Geo Clusters Class
##############################################################################################################################
class convertData():
    def __init__(self):
        self.dtype     = None
        self.datadf    = None
        self.datalist  = None
        self.dataarray = None
        self.dim       = None
        self.colnames  = None
        
        
    def setColumns(self, colnames):
        self.dim = len(colnames)
        self.colnames = colnames
        
        
    def setData(self, data):
        if data is None:
            raise ValueError("Data is None in convert Data!")
            

        if isDataFrame(data):
            self.dtype  = DataFrame
            self.datadf = data
        elif isArray(data):
            self.dtype = ndarray
            self.dataarray = data
        elif isinstance(data, list):
            self.dtype = list
            if len(data) == 2:
                self.datalist = list(zip(data[0], data[1]))
            else:
                self.datalist = data
        else:
            raise ValueError("Did not recognize data format {0}".format(type(data)))
            
            
    def dType(self):
        print(self.dtype)
        
    def getData(self):
        print("Data is type {0}".format(self.dtype))
        if self.datadf is not None:
            return self.datadf
        elif self.datalist is not None:
            return self.datalist
        elif self.dataarray is not None:
            return self.dataarray
        else:
            print("Data is not set!")
        
        
    
    #########################################################################################################
    # Return Types
    #########################################################################################################
    def getDataFrame(self):
        if self.datadf is None:
            self.datadf = self.convertToDataFrame()
        return self.datadf
    
    def getList(self):
        if self.datalist is None:
            self.datalist = self.convertToList()
        return self.datalist
    
    def getArray(self):
        if self.dataarray is None:
            self.dataarray = self.convertToArray()
        return self.dataarray
    
        
        
    #########################################################################################################
    # Convert Data to Correct Format
    #########################################################################################################
    def convertToDataFrame(self, debug=False):
        if self.dtype is None:
            raise ValueError("Data type is unknown!")
            
            
        ## Test array data first
        datadf = None
        if self.dataarray is not None:
            try:
                datadf = DataFrame(self.dataarray)
            except:
                print("Could not convert NP array to DataFrame!")

        if isinstance(datadf, DataFrame):
            return datadf
        
        
        ## Test list data second
        if self.datalist is not None:
            try:
                datadf = DataFrame(self.datalist)
            except:
                print("Could not convert list to DataFrame!")
                
        if isinstance(datadf, DataFrame):
            return datadf
            
            
        print("Could not convert data of type {0} to NP array".format(self.dtype))
        return None
            
            

    def convertToArray(self, debug=False):
        if self.dtype is None:
            raise ValueError("Data type is unknown!")
            
        ## Test list data first
        dataarray = None
        if self.datalist is not None:
            dataarray = asarray(self.datalist)
            try:
                dataarray = asarray(self.datalist)
            except:
                print("Could not convert list data to NP array!")

        if isinstance(dataarray, ndarray):
            return dataarray
        
        
        ## Test pandas data second
        if self.datadf is not None:
            try:
                dataarray = self.datadf.values
            except:
                print("Could not convert DataFrame data to NP array!")
                
        if isinstance(dataarray, ndarray):
            return dataarray
            
            
        print("Could not convert data of type {0} to NP array".format(self.dtype))
        return None
    

    def convertToList(self, debug=False):
        if self.dtype is None:
            raise ValueError("Data type is unknown!")
            
            
        ## Test list data first
        datalist = None
        if self.dataarray is not None:
            try:
                datalist = self.dataarray.tolist()
            except:
                print("Could not convert NP array data to list!")

        if isinstance(datalist, list):
            return datalist
        
        
        ## Test pandas data second
        if self.datadf is not None:
            try:
                dataarray = self.datadf.values.tolist()
            except:
                print("Could not convert DataFrame data to NP array!")
                
        if isinstance(datalist, list):
            return datalist
            
            
        print("Could not convert data of type {0} to NP array".format(self.dtype))
        return None