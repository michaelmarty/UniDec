"""
     Waters
    MassLynx Python SDK
"""

import ctypes
from ctypes import*

from unidec.UniDecImporter.Waters.MassLynxRawReader import MassLynxRawReader
from unidec.UniDecImporter.Waters.MassLynxRawReader import MassLynxBaseType


class MassLynxRawInfoReader(MassLynxRawReader):

    def __init__(self, source ):
        super().__init__(source, MassLynxBaseType.INFO)


    def GetNumberofFunctions( self ):
        size = c_int(0)
        getFunctionCount =  MassLynxRawReader.massLynxDll.getFunctionCount
        getFunctionCount.argtypes = [c_void_p, POINTER(c_int)]
        super().CheckReturnCode( getFunctionCount(self._getReader(),size) )

        return size.value

    def GetScansInFunction( self, whichFunction ):
        size = c_int(0)   
        getScanCount = MassLynxRawReader.massLynxDll.getScanCount
        getScanCount.argtypes = [c_void_p, c_int, POINTER(c_int)]
        super().CheckReturnCode( getScanCount(self._getReader(),whichFunction,size) )
 
        return size.value

    def GetAcquisitionMassRange( self, whichFunction ):
        lowMass = c_float(0)
        highMass = c_float(0)
        getAcquisitionMassRange = MassLynxRawReader.massLynxDll.getAcquisitionMassRange
        getAcquisitionMassRange.argtypes = [c_void_p, c_int, c_int, POINTER(c_float), POINTER(c_float)]
        super().CheckReturnCode( getAcquisitionMassRange(self._getReader(),whichFunction, 0,lowMass,highMass) )
 
        return lowMass.value, highMass.value

    def GetAcquisitionTimeRange( self, whichFunction ):
        startTime = c_float(0)
        endTime = c_float(0)
        getAcquisitionTimeRange = MassLynxRawReader.massLynxDll.getAcquisitionTimeRange
        getAcquisitionTimeRange.argtypes = [c_void_p, c_int, POINTER(c_float), POINTER(c_float)]
        super().CheckReturnCode( getAcquisitionTimeRange(self._getReader(),whichFunction,startTime,endTime) )
 
        return startTime.value, endTime.value

    def GetFunctionType( self, whichFunction ):
        functionType = c_int(0)
        getFunctionType = MassLynxRawReader.massLynxDll.getFunctionType
        getFunctionType.argtypes = [c_void_p, c_int, POINTER(c_int)]
        super().CheckReturnCode( getFunctionType(self._getReader(),whichFunction, functionType) )
 
        return functionType.value

    def GetFunctionTypeString( self, functionType ):
        functionTypeString = c_char_p()
        temp = c_char_p()
        getFunctionTypeString = MassLynxRawReader.massLynxDll.getFunctionTypeString
        getFunctionTypeString.argtypes = [c_void_p, c_int, POINTER(c_char_p)]
        self.CheckReturnCode( getFunctionTypeString(self._getReader(),functionType, temp) )

        return self.ToString( temp )

    def IsContinuum( self, whichFunction ):
        continuum = c_bool(0)
        isContinuum = MassLynxRawReader.massLynxDll.isContinuum
        isContinuum.argtypes = [c_void_p, c_int, POINTER(c_bool)]
        super().CheckReturnCode( isContinuum(self._getReader(),whichFunction, continuum) )
 
        return continuum.value

    def GetIonMode( self, whichFunction ):
        ionMode = c_int()
        getIonMode = MassLynxRawReader.massLynxDll.getIonMode
        getIonMode.argtypes = [c_void_p, c_int, POINTER(c_int)]
        super().CheckReturnCode(getIonMode(self._getReader(),whichFunction, ionMode ))
            
        return ionMode.value

    def GetIonModeString( self, ionMode ):
        ionModeString = c_char_p()
        temp = c_char_p()
        getIonModeString = MassLynxRawReader.massLynxDll.getIonModeString
        getIonModeString.argtypes = [c_void_p, c_int, POINTER(c_char_p)]
        super().CheckReturnCode(getIonModeString(self._getReader(),ionMode,  temp ))
        return super().ToString(temp)

    def GetHeaderItems( self, whichItems ):
        itemString = c_char_p()
        nItems = len(whichItems )

        temp = c_char_p()
        items = (c_int * nItems)(*whichItems)
        delimiter = ctypes.create_string_buffer(1)

        getHeaderItem = MassLynxRawReader.massLynxDll.getHeaderItems
        getHeaderItem.argtypes = [c_void_p, POINTER(c_int), POINTER(c_char_p), c_int, POINTER(c_char)]
        super().CheckReturnCode( getHeaderItem( self._getReader(), items, temp, nItems, delimiter))

        itemString = super().ToString(temp) #temp.value.decode()
#        MassLynxRawReader.ReleaseMemory( temp)

        delim = delimiter.value.decode()
        return itemString.split(delim)

    def GetHeaderItem( self, whichItem ):
        whichItems =  list()
        whichItems.append( whichItem )
        values =  self.GetHeaderItems( whichItems)
        
        return values[0] 

    # scan stats
    def GetScanItem( self, whichFunction, whichScan, whichItem ):
        whichItems = list()
        whichItems.append( whichItem )
        values =  self.GetScanItems(whichFunction, whichScan, whichItems)

        return values[0] 

    def GetScanItems(self, whichFunction, whichScan, whichItems):
        nItems = len(whichItems)
        temp = c_char_p()
        items = (c_int * nItems)(*whichItems)
        delimiter = ctypes.create_string_buffer(1)

        getScanItem = MassLynxRawReader.massLynxDll.getScanItems
        getScanItem.argtypes = [c_void_p, c_int, c_int, POINTER(c_int), POINTER(c_char_p), c_int, POINTER(c_char)]
        super().CheckReturnCode( getScanItem( self._getReader(), whichFunction, whichScan, items, temp, nItems, delimiter))

        itemString = super().ToString(temp)
 
        delim = delimiter.value.decode()
        return itemString.split(delim)

    def GetScanItemString(self, whichItems):
         # get the array of items
        nItems = len(whichItems)
        temp = c_char_p()
        items = (c_int * nItems)(*whichItems)
        delimiter = ctypes.create_string_buffer(1)

        getScanItemNames = MassLynxRawReader.massLynxDll.getScanItemNames
        getScanItemNames.argtypes = [c_void_p, POINTER(c_int), POINTER(c_char_p), c_int, POINTER(c_char)]
        super().CheckReturnCode( getScanItemNames( self._getReader(), items, temp, nItems, delimiter))

        itemString = super().ToString(temp)
 
        delim = delimiter.value.decode()
        return itemString.split(delim)

    def GetItemsInFunction( self, whichFunction, nWhichScan ):
        size = c_int(0)
        pItems = c_void_p()
        getItemsInFunction = MassLynxRawReader.massLynxDll.getItemsInFunction
        getItemsInFunction.argtypes = [c_void_p, c_int, POINTER(c_void_p), POINTER(c_int)]
        super().CheckReturnCode( getItemsInFunction(self._getReader(),whichFunction,pItems,size) )

        # fill the array
        pI = cast(pItems,POINTER(c_int))
        items = pI[0:size.value]
       
        # dealocate memory
        MassLynxRawReader.ReleaseMemory( pItems)

        return items

    def GetRetentionTime( self, whichFunction, nWhichScan ):
        retentionTime = c_float(0)
        getRetentionTime = MassLynxRawReader.massLynxDll.getRetentionTime
        getRetentionTime.argtypes = [c_void_p, c_int, c_int, POINTER(c_float)]
        super().CheckReturnCode( getRetentionTime(self._getReader(),whichFunction,nWhichScan,retentionTime) )
 
        return retentionTime.value

    def GetDriftTime( self, whichFunction, nWhichDrift ):
        driftTime = c_float(0)
        getDriftTime = MassLynxRawReader.massLynxDll.getDriftTime
        getDriftTime.argtypes = [c_void_p, c_int, c_int, POINTER(c_float)]
        super().CheckReturnCode( getDriftTime(self._getReader(),whichFunction,nWhichDrift,driftTime) )
 
        return driftTime.value

    def CanLockMassCorrect( self ):
        canApply = c_bool(0)
        canLockMassCorrect =  MassLynxRawReader.massLynxDll.canLockMassCorrect
        canLockMassCorrect.argtypes = [c_void_p, POINTER(c_bool)]
        super().CheckReturnCode( canLockMassCorrect(self._getReader(), canApply) )

        return canApply.value

    
    def IsLockMassCorrected( self ):
        corrected = c_bool(0)
        mlMethod =  MassLynxRawReader.massLynxDll.isLockMassCorrected
        mlMethod.argtypes = [c_void_p, POINTER(c_bool)]
        super().CheckReturnCode( mlMethod(self._getReader(), corrected) )

        return corrected.value


