''' Waters 
    MassLynx Python Chromatogram reader SDK
'''
import copy

import ctypes
from ctypes import*

from array import *

from unidec_modules.waters_importer.MassLynxRawReader import *
#from MassLynxInfoReader import *

class MassLynxRawChromatogramReader(MassLynxRawReader):
    """Read masslynx chromatogram data"""
    def __init__(self, source ):
        super().__init__(source, MassLynxBaseType.CHROM)
   
    #@classmethod
    #def CreateFromPath( cls, path ):                      # alternative constructor - pass class to constructor
    #    return cls(MassLynxRawReader.fromPath( path, 3 ))                     # initalise with reader

    #@classmethod                 
    #def CreateFromReader( cls, sourceReader ):                  # alternative constructor - pass class to constructor
    #     return cls(MassLynxRawReader.fromReader( sourceReader, 3 ))                     # initalise with reader

    def ReadTIC( self, whichFunction ):
        times = []
        intensities = []

        # create the retrun values
        size = c_int(0)
        pTimes = c_void_p()
        pIntensities = c_void_p()
            
        # read tic
        readTIC = MassLynxRawReader.massLynxDll.readTICChromatogram
        readTIC.argtypes = [c_void_p, c_int, POINTER(c_void_p), POINTER(c_void_p), POINTER(c_int)]
        super().CheckReturnCode( readTIC(self._getReader(),whichFunction, pTimes, pIntensities, size) )

        # fill the array
        pT = cast(pTimes,POINTER(c_float))
        pI = cast(pIntensities,POINTER(c_float))

        times = pT[0:size.value]
        intensities = pI[0:size.value]

        # dealocate memory
        MassLynxRawReader.ReleaseMemory( pTimes)
        MassLynxRawReader.ReleaseMemory( pIntensities)

        return times, intensities

    def ReadBPI( self, whichFunction ):
        #times = []
        #intensities = []
        #try:
        # create the retrun values
        size = c_int(0)
        pTimes = c_void_p()
        pIntensities = c_void_p()
            
        # read tic
        readBPI = MassLynxRawReader.massLynxDll.readTICChromatogram
        readBPI.argtypes = [c_void_p, c_int, POINTER(c_void_p), POINTER(c_void_p), POINTER(c_int)]
        super().CheckReturnCode( readBPI(self._getReader(),whichFunction, pTimes, pIntensities, size) )

        # fill the array
        pT = cast(pTimes,POINTER(c_float))
        pI = cast(pIntensities,POINTER(c_float))

        times = pT[0:size.value]
        intensities = pI[0:size.value]

        # dealocate memory
        MassLynxRawReader.ReleaseMemory( pTimes)
        MassLynxRawReader.ReleaseMemory( pIntensities)

        #except MassLynxException as e:
        #    e.Handler()
        #    return [], []

        return times, intensities

    def ReadMassChromatogram( self, whichFunction, whichMass, massTollerance, daughters ):
        # just call multiple mass with list of 1
        whichMasses = [whichMass]
        times, intensities =  self.ReadMassChromatograms( whichFunction, whichMasses, massTollerance, daughters )
    
        return times, intensities[0]

    def ReadMassChromatograms( self, whichFunction, whichMasses, massTollerance, daughters ):
        times = []
        intensities = []

        # get the array of masses
        numMasses = len(whichMasses)
        masses = (c_float * numMasses)(*whichMasses)

        # create the retrun values
        size = c_int(0)
        pTimes = c_void_p()

        # create array of pointers to hold return intensities
        pIntensities = c_void_p()

        readMassChroms = MassLynxRawReader.massLynxDll.readMassChromatograms
        readMassChroms.argtypes = [c_void_p, c_int, POINTER(c_float), c_int, POINTER(c_void_p), POINTER(c_void_p), c_float, c_bool, POINTER(c_int)]
        super().CheckReturnCode( readMassChroms( self._getReader(), whichFunction, masses, numMasses, pTimes, pIntensities, massTollerance, daughters, size))

        # fill the array and free memory
        pT = cast(pTimes,POINTER(c_float))
        times = pT[0:size.value]
        MassLynxRawReader.ReleaseMemory( pTimes)

        # fill in the mass chroms and free memory
        pI = cast(pIntensities ,POINTER(c_float))
        for index in range(0, numMasses ):  
            intensities.append( pI[index * size.value :(index + 1)* size.value ])
#            MassLynxRawReader.ReleaseMemory( ppIntensities[ index] )
        MassLynxRawReader.ReleaseMemory( pIntensities )


        return times, intensities

   # def getMRMsinFunction( self, whichFunction ):
   #     try:
   #         # get the number of MRM transitions
   #         numberMRMs = 0
   #         getMRMsInFunction = RawReader.massLynxDll.getMRMsInFunction
   #         getMRMsInFunction.argtypes = [c_void_p, c_int, POINTER(c_int)]
   #         RawReader.CheckReturnCode(getMRMsInFunction(RawReader.getReader(self), whichFunction, numberMRMs))

   #     except MassLynxException as e:
   #         e.Handler()
   #         return 0

   #     return numberMRMs;

   # def readMRM( self, whichFunction, whichMRM ):
   #     try:
   #         # check we are requesting a valid MRM index - no...
   #         # client  code should not know about data
   ##         if ( self.getMRMsinFunction( whichFunction ) > whichMRM ):

   #         # create the float array
   #         scans = self.getScansInFunction(whichFunction)
   #         times = (c_float*scans)()
   #         intensities = (c_float*scans)()

   #         # check we have some scans ?

   #         readMRM = RawReader.massLynxDll.readMRMChromatogram
   #         readMRM.argtypes = [c_void_p, c_int, c_int, POINTER(c_float), POINTER(POINTER(c_float))]
   #         RawReader.CheckReturnCode( readMRM( RawReader.getReader(self), whichFunction, whichMRM, times, intensities)) # test with invalid MRM
        
   #     except MassLynxException as e:
   #         e.Handler()
   #         return [], []
    
   #     return list(times), list(intensities)


   # def readMRMs( self, whichFunction ):
   #     try:
   #         numberMRMs = self.getMRMsinFunction( whichFunction )

   #         # create the float arrays
   #         scans = self.getScansInFunction(whichFunction)
   #         times = (c_float*scans)()

   #         # create array of pointers to hold return intensities
   #         intensities = (POINTER(c_float) * numberMRMs)()
   #         for index in range(0, numberMRMs ):
   #             intensities[ index ] = (c_float * scans)()
     
   #         readMRMs = RawReader.massLynxDll.readMRMChromatograms
   #         readMRMs.argtypes = [c_void_p, c_int, c_int, POINTER(c_float), POINTER(POINTER(c_float))]
   #         RawReader.CheckReturnCode( readMRMs( RawReader.getReader(self), whichFunction, numberMRMs, times, intensities,))

   #     except MassLynxException as e:
   #         e.Handler()
   #         return [], []

   #     # create the return types
   #     intensities_lists = []
   #     for index in range(0, numberMRMs ):
   #         intensities_lists.append( intensities[index][0:scans] )

   #     return list(times), intensities_lists
