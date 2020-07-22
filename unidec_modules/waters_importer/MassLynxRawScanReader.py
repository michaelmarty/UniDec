''' Waters 
    MassLynx Python SDK
'''
import copy
import ctypes
from ctypes import*
from array import *
from unidec_modules.waters_importer.MassLynxRawReader import *

class MassLynxRawScanReader(MassLynxRawReader):
    """Read masslynx scan data"""
   
    def __init__(self, source ):
        super().__init__(source, MassLynxBaseType.SCAN)


    def ReadScan( self, whichFunction, whichScan ):
        # create the retrun values
        size = c_int(0)
        pMasses = c_void_p()
        pIntensities = c_void_p()

        # read scan
        readScan = MassLynxRawReader.massLynxDll.readScan
        readScan.argtypes = [c_void_p, c_int, c_int, POINTER(c_void_p), POINTER(c_void_p), POINTER(c_int)]
        super().CheckReturnCode( readScan(self._getReader(),whichFunction,whichScan,pMasses,pIntensities,size) )

        # fill the array
        pM = cast(pMasses,POINTER(c_float))
        pI = cast(pIntensities,POINTER(c_float))

        masses = pM[0:size.value]
        intensities = pI[0:size.value]

        # dealocate memory
        #MassLynxRawReader.ReleaseMemory( pMasses)
        #MassLynxRawReader.ReleaseMemory( pIntensities)

        return masses, intensities

    def CombineScan(self, whichFunction, scans):
        # create the retrun values
        size = c_int(0)
        pMasses = c_void_p()
        pIntensities = c_void_p()

        #Create the input data
        nScans = len(scans)
        pScans = (c_int * len(scans))(*scans)

        self.processor = MassLynxRawProcessor(self)
        # read create the function and arguments
        combineScan = MassLynxRawReader.massLynxDll.combineScan
        combineScan.argtypes = [c_void_p, c_int, POINTER(c_int), c_int]
        # run combinescan and check errors
        out = combineScan(self.processor._getProcessor(), whichFunction, pScans, nScans)
        self.processor._codeHandler.CheckReturnCode(out)


        #combineScan.argtypes = [c_void_p, c_int, c_int]
        #out2 = combineScan(self.processor._getProcessor(), whichFunction, 1)
        #print(out2)
        #self.processor._codeHandler.CheckReturnCode(out2)

        # Get scan from the combined
        getScan = MassLynxRawReader.massLynxDll.getScan
        getScan.argtypes = [c_void_p, POINTER(c_void_p), POINTER(c_void_p), POINTER(c_int)]
        out3 = getScan(self.processor._getProcessor(), pMasses, pIntensities, size)
        self.processor._codeHandler.CheckReturnCode(out3)


        # fill the array
        pM = cast(pMasses, POINTER(c_float))
        pI = cast(pIntensities, POINTER(c_float))

        masses = pM[0:size.value]
        intensities = pI[0:size.value]

        # dealocate memory
        # MassLynxRawReader.ReleaseMemory( pMasses)
        # MassLynxRawReader.ReleaseMemory( pIntensities)

        return masses, intensities


    def ReadScanFlags( self, whichFunction, whichScan ):      
        masses = []
        intensities = []
        flags = []
           
        # create the retrun values
        size = c_int(0)
        pMasses = c_void_p()
        pIntensities = c_void_p()
        pFlags = c_void_p()

        # read scan
        readScanFlags = MassLynxRawReader.massLynxDll.readScanFlags
        readScanFlags.argtypes = [c_void_p, c_int, c_int, POINTER(c_void_p), POINTER(c_void_p), POINTER(c_void_p), POINTER(c_int)]
        super().CheckReturnCode( readScanFlags(self._getReader(),whichFunction,whichScan,pMasses,pIntensities,pFlags,size) )

        # fill the array
        pM = cast(pMasses,POINTER(c_float))
        pI = cast(pIntensities,POINTER(c_float))

        masses = pM[0:size.value]
        intensities = pI[0:size.value]

        # check for flags
        if None != pFlags.value:
            pF = cast(pFlags,POINTER(c_byte))
            flags = pF[0:size.value]


        return masses, intensities, flags

    def ReadDriftScan( self, whichFunction, whichScan, whichDrift ):
        masses = []
        intensities = []
               
        # create the retrun values
        size = c_int(0)
        pMasses = c_void_p()
        pIntensities = c_void_p()

        # read scan
        readDriftScan = MassLynxRawReader.massLynxDll.readDriftScan
        readDriftScan.argtypes = [c_void_p, c_int, c_int, c_int, POINTER(c_void_p), POINTER(c_void_p), POINTER(c_int)]
        super().CheckReturnCode( readDriftScan(self._getReader(),whichFunction, whichScan, whichDrift, pMasses,pIntensities,size) )

        # fill the array
        pM = cast(pMasses,POINTER(c_float))
        pI = cast(pIntensities,POINTER(c_float))

        masses = pM[0:size.value]
        intensities = pI[0:size.value]

        # dealocate memory
        #MassLynxRawReader.ReleaseMemory( pMasses)
        #MassLynxRawReader.ReleaseMemory( pIntensities)

        return masses, intensities

    #def readDaughterScan( self, whichFunction, whichScan ):
    #    try:
    #        size = self.getScanSize( whichFunction,whichScan )
            
    #    	# get the daughter scan size
    #        daughtersize = c_int(0)
            
    #        # get daughter size         
    #        getDaughterScanSize =  RawReader.massLynxDll.getDaughterScanSize
    #        getDaughterScanSize.argtypes = [c_void_p, c_int, c_int, POINTER(c_int)]
    #        RawReader.CheckReturnCode( getDaughterScanSize(RawReader.getReader(self),whichFunction,whichScan, daughtersize) )

	   #      # create the float arrays
    #        masses = (c_float*size)()
    #        intensities = (c_float*size)()
    #        daughters = (c_float*daughtersize.value)()
             
    #        # read daughter size
    #        readSpectrumDaughters = RawReader.massLynxDll.readSpectrumDaughters
    #        readSpectrumDaughters.argtypes = [c_void_p, c_int, c_int,  POINTER(c_float), POINTER(c_float), POINTER(c_float)]
    #        RawReader.CheckReturnCode( readSpectrumDaughters(RawReader.getReader(self), whichFunction, whichScan, masses, intensities, daughters) )

    #    except RawReaderException as e:
    #        e.Handler()
    #        return [], [], []

    #    return list(masses), list(intensities), list(daughters)
