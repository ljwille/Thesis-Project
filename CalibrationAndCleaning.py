#!/usr/bin/env python

#---------------------------------------------------------------
# Do wave calibrations and wave splitting. Then feature extraction
# and SRT and TW cleaning.
#---------------------------------------------------------------

import os
from I3Tray import *
from os.path import *
import sys

from icecube import tableio
from icecube.tableio import I3TableWriter
#from icecube.rootwriter import I3ROOTTableService
from icecube import icetray, dataclasses, dataio, WaveCalibrator
from icecube.STTools.seededRT.configuration_services import I3DOMLinkSeededRTConfigurationService

@icetray.traysegment
def CalibrationAndCleaning(tray, name): #Adapted from Chang Hyon's EHE L2
    def CheckRawData(frame):
        if frame.Has("InIceRawData"):
            #print "Found InIceRawData for This Event.."
            return True
        else:
            #print "fail"
            return False
    tray.AddModule(CheckRawData, "check-raw-data", Streams = [icetray.I3Frame.Physics, icetray.I3Frame.DAQ])
    #tray.Add(CheckRawData)
    
    tray.AddModule("I3WaveCalibrator", "calibrator")(
        #("Launches", "CleanInIceRawData"),  # i.e. Launches
        ("Launches", "InIceRawData"),  # EHE burn sample IC86
        ("Waveforms", "CalibratedWaveforms"),
        ("ATWDSaturationMargin",123), # 1023-900 == 123
        #("DOMsimulatorWorkArounds", False), # i.e. correctForDOMSimulator
        ("FADCSaturationMargin",  0), # i.e. FADCSaturationMargin
        ("Errata", "OfflineInIceCalibrationErrata"), # SAVE THIS IN L2 OUTPUT?
        )
    tray.AddModule("I3WaveformSplitter", "waveformsplit")(
        ("Input","CalibratedWaveforms"),
        ("HLC_ATWD","CalibratedWaveformsHLCATWD"),
        ("HLC_FADC","CalibratedWaveformsHLCFADC"),
        ("SLC","CalibratedWaveformsSLC"),
        ("Force",True),
        )

        
    return
