#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/icetray-start
#METAPROJECT icerec/V05-02-00

import icecube
from icecube import dataclasses, dataio#, weighting
from I3Tray import *
from icecube.tableio import I3TableWriter
from icecube.dataclasses import I3VectorDouble
from icecube import hdfwriter

import os
import sys
import glob 
import stat
from optparse import OptionParser
from os.path import expandvars
import time
import math
import numpy as np

import PolygonContainment
import Reconstruction
import EventType
import CalibrationAndCleaning
import TauDP_DifferentialWaveform_reOP
import QTot

#from icecube.AtmCscdEnergyReco import AtmCscdEnergyRecoParams
from icecube.recclasses import I3DipoleFitParams
from icecube.recclasses import I3LineFitParams
from icecube.recclasses import I3CLastFitParams
from icecube.recclasses import I3CscdLlhFitParams
from icecube.recclasses import I3TensorOfInertiaFitParams
from icecube.recclasses import I3FillRatioInfo
from icecube.recclasses import CramerRaoParams
from icecube.recclasses import I3StartStopParams
from icecube.gulliver import I3LogLikelihoodFitParams
from icecube.millipede import MillipedeFitParams


def list_callback(option, opt, value, parser):
    noWhiteSpace = value.replace(' ', '')
    noPyAppend   = noWhiteSpace.replace('.py', '')
    cleanList    = noPyAppend.split(',')
    setattr(parser.values, option.dest, cleanList)
    
def string_callback(option, opt, value, parser):
    lower_case   = value.lower()
    no_hyphens   = lower_case.replace('-', '')
    no_underscore= no_hyphens.replace('_', '')
    setattr(parser.values, option.dest, no_underscore)

usage = "%prog [options] <inputfiles>"

parser = OptionParser(usage=usage)

parser.add_option(
    "-i", "--inputfile",
    type      = "string",
    action    = "callback",
    callback  = list_callback,
    default   = "",
    metavar   = "<input file name>",
    help      = "List of the input files",
    )

parser.add_option(
    "-g", "--gcdfile",
    type      = "string",
    action    = "store",
    default   = "",
    metavar   = "<geo file>",
    help      = "Name of GCD file",
    )

parser.add_option(
    "-m", "--MC",
    type      = "int",
    action    = "store",
    default   = "1",
    metavar   = "<0 if data, 1 if neutrino MC, 2 if corsika>",
    help      = "Input data file type, 0 if real data, 1 if neutrino MC, 2 if corsika",
    )

parser.add_option(
    "-y", "--year",
    type      = "int",
    action    = "store",
    default   = "2012",
    metavar   = "<year of files>",
    help      = "Year the files were produced for",
    )

parser.add_option(
    "-o", "--outputfile",
    type      = "string",
    action    = "store",
    default   = "TEST",
    metavar   = "<output file(s) name>",
    help      = "Name of the output file(s), i.e. .root and .i3.gz names",
    )

parser.add_option(
    "-p", "--outputpath",
    type      = "string",
    action    = "store",
    default   = "",
    metavar   = "<output path name>",
    help      = "Path to the output file(s), i.e. .root and .i3.gz names",
    )

parser.add_option(
    "-c", "--condense",
    action="store_true",
    default   = False,
    metavar   = "<Condense output files>",
    help      = "Removes empty files and stores GCD if true",
    )

(options, args) = parser.parse_args()

tray = I3Tray()

files = options.inputfile

if options.gcdfile != "":
    files = [options.gcdfile] + files

print "Input files: ", files

year = options.year

tray.AddModule('I3Reader', 'reader',
               SkipKeys=["CalibrationErrata", "I3ModuleGeoMap", "I3OMGeoMap", "I3StationGeoMap", "Subdetectors", "OfflineInIceCalibrationErrata", "CalibrationWaveformRange", "CalibratedWaveformRange",'OfflinePulsesTimeRange'], FilenameList = files)

#tray.AddSegment(Filter.Filter13, "ehe-filtering")


# Waveform calibration

def EHEfilter(frame):
    if frame["FilterMask"].has_key("EHEFilter_10"): #IC79
        if frame["FilterMask"]["EHEFilter_10"].condition_passed:
            return True
        else:
            return False

if year == 2010 or year == 2011:
    tray.AddSegment(QTot.CalQTot, "selfveto-qtot", pulses='OfflinePulses')
    
if year >=2012:
    tray.AddModule(lambda frame: frame["I3EventHeader"].sub_event_stream == 'InIceSplit' , "Split",Streams = [icetray.I3Frame.Physics])
    tray.AddSegment(QTot.CalQTot, "selfveto-qtot", pulses='SplitInIcePulses') 
    
tray.AddSegment(CalibrationAndCleaning.CalibrationAndCleaning, "cal-and-clean")        

# Raise charge cut to log10(QTot)>3.3, changed on July 17, 2014
tray.AddModule(lambda frame: frame["CausalQTot"].value>1995.26231497, "qtotcut")

# Double pulse waveform identification. DISCARD events don't have double pulse waveforms.
tray.AddSegment( TauDP_DifferentialWaveform_reOP.CalWaveformDerivatives, "dp_cutsLC",
                           DiscardEvents=True, 
                           WfQtot=10000,
                           bins_p1=1,
                           bins_trailing=2,
                           bins_p2=3,
                           Amp1LC=1.,
                           Amp2LC=12.,
                           AmpTrailingLC=-0.5,
                           Amp1SD = 10.,
                           Amp2SD = 18.,
                           AmpTrailingSD = -17.,
                           der_step = 4)

    
tray.AddSegment(Reconstruction.OfflineCascadeReco, "CscdReco", suffix="_DP", Pulses='HLCPulses')
tray.AddSegment(Reconstruction.MuonReco, "MuonReco", Pulses='HLCPulses')


# rlogl cut
tray.Add(lambda fr: (fr['SPEFit32_DPFitParams'].rlogl - fr['CascadeLlhVertexFit_DPParams'].ReducedLlh)>-0.5, 'rlogl_cut')
tray.Add(lambda fr: (fr['depthFirstHit'].value)<475, 'z_cut')


def rlogl(fr):
    fr['rlogl'] = dataclasses.I3Double(fr['SPEFit32_DPFitParams'].rlogl - fr['CascadeLlhVertexFit_DPParams'].ReducedLlh)
    return 1
tray.AddModule(rlogl,'rlogl')

if year == 2010:
    geo ='ic79'
elif year >= 2011:
    geo = 'ic86'
    
pulses='HLCPulses'
    
tray.AddSegment(Reconstruction.OfflineCascadeReco_noDC, "CscdReco_noDC", suffix="_noDC_DP", Pulses=pulses)
tray.AddSegment(Reconstruction.MuonReco_noDC, "MuonReco_noDC", Pulses=pulses)

tray.AddSegment(PolygonContainment.PolygonContainment, 'polyfit', geometry = geo,RecoVertex='VHESelfVetoVertexPos',outputname='_Veto')

tray.Add(lambda fr: (fr['SPEFit32_noDC_DPFitParams'].rlogl - fr['CascadeLlhVertexFit_noDC_DPParams'].ReducedLlh)>-0.15, 'noDC_rlogl_cut')

contain = 10
b = 400
m = - 1./3.
bottomz = -200
bottome = 75

tray.Add(lambda fr: (fr["LeastDistanceToPolygon_Veto"].value)>contain, 'contain')
tray.Add(lambda fr: (fr["LeastDistanceToPolygon_Veto"].value)>bottome  or fr['VHESelfVetoVertexPos'].z>bottomz, 'bottom corner')
tray.Add(lambda fr: ((m * fr["LeastDistanceToPolygon_Veto"].value + fr['VHESelfVetoVertexPos'].z) < b), 'top corner')

pframes = 0
def count_physics(fr):
    global pframes
    pframes += 1
    print fr['I3EventHeader']

def fixDST(fr, pulseMaskName, newPulseMapName):
    if pulseMaskName in fr:
        pulsemap = fr[pulseMaskName].apply(fr)
        newMap = dataclasses.I3RecoPulseSeriesMap()
        for (omkey, pulses) in pulsemap:
            series = dataclasses.I3RecoPulseSeries()
            for p in pulses:
                if p.width <= 0.:
                    p.time -= 0.51 * I3Units.ns
                    p.width = 0.5 * I3Units.ns
                series.append(p)
            newMap[omkey] = series
        fr[newPulseMapName] = newMap
        
doFixDST = lambda fr: fixDST(fr, pulses, 'Fixed'+pulses)
tray.Add(doFixDST, Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])
    
tray.AddModule(count_physics, "count",Streams = [icetray.I3Frame.Physics])
tray.AddSegment(Reconstruction.Monopod, 'monopod', Pulses='Fixed'+pulses)
#tray.AddSegment(Reconstruction.Monopod, 'monopod', Pulses=pulses)


def CorsEventType (frame):
    frame["EventType"] = icetray.I3Int(8)
    return

if options.MC == 2:
    tray.AddModule( CorsEventType, "cors", Streams = [icetray.I3Frame.Physics])
    tray.Add(weighting.get_weighted_primary, 'pr')

if options.MC == 1:
    tray.AddSegment(EventType.EventType, "nufinder")

if options.MC == 0 and options.condense == True:
    outputstreams = [icetray.I3Frame.Geometry, icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus,icetray.I3Frame.Physics, icetray.I3Frame.DAQ]
else:
    outputstreams = [icetray.I3Frame.Physics, icetray.I3Frame.DAQ]

if options.condense == True:
    outputfile = os.getenv("_CONDOR_SCRATCH_DIR") + '/lwille_tmp'
else:
    outputfile = options.outputpath + options.outputfile

   
    
FrameObjectsToKeep = ['Amp1',
                      'Amp2',
                      'AmpTrailing',
                      'BinsToT1',
                      'BinsTbT',
                      'BinsToT2',
                      'DPWaveformQTot',
                      'DP_OMs',
                      'EventType',
                      'DPFlagSD',
                      'DPFlagLC',
                      'QTot', 
                      'I3MCWeightDict',
                      'IntPos',
                      'MCPrimary',
                      'TauLength', 
                      'TauDecay',
                      'FilterMask',
                      'I3EventHeader',
                      'SPEFit2',
                      'SPEFit2FitParams',
                      'SPEFit32_DP',
                      'SPEFit32_DPFitParams',
                      'SPEFit32_DP_Characteristics',
                      'SPEFit32_noDC_DP',
                      'SPEFit32_noDC_DPFitParams',
                      'MPEFit',
                      'MPEFitFitParams',
                      'MPEFitCramerRaoParams',
                      'MPEFitCharacteristics',
                      'MPEFitMuEX',
                      'LineFit',
                      'LineFitParams',
                      'CascadeContainmentTagging_DP',
                      'CascadeDipoleFit_DP',
                      'CascadeDipoleFit_DPParams',
                      'CascadeFillRatio_DP',
                      'CascadeFillRatio_noDC_DP',
                      'CascadeLast_DP',
                      'CascadeLast_DPParams',
                      'CascadeLineFit_DP',
                      'CascadeLineFit_DPParams',
                      'CascadeLlhVertexFit_DP',
                      'CascadeLlhVertexFit_DPParams',
                      'CascadeLast_noDC_DP',
                      'CascadeLast_noDC_DPParams',
                      'CascadeLineFit_noDC_DP',
                      'CascadeLineFit_noDC_DPParams',
                      'CascadeLlhVertexFit_noDC_DP',
                      'CascadeLlhVertexFit_noDC_DPParams',
                      'L4VetoLayer0',
                      'L4VetoLayer1',
                      'L4VetoLayerQTot',
                      'L4VetoLayerVertexTime',
                      'VHESelfVeto',
                      'VHESelfVetoVertexPos',
                      'VHESelfVetoVertexTime',
                      'CausalQTot',
                      'depthFirstHit',
                      'earliestLayer',
                      'maxDomChargeLayer',
                      "LeastDistanceToPolygon",
                      "LeastDistanceToCornerString",
                      'MonopodAngular', 
                      'MonopodAngularFitParams', 
                      'MonopodAngularGuess', 
                      'MonopodAngularGuessFitParams',
                      'cog_noDC_10_bpercent',
                      'cog_noDC_10_bpercent_qweighted',
                      'cog_noDC_20_bpercent',
                      'cog_noDC_20_bpercent_qweighted',
                      'cog_noDC_50_bpercent',
                      'cog_noDC_50_bpercent_qweighted',
                      'cog_noDC_50_epercent',
                      'cog_noDC_50_epercent_qweighted',
                      'cog_noDC_80_epercent',
                      'cog_noDC_80_epercent_qweighted',
                      'cog_noDC_90_epercent',
                      'cog_noDC_90_epercent_qweighted',
                      'MuonWeight',
                      'SplineMPErecommended',
                      'SplineMPErecommendedFitParams',
                      "LeastDistanceToPolygon_Veto"
                      ]

tray.AddSegment( hdfwriter.I3HDFWriter, 'hdf5writer',
                Output = outputfile +'.h5',
                Keys = FrameObjectsToKeep,
                Types=[],
                SubEventStreams=['InIceSplit','nullsplit','inice','in_ice']
               )

tray.AddModule( 'TrashCan' , 'Done' )

tray.Execute()

tray.Finish()

if pframes > 0 and options.condense == True:
    tray2 = I3Tray()
    tray2.AddModule('I3Reader', FilenameList=[outputfile+ ".i3.gz"])
    tray2.AddModule( 'I3Writer', 'EventWriter',
                CompressionLevel = 9,
                Filename          = options.outputpath + options.outputfile + ".i3.gz",
                Streams           = outputstreams,
                SkipKeys          = ["ATWDPortiaPulse", "FADCPortiaPulse"]
                )
    tray2.AddSegment( hdfwriter.I3HDFWriter, 'hdf5writer',
                Output = options.outputpath + options.outputfile + '.h5',
                Keys = FrameObjectsToKeep,
                Types=[],
                SubEventStreams=['InIceSplit','nullsplit','inice','in_ice']
               )
                
    tray2.Execute()
