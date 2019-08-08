#---------------------------------------------------------------
# Run some reconstructions to separate muons and nutau
#---------------------------------------------------------------

#!/usr/bin/env python

import os
import math
import time
import sys

from I3Tray import *
from icecube import icetray, dataclasses, dataio, tableio, CascadeVariables, millipede, photonics_service, wavedeform
from icecube.icetray import I3Units
from icecube.dataclasses import I3RecoPulse, I3Particle, I3Waveform, I3Double
from icecube import linefit, dipolefit, clast, cscd_llh, fill_ratio, tensor_of_inertia
from icecube import truncated_energy
from icecube import spline_reco
import icecube


load("libgulliver")                 # Gulliver
load("libgulliver-modules")         # Gulliver Module
load("liblilliput")                 # 
load("libcscd-llh")
load("libphotonics-service")
load("libDomTools")    

@icetray.traysegment
def OfflineCascadeReco( tray, name, If = lambda f: True, suffix = '',
                        Pulses = '', 
                        CascadeLineFit = 'CascadeLineFit',
                        CascadeDipoleFit = 'CascadeDipoleFit',
                        CascadeLast = 'CascadeLast',
                        CascadeLlhVertexFit = 'CascadeLlhVertexFit',
                        #CascadeLlhVertexFitSplit = 'CascadeLlhVertexFitSplit',
                        BadDOMListName = 'BadDomsList',
                        PhotonicsServiceName = '',
                        #AtmCscdEnergyReco = 'AtmCscdEnergyReco',
                        CascadeFillRatio = 'CascadeFillRatio',
                        CascadeFillRatio_noDC = 'CascadeFillRatio_noDC',
                        #CascadeSplitPulses = 'CascadeSplitPulses',
                        #CascadeLineFitSplit = 'CascadeLineFitSplit',
                        #CascadeToISplit = 'CascadeToISplit',
                        #CascadeImprovedLineFit = 'CascadeImprovedLineFit',
                        CascadeContainmentTagging = 'CascadeContainmentTagging',
                        ):

    tray.AddModule( 'I3LineFit', name + '_CascadeLinefit' + suffix,
                    Name = CascadeLineFit + suffix, # ! Name of fit
                    InputRecoPulses = Pulses, 
                    LeadingEdge = 'FLE', # ! Use only first leading edge, for Cascades especially
                    If = If,
                    )

    tray.AddModule( 'I3DipoleFit', name + '_CascadeDipolefit' + suffix,
                    AmpWeightPower = 0,
                    DipoleStep = 0,
                    InputRecoPulses = Pulses, 
                    MinHits =  5,
                    Name = CascadeDipoleFit + suffix,
                    If = If,
                    )
    
    tray.AddModule('I3CLastModule', name + '_CascadeLast' + suffix, 
                   Name = CascadeLast + suffix,
                   InputReadout = Pulses,
                   If = If,
                   )
    
    tray.AddModule( 'I3CscdLlhModule', name + '_CascadeLlh' + suffix,
                    InputType = 'RecoPulse', # ! Use reco pulses
                    RecoSeries = Pulses, # ! Name of input pulse series
                    FirstLE = True, # Default
                    SeedWithOrigin = False, # Default
                    SeedKey = CascadeLast + suffix, # ! Seed fit - CLast reco
                    MinHits = 8, # ! Require 8 hits
                    AmpWeightPower = 0.0, # Default
                    ResultName = CascadeLlhVertexFit + suffix, # ! Name of fit result
                    Minimizer = 'Powell', # ! Set the minimizer to use
                    PDF = 'UPandel', # ! Set the pdf to use
                    ParamT = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    ParamX = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    ParamY = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    ParamZ = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    If = If,
                    )
    # fill ratio looks in D frame for BadDomsList. If you want others, you must add them with the BadOMs paramters
    tray.AddModule('I3FillRatioModule', name + '_CascadeFillRatio' + suffix,
                   AmplitudeWeightingPower = 0,
                   #DetectorConfig = 1, #! Full Detector
                   BadDOMList = BadDOMListName, #! Get the list of Bad DOMs from the frame
                   RecoPulseName = Pulses,
                   ResultName = CascadeFillRatio + suffix,
                   SphericalRadiusMean = 1.0,
                   SphericalRadiusRMS = 1.0,
                   SphericalRadiusMeanPlusRMS = 1.0,
                   SphericalRadiusNCh = 1.0,
                   VertexName = CascadeLlhVertexFit + suffix,
                   If = If,
                   )
    #make DeepCore OMkeys
    deepcore=[]
    for dstring in range(79,87):
        for ddom in range(1,61):
            deepcore.append(icetray.OMKey(dstring,ddom))
    #excluding DeepCore for fill ratio calculation
    tray.AddModule('I3FillRatioModule', name + '_CascadeFillRatio_noDC' + suffix,
                   AmplitudeWeightingPower = 0,
                   #DetectorConfig = 1, #! Full Detector                       
                   BadDOMList = BadDOMListName, #! Get the list of Bad DOMs from the frame             
                   BadOMs = deepcore,
                   RecoPulseName = Pulses,
                   ResultName = CascadeFillRatio_noDC + suffix,
                   SphericalRadiusMean = 1.0,
                   SphericalRadiusRMS = 1.0,
                   SphericalRadiusMeanPlusRMS = 1.0,
                   SphericalRadiusNCh = 1.0,
                   VertexName = CascadeLlhVertexFit + suffix,
                   If = If,
                   )

    #tray.AddSegment( improvedLinefit.simple, CascadeImprovedLineFit+suffix, inputResponse = Pulses, fitName = CascadeImprovedLineFit+suffix, If = If )

    # Containment Tagging
    tray.AddModule('I3VetoModule', 'ContainmentTag'+suffix, 
                   HitmapName=Pulses,
                   OutputName=CascadeContainmentTagging + suffix,
                   DetectorGeometry=86, 
                   useAMANDA=False,
                   FullOutput=True,
                   If = If,
                   )
    def makedouble(frame):
        if frame.Has(CascadeContainmentTagging + suffix):
            contmap=frame[CascadeContainmentTagging + suffix]        
            frame['earliestLayer']=dataclasses.I3Double(contmap.earliestLayer)
            frame['maxDomChargeLayer']=dataclasses.I3Double(contmap.maxDomChargeLayer)
            frame['depthFirstHit']=dataclasses.I3Double(contmap.depthFirstHit)

    tray.AddModule(makedouble,'makedoubles')

    '''
    tray.AddModule('AtmCscdEnergyReco', name + '_CascadeEnergy' + suffix,
                   InputRecoPulses = Pulses,
                   Output = AtmCscdEnergyReco + suffix,
                   BadDOMListName = BadDOMListName,
                   CascadeVertex = CascadeLlhVertexFit + suffix,  # ! Seed fit - CascadeLlh reco
                   PhotonicsServiceName = PhotonicsServiceName,
                   NoiseRate = 450 * I3Units.hertz, 
                   If = If,
                   )
                   '''
    
@icetray.traysegment
def OfflineCascadeReco_noDC( tray, name, If = lambda f: True, suffix = '_noDC_DP',
                        Pulses = '', 
                        CascadeLineFit = 'CascadeLineFit',
                        #CascadeDipoleFit = 'CascadeDipoleFit',
                        CascadeLast = 'CascadeLast',
                        CascadeLlhVertexFit = 'CascadeLlhVertexFit',
                        #CascadeLlhVertexFitSplit = 'CascadeLlhVertexFitSplit',
                        BadDOMListName = 'BadDomsList',
                        PhotonicsServiceName = '',
                        #AtmCscdEnergyReco = 'AtmCscdEnergyReco',
                        #CascadeFillRatio = 'CascadeFillRatio',
                        #CascadeFillRatio_noDC = 'CascadeFillRatio_noDC',
                        CascadeSplitPulses = 'CascadeSplitPulses',
                        #CascadeLineFitSplit = 'CascadeLineFitSplit',
                        #CascadeToISplit = 'CascadeToISplit',
                        #CascadeImprovedLineFit = 'CascadeImprovedLineFit',
                        #CascadeContainmentTagging = 'CascadeContainmentTagging',
                        ):
    
    #make a pulse map without DeepCore
    def pulse_noDC(frame):
        frame[Pulses + '_noDC'] = dataclasses.I3RecoPulseSeriesMapMask(frame, Pulses, lambda omkey, index, pulse: omkey.string < 79)
    tray.AddModule(pulse_noDC)

    tray.AddModule( 'I3LineFit', name + '_CascadeLinefit' + suffix,
                    Name = CascadeLineFit + suffix, # ! Name of fit
                    InputRecoPulses = Pulses + '_noDC', 
                    LeadingEdge = 'FLE', # ! Use only first leading edge, for Cascades especially
                    If = If,
                    )

    tray.AddModule('I3CLastModule', name + '_CascadeLast' + suffix, 
                   Name = CascadeLast + suffix,
                   InputReadout = Pulses + '_noDC',
                   If = If,
                   )
    
    tray.AddModule( 'I3CscdLlhModule', name + '_CascadeLlh' + suffix,
                    InputType = 'RecoPulse', # ! Use reco pulses
                    RecoSeries = Pulses + '_noDC', # ! Name of input pulse series
                    FirstLE = True, # Default
                    SeedWithOrigin = False, # Default
                    SeedKey = CascadeLast + suffix, # ! Seed fit - CLast reco
                    #ExcludedOMs = deepcore ### Excluding DeepCore!!!
                    MinHits = 8, # ! Require 8 hits
                    AmpWeightPower = 0.0, # Default
                    ResultName = CascadeLlhVertexFit + suffix, # ! Name of fit result
                    Minimizer = 'Powell', # ! Set the minimizer to use
                    PDF = 'UPandel', # ! Set the pdf to use
                    ParamT = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    ParamX = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    ParamY = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    ParamZ = '1.0, 0.0, 0.0, false',   # ! Setup parameters
                    If = If,
                    )
    
@icetray.traysegment
def MuonReco(tray, name, Pulses='', If=lambda f: True):
##-----------This part was adapted from Chang Hyon's Big Bird Reonstruction script------------
## SPE 32 reconstruction
##############################################################################################
    # Services to do Gulliver reconstruction
    tray.AddService("I3SimpleParametrizationFactory","SimpleTrack",
                    StepX=20*I3Units.m,                                     # ! 20m step size
                    StepY=20*I3Units.m,                                     # ! 20m step size
                    StepZ=20*I3Units.m,                                     # ! 20m step size
                    StepZenith=0.1 * I3Units.radian,                        # ! 0.1 radian step size in zenith
                    StepAzimuth=0.2 * I3Units.radian,                       # ! 0.2 radian step size in azimuth
                    StepT= 0. ,                                             # Default
                    StepLinE= 0. ,                                          # Default
                    StepLogE= 0. ,                                          # Default
                    StepLinL= 0. ,                                             # Default
                    BoundsX= [ -2000 * I3Units.m,2000 * I3Units.m ] ,       # ! Set bounds to +-2000m
                    BoundsY= [ -2000 * I3Units.m,2000 * I3Units.m] ,        # ! Set bounds to +-2000m
                    BoundsZ= [ -2000 * I3Units.m,2000 * I3Units.m ] ,       # ! Set bounds to +-2000m
                    BoundsZenith= [ 0., 0. ] ,                              # Default
                    BoundsAzimuth= [ 0., 0. ] ,                             # Default
                    BoundsT= [ 0., 0. ] ,                                   # Default
                    BoundsLinE= [ 0., 0. ] ,                                   # Default
                    BoundsLinL= [ 0., 0. ] ,                                   # Default
                    )
    
    # Define the gulliver minimization sevice to use
    tray.AddService( "I3GulliverMinuitFactory", "Minuit",
                     Algorithm= "SIMPLEX",                                  # Default
                     Tolerance= 0.01,                                       # ! change to 0.01
                     MaxIterations= 10000,                                  # Default
                     MinuitPrintLevel= -2,                                  # Default
                     MinuitStrategy= 2,                                     # Default
                     FlatnessCheck= True,                                   # Default
                     )
    
    # Use convoluted pandel as the PDF for the likelihood
    tray.AddService("I3GulliverIPDFPandelFactory","Pandel",
                    InputReadout= Pulses,                               # ! Name of pulses to use
                    Likelihood= "SPE1st",                                  # Default
                    PEProb= "GaussConvoluted",                                  # Default
                    IceModel= 2,                                           # Default
                    IceFile= "",                                           # Default
                    AbsorptionLength= 98.0 * I3Units.m,                    # Default
                    JitterTime= 15.0 * I3Units.ns,                         # Default
                    NoiseProbability= 1.0*I3Units.hertz * 10.0*I3Units.ns,  # ! Added a little noise term
                    )

    # linefit seed service
    tray.AddService( "I3BasicSeedServiceFactory", "LFSeed",
                     FirstGuesses= ["CascadeLineFit_DP"],                   # ! Use spe2 and LF (in case spe2 fails) as seeds
                     InputReadout= Pulses,                                  # ! Use pulses for vertex correction
                     TimeShiftType= "TFirst",                               # ! Use TFirst for vertex correction
                     SpeedPolice= True,                                     # Default
                     MaxMeanTimeResidual= 1000.0 * I3Units.ns,              # Default
                     )

    # track fit
    tray.AddModule( "I3IterativeFitter", "SPEFit32_DP",
                    RandomService= "SOBOL",                                # Default
                    NIterations= 32,                                       # ! Nunmber of iterations
                    SeedService= "LFSeed",                               # ! Name of seed service
                    Parametrization= "SimpleTrack",                        # ! Name of track parametrization service
                    LogLikelihood= "Pandel",                               # ! Name of likelihood service
                    CosZenithRange= [ -1, 1 ],                             # Default
                    Minimizer= "Minuit",                                   # ! Name of minimizer service
                    #If = InIceCscd,
                    If = If
                    )
    
    return

@icetray.traysegment
def MuonReco_noDC(tray, name, Pulses='', If=lambda f: True):
##-----------This part was adapted from Chang Hyon's Big Bird Reonstruction script------------
## SPE 32 reconstruction
##############################################################################################
    # Services to do Gulliver reconstruction
    tray.AddService("I3SimpleParametrizationFactory","SimpleTrack_noDC",
                    StepX=20*I3Units.m,                                     # ! 20m step size
                    StepY=20*I3Units.m,                                     # ! 20m step size
                    StepZ=20*I3Units.m,                                     # ! 20m step size
                    StepZenith=0.1 * I3Units.radian,                        # ! 0.1 radian step size in zenith
                    StepAzimuth=0.2 * I3Units.radian,                       # ! 0.2 radian step size in azimuth
                    StepT= 0. ,                                             # Default
                    StepLinE= 0. ,                                          # Default
                    StepLogE= 0. ,                                          # Default
                    StepLinL= 0. ,                                             # Default
                    BoundsX= [ -2000 * I3Units.m,2000 * I3Units.m ] ,       # ! Set bounds to +-2000m
                    BoundsY= [ -2000 * I3Units.m,2000 * I3Units.m] ,        # ! Set bounds to +-2000m
                    BoundsZ= [ -2000 * I3Units.m,2000 * I3Units.m ] ,       # ! Set bounds to +-2000m
                    BoundsZenith= [ 0., 0. ] ,                              # Default
                    BoundsAzimuth= [ 0., 0. ] ,                             # Default
                    BoundsT= [ 0., 0. ] ,                                   # Default
                    BoundsLinE= [ 0., 0. ] ,                                   # Default
                    BoundsLinL= [ 0., 0. ] ,                                   # Default
                    )
    
    # Define the gulliver minimization sevice to use
    tray.AddService( "I3GulliverMinuitFactory", "Minuit_noDC",
                     Algorithm= "SIMPLEX",                                  # Default
                     Tolerance= 0.01,                                       # ! change to 0.01
                     MaxIterations= 10000,                                  # Default
                     MinuitPrintLevel= -2,                                  # Default
                     MinuitStrategy= 2,                                     # Default
                     FlatnessCheck= True,                                   # Default
                     )
    
    # Use convoluted pandel as the PDF for the likelihood
    tray.AddService("I3GulliverIPDFPandelFactory","Pandel_noDC",
                    InputReadout= Pulses + '_noDC',                               # ! Name of pulses to use
                    Likelihood= "SPE1st",                                  # Default
                    PEProb= "GaussConvoluted",                                  # Default
                    IceModel= 2,                                           # Default
                    IceFile= "",                                           # Default
                    AbsorptionLength= 98.0 * I3Units.m,                    # Default
                    JitterTime= 15.0 * I3Units.ns,                         # Default
                    NoiseProbability= 1.0*I3Units.hertz * 10.0*I3Units.ns,  # ! Added a little noise term
                    )

    # linefit seed service
    tray.AddService( "I3BasicSeedServiceFactory", "LFSeed_noDC",
                     FirstGuesses= ["CascadeLineFit_noDC_DP"],                   # ! Use spe2 and LF (in case spe2 fails) as seeds
                     InputReadout= Pulses + '_noDC',                                  # ! Use pulses for vertex correction
                     TimeShiftType= "TFirst",                               # ! Use TFirst for vertex correction
                     SpeedPolice= True,                                     # Default
                     MaxMeanTimeResidual= 1000.0 * I3Units.ns,              # Default
                     )
    # track fit
    tray.AddModule( "I3IterativeFitter", "SPEFit32_noDC_DP",
                    RandomService= "SOBOL",                                # Default
                    NIterations= 32,                                       # ! Nunmber of iterations
                    SeedService= "LFSeed_noDC",                               # ! Name of seed service
                    Parametrization= "SimpleTrack_noDC",                        # ! Name of track parametrization service
                    LogLikelihood= "Pandel_noDC",                               # ! Name of likelihood service
                    CosZenithRange= [ -1, 1 ],                             # Default
                    Minimizer= "Minuit_noDC",                                   # ! Name of minimizer service
                    OutputName = 'SPEFit32_noDC_DP',
                    #If = InIceCscd,
                    If = If
                    )

@icetray.traysegment
def SplineReco(tray, name, Pulses='SplitInIcePulses', If=lambda f: True):  

    TruncatedName = "SplineMPE" + "TruncatedEnergy" + "_SPICEMie"

    tray.AddService( "I3PhotonicsServiceFactory", "PhotonicsServiceMu",
        PhotonicsTopLevelDirectory  ="/cvmfs/icecube.opensciencegrid.org/data/photon-tables/SPICEMie/",
        DriverFileDirectory         ="/cvmfs/icecube.opensciencegrid.org/data/photon-tables/SPICEMie/driverfiles",
        PhotonicsLevel2DriverFile   = "mu_photorec.list", 
        PhotonicsTableSelection     = 2,
        ServiceName                 = "PhotonicsServiceMu" )
    
    tray.AddModule( "I3TruncatedEnergy", 
        RecoPulsesName          = Pulses,           
        RecoParticleName        = "SPEFit32_DP",
        ResultParticleName      = TruncatedName,   
        I3PhotonicsServiceName  = "PhotonicsServiceMu", 
        UseRDE                  = True,             
        If                      = If )
   

    timingSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_prob_z20a10_V2.fits'
    amplitudeSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_abs_z20a10_V2.fits'
    stochTimingSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfHighEStoch_mie_prob_z20a10.fits'
    stochAmplitudeSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfHighEStoch_mie_abs_z20a10.fits'
    EnEstis = ["SplineMPETruncatedEnergy_SPICEMie_AllDOMS_Muon",
        "SplineMPETruncatedEnergy_SPICEMie_DOMS_Muon",
        "SplineMPETruncatedEnergy_SPICEMie_AllBINS_Muon",
        "SplineMPETruncatedEnergy_SPICEMie_BINS_Muon"]
    
    tray.AddSegment(spline_reco.SplineMPE, "SplineMPErecommended",
        configuration="recommended", PulsesName=Pulses,
        #TrackSeedList=["MPEFit_TT"], BareMuTimingSpline=timingSplinePath, # we don't have in test data
        TrackSeedList=["SPEFit32_DP"], BareMuTimingSpline=timingSplinePath,
        BareMuAmplitudeSpline=amplitudeSplinePath, EnergyEstimators=EnEstis)

@icetray.traysegment
def Monopod(tray, name, Pulses='OfflinePulses', If=lambda f: True):
    
    def forgeseed(fr):
        seed = dataclasses.I3Particle()
        seed.pos=fr['CascadeLlhVertexFit_DP'].pos
        seed.time=fr['CascadeLlhVertexFit_DP'].time
        seed.dir=fr['SPEFit32_DP'].dir
        seed.shape=dataclasses.I3Particle.Cascade
        seed.type=dataclasses.I3Particle.EMinus
        seed.fit_status=dataclasses.I3Particle.OK
        fr['CascadeSeed']=seed
    tray.AddModule(forgeseed, 'cscd-seed')

    cascade_service_spicemie = photonics_service.I3PhotoSplineService('/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.abs.fits', '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.prob.fits', 0)
    
    #add missing time window

    tray.AddModule(icecube.wavedeform.AddMissingTimeWindow, "add-timewindow", Pulses=Pulses, WaveformTimeRange="CalibratedWaveformRange")
    
    #def time_window(frame):
    #    print "Pulse start time: ", frame['OfflinePulsesTimeRange'].start
    #    print "Pulse stop time: ", frame['OfflinePulsesTimeRange'].stop
    #tray.AddModule(time_window, 'time-win')


    tray.AddSegment(millipede.MonopodFit, 'MonopodAngularGuess', 
                    #Pulses='SplitOfflinePulses',
                    Pulses=Pulses,
                    Iterations=1, 
                    Parametrization="SIMPLE",
                    #Seed="CascadeLLH",
                    Seed="CascadeSeed",
                    CascadePhotonicsService=cascade_service_spicemie,
                    PhotonsPerBin=15, 
                    #ExcludedDOMs=['CalibrationErrata', 'SaturatedDOMs', 'BadDomsListSLC', 'InIceErrata'],
                    ExcludedDOMs=['OfflineInIceCalibrationErrata', 'SaturatedDOMs', 'BadDomsList', 'BadDomsListSLC', 'InIceErrata'],
                    PartialExclusion=False)
    print 'first attempt complete'

    tray.AddSegment(millipede.MonopodFit, 'MonopodAngular',
                    #Pulses='SplitOfflinePulses',
                    Pulses=Pulses,
                    #Iterations=4,
                    Iterations=32,
                    Parametrization="HalfSphere",
                    Seed="MonopodAngularGuess",
                    CascadePhotonicsService=cascade_service_spicemie,
                    PhotonsPerBin=15,
                    #ExcludedDOMs=['CalibrationErrata', 'SaturatedDOMs', 'BadDomsListSLC', 'InIceErrata'],
                    ExcludedDOMs=['OfflineInIceCalibrationErrata', 'SaturatedDOMs', 'BadDomsList', 'BadDomsListSLC', 'InIceErrata'],
                    PartialExclusion=False)
    
    return
