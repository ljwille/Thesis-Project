#!/usr/bin/env python
from I3Tray import *
from icecube import icetray, dataio, dataclasses, DomTools
#import sys, numpy

load('VHESelfVeto')

#Find charge deposited using IceTray tools.
@icetray.traysegment
def CalQTot(tray, name, pulses='', If=True):
	
	tray.AddModule('HomogenizedQTot', 'qtot_total', Pulses=pulses)
	#tray.AddModule(lambda fr: fr['QTot'].value > 1000, 'qtotcut')
	tray.AddModule('I3LCPulseCleaning', 'cleaning', OutputHLC='HLCPulses',
		       OutputSLC='', Input=pulses, If=lambda frame: 'HLCPulses' not in frame)
	tray.AddModule('VHESelfVeto', 'selfveto', TimeWindow=3000, VertexThreshold=250, DustLayer=-160, DustLayerWidth=60, VetoThreshold=3, Pulses='HLCPulses')
	tray.AddModule('HomogenizedQTot', 'qtot_causal', Pulses=pulses,
		       Output='CausalQTot', VertexTime='VHESelfVetoVertexTime')
	
	return
