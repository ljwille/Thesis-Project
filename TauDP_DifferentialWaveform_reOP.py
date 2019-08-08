#---------------------------------------------------
# Calculate the derivatives of ATWD waveforms. Those
# derivatives should show more explicitly the trend
# a waveform goes--ascending or descending.
# Author: Logan Wille and Donglian Xu      March 2017
#---------------------------------------------------
from I3Tray import *
from icecube import icetray, dataclasses, dataio, tableio
from icecube.icetray import I3Units
import numpy as np
import math

#@icetray.traysegment
def CalWaveformDerivatives(tray, name, 
                           WfQtot=10000,
                           bins_p1=1,
                           bins_trailing=2,
                           bins_p2=3,
                           Amp1LC=1,
                           Amp2LC=12,
                           AmpTrailingLC=-1,
                           Amp1SD = 10,
                           Amp2SD = 18,
                           AmpTrailingSD = -17,
                           der_step = 4,
                           DiscardEvents=False):
    
    def WaveformDerivatives(frame):
        #frame["DifferentialWaveforms"] = dataclasses.I3MapKeyVectorDouble()
        frame["DoublePulseWaveforms"] = dataclasses.I3WaveformSeriesMap()
        #frame["DoublePulseSeriesMap"] = dataclasses.I3RecoPulseSeriesMap()
        frame["BinsToT1"] = dataclasses.I3VectorDouble()# bins with time over threshold (ToT) for the first pulse of a derivative waveform
        frame["BinsToT2"] = dataclasses.I3VectorDouble()# bins with time over threshold (ToT) for the second pulse of a derivative waveform
        frame["BinsTbT"] = dataclasses.I3VectorDouble()# bins with time below threshold (TbT) for the first trailing edge of a derivative waveform
        frame["Amp1"] = dataclasses.I3VectorDouble()# cumulative derivative amplitude of time over threshold (ToT) for the first pulse
        frame["Amp2"] = dataclasses.I3VectorDouble()# cumulative derivative amplitude of time over threshold (ToT) for the second pulse
        frame["AmpTrailing"] = dataclasses.I3VectorDouble()# cumulative derivative amplitude of time below threshold (ToT) for the trailing edge
        frame["DPWaveformQTot"] = dataclasses.I3VectorDouble()
        frame["DP_T1"] = dataclasses.I3VectorDouble()
        frame["DP_T2"] = dataclasses.I3VectorDouble()
        frame["DPWaveformPulse1Amp"] = dataclasses.I3VectorDouble()
        frame["DPWaveformPulse2Amp"] = dataclasses.I3VectorDouble()
        frame["DP_OMs"] = dataclasses.I3VectorOMKey()

        SDtag = False
        LCtag = False

        if frame.Has("CalibratedWaveformsHLCATWD"):

            wf_map=frame["CalibratedWaveformsHLCATWD"]
            dp_count=0 #count number of dp waveforms from this event
            dp_om_keys=[]
            for om, wf_series in wf_map:
                for wf in wf_series:
                    wf_vect=wf.waveform
                    #print "sum of wf_vect:", sum(wf_vect)
                    #print "sum of wf_vect [mV]:", sum(wf_vect)/I3Units.mV
                    wf_status=wf.status
                    wf_binwidth=wf.bin_width
                    start_time=wf.time
                    if wf_status == 0: #0:VIRGINAL, 2:COMBINED, 4:SATURATED, 8:UNDERSHOT
                        wf_qtot=0.0
                        #print wf_series
                        for i in range(128): 
                            wf_vect[i] = wf_vect[i]/I3Units.mV #transform wf unit to mV
                            wf_qtot += wf_vect[i] #individual wf integrated charge in mV  
                        #Calculate Sliding Time Window (STW) derivatives for waveforms
                        stw_vect = []
                        stepsize=2
                        for i in range(1, 127):
                            stw_vect += [(wf_vect[i+1]-wf_vect[i-1])/(wf_binwidth*stepsize)]
                        der_wf_bins=len(stw_vect)
                        #print "STW_wf_vector length: ", der_wf_bins

                        der_flag=0
                        for i in range(1, der_wf_bins-6): # 3, 4, 5, 6 bins tried. With 6 bins required, all the manually-selected nice double pulse waveforms picked up by the algorithm.   
                            if (stw_vect[i-1]>0 and stw_vect[i]>0 and stw_vect[i+1]>0 and stw_vect[i+2]>0 and stw_vect[i+3]>0 and stw_vect[i+4]>0):
                                der_flag=i-1
                                #print "waveform is constantly rising from bin %i"%(der_flag)
                                break
                        #sum up the amplitude of the rising edge
                        amp_rising=0.
                        for i1 in range(der_flag, der_flag+6):
                            amp_rising += stw_vect[i1]

                        wf_diff=[]
                        for i in range(der_flag, 128-der_step, der_step):
                            wf_diff += [(wf_vect[i+der_step] - wf_vect[i])/(wf_binwidth*der_step)]
                            #print "Derivative wf: ", wf_diff
                        diff_wf_bins=len(wf_diff)
                        #print "Derivative waveform vector length: ", diff_wf_bins

                        #search for the first and second pulses!
                        pulse1=False
                        pulse2=False
                        trailing=False
                        i_flag=0
                        j_flag=0
                        k_flag=0
                        amp_p1=0   #amplitude of first rising edge in the derivative waveform
                        amp_p2=0   #amplitude of second rising edge in the derivative waveform
                        amp_trailing=0  #amplitude of first trailing edge in the derivative waveform

                        binsToT1=0
                        binsToT2=0
                        binsTbT=0
                        j1_flag=0
                        j2_flag=0
                        j3_flag=0
                        bins_p21=0
                        bins_p22=0
                        dp_t1=0                      
                        dp_t2=0
                        dp_amp1=0 # amplitude of first pulse
                        dp_amp2=0 # amplitude of second pulse
                        #first pulse search 
                        for i in range(diff_wf_bins - bins_p2 - bins_trailing):
                            #if wf_diff[i]>0: #and wf_diff[i+1]>0:
                            if np.prod([True if wf_diff[n]>0.0 else False for n in range(i,i+bins_p1)]):
                                pulse1=True
                                i_flag=i
                                #print "First pulse starts at bin %i"%(i_flag)
                                #print "wf_diff[%i]: "%(i_flag), wf_diff[i_flag]
                                break

                        #trailing edge search
                        if pulse1:
                            # firstly move the index out of the positive region
                            for k in range(i_flag+bins_p1, diff_wf_bins-bins_trailing):#
                                if np.prod([True if wf_diff[n] < 0.0 else False for n in range(k,k+bins_trailing)]):
                                    k_flag=k
                                    trailing=True
                                    l_flag=k_flag+bins_trailing
                                    #print "First pulse ends at bin %i "%(k_flag)                SHOULDN'T THIS SECTION BE THE ONE THE DETERMINES TRAILING EDGE?
                                    #print "wf_diff[%i]: "%(k_flag), wf_diff[k_flag]
                                    break
                                    dp_t1=start_time+der_flag*422.4/128 #beginning time of the first pulse

                        #second pulse search
                        if pulse1 and trailing:
                            for j in range(l_flag, (diff_wf_bins - bins_p2)):
                                #if wf_diff[j]>0: # and wf_diff[j+1]>0 and wf_diff[j+2]>0:
                                if np.prod([True if wf_diff[n]>0.0 else False for n in range(j,j+bins_p2)]):
                                    pulse2=True
                                    j_flag=j
                                    #print "Second pulse starts at bin %i"%(j_flag)
                                    #print "wf_diff[%i]: "%(j_flag), wf_diff[j_flag]
                                    break

                            for q in range(der_flag, der_flag+j_flag*der_step): # transform the derivative waveform index back to the waveform vector basis
                                dp_amp1 += wf_vect[q] #calculate amplitude of first pulse

                            #move the index out of the positive region, and do another round of big pulse search.
                            #If the later pulse is bigger than the earlier one, point the index to the bigger one.
                            #Having more than (including) two pulses with ToT of 3 derivetive waveform bins, which is
                            #duration of 3*4*3.3 ns, is quite rare. Therefore, another round of searching should be sufficient enough.
                            if j_flag and j_flag<(diff_wf_bins-bins_p2):
                                for j1 in range(j_flag, (diff_wf_bins - bins_p2)):
                                    if wf_diff[j1]<0:
                                        j1_flag=j1
                                        break
                            if j1_flag and j1_flag<(diff_wf_bins-bins_p2):
                                for j2 in range(j1_flag, (diff_wf_bins - bins_p2)):
                                    if np.prod([True if wf_diff[n]>0.0 else False for n in range(j2,j2+bins_p2)]):
                                        j2_flag=j2
                                        break
                                #move the index out of the positive region
                            if j2_flag and j2_flag<(diff_wf_bins-bins_p2):
                                j3_flag=0
                                for g in range(j2_flag+bins_p2, diff_wf_bins):
                                    if wf_diff[g]<0:
                                        j3_flag=g
                                        #print "Second pulse ends at bin %i "%(h_flag)
                                        break
                                if j3_flag == 0:
                                    j3_flag=diff_wf_bins
                                bins_p21=j1_flag-j_flag
                                bins_p22=j3_flag-j2_flag
                                pulse3=True
                                if bins_p22>bins_p21:
                                    j_flag=j2_flag

                            dp_t2=start_time+(der_flag+j_flag*der_step)*422.4/128 #starting time of second pulse
                            for q in range(der_flag+j_flag*der_step, 128): # transform the derivative waveform index back to the waveform vector basis
                                dp_amp2 += wf_vect[q] #calculate amplitude of second pulse 

                        if pulse1 and trailing and pulse2:
                            # calculate the cumulative derivative amplitude for potential first pulse
                            for s in range(i_flag, k_flag):
                                amp_p1 += wf_diff[s]
                                #print "Derivative amplitude p1: ", amp_p1
                            binsToT1=k_flag-i_flag
                            #print "p1 ToT bins: ", binsToT1

                            #move the index out of the first trailing edge region
                            r_flag=0
                            for r in range(l_flag, (diff_wf_bins - bins_p2)):
                                if wf_diff[r]>0:
                                    r_flag=r
                                    #print "The first trailing edge ends at bin %i"%(r_flag)
                                    break
                            for q in range(k_flag, r_flag):
                                amp_trailing += wf_diff[q]
                                #print "Derivative amplitude trailing edge: ", amp_trailing
                            binsTbT=r_flag-k_flag
                            #print "Trailing edge TbT bins: ", binsTbT

                            # move the index out of the second positive region
                            h_flag=0
                            for h in range(j_flag+bins_p2, diff_wf_bins):
                                if wf_diff[h]<0:
                                    h_flag=h
                                    #print "Second pulse ends at bin %i "%(h_flag)
                                    break
                            if h_flag == 0:
                                h_flag=diff_wf_bins
                                 #print "Second pulse continue till end of waveform..."
                            # calculate the cumulative derivative amplitude for second pulse
                            for t in range(j_flag, h_flag):
                                amp_p2 += wf_diff[t]
                                #print "Derivative amplitude p2: ", amp_p2
                            binsToT2=h_flag-j_flag

                            #print "wf_qtot: ", wf_qtot
                            #print "BinsToT1: ", binsToT1
                            #print "BinsTbT: ", binsTbT
                            #print "BinsToT2: ", binsToT2
                            #print "Amp1: ", amp_p1
                            #print "AmpTrailing: ", amp_trailing
                            #print "Amp2: ", amp_p2
                            #print "============="



                        #if wf_qtot>WfQtot and binsToT1>=BinsToT1 and binsTbT>=BinsTbT and binsToT2>=BinsToT2 and amp_p1>Amp1 and amp_p2>Amp2 and amp_trailing<AmpTrailing: # if this waveform has double peak feature.
                        if wf_qtot>WfQtot and pulse1 and trailing and pulse2 and amp_p1>Amp1LC and amp_p2>Amp2LC and amp_trailing<AmpTrailingLC: # if this waveform has double peak feature.
                            if amp_p1>Amp1SD and amp_p2>Amp2SD and amp_trailing<AmpTrailingSD and binsToT1 > 1: #if this waveform passes single dom double peak:
                                SDtag = True
                                #print "Found Single DoM double pulse

                            #print "caught one double pulse waveform! Hooray------------------------------------------------------------------"
                            dp_count += 1
                            dp_om_keys.append(om)
                            
                            frame["Amp1"].append(amp_p1)
                            frame["Amp2"].append(amp_p2)
                            frame["AmpTrailing"].append(amp_trailing)
                            frame["BinsToT1"].append(binsToT1)
                            frame["BinsTbT"].append(binsTbT)
                            frame["BinsToT2"].append(binsToT2)
                            frame["DPWaveformQTot"].append(wf_qtot)
                            frame["DP_OMs"].append(om)
                            frame["DoublePulseWaveforms"].update({om:wf_series})
                            frame["DPWaveformPulse1Amp"].append(dp_amp1)
                            frame["DPWaveformPulse2Amp"].append(dp_amp2)
                            frame["DP_T1"].append(dp_t1)
                            frame["DP_T2"].append(dp_t2)

                    #end if
                #end for 
            #end for

           
            for i in dp_om_keys:
                for j in dp_om_keys:
                    if i[0] == j[0]:
                        if 0 < np.abs(np.int(j[1]) - np.int(i[1])) < 3:
                            LCtag = True
                            #print "Found local coincident DPs"
                            break
                if LCtag:
                    break
                    
            frame["DPFlagSD"] = icetray.I3Int(SDtag)
            frame["DPFlagLC"] = icetray.I3Int(LCtag)


        else:
            #print "no CalibratedWaveformsHLCATWD in frame"
            return False



        #Keeping the double pulse events only
        if DiscardEvents:
            if LCtag or SDtag:
                return True
            return False

        #    tray.AddModule( DP_bg, "background_dp", Streams = [icetray.I3Frame.Physics])
    tray.AddModule( WaveformDerivatives, "differential_wf", Streams = [icetray.I3Frame.Physics])


    return
    
