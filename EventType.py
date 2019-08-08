#------------------------------------------------------------------------------------
#Determine Event Type and neutrino information for studying classification performance.
#-------------------------------------------------------------------------------------
from I3Tray import *
from icecube import icetray, dataclasses, dataio

@icetray.traysegment
def EventType(tray, name):
    def IsNu(frame):
        #Finding neutrino that interacted in simulated data.
        if frame.Has("I3MCTree"):
            mct = frame["I3MCTree"]
            w_map=frame["I3MCWeightDict"]
            intT = w_map["InteractionType"]
            
            for pt in mct:
                if (pt.is_neutrino and pt.location_type == pt.InIce and len(mct.children(pt)) > 0):
                    p = pt
                    break
            try:
                p
            except NameError:
                return False
            frame['IntPos'] = mct.children(p)[0].pos

                        
            if p.pdg_encoding == 16 or p.pdg_encoding == -16:
                #Tagging tau neutrino interaction type for sub-classification and study of classification algorithm performance.
                if intT == 1:
                    frame["EventType"] = icetray.I3Int(1)
                    for pc in mct.children(p):
                        if pc.pdg_encoding == 15 or pc.pdg_encoding == -15:
                            frame["TauLength"] = dataclasses.I3Double(pc.length)
                            TauDecay = 0
                            for pgc in mct.children(pc):
                                if pgc.pdg_encoding == 11 or pgc.pdg_encoding == -11:
                                    TauDecay = 1
                                if pgc.type_string == 'Hadrons':
                                    TauDecay = 1
                                if pgc.pdg_encoding == 13 or pgc.pdg_encoding == -13:
                                    TauDecay = 0
                            frame["TauDecay"] = dataclasses.I3Double(TauDecay)
                elif intT == 2:
                    frame["EventType"] = icetray.I3Int(2)
                else:
                    print "Unknown Event Type"
                    print intT
                    return False
                return True
                #print "energy", p.energy
            elif p.pdg_encoding == 12 or p.pdg_encoding == -12:
                #Tagging background neutrino event.
                if intT == 1:
                    frame["EventType"] = icetray.I3Int(3)
                elif intT == 2:
                    frame["EventType"] = icetray.I3Int(4)
                elif intT == 3:
                    frame["EventType"] = icetray.I3Int(7)
                else:
                    print "Unknown Event Type"
                    print intT
                    return False
                return True
            elif p.pdg_encoding == 14 or p.pdg_encoding == -14:
                #Tagging background neutrino event.
                if intT == 1:
                    frame["EventType"] = icetray.I3Int(5)
                elif intT == 2:
                    frame["EventType"] = icetray.I3Int(6)
                else:
                    print "Unknown Event Type"
                    print intT
                    return False
                return True
            else:
                return False
        else:
            print "There is NO MCTree"
            return False
            
    tray.AddModule( IsNu, "TauCC_Check", Streams = [icetray.I3Frame.Physics, icetray.I3Frame.DAQ])
    
    return