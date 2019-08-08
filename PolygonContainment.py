#-------------------------------------------------------------------------------------------------
# Calculate perpendicular distances from an event x-y vertex to the detector polygon side surfaces.
# The shortest distance to one of the surfaces can be used as a containment cut variable
# Author: Donglian Xu   Date: 04/06/2014
# ------------------------------------------------------------------------------------------------

import icecube
from I3Tray import *
from icecube import dataclasses, icetray, dataio, tableio
from icecube.dataclasses import I3Particle, I3Double, I3OMGeo, I3Position
from icecube.icetray import OMKey

import os, sys, math
import numpy as np

@icetray.traysegment
def PolygonContainment(tray, name, RecoVertex="MonopodAngular", geometry='ic86',outputname = ''):
    
    def point_in_polygon(gx, gy, rx, ry):
        # this utility is borrowed from Jakob VS' IC79 level3 cuts:
        # http://code.icecube.wisc.edu/projects/icecube/browser/sandbox/jvansanten/CascadeL3_IC79/candidates/V12-04-00_rc2/python/level3/polygon.py
        
        rx = np.asarray(rx)
        ry = np.asarray(ry)
        gx = np.asarray(gx); gx=gx.reshape((gx.size, 1))
        gy = np.asarray(gy); gy=gy.reshape((gy.size, 1))
        
        # Coordindates of the "next" vertex. Since the polygon is closed, the
        # last vertex is next to the first.
        nx = np.roll(gx, 1, axis=0)
        ny = np.roll(gy, 1, axis=0)
        
        # Draw a horizontal line at y. Which edges does it intersect?
        crossings = ((gy <= ry)&(ry < ny))|((ny <= ry)&(ry < gy))
        
        # Now, cast a ray to the right. Which edges does it intersect?
        crossings &= (rx < (nx-gx)*(ry-gy)/(ny-gy) + gx)
        
        # Count the number of crossings.
        inside = (crossings.sum(axis=0) % 2) != 0
        if inside.size == 1:
            return inside[0]
        else:
            return inside

    def cal_distance(x1, y1, x2, y2):
        return math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))

    def least_distance_to_polygon(fr):
        if fr.Has("I3Geometry"):
            geo_map=fr["I3Geometry"]
        else:
            print "ERROR: No I3Geometry in frame!"
            exit(1)
        
        #String numbers that defines the polygon. String 72 is repeated so that the polygon is closed when connecting those strings in order.
        if geometry=='ic86':
            poly_strings=[72, 74, 50, 6, 1, 31, 75, 78]
        if geometry=='ic79': 
            poly_strings=[72, 74, 50, 6, 2, 41, 75, 78]
        # read in x,y coordinates from the gcd file for the polygon strings
        om=OMKey()
        gx=[]
        gy=[]
        for i in range(len(poly_strings)):
            om.string=poly_strings[i]
            om.om=60
            gx += [geo_map.omgeo[om].position.x]
            gy += [geo_map.omgeo[om].position.y]
        
        #print "polygon x coordinates: ", gx
        #print "polygon y coordinates: ", gy
        
        # pull the x,y coordinates of string 78, 72, and 74 for later usage
        x72=gx[0]
        y72=gy[0]
        x74=gx[1]
        y74=gy[1]
        x78=gx[7]
        y78=gy[7]
            
        # x,y coordinates for event vertex 
        rx=fr[RecoVertex].x
        ry=fr[RecoVertex].y
                 
        inside=point_in_polygon(gx, gy, rx, ry)
        
        if inside:
            print "Event vertex is inside the IC86 polygon. Now calculate the perpendicular distances to the polygon edges.." 
            str_dist=[]
            event_dist=[]
            for i in range(1, len(poly_strings)):
                str_dist += [cal_distance(gx[i-1], gy[i-1], gx[i], gy[i])]
                event_dist += [cal_distance(rx, ry, gx[i-1], gy[i-1])]
            
            str_dist += [cal_distance(gx[-1], gy[-1], gx[0], gy[0])] # add the last polygon edge to close the polygon
            event_dist += [cal_distance(rx, ry, gx[-1], gy[-1])] # add the last polygon edge to close the polygon
            #print "str_dist: ", str_dist
            #print "event_dist: ", event_dist
            #event_to polystring_dist.pop() #remove the last repeated distances between event vertex and string 31
           
            """
             The perpendicular distances from an event vertex to one polygon edge is the height of the triangle formed by the event vertex and
             the two string which defines the polygon edge. According to Heron's theorem: area=sqrt(s(s-a)(s-b)(s-c)), semiperimeter s=1/2(a+b+c)
             Also, area=1/2*base*height. Hence, height=2*area/base=2*sqrt(s(s-a)(s-b)(s-c))/base. Here, the bases are the polygon edges.
            """
            
            triad_heights=[]
            semiperimeters=[]
            for i in range(1, len(str_dist)):
                semiperimeters += [0.5*(str_dist[i-1]+event_dist[i-1]+event_dist[i])]
            semiperimeters += [0.5*(str_dist[-1]+event_dist[-1]+event_dist[0])] # add the last polygon edge element to close the polygon
            #print "semiperimeter: ", semiperimeters

            for i in range(1, len(str_dist)):
                if semiperimeters[i-1]*(semiperimeters[i-1]-event_dist[i-1])*(semiperimeters[i-1]-str_dist[i-1])*(semiperimeters[i-1]-event_dist[i]) < 0:
                    triad_heights += [0]
                    continue
                triad_heights += [2*np.sqrt(semiperimeters[i-1]*(semiperimeters[i-1]-event_dist[i-1])*(semiperimeters[i-1]-str_dist[i-1])*(semiperimeters[i-1]-event_dist[i]))/str_dist[i-1]]
            if semiperimeters[-1]*(semiperimeters[-1]-event_dist[-1])*(semiperimeters[-1]-str_dist[-1])*(semiperimeters[-1]-event_dist[0]) > 0:
                triad_heights += [2*np.sqrt(semiperimeters[-1]*(semiperimeters[-1]-event_dist[-1])*(semiperimeters[-1]-str_dist[-1])*(semiperimeters[-1]-event_dist[0]))/str_dist[-1]] # add the last polygon edge element to close the polygon 
            else: 
                triad_heights += [0]
           
            """
            Now divide the polygon into three regions: region A: perpendicular distance to edge 78-72 is inside of the polygon; region B: perpendicular distances to both edges 78-72 and 72-74 are inside of the polygon; 
            region C: perpendicular distance to edge 72-74 is inside of the polygon. These three regions are formed by two perpendicular lines through string 72, w.r.t. line 78-72 and line 72-74 respectively.
            The two lines deviding these regions are defined as: l_{AB} := -(x78-x72)/(y78-y72)(x-x72)+y72; l_{BC} := -(x74-x72)/(y74-y72)(x-x72)+y72
            For points fall within region A, distance to edge 72-74 should not be considered as minimum distance to polygon edges. Likewise, points that fall within region B, both distances to edges 78-72 and 72-74 should not
            be considered, and for region C, distance to edge 78-72 should not be considered. 
            
            """
            # CAUTION: the following manipulation is dependent on the ordering of the triad_heights elements: distance to edge 72-74 is first element, while distance to edge 72-78 is last element.
            if (ry-y72+(x78-x72)/(y78-y72)*(rx-x72))>=0: # point in region A
                #print "point in region A!"
                del triad_heights[0] #remove distance to edge 72-74
                
            elif ((ry-y72+(x78-x72)/(y78-y72)*(rx-x72))<0) and ((ry-y72+(x74-x72)/(y74-y72)*(rx-x72))<0): # point in region B 
                #print "point in region B!"
                del triad_heights[0]
                del triad_heights[-1]
                
            elif (ry-y72+(x74-x72)/(y74-y72)*(rx-x72))>=0: # point in region C
                #print "point in region C!"
                del triad_heights[-1]
                
            #print "Perpendicular distances to polygon edges: ", triad_heights
            print "Minimum perpendicular distance to polygon edges: ", min(triad_heights)
            print "Minimum distance to polygon corner strings: ", min(event_dist)
            
            least_dist_to_polygon=min(triad_heights)
            #if min(event_dist)<min(triad_heights):
            #    least_dist_to_polygon=min(event_dist)

            fr["LeastDistanceToPolygon"+outputname] = dataclasses.I3Double(least_dist_to_polygon)
            fr["LeastDistanceToCornerString"+outputname] = dataclasses.I3Double(min(event_dist))
 
        else:
            print "Event vertex is outside of the IC86 polygon.."
            fr["LeastDistanceToPolygon"+outputname]=dataclasses.I3Double(-1)
            fr["LeastDistanceToCornerString"+outputname] = dataclasses.I3Double(-1)
            
    tray.AddModule(least_distance_to_polygon, "cal-dist")
    
    return
