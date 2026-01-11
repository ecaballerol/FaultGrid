# -*- coding: utf-8 -*-
'''
This script:
- Load indivual faults from the CFM New Zealand
- Compare with the GCMT location and select the faults closer to a thhreshold
- Use the a priori GCMT and write a GCMT file for each file
- In option : writes the final preseismic offset for data and synthetics.
'''

# Import externals
from shutil import copyfile
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from obspy import UTCDateTime
import copy
import pyproj as pp
import glob
import time
from os.path import join
import pickle
from multiprocessing import Pool

# Import csi
import csi.seismiclocations as sl
import csi.TriangularPatches as triangleFault
import csi.geodeticplot as geoplt
import cmt
from Arguments import *
from calcFaultDist import calcFaultSel,wrapperfaultSel

Vp = 5.7
Vs = 3.2

#----------------------------------------------------
# Build SeismicLocation objects From WCMT Solutions
#----------------------------------------------------
# Initialize SeismicLocation Object
seismo = sl('Earthquakes',utmzone=None,lon0=173.0,lat0=0,ellps="GRS80")
# Correct for New Zealand NZTM2000
seismo.utm = pp.CRS("EPSG:2193")
seismo.proj2utm = pp.Transformer.from_crs(seismo.wgs, seismo.utm, always_xy=True)
seismo.proj2wgs = pp.Transformer.from_crs(seismo.utm, seismo.wgs, always_xy=True)
# Read the correspondig CMT file
seismo.read_CMTSolutions(i_cmt_file)

Mw = seismo.mag
RLD = 10 ** (-2.44 + 0.59 * Mw)
dist_fault = RLD / 2
dist_fault = 1.1 *dist_fault
# dist_fault = 50
Est_Area = 10 ** (-3.49 + 0.91 * Mw)
Thrs_Area= Est_Area*0.2
#----------------------------------------------------
# Build Faults list from the NW CommunityFaultModel
#----------------------------------------------------
start = time.time()

faultnames = sorted(glob.glob('Tsurf_individual_v2/*.ts'))
seismotmp = copy.deepcopy(seismo)

if __name__=='__main__':
        
    inputs = []
    for ifault in faultnames[-350:]:
        inputs.append((ifault,seismotmp,dist_fault))

    ncores = 5
    pool = Pool(ncores)
    o = pool.map(wrapperfaultSel ,inputs)

    faults = []
    for ifault in o:
        if len (ifault) != 0:
            faults.append(ifault[0])

    seismo.distance2fault(faults,dist_fault)
    print('Faults selected')
    print(seismo.distancetofaults)

    end = time.time()
    print("Total running time: " + str(end - start))

    #Initiating CMT object
    filecmt = cmt.cmt()
    filecmt.read(i_cmt_file)
    filecmt.set_hypofromPDE()
    nodal1, nodal2 = filecmt.nodalplanes()

    #+++++++++++++++++++++++++++++++++++++++
    #Selecting faults by area and orientation
    selec_fault = {}
    numfaults = 0
    for ifault in faults:
        center_tmp = np.array(ifault.getcenters())
        flt_centroid = np.mean(center_tmp,axis=0)
        distances = np.linalg.norm(center_tmp - flt_centroid, axis=1)
        best_idx = np.argmin(distances)
        closest_patch = ifault.patch[best_idx]
        xy = ifault.getcenter(int(best_idx))
        lon,lat = ifault.xy2ll(xy[0],xy[1])
        *_, strike, dip = ifault.getpatchgeometry(int(best_idx))
        xhy,yhy = ifault.ll2xy(filecmt.hypo_lon,filecmt.hypo_lat)
        dist_hypo = np.sqrt( (xhy - xy[0])**2\
                        +(yhy - xy[1])**2\
                        +(filecmt.hypo_dep - xy[2])**2)
        ifault.computeArea()
        area = np.array(ifault.area).sum()
        strike = np.rad2deg(strike)
        dip = np.rad2deg(dip)
         #To be in accordance with Chamberlain:
        if ifault.name =='The_Humps':
            strike +=180
        if area > Thrs_Area[0]:
            if (strike > nodal1[0]-105) and (strike < nodal1[0]+105):
                selec_fault[ifault.name] = {}
                selec_fault[ifault.name]['center'] = [lon,lat,xy[2]]
                
                selec_fault[ifault.name]['sd'] = [strike,dip]
                selec_fault[ifault.name]['Area'] = area
                # selec_fault[ifault.name]['time'] = [dist_hypo/Vs-Vs*2,dist_hypo/Vs+Vs*2]
                selec_fault[ifault.name]['time'] = [np.round(dist_hypo/(Vs*5/4)),np.round(dist_hypo/(Vs*4/5)+Vs)]
                numfaults += 1

    for idx, (ifault, val) in enumerate(sorted(selec_fault.items(),key=lambda item: min(abs(t) for t in item[1]['time']))):
    # for idx, ifault in enumerate(selec_fault):
        filetmp = copy.deepcopy(filecmt)
        lon,lat,dep = selec_fault[ifault]['center']
        filetmp.lon, filetmp.lat, filetmp.dep = lon,lat,dep
        if selec_fault[ifault]['sd'][0]<180:
            if selec_fault[ifault]['sd'][1]<50:
                selec_fault[ifault]['rake'] = 110
            else:
                selec_fault[ifault]['rake'] = 135
        elif selec_fault[ifault]['sd'][0]>180:
            selec_fault[ifault]['rake'] = 135

        filetmp.sdr2MT(selec_fault[ifault]['sd'][0],selec_fault[ifault]['sd'][1],selec_fault[ifault]['rake'])
        filetmp.ts = 0 #Cause we will change it directly in the inversion code
        filetmp.pdeline = filetmp.pdeline.replace('SOUTH ISLAND',ifault)
        idxfmt = f"{idx:02d}"
        filetmp.write(join(GCMT_dir,i_cmt_file+'_'+idxfmt))

    with open('saved_dictionary.pkl', 'wb') as f:
        pickle.dump(selec_fault, f)
        print('Saving selected fault parameters')