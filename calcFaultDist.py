#!/usr/bin/env python

# Import externals

# Import Wavemod
import numpy as np
import wavemod as wm
from copy import deepcopy

import csi.seismiclocations as sl
import csi.TriangularPatches as triangleFault
import csi.geodeticplot as geoplt
import pyproj as pp

import cmt
from Arguments import *

def calcFaultSel(ifault,seisloc,dist_fault):

    # Assign fault parameters
    seismotmp = deepcopy(seisloc)
    fname = ifault.removesuffix('.ts').removeprefix('Tsurf_individual_v2/')
    # Initialize a fault from data
    fault = triangleFault(fname,utmzone=None,lon0=173.0,lat0=0,ellps="GRS80")
    # Correct for New Zealand NZTM2000
    fault.utm = pp.CRS("EPSG:2193")
    fault.proj2utm = pp.Transformer.from_crs(fault.wgs, fault.utm, always_xy=True)
    fault.proj2wgs = pp.Transformer.from_crs(fault.utm, fault.wgs, always_xy=True)

    # Get the fault patches from GoCAD
    fault.readGocadPatches(ifault,neg_depth=True,utm=True,factor_xy=1e-3,factor_depth=1e-3)
    fault.initializeslip(values='depth')
    fault.setTrace(delta_depth=1)

    try:
        seismotmp.distance2fault([fault],dist_fault)
        return fault
    except Exception as e:
        print(f"fault and earthquake too far ")
        return None


def wrapperfaultSel(input_args):
    ifault,seisloc,dist_fault = input_args
    ## Select model
    fault = calcFaultSel(ifault,seisloc,dist_fault)
    if fault is not None:
        return [fault]
    else:
        return []

