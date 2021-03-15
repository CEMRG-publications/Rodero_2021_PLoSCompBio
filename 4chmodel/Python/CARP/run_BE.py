#!/usr/bin/env python2

import math, os

import numpy as np
import pandas as pd
from carputils import tools

import sys


# def parser():
#     parser = tools.standard_parser()
#     parser.add_argument('--biv_path',
#                         default='',
#                         help='Path to the biventricular mesh.')
#     parser.add_argument('--UVC_path',
#                         default='',
#                         help='Path to the UVC files.')
#     return parser
# def jobID(args):
#     """
#     Generate name of top level output directory.
#     """

#     return ''
def run(current_case="0"):

    # args = parser()

    meshbase = "/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case" + current_case + "/meshing/1000um/BiV/BiV"
    lap_apba = "/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case" + current_case + "/meshing/1000um/BiV/UVC/COORDS_Z.dat"
    lap_circ = "/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case" + current_case + "/meshing/1000um/BiV/UVC/COORDS_PHI.dat"
    # meshbase = "/home/crg17/Desktop/scripts/4chmodel/testUVC/mesh/its_bad/BiV"
    # lap_apba = "/home/crg17/Desktop/scripts/4chmodel/testUVC/UVC_its_bad_minimal/UVC/COORDS_Z.dat"
    # lap_circ = "/home/crg17/Desktop/scripts/4chmodel/testUVC/UVC_its_bad_minimal/UVC/COORDS_PHI.dat"

  # read mesh elements and points
    print('Loading mesh...\n')
    if os.path.isfile(meshbase+'.elem'):
        # print("hey")
        elems = pd.read_csv(meshbase+'.elem', skiprows=1, delim_whitespace=True, usecols=(1, 2, 3, 4, 5), dtype=int, header=None).values.squeeze()
        pnts = pd.read_csv(meshbase+'.pts', skiprows=1, delim_whitespace=True, dtype=float, header=None).values.squeeze()
    else: 
        raiseIOError('Could not load the mesh!')
    print('Reading data....\n')
  # read laplace solutions in apical-basal-direction and phi-direction
    apba = pd.read_csv(lap_apba, delim_whitespace=True, dtype=float, header=None).values.squeeze()
    circ = pd.read_csv(lap_circ, delim_whitespace=True, dtype=float, header=None).values.squeeze()

  # normalize phi-directional laplace solution
    circmin, circmax = np.amin(circ), np.amax(circ)
    circ = (circ - circmin)/(circmax - circmin)

  # define bullseye sections and indices
    apex   = {(-1.000, 2.000): 17}
    apical = {(-1.000, 0.125): 16, (0.125, 0.375): 15, (0.375, 0.635): 14, (0.635, 0.875): 13, (0.875, 2.000): 16}
    midcav = {(-1.000, 0.166): 11, (0.166, 0.333): 10, (0.333, 0.500):  9, (0.500, 0.666):  8, (0.666, 0.833):  7, (0.833, 2.000): 12}
    basal  = {(-1.000, 0.166):  5, (0.166, 0.333):  4, (0.333, 0.500):  3, (0.500, 0.666):  2, (0.666, 0.833):  1, (0.833, 2.000):  6}
    aorta  = {(-1.000, 2.000): -1}
    idxmap = {0:apex, 1: apex, 2: apical, 3: midcav, 4: basal, 5: aorta}
  # retag the elements

    # apba = (apba / 0.75)*(apba <= 0.75) + (apba > 0.75)*1.
    
    print('Tagging...\n')
    tagidx = []  
    for i, elem in enumerate(elems):
        # print(elem)
        if (elem[4] == 1) and (np.mean(np.take(apba, elem[0:4]))>=0):
            apbaval = int(math.ceil(np.mean(np.take(apba, elem[0:4]))*5.0))
            circvals = np.take(circ, elem[0:4])
            circval = np.mean(circvals)
            if (np.amin(circvals) < 0.1) and (np.amax(circvals) > 0.9):
                circval += 0.5*(-1.0 if circval < 0.5 else 1.0)    
            sections = idxmap.get(apbaval, None)
            if sections is None:
                tagidx.append(-1)        
            else:
                for rng, idx in sections.iteritems():
                    if rng[0] < circval <= rng[1]:
                        tagidx.append(idx)
                        break
                else:
                    tagidx.append(-1)        
        else:
            tagidx.append(-1)

    tagidx = np.array(tagidx, dtype=int)
    print('Writing...\n')
  # write new element data
    with open(meshbase+'.retag.elem', 'w') as fp:
        fp.write('{}\n'.format(elems.shape[0]))
        for i, elem in enumerate(elems):
            fp.write('Tt {0[0]} {0[1]} {0[2]} {0[3]} {1}\n'.format(elem, tagidx[i]))

def scripthelp():
    print('Script to create a bullseye map (AHA) on a left ventricle, once the UVC has been created. The argument is the case number with two digits.')


if __name__ == '__main__':
    current_case = sys.argv[1]
    if current_case == '-h':
        scripthelp()
    else:
       run(current_case=current_case)
