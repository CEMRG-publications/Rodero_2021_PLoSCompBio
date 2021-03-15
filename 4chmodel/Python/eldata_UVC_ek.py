#!/usr/bin/env python2

import os
import sys
import math
import subprocess
import copy
#import matplotlib as mpl
#import matplotlib.pyplot as mplplot
import numpy as np
import sets
import scipy.linalg as spla
import scipy.signal as spsig

import argparse



def read_pnts(filename):
    return np.loadtxt(filename, dtype=float, skiprows=1)


def write_pnts(filename, pts):
    assert len(pts.shape) == 2 and pts.shape[1] == 3
    with open(filename, 'w') as fp:
        fp.write('{}\n'.format(pts.shape[0]))
        for pnt in pts:
            fp.write('{0[0]}\t{0[1]}\t{0[2]}\n'.format(pnt))


def read_surface(filename):
    return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(1,2,3))


def write_surface(filename, surfs):
    assert len(surfs.shape) == 2 and surfs.shape[1] == 3
    with open(filename, 'w') as fp:
        fp.write('{}\n'.format(surfs.shape[0]))
        for tri in surfs:
            fp.write('Tr {0[0]}\t{0[1]}\t{0[2]}\n'.format(tri))


def read_neubc(filename):
   return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(0,1,2,3,4,5)) 


def read_elems(filename):
  return np.loadtxt(filename, dtype=int, skiprows=1, usecols=(1,2,3,4,5))


def vector_cprod(vec1, vec2):
  return np.array([vec1[1]*vec2[2]-vec1[2]*vec2[1],
                   vec1[2]*vec2[0]-vec1[0]*vec2[2],
                   vec1[0]*vec2[1]-vec1[1]*vec2[0]])

def read_dat(filename):
    return np.loadtxt(filename, dtype=float, skiprows=0)

def read_vtx(filename):
    return np.loadtxt(filename, dtype=int, skiprows=2)

def vector_sprod(vec1, vec2):
  return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]


def create_csys(vec):
    vec0 = None
    vec1 = None
    if (vec[0] < 0.5) and (vec[1] < 0.5):
        tmp = math.sqrt(vec[1]*vec[1]+vec[2]*vec[2])
        vec1 = np.array([0.0, -vec[2]/tmp, vec[1]/tmp])
        vec0 = vector_cprod(vec, vec1)
    else:
        tmp = math.sqrt(vec[0]*vec[0]+vec[1]*vec[1])
        vec1 = np.array([vec[1]/tmp, -vec[0]/tmp, 0.0])
        vec0 = vector_cprod(vec, vec1)
    return [vec0, vec1, vec]

def main(args):

    # subprocess.call('clear')

    outfile='PM_' + args.heart + '.dat'

    basedir = '/home/crg17/Desktop/Seg3DProjects/healthy_CT_MRI/h_case' + args.heart + '/meshing/1000um'
    basename = os.path.join(basedir, 'h_case' + args.heart)
    print('Reading points...\n')
    pnts = read_pnts(basename+'.pts')
    print('Reading elements...\n')    
    elem = read_elems(basename+'.elem')    

    datdir = basedir+'/BiV/UVC_PM/UVC'
    datname = os.path.join(datdir, 'COORDS_Z')
    print('Reading UVC...\n')
    UVCpts = read_dat(datname+'.dat')

    vtxdir = basedir+'/BiV'
    vtxname = os.path.join(vtxdir, 'BiV')
    print('Reading vtx...\n')
    vtxsubmsh = read_vtx(vtxname+'.vtx')   
    vtxsubmsh = np.array(vtxsubmsh, dtype=int)

    pcatags = [1, 2, 25, 26]
    UVCptsmsh = np.zeros(len(pnts[:,0]))
    print('Creating auxiliar variables...\n')
    for i, ind in enumerate(vtxsubmsh):
        UVCptsmsh[ind] = UVCpts[i]

    UVCelem = []
    print('Creating UVC elementwise...\n')
    for i, elm in enumerate(elem):
        if ( elem[i,4] in pcatags ):
            UVCelem.append((UVCptsmsh[elm[0]]+UVCptsmsh[elm[1]]+UVCptsmsh[elm[2]]+UVCptsmsh[elm[3]])*0.25)
        else:
            UVCelem.append(0.0)
    print('Saving file...\n')
    np.savetxt(os.path.join(basedir, 'UVC_elem.dat'), UVCelem, fmt='%.8f')

    # compute the data on the elements
    p1 = 1.5266
    p2 = -0.37
    p3 = 0.4964
    p4 = 0

    th = 0.82
    print('Creating penalty map...\n')
    elemdat = []
    for i, l in enumerate(UVCelem):
        if ( elem[i,4] in pcatags ):
            if (UVCelem[i] >= th):
                elemdat.append(0.0) 
            else: 
                x = UVCelem[i]
                x_m = th-x
                elemdat.append(p1*x_m**3 + p2*x_m**2 + p3*x_m + p4)
        else:
            elemdat.append(0.0)    
    print('Saving file...\n')
    np.savetxt(os.path.join(basedir,outfile), elemdat, fmt='%.8f') 

    print 'DONE'

    cmd = '/home/common/bin/GlVTKConvert -m '+basename+' -e '+os.path.join(basedir, outfile)+' -e '+os.path.join(basedir, 'UVC_elem.dat')+' -F bin -o '+basename+'_elem_dat_UVC'
    os.system(str(cmd))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Optional app description')
    parser.add_argument('heart',type=str)
    args = parser.parse_args()
    print('This is case ')
    print(args.heart)
    print('\n')
    main(args)
