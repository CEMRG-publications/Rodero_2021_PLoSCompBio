#!/usr/bin/python

import os
from datetime import date
from glob import glob

from carputils import settings
from carputils import tools
from carputils import ep
from carputils import model
from carputils.resources import petsc_block_options

import sys


def scripthelp():
    print('Script to solve a laplacian solution from a vertex to a surface.')
    print('Mandatory arguments: --path, --imsh, --source, --ground, --ID.')
    print('Use the help of each argument for more information.')

def parser():
    parser = tools.standard_parser()
    parser.add_argument('--path',
                        default='',
                        help='Path to the mesh.')
    parser.add_argument('--imsh',
                    default='',
                    help='Name of the mesh.')
    parser.add_argument('--source',
                        default='',
                        help='Vtx file from where the activation starts.')
    parser.add_argument('--ground',
                        default='',
                        help='Vtx file where simulation finishes.')
    parser.add_argument('--simulation_name',
                        default='',
                        help='Name of the simulation.')
    return parser


def jobID(args):
    """
    Generate name of top level output directory.
    """
   

    return '{}/{}'.format(args.path,args.simulation_name)

@tools.carpexample(parser, jobID)
def run(args, job):

    meshname   = '{}/{}'.format(args.path,args.imsh)

        # Get basic command line, including solver options
    cmd = tools.carp_cmd()

    cmd += ['-experiment',   2, # perform Laplace solve only
           '-bidomain',     1,
           '-simID',        job.ID,
           '-meshname',     '{}/{}'.format(args.path,args.imsh)]

    stimuli = ['-num_stim', 2,
               '-stimulus[0].stimtype', 3, # extracellular ground
               '-stimulus[0].vtx_file','{}/{}'.format(args.path,args.ground),
               '-stimulus[1].stimtype', 2, #extracellular voltage
               '-stimulus[1].vtx_file','{}/{}'.format(args.path,args.source),
               '-stimulus[1].start',   0,
               '-stimulus[1].duration', 1,
               '-stimulus[1].strength', 1]

    # set isotropic conductivities everywhere
    cond = ['-num_gregions'   , 1,
            '-gregion[0].g_il', 1,
            '-gregion[0].g_it', 1,
            '-gregion[0].g_el', 1,
            '-gregion[0].g_et', 1]

    cmd += stimuli
    cmd += cond

    # Run simulation
    job.carp(cmd)


    # Create .dat file
    phie = os.path.join(args.path, args.simulation_name, 'phie.igb')
    dat  = os.path.join(args.path, args.simulation_name, 'phie.dat')
    cmd  = [settings.execs.IGBEXTRACT, phie,
           '-o', 'ascii_1pL',
           '-O', dat]

    job.bash(cmd,None)

    cmd = [settings.execs.GLVTKCONVERT,
               '-m', '{}/{}'.format(args.path,args.imsh),
               '-n', dat,
               '-o', '{}/{}_{}'.format(args.path,args.imsh,args.simulation_name),
               '-F', 'bin']


    job.bash(cmd, None)

if __name__ == '__main__':
    if sys.argv[1] == '-h':
        scripthelp()
    else:
        run()

