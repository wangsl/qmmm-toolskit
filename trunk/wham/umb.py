#!/usr/bin/python

import string
from os import unlink, path, getcwd, system, listdir
from glob import glob
from sys import exit

def delete_distance_files() :
    for f in glob("distance.*.dat") :
        unlink(f)
    return

def delete_rho_files() :
    for d in glob("Window.*.rho") :
        unlink(d)
    return

def plot_reaction_coordinates() :
    dat_files = glob("distance.*.dat")
    dat_files.sort()
    cmd = "xmgr "
    for f in dat_files :
        cmd += "-block " + f + " -bxy 1:4 "
    cmd += " &"
    system(cmd)
    return



    

