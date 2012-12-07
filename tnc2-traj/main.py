#!/bin/env python

# $Id$ 

from sys import exit, argv
from utility import *
from TNC2 import MDTNC2
from TinkerXYZ import TinkerXYZ

tnc_file = argv[1]
check_file_exists(tnc_file)

tinker_xyz_file = argv[2]
check_file_exists(tinker_xyz_file)

tnc2 = MDTNC2(tnc_file)
tinker_xyz = TinkerXYZ(tinker_xyz_file)

i = 0
while i < tnc2.n_snapshots :
    print " Snapshot index: %i" % i
    index = "00000" + str(i)
    index = index[-5:]

    xyz_coordinates = tnc2.coordinates(i)

    xyz_file_name = index + ".xyz"
    xyz_out = open(xyz_file_name, "w")
    tinker_xyz.print_Tinker_xyz_with_coordinates(xyz_coordinates, xyz_out)
    xyz_out.close()

    if tnc2.qm_esp_charges_saved() :
        fout = open(index + ".chg", "w")
        tnc2.write_qm_esp_charges(fout)
        fout.close()

    if tnc2.active_atom_gradients_saved() :
        fout = open(index + ".grad", "w")
        tnc2.write_active_gradients(fout)
        fout.close()

    if tnc2.active_atom_velocities_saved() :
        fout = open(index + ".vel", "w")
        tnc2.write_active_velocities(fout)
        fout.close()
        
    i += 1    

tnc2.delete()

    
