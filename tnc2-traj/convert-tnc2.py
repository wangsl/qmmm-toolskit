#!/usr/bin/python

from sys import argv, exit, stdout
from utility import check_file_exists
from TNC2TinkerXYZ import MDTNC2TinkerXYZ

def convert_tnc2_to_Tinker_xyz(Tinker_xyz_file, tnc2_file,
                               output_file = None) :
    check_file_exists(Tinker_xyz_file)
    check_file_exists(tnc2_file)

    tnc2 = MDTNC2TinkerXYZ(Tinker_xyz_file = Tinker_xyz_file,
                           TNC2_file = tnc2_file)

    output = "tnc2-Tinker"
    if output_file :
        output = output_file

    isnapshot = 0
    while isnapshot < tnc2.n_snapshots :
        print " snapshot index: %i" % isnapshot

        fout_xyz = open("%s-%d.xyz" % (output, isnapshot+1), "w")
        tnc2.read_snapshot_and_write_as_Tinker_xyz(isnapshot, fout_xyz)
        fout_xyz.close()

        if tnc2.TNC2.qm_esp_charges_saved() :
            fout_esp = open("%s-%d.esp" % (output, isnapshot+1), "w")
            tnc2.TNC2.write_qm_esp_charges(fout_esp)
            fout_esp.close()

        if tnc2.TNC2.active_atom_gradients_saved() :
            fout_grads = open("%s-%d.grads" % (output, isnapshot+1), "w")
            tnc2.TNC2.write_active_gradients(fout_grads)
            fout_grads.close()

        isnapshot += 1

    return

if __name__ == "__main__" :

    Tinker_xyz_file = None
    tnc2_file = None
    output_file = None

    if len(argv) == 3 :
        Tinker_xyz_file = argv[1]
        tnc2_file = argv[2]
    elif len(argv) == 4 :
        Tinker_xyz_file = argv[1]
        tnc2_file = argv[2]
        output_file = argv[3]
    else :
        print "usage: %s Tiner_xyz_template tnc2_file output_file" % argv[0]
        exit()

    if output_file and output_file[-4:] == ".xyz" :
        output_file = output_file[:-4]
        
    convert_tnc2_to_Tinker_xyz(Tinker_xyz_file = Tinker_xyz_file,
                               tnc2_file = tnc2_file,
                               output_file = output_file)

    exit()
