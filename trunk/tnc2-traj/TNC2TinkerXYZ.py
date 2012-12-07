#!/bin/env python

# $Id$

from sys import stdout
from utility import check_file_exists, die
from TinkerXYZ import TinkerXYZ
from TNC2 import MDTNC2

class MDTNC2TinkerXYZ :
    
    def __init__(self, Tinker_xyz_file, TNC2_file) :
        check_file_exists(Tinker_xyz_file)
        self.Tinker_xyz = TinkerXYZ(Tinker_xyz_file)
        check_file_exists(TNC2_file)
        self.TNC2 = MDTNC2(TNC2_file)

        self.__check_atom_names()
        self.__check_atom_types()

        self.n_atoms = self.Tinker_xyz.n_atoms
        self.n_snapshots = self.TNC2.n_snapshots
        
        return

    def __check_atom_names(self) :
        assert self.Tinker_xyz.n_atoms == self.TNC2.n_atoms
        iatom = 0
        while iatom < self.Tinker_xyz.n_atoms :
            if self.Tinker_xyz.atoms[iatom].name != self.TNC2.atom_names[iatom] :
                die(" atom name error: %i %s %s" % (iatom, self.Tinker_xyz.atoms[iatom].name,
                                                    self.TNC2.atom_names[iatom]))
            iatom += 1
        return

    def __check_atom_types(self) :
        assert self.Tinker_xyz.n_atoms == self.TNC2.n_atoms
        iatom = 0
        while iatom < self.Tinker_xyz.n_atoms :
            if self.Tinker_xyz.atoms[iatom].type != self.TNC2.atom_types[iatom] :
                die(" atom name error: %i %i %i" % (iatom, self.Tinker_xyz.atoms[iatom].type,
                                                    self.TNC2.atom_types[iatom]))
            iatom += 1
        return

    def read_snapshot_and_write_as_Tinker_xyz(self, isnapshot, fout = stdout) :
        assert isnapshot < self.n_snapshots
        xyz = self.TNC2.coordinates(isnapshot)
        self.Tinker_xyz.print_Tinker_xyz_with_coordinates(xyz, fout)
        return
        

        

    
