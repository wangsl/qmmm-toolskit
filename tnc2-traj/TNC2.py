#!/bin/env python

from sys import stdout
from utility import check_file_exists
from QMMMMDtnc2 import *

class MDTNC2 :

    def __init__(self, file_name) :
        check_file_exists(file_name)
        self.n_atoms = 0
        self.n_snapshots = 0
        self.atom_types = []
        self.atom_names = []
        create_read_tnc2_wrapper(file_name)
        self.__set_n_atoms()
        self.__set_n_snapshots()
        self.__set_atom_types()
        self.__set_atom_names()
        return

    def delete(self) :
        delete_read_tnc2_wrapper()
        return

    def __set_n_atoms(self) :
        self.n_atoms = n_atoms_wrapper()
        return

    def __set_n_snapshots(self) :
        self.n_snapshots = snapshot_size_wrapper()
        return

    def __set_atom_types(self) :
        self.atom_types = atom_types_wrapper()
        return

    def __set_atom_names(self) :
        self.atom_names = []
        names_from_tnc2 = atom_names_wrapper()
        iatom = 0
        while iatom < self.n_atoms :
            name_length = len(names_from_tnc2[iatom])
            name = ""
            j = 0
            while j < name_length :
                if names_from_tnc2[iatom][j] != ' ' :
                    name += names_from_tnc2[iatom][j]
                    j += 1
                else :
                    j = name_length
                    
            self.atom_names.append(name)
            iatom += 1
            
        assert len(self.atom_names) == self.n_atoms
        return
    
    def coordinates(self, index) :
        return snapshot_wrapper(index)

    def active_atom_list(self) :
        return active_atom_list_wrapper()

    def qm_esp_charges_saved(self) :
        return qm_esp_charges_saved_wrapper()

    def active_atom_gradients_saved(self) :
        return gradients_saved_wrapper()

    def active_atom_velocities_saved(self) :
        return velocities_saved_wrapper()

    def qm_esp_charges(self) :
        return qm_esp_charges_wrapper()

    def active_atom_gradients(self) :
        return gradients_wrapper()

    def active_atom_velocities(self) :
        return velocities_wrapper()

    def write_qm_esp_charges(self, fout = stdout) :
        fout.write(" %6d\n" % len(self.qm_esp_charges()))
        i = 1
        for e in self.qm_esp_charges() :
            fout.write(" %6d%16.10f\n" % (i, e))
            i += 1
        return

    def write_active_gradients(self, fout = stdout) :
        use = self.active_atom_list()
        g = self.active_atom_gradients()
        assert len(use) == len(g)
        fout.write(" %6d\n" % len(use))
        i = 0
        while i < len(use) :
            fout.write(" %6d%16.10f%16.10f%16.10f\n" % (use[i], g[i][0], g[i][1], g[i][2]))
            i += 1
        return

    def write_active_velocities(self, fout = stdout) :
        use = self.active_atom_list()
        v = self.active_atom_velocities()
        assert len(use) == len(v)
        fout.write(" %6d\n" % len(use))
        i = 0
        while i < len(use) :
            fout.write(" %6d%16.10f%16.10f%16.10f\n" % (use[i], v[i][0], v[i][1], v[i][2]))
            i += 1
        return



