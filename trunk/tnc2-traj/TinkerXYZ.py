#!/bin/env python

from string import atoi, atof
from sys import stdout
from utility import check_file_exists, die
from TinkerAtom import TinkerAtom
from Cartesian import Cartesian

class TinkerXYZ :

    def __init__(self, file_name) :
        self.__file_name = file_name
        check_file_exists(file_name)
        self.n_atoms = 0
        self.n_qm_atoms = 0
        self.n_pseudo_atoms = 0
        self.n_zero_charge_atoms = 0
        self.n_mm_atoms = 0
        self.charge = 0
        self.spin = 1
        self.atoms = []
        self.__read_from_Tinker_xyz_file()
        return

    def __check_atom_numbers(self) :
        assert self.n_atoms == (self.n_qm_atoms + 
                                self.n_pseudo_atoms +
                                self.n_zero_charge_atoms +
                                self.n_mm_atoms)
        return

    def __read_from_Tinker_xyz_file(self) :
        self.n_atoms = 0
        self.n_qm_atoms = 0
        self.n_pseudo_atoms = 0
        self.n_zero_charge_atoms = 0
        self.n_mm_atoms = 0
        self.charge = 0
        self.spin = 1
        self.atoms = []

        fin = open(self.__file_name, "r")

        line = fin.readline()
        assert line
        tmp = line.split()

        if len(tmp) >= 7 :
            self.n_atoms = atoi(tmp[0])
            self.n_qm_atoms = atoi(tmp[1])
            self.n_pseudo_atoms = atoi(tmp[2])
            self.n_zero_charge_atoms = atoi(tmp[3])
            self.n_mm_atoms = atoi(tmp[4])
            self.charge = atoi(tmp[5])
            self.spin = atoi(tmp[6])
            self.__check_atom_numbers()
        elif len(tmp) >= 5 :
            self.n_atoms = atoi(tmp[0])
            self.n_qm_atoms = atoi(tmp[1])
            self.n_pseudo_atoms = atoi(tmp[2])
            self.n_zero_charge_atoms = atoi(tmp[3])
            self.n_mm_atoms = atoi(tmp[4])
            self.__check_atom_numbers()
        elif len(tmp) >= 1 :
            self.n_atoms = atoi(tmp[0])
            self.n_mm_atoms = self.n_atoms
            self.__check_atom_numbers()
        else :
            die(self.__file_name + ": first line error")

        i = 1
        while i <= self.n_atoms :
            line = fin.readline()
            assert line
            tmp = line.split()
            assert len(tmp) >= 6
            index = atoi(tmp[0])
            assert index == i
            name = tmp[1]
            x = atof(tmp[2])
            y = atof(tmp[3])
            z = atof(tmp[4])
            type = atoi(tmp[5])
            atom = TinkerAtom(index = index, name = name, coordinate = Cartesian(x, y, z),
                              type = type, bond_list = map(int, tmp[6:]))
            self.atoms.append(atom)
            i += 1
        fin.close()

        assert len(self.atoms) == self.n_atoms
        return

    def __print_atom_numbers(self, fout = stdout) :
        fout.write("%6i%6i%6i%6i%6i%6i%6i\n" % (self.n_atoms, self.n_qm_atoms,
                                                self.n_pseudo_atoms, self.n_zero_charge_atoms,
                                                self.n_mm_atoms, self.charge, self.spin))
        return
    
    def print_Tinker_xyz(self, fout = stdout) :
        assert len(self.atoms) == self.n_atoms
        self.__print_atom_numbers(fout)

        iatom = 0
        while iatom < self.n_atoms :
            fout.write("%s\n" % self.atoms[iatom])
            iatom += 1
        return

    def print_Tinker_xyz_with_coordinates(self, coordinates, fout = stdout) :
        assert len(coordinates) == self.n_atoms

        assert len(self.atoms) == self.n_atoms
        self.__print_atom_numbers(fout)

        iatom = 0
        while iatom < self.n_atoms :
            self.atoms[iatom].set_coordinate(Cartesian(coordinates[iatom][0],
                                                       coordinates[iatom][1],
                                                       coordinates[iatom][2]))
            fout.write("%s\n" % self.atoms[iatom])
            iatom += 1
        return
        
