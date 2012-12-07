#!/bin/env python

# $Id$

from Cartesian import Cartesian

class TinkerAtom :

    def __init__(self,
                 index = None, name = None, coordinate = None,
                 type = None, bond_list = []) :
        self.index = index
        self.name = name
        self.coordinate = coordinate
        self.type = type
        self.bond_list = bond_list
        return

    def __repr__(self) :
        atom = ""
        assert self.index
        atom += "%6d" % (self.index)
        assert self.name
        atom += "  %-3s" % self.name
        assert self.coordinate and len(self.coordinate) == 3
        atom += "  %16.10f%16.10f%16.10f" % (self.coordinate.x(),
                                             self.coordinate.y(),
                                             self.coordinate.z())
        assert self.type
        atom += "  %6i" % (self.type)
        if self.bond_list and len(self.bond_list) :
            for b in self.bond_list :
                atom += " %6i" % (b)
        return atom

    def set_coordinate(self, coordinate) :
        assert isinstance(coordinate, Cartesian)
        self.coordinate = coordinate
        return
    
if __name__ == "__main__" :

    atom = TinkerAtom(index = 123, name = "Cl-", coordinate = Cartesian(0.1223232, 0.22222, 0.2),
                      type = 960, bond_list = [1, 2, 3, 404])

    print atom

    atom.set_coordinate(Cartesian(-1.22, 2.2, 3.4))
    print atom
                 
