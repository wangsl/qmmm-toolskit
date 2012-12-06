#!/usr/bin/python

from time import localtime, asctime
import sys, re, os
from os import path

element_symbol = ["X", "H", "He",
                  "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                  "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
                  "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",
                  "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
                  "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
                  "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
                  "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
                  "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
                  "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
                  "Am", "Cm", "Bk", "Cf", "Es", "Fm"]

blank_line_pattern = re.compile(r'^\s?$')

debye = 4.8033324

def die(words) :
    print words
    sys.exit()
    
def is_blank_line(line) :
    if blank_line_pattern.match(line) :
        return True
    else :
        return None

def remove_file(f) :
    if os.path.exists(f) :
        os.remove(f)
    return

def current_time() :
    return asctime(localtime())

def check_file_exists(file_name) :
    if(not path.exists(file_name)) :
        die("File '" + file_name + "' does not exit")
    return

def check_directory_exists(dir_name) :
    if(not path.exists(dir_name)) :
        die("Directory '" + dir_name + "' does not exit")
    return




