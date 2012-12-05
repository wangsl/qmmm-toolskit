#!/usr/bin/python

"""
Note WHAM iteration by Shenglong Wang Oct. 6, 2005
Department of Chemistry, New York University
Reference: 
(1) Marc Souaille and Benoit Roux  Comput. Phys. Commun. 135, 40 (2001)
(2) Benoit Roux  Comput. Phys. Commun. 91, 275 (1995)
(3) Erik M. Boczko and Charles L. Brooks, III  J. Phys. Chem. 97, 4509 (1993)
(4) Alan  Grossfield  http://dasher.wustl.edu/alan/wham/
(5) Shankar Kumar et.al  J. Comput. Chem. 13, 1011 (1992)
"""

import string, os
from math import exp, pi, log
from utility import die, kB

class Window :
    def __init__(self, file_name,
                 data_begin = 0, data_end = 2**31,
                 temperature = 298.15,
                 bin_minimum = -200.0, bin_maximum = -100.0,
                 number_of_bins = 100) :
        if not os.path.exists(file_name) :
            die("Data file '" + file_name + "' does not exist")
        self.__file_name = file_name
        self.__data_begin = data_begin
        self.__data_end = data_end
        if(self.__data_end <= self.__data_begin) :
            die("For data file '%s', data_begin = %d, data_end = %d, *** error ***" %
                (self.__file_name, self.__data_begin, self.__data_end))
        self.__temperature = temperature
        self.__bin_minimum = bin_minimum
        self.__bin_maximum = bin_maximum
        self.__number_of_bins = number_of_bins
        self.beta = 1.0/(kB*temperature)
        self.free_energy_shift = 0.0
        self.data = []
        self.biased_distribution = []
        self.unbiased_distribution = []
        self.unbiased_potential_of_mean_force = []
        return

    def number_of_data(self) :
        return len(self.data)

    def file_name(self) :
        return self.__file_name

    def temperature(self) :
        return self.__temperature

    def bin_minimum(self) :
        return self.__bin_minimum

    def bin_maximum(self) :
        return self.__bin_maximum

    def number_of_bins(self) :
        return self.__number_of_bins
    
    def set_temperature(self, T) :
        self.__temperature = T
        self.beta = 1.0/(kB*T)
        return

    def __repr__(self) :
        s = "file_name: %s\n" % (self.__file_name)
        s += "data_begin: %d, " % (self.__data_begin)
        if(self.__data_end >= 2**30) :
            s += "data_end: None"
        else :
            s += "data_end: %d" % (self.__data_end)
        s += "\n"
        s += "number_of_data: %d\n" % (self.number_of_data())
        s += "bin_minimum: %.4f, bin_maximum: %.4f\n" % \
             (self.__bin_minimum, self.__bin_maximum)
        s += "number_of_bins: %d\n" % (self.__number_of_bins)
        s += "temperature: %.4f K, beta: %.4f (kcal/mol)^-1\n" % \
             (self.__temperature, self.beta)
        s += "free_energy_shift: %.5f kcal/mol\n" % (self.free_energy_shift)
        return s

    def read_data_from_file(self, every = 1) :
        fin = open(self.__file_name, "r")
        n = 0
        while 1 :
            n = n+1
            if n > self.__data_begin : break
            line = fin.readline()
            if not line : break

        if n < self.__data_begin :
            die("The total lines in file '" + self.__file_name +
                "' is less than " + str(self.__data_begin))

        #n = 0
        while 1 and n <= self.__data_end :
            line = fin.readline()
            if not line : break
            if not n%every :
                if line[-1] != '\n' : continue
                tmp = line.split()
                x = None
                if len(tmp) == 0 :
                    print "The data format in file '" + self.__file_name + \
                          "' is not correct, which should be 1 or 2 or 4 columns"
                    pass
                else :
                    x = string.atof(tmp[-1])
                if x > self.__bin_minimum and x < self.__bin_maximum :
                    self.data.append(x)
            n = n+1
        return

    def exp_beta_free_energy_shift(self) :
        beta = self.beta
        f = self.free_energy_shift
        return exp(beta*f)

    def umbrella_potential(self, coordinate) :
        print "\n*** umbrella_potential To be implemented in other class ***\n"
        return None

    def exp_beta_umb_potential(self, coordinate) :
        beta = self.beta
        p = self.umbrella_potential(coordinate)
        return exp(-beta*p)

    def update_free_energy_shift(self, f = 0.0) :
        self.free_energy_shift = f
        return

    def calculate_biased_distribution(self) :
        assert(self.__bin_maximum > self.__bin_minimum)
        if len(self.biased_distribution) : return
        dx = (self.__bin_maximum - self.__bin_minimum)/self.__number_of_bins
        for i in range(self.__number_of_bins) :
            x = (i+0.5)*dx + self.__bin_minimum
            self.biased_distribution.append([x, 0])

        ntot = 0
        for d in self.data :
            if d >= self.__bin_minimum and d <= self.__bin_maximum :
                n = int((d-self.__bin_minimum)/dx)
                self.biased_distribution[n][1] += 1
                ntot += 1
        assert(ntot == self.number_of_data())

        factor = dx*self.number_of_data()
        for d in self.biased_distribution :
            d[1] /= factor
        return
            
    def show_biased_distribution(self) :
        self.calculate_biased_distribution()
        for x, p in self.biased_distribution :
            print("%10.4f%14.6f" % (x, p))
        return

    def write_biased_distribution_to_file(self, fname) :
        fout = open(fname, "w")
        self.calculate_biased_distribution()
        for x, p in self.biased_distribution :
            fout.write("%10.4f%14.6f\n" % (x, p))
        fout.close()
        return

    def calculate_unbiased_distribution(self) :
        assert(self.__bin_maximum > self.__bin_minimum)
        if len(self.unbiased_distribution) : return
        self.calculate_biased_distribution()
        beta = self.beta
        factor = exp(-beta*self.free_energy_shift)
        for x, p in self.biased_distribution :
            V = self.umbrella_potential(x)
            p *= exp(beta*V)*factor
            self.unbiased_distribution.append([x, p])
        return

    def show_unbiased_distribution(self) :
        self.calculate_unbiased_distribution()
        for x, p in self.unbiased_distribution :
            print("%10.4f%14.6f" % (x, p))
        return

    def write_unbiased_distribution_to_file(self, fname) :
        fout = open(fname, "w")
        self.calculate_unbiased_distribution()
        for x, p in self.unbiased_distribution :
            fout.write("%10.4f%14.6f\n" % (x, p))
        fout.close()
        return

    def calculate_unbiased_potential_of_mean_force(self) :
        if len(self.unbiased_potential_of_mean_force) : return
        self.calculate_biased_distribution()
        kBT = 1.0/self.beta
        for x, p in self.biased_distribution :
            if p <= 0.0 :
                w = 100.0 - self.umbrella_potential(x)
            else :
                w = -kBT*log(p) - self.umbrella_potential(x)
            w += self.free_energy_shift
            self.unbiased_potential_of_mean_force.append([x, w])
        return
    
    def show_unbiased_potential_of_mean_force(self) :
        self.calculate_unbiased_potential_of_mean_force()
        for x, w in self.unbiased_potential_of_mean_force :
            print("%10.4f%14.6f" % (x, w))
        return

    def set_data_begin(self, data_begin) :
        assert(data_begin > 0)
        self.__data_begin = data_begin
        assert(self.__data_begin < self.__data_end)
        return

    def set_data_end(self, data_end) :
        assert(data_end > 0)
        self.__data_end = data_end
        assert(self.__data_begin < self.__data_end)
        return
    
    
