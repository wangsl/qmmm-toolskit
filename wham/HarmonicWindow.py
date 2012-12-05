#!/usr/bin/python

import sys
from Window import Window
from math import pi

class HarmonicWindow(Window) :
    def __init__(self, file_name,
                 data_begin = 0,
                 data_end = 2**31,
                 temperature = 298.15,
                 bin_minimum = -200.0, bin_maximum = -100.0,
                 number_of_bins = 100,
                 force_constant = 0.0, minimum = 0.0) :
        Window.__init__(self, file_name = file_name,
                        data_begin = data_begin, data_end = data_end,
                        temperature = temperature,
                        bin_minimum = bin_minimum, bin_maximum = bin_maximum,
                        number_of_bins = number_of_bins)
        self.__force_constant = force_constant
        self.__minimum = minimum
        return

    def umbrella_potential(self, coordinate) :
        dx = coordinate - self.__minimum
        w = self.__force_constant * dx*dx
        return w

    """
    def calculate_biased_distribution(self) :
        bin_maximum = self.bin_maximum()
        bin_minimum = self.bin_minimum()
        number_of_bins = self.number_of_bins()
        assert(bin_maximum > bin_minimum)
        if len(self.biased_distribution) : return
        dx = (bin_maximum - bin_minimum)/number_of_bins
        for i in range(number_of_bins) :
            x = (i+0.5)*dx + bin_minimum
            self.biased_distribution.append([x, 0])

        ntot = 0
        for d in self.data :
            if d > bin_minimum and d < bin_maximum :
                n = int((d-bin_minimum)/dx)
                self.biased_distribution[n][1] += 1
                ntot += 1
        assert(ntot == self.number_of_data())

        factor = 4*pi*dx*self.number_of_data()
        for d in self.biased_distribution :
            d[1] /= factor*d[0]*d[0]
    """


        
        
                 
                 
