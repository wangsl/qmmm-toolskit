#!/usr/bin/python

import os, string
from HarmonicWindow import HarmonicWindow
from utility import die, is_blank_line

class HarmonicUmbrellaWindows:
    def __init__(self, file_name,
                 data_begin = 0, data_end = 2**31,
                 temperature = 298.15,
                 bin_minimum = -200.0, bin_maximum = -100.0, number_of_bins = 100) :
        if not os.path.exists(file_name) :
            die("Data file '" + file_name + "' does not exist")
        self.windows = []
        fin = open(file_name, "r")
        while 1 :
            line = fin.readline()
            if not line : break
            if not is_blank_line(line) :
                tmp = line.split()
                fname = tmp[0]
                minimum = string.atof(tmp[1])
                force_constant = string.atof(tmp[2])
                self.windows.append(HarmonicWindow(file_name = fname,
                                                   data_begin = data_begin,
                                                   data_end = data_end,
                                                   temperature = temperature,
                                                   bin_minimum = bin_minimum,
                                                   bin_maximum = bin_maximum,
                                                   number_of_bins = number_of_bins,
                                                   force_constant = force_constant,
                                                   minimum = minimum))
        return

    def number_of_windows(self) :
        return len(self.windows)

    def read_window_data(self, every = 1) :
        for w in self.windows :
            w.read_data_from_file(every = every)
        return

    def show_windows(self) :
        dashes = '-'
        for i in range(55) :
            dashes += '-'
        n = 0
        for w in self.windows :
            print dashes
            print("Window number: %4i" % (n))
            n += 1
            print w
        print dashes
        print dashes
        return

    def calculate_biased_distribution(self) :
        for w in self.windows :
            w.calculate_biased_distribution()
        return

    def show_biased_distribution(self) :
        for w in self.windows :
            w.show_biased_distribution()
        return

    def update_free_energy_shift(self, fshift) :
        assert(len(fshift) == self.number_of_windows())
        for i in range(self.number_of_windows()) :
            self.windows[i].update_free_energy_shift(fshift[i])
        return

    def write_biased_distributions_to_files(self) :
        n = 0
        for w in self.windows :
            name = "000" + str(n)
            n = n+1
            end = len(name)
            begin = end-3
            fname = "Window." + name[begin:end] + ".rho"
            w.write_biased_distribution_to_file(fname)
        return

    def set_data_begin_to_window(self, window_index, data_begin) :
        assert(window_index >= 0 and window_index < len(self.windows))
        assert(data_begin > 0)
        self.windows[window_index].set_data_begin(data_begin)
        return

    def set_data_end_to_window(self, window_index, data_end) :
        assert(window_index >= 0 and window_index < len(self.windows))
        assert(data_end > 0)
        self.windows[window_index].set_data_end(data_end)
        return
    
               

    
    
    

        
    
                              
