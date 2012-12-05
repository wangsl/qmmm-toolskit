#!/usr/bin/python

"""
Note WHAM iteration by Shenglong Wang Oct. 6, 2005
Department of Chemistry, New York University
Reference: 
(1) Marc Souaille and Benoit Roux  Comput. Phys. Commun. 135, 40 (2001)
(2) Benoit Roux  Comput. Phys. Commun. 91, 275 (1995)
(3) Erik M. Boczko and Charles L. Brooks, III  J. Phys. Chem. 97, 4509 (1993)
(4) Alan Grossfield  http://dasher.wustl.edu/alan/wham/
(5) Shankar Kumar et.al  J. Comput. Chem. 13, 1011 (1992)
"""

# I8 for Fortran compiler and -DI64 for C compiler
# from numpy import int64 as intpy
# No I8 for Fortran compiler and -DI64 for C compiler
from numpy import int32 as intpy

from HarmonicUmbrellaWindows import HarmonicUmbrellaWindows
from numpy import zeros, float64, resize
from whamit import wham_iteration as wham_iteration_fortran
from whamit import cpu_numbers
from math import log
from utility import kB, is_blank_line
import string, sys, os

class WHAM :
    def __init__(self, file_name,
                 data_begin = 0, data_end = 2**31,
                 every = 1,
                 temperature = 298.15,
                 bin_minimum = -200.0, bin_maximum = -100.0,
                 number_of_bins = 100) :

        assert(type(file_name) == type(""))
        assert(type(data_begin) == type(0))
        assert(type(data_end) == type(0))
        assert(type(every) == type(0))
        assert(type(temperature) == type(0.0))
        assert(type(bin_minimum) == type(0.0))
        assert(type(bin_maximum) == type(0.0))
        assert(type(number_of_bins) == type(0))        
        
        self.every = every
        self.__have_read_data = None
        self.temperature = temperature
        self.bin_minimum = bin_minimum
        self.bin_maximum = bin_maximum
        self.number_of_bins = number_of_bins
        self.maximum_iterations = 2000
        self.tolerance = 1.0e-5
        self.have_combined_biased_probability = None
        self.iteration_results_output_file_name = "ExpBetaF.out"
        self.umbrella_windows = HarmonicUmbrellaWindows(file_name = file_name,
                                                        data_begin = data_begin,
                                                        data_end = data_end,
                                                        temperature = temperature,
                                                        bin_minimum = bin_minimum,
                                                        bin_maximum = bin_maximum,
                                                        number_of_bins = number_of_bins)
        return
        
    def set_maximum_iterations(self, max_iter) :
        self.maximum_iterations = max_iter
        return

    def set_iteration_output_filename(self, filename) :
        self.iteration_results_output_file_name = filename
        return

    def set_tolerance(self, tol) :
        self.tolerance = tol
        return

    def windows(self) :
        return self.umbrella_windows.windows

    def number_of_windows(self) :
        return self.umbrella_windows.number_of_windows()

    def read_window_data(self) :
        if not self.__have_read_data :
            print "*** WHAM is reading window data ***"
            self.umbrella_windows.read_window_data(every = self.every)
            self.__have_read_data = True
            self.check_window_data()
        return

    def check_window_data(self) :
        n = 0
        for w in self.windows() :
            if not w.number_of_data() :
                print "\n\n ***** Window " + str(n) + " has 0 data *****\n"
                sys.exit()
            n += 1
        return

    def show_windows(self) :
        self.umbrella_windows.show_windows()
        return

    def set_free_energy_shift(self, fs) :
        assert(len(fs) == self.number_of_windows())
        self.umbrella_windows.update_free_energy_shift(fs)
        self.have_combined_biased_probability = None
        return

    def set_exp_beta_free_energy_shift(self, ebf) :
        assert(len(ebf) == self.number_of_windows())
        kT = kB*self.temperature
        free_energy_shift = []
        for i in range(self.number_of_windows()) :
            free_energy_shift.append(kT*log(ebf[i]))
        self.set_free_energy_shift(free_energy_shift)
        return

    def set_exp_beta_free_energy_shift_from_file(self, fname) :
        fin = open(fname, "r")
        ebf = []
        while 1 :
            line = fin.readline()
            if not line : break
            if not is_blank_line(line) :
                ebf.append(string.atof(line))
        fin.close()
        for i in range(len(ebf), self.number_of_windows()+1) :
            ebf.append(1.0)
        assert(len(ebf) >= self.number_of_windows())
        self.set_exp_beta_free_energy_shift(ebf[:self.number_of_windows()])
        return

    def combine_biased_probability_from_all_windows(self) :
        if self.have_combined_biased_probability : return
        print "WHAM is combining the biased probability ... \n"
        
        windows = self.windows()

        factor = []
        for w in windows :
            n = w.number_of_data()
            ebf = w.exp_beta_free_energy_shift()
            factor.append(n*ebf)

        self.umbrella_windows.calculate_biased_distribution()
        
        self.rho = []
        dx = (self.bin_maximum-self.bin_minimum)/self.number_of_bins
        for k in range(self.number_of_bins) :
            x = (k+0.5)*dx + self.bin_minimum
            bottom = 0.0
            for j in range(self.number_of_windows()) :
                bottom += factor[j] * windows[j].exp_beta_umb_potential(x)
            top = 0.0
            for i in range(self.number_of_windows()) :
                d = windows[i].biased_distribution[k]
                assert(x == d[0])
                top += windows[i].number_of_data() * d[1]
            self.rho.append([x, top/bottom])

        self.have_combined_biased_probability = True
        return

    def free_energy(self) :
        self.combine_biased_probability_from_all_windows()
        kT = kB*self.temperature
        pmf = []
        for x, p in self.rho :
            if p <= 0.0 :
                w = 100.0
            else :
                w = -kT*log(p)
            pmf.append([x, w])
        return pmf

    def potential_of_mean_force(self) :
        return self.free_energy()

    """
    def set_OMP_NUM_THREADS(self, OMP_NUM_THREADS = 2) :
        os.environ['OMP_NUM_THREADS'] = str(OMP_NUM_THREADS)
        return
    """
    
    def set_threads_number(self, threads_number = 2) :
        self.set_OMP_NUM_THREADS(threads_number)
        return

    def set_OMP_NUM_THREADS(self, OMP_NUM_THREADS = 0) :
        ncpu = OMP_NUM_THREADS
        if OMP_NUM_THREADS == 0 :
            ncpu = cpu_numbers()
        os.environ['OMP_NUM_THREADS'] = str(ncpu)
        return
    
    def run_WHAM_iterations_from_fortran(self) :
        self.read_window_data()
        Nwind = self.number_of_windows()
        nt = zeros(Nwind, intpy)
        kBT = zeros(Nwind, float64)
        ebf = zeros(Nwind, float64)

        i = 0
        for w in self.windows() :
            nt[i] = w.number_of_data()
            w_ebf = w.exp_beta_free_energy_shift()
            ebf[i] = w_ebf
            kBT[i] = 1.0/w.beta
            i += 1

        print "\n*** WHAM is preparing the data for iteration ***"

        if not 'OMP_NUM_THREADS' in os.environ.keys() :
            self.set_OMP_NUM_THREADS()

        windows = self.windows()
        max_nt = max(nt)
        ebw = zeros([Nwind, max_nt, Nwind], float64)
        for i in range(Nwind) :
            for l in range(nt[i]) :
                ril = windows[i].data[l]
                for k in range(Nwind) :
                    ebw[i][l][k] = windows[k].exp_beta_umb_potential(ril)

        print "\n*** WHAM begins to run the iteration to get\n" \
              "the free energy shift for every umbrella window ***\n"

        ebf2 = wham_iteration_fortran(ebw, self.maximum_iterations, self.tolerance,
                                      nt, kBT, ebf, Nwind, max_nt,
                                      self.iteration_results_output_file_name)
        
        resize(ebw, (1,1,1))
        assert(len(ebf2) == self.number_of_windows())
        
        free_energy_shift = []
        for i in range(self.number_of_windows()) :
            free_energy_shift.append(-kBT[i]*log(ebf2[i]))
        self.umbrella_windows.update_free_energy_shift(free_energy_shift)

        print "*** Free energy shift after WHAM iterations ***"
        i = 0
        for f in free_energy_shift :
            print("Window index: %4d   Free energy shift: %14.8f kcal/mol" % (i, f))
            i += 1
        print
        return

    def run_WHAM_iterations(self) :
        self.run_WHAM_iterations_from_fortran()
        return

    def write_biased_distributions_to_files(self) :
        self.umbrella_windows.write_biased_distributions_to_files()
        return

    def set_data_begin_to_window(self, window_index, data_begin) :
        self.umbrella_windows.set_data_begin_to_window(window_index, data_begin)
        return

    def set_data_end_to_window(self, window_index, data_end) :
        self.umbrella_windows.set_data_end_to_window(window, data_end)
        return
    
    
