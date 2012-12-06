#!/bin/csh -f

if(`uname -m` == "ppc64") then
    make -j4 -f Makefile-OpenMP-PPC970 all
    exit
endif

if(`uname -m` == "x86_64") then
    set job = "make -j8 -f Makefile-Tinker-Intel"
    setenv SRC `pwd`
    setenv  B  "$SRC/../Bin-Serial"
    setenv OpenMP "" #-DI64"
    $job

    setenv B  "$SRC/../Bin-OpenMP"
    setenv OpenMP  "-openmp"
    $job
    
    exit
endif
