#!/usr/bin/python

from numpy import zeros, float64

n = 10000

a = zeros(n, float64)

i = 0
while i < n :
    a[i] = i
    i += 1

print a

i = 0
while i < len(a) :
    print a[i]
    i += 1
    
