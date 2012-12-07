#!/bin/env python

# $Id$

class Cartesian :
    def __init__(self, x = 0.0, y = 0.0, z = 0.0) :
        self.__elements = [x, y, z]
        return
    
    def __repr__(self) :
        return "%15.8f%15.8f%15.8f" % (self.__elements[0], 
                                       self.__elements[1],
                                       self.__elements[2])
    
    def __iadd__(self, that) :
        for i in range(3) :
            self.__elements[i] += that.__elements[i]
        return self
        
    def __isub__(self, that) :
        for i in range(3) :
            self.__elements[i] -= that.__elements[i]
        return self

    def __imul__(self, s) :
        for i in range(3) :
            self.__elements[i] *= s
        return self

    def __add__(self, that) :
        return Cartesian(self.x() + that.x(),
                         self.y() + that.y(),
                         self.z() + that.z())

    def __sub__(self, that) :
        return Cartesian(self.x() - that.x(),
                         self.y() - that.y(),
                         self.z() - that.z())

    def __mul__(self, s) :
        return Cartesian(self.x()*s, self.y()*s, self.z()*s)

    def __rmul__(self, s) :
        return Cartesian(self.x()*s, self.y()*s, self.z()*s)
        
    def __div__(self, s) :
        return Cartesian(self.x()/s, self.y()/s, self.z()/s)

    def __neg__(self) :
        return Cartesian(-self.x(), -self.y(), -self.z())

    def __abs__(self) :
        return self.magnitude()

    def __getitem__(self, i) :
        assert(0 <= i and i < 3)
        return self.__elements[i]

    def __len__(self) :
        return len(self.__elements)

    def dot(self, that) :
        return self.x()*that.x() + self.y()*that.y() + self.z()*that.z()
    
    def cross(self, that) :
        return Cartesian(self.y()*that.z() - self.z()*that.y(),
                         self.z()*that.x() - self.x()*that.z(),
                         self.x()*that.y() - self.y()*that.x())

    def magnitude(self) :
        from math import sqrt
        return sqrt(self.dot(self))

    def magnitude_sq(self) :
        return self.dot(self)

    def distance(self, that) :
        return (self-that).magnitude()

    def distance_sq(self, that) :
        return (self-that).magnitude_sq()
    
    def apply_function(self, func) :
        for i in range(3) :
            self.__elements[i] = func(self.__elements[i])
        return

    def map_function(self, func) :
        return Cartesian(func(self.__elements[0]),
                         func(self.__elements[1]),
                         func(self.__elements[1]))

    def scale(self, s) :
        for i in range(3) :
            self.__elements[i] *= s
        return self

    def copy(self) :
        return Cartesian(self.x(), self.y(), self.z())

    def scale_to_unit_vector(self) :
        return self.scale(1/self.magnitude())

    def as_unit_vector(self) :
        r = self.copy()
        r.scale(1/r.magnitude())
        return r

    def tuple(self) :
        return tuple(self.__elements)
    
    def x(self) : return self.__elements[0]
    def y(self) : return self.__elements[1]
    def z(self) : return self.__elements[2]

    def set_x(self, x) :
        self.__elements[0] = x;
        return

    def set_y(self, y) :
        self.__elements[1] = y;
        return

    def set_z(self, z) :
        self.__elements[2] = z;
        return
    
    def elements(self) :
        return self.__elements

    def set(self, that) :
        if isinstance(that, Cartesian) :
            self.__elements = [that.x(), that.y(), that.z()]
            return
        else :
            raise ValueError, "Cartesian::set Expected a Cartesian coordinates"
    
if __name__ == "__main__" :

    from math import pi
    r1 = Cartesian(1.1, -2.3, 5.6);
    print r1

    print len(r1)
    print r1[0], r1[2]
    print r1.magnitude()
    print abs(r1)
