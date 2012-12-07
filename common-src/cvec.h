/* $Id: cvec.h 2 2007-10-20 14:59:10Z wangsl $ */

#ifndef CVEC_H
#define CVEC_H

#include <string.h>
#include "vec.h"
#include "coord.h"

class CVec : public Vec<Cartesian>
{
public:
  CVec(int sz = 0) : Vec<Cartesian>(sz) { } 
  CVec(const CVec &v) : Vec<Cartesian>(v) { }
  CVec(int sz, const Cartesian &t) : Vec<Cartesian>(sz,t) { }

  // By Shenglong Wang
  CVec(int sz, Cartesian *q) : Vec<Cartesian>(sz, q) { }
  CVec(int sz, double *q) : Vec<Cartesian>(sz, (Cartesian *) q) { }

  CVec & operator=(const CVec &v) 
  { return (CVec &) Vec<Cartesian>::operator=(v); }
  CVec & apply(double (*f)(double))
  { 
    for (int i = 0; i < size(); i++)
	(*this)[i].apply(f);
    return *this;
  }
  CVec map(Cartesian (Cartesian::*f)() const) const
  {
    CVec v(size());
    for (int i = 0; i < size(); i++)
      v[i] = ((*this)[i].*f)();
    return v;
  }
  CVec as_unit_vector() const
  { return map(&Cartesian::as_unit_vector); }
  CVec & copy(const Cartesian *v)
  {
    memcpy((Cartesian *) *this, v, size()*sizeof(Cartesian));
    return *this;
  }
  CVec copy() const
  {
    CVec v(size());
    v.copy(*this);
    return v;
  }
  void copy_from(const double *x, const double *y, const double *z)
  {
    Cartesian *c = (*this);
    int n = size();
    for (int i = 0; i < n; i++) {
      c[i].x = x[i];
      c[i].y = y[i];
      c[i].z = z[i];
    }
  }
  void copy_to(double *x, double *y, double *z) const
  {
    const Cartesian *c = (*this);
    int n = size();
    for (int i = 0; i < n; i++) {
      x[i] = c[i].x;
      y[i] = c[i].y;
      z[i] = c[i].z;
    }
  }
  void copy_from(const double *x)
  { memcpy((Cartesian *) *this, (const Cartesian *) x, size()*sizeof(Cartesian)); }
  void copy_to(double *x) const
  { memcpy((Cartesian *) x, (const Cartesian *) *this, size()*sizeof(Cartesian)); }
  void copy_to(Cartesian *x) const
  { memcpy(x, (const Cartesian *) *this, size()*sizeof(Cartesian)); }
  void copy_from(const double xyz[][3])
  { memcpy((Cartesian *) *this, (const Cartesian *) xyz, size()*sizeof(Cartesian)); }
  void copy_to(double xyz[][3]) const
  { memcpy((Cartesian *) xyz, (const Cartesian *) *this, size()*sizeof(Cartesian)); }
  void add_to(double xyz[][3]) const
  {
    for (int i = 0; i < size(); i++)
      *((Cartesian *) xyz[i]) += (*this)[i];
  }
  void zero()
  { memset(*this, 0, size()*sizeof(Cartesian)); }
  friend double operator*(const CVec a, const CVec b)
  {
    double c = 0.0;
    assert(a.is_conformable_with(b));
    for (int i = 0; i < a.size(); i++)
      c += a[i]*b[i];
    return c;
  }
  double sq() const { return (*this)*(*this); }
  CVec & operator+=(const CVec a)
  {
    assert(is_conformable_with(a));
    for (int i = 0; i < size(); i++)
      (*this)[i] += a[i];
    return *this;
  }
  CVec & operator-=(const CVec a)
  {
    assert(is_conformable_with(a));
    for (int i = 0; i < size(); i++)
      (*this)[i] -= a[i];
    return *this;
  }
  CVec & operator*=(double x)
  {
    for (int i = 0; i < size(); i++)
      (*this)[i] *= x;
    return *this;
  }
  Cartesian sum() const
  {
    Cartesian c(0,0,0);
    for (int i = 0; i < size(); i++)
      c += (*this)[i];
    return c;
  }
  Cartesian average() const
  { return sum()/size(); }
  Tensor outer_product(const CVec bv) const
  {
    assert(is_conformable_with(bv));
    Tensor t(0);
    const int n = size();
    const Cartesian *a = *this;
    const Cartesian *b = bv;
    for (int i = 0; i < n; i++) {
      t.xx += a[i].x * b[i].x;
      t.xy += a[i].x * b[i].y;
      t.xz += a[i].x * b[i].z;
      t.yx += a[i].y * b[i].x;
      t.yy += a[i].y * b[i].y;
      t.yz += a[i].y * b[i].z;
      t.zx += a[i].z * b[i].x;
      t.zy += a[i].z * b[i].y;
      t.zz += a[i].z * b[i].z;
    }
    return t;
  }
};
  
#endif /* CVEC_H */
