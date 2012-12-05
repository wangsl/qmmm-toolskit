
/* $Id: table.h 148 2009-03-21 15:52:57Z wangsl $ */

#ifndef TABLE_H
#define TABLE_H

#include <cassert>
#include <cstdio>
#include <iostream>
using namespace std;
#include <cstring>
#include "indent.h"

#include "die.h"

/***
 *
 *  Table of type T indexed by char *
 *
 ***/

template<class T> class TableIterator;
template<class T> class ConstTableIterator;

template<class T> class Table {
  friend class TableIterator<T>;
  friend class ConstTableIterator<T>;
  struct item {
    char *key;
    T val;
    item *next;
    item(const char *k, const T &ob, item *n) : key(new char[strlen(k)+1]), 
						val(ob), next(n) 
    { strcpy(key,k); }
    ~item() 
    { delete[] key; delete next; key = 0; next = 0; }
  };
  
  item **items;
  int n;
  unsigned hash(const char *) const;
  
public:
  T default_value;
  inline const T& operator[](const char *) const;
  inline T& operator[](const char *);
  
  Table(int sz = 101) : n(sz) 
  {
    items = new item *[n];
    for (int i = 0; i < n; i++)
      items[i] = 0;
  }
  
  Table(int sz, const T &def) : n(sz), default_value(def)
  {
    items = new item *[n];
    for (int i = 0; i < n; i++)
      items[i] = 0;
  }
  
  Table(const Table<T> &t) : n(t.n)
  {
    int i;
    items = new item *[n];
    for(i = 0; i < n; i++)
      items[i] = 0;
    for (i = 0; i < n; i++)
      for (item *j = t.items[i]; j; j = j->next)
	(*this)[j->key] = j->val;
    default_value = t.default_value;
  }
  
  Table<T> & operator=(const Table<T> &t)
  {
    int i;
    remove_all_entries();
    delete[] items;
    n = t.n;
    items = new item *[n];
    for (i = 0; i < n; i++)
      items[i] = 0;
    for (i = 0; i < n; i++)
      for (item *j = t.items[i]; j; j = j->next)
	(*this)[j->key] = j->val;
    default_value = t.default_value;
    return *this;
  }
  
  ~Table();
  int entries() const;
  int size() const
  { return entries(); }
  int exists(const char *) const;
  const char *key(const char *) const;
  void remove_all_entries();
  int is_empty() const;
  
  // By Shenglong Wang 11-24-2004
  void combine(Table<T> b);
};

template<class T> inline Table<T>::~Table()
{
  remove_all_entries();
  delete[] items;
  items = 0;
}

template<class T> inline T& Table<T>::operator[](const char *s)
{
  int k = hash(s) % n;
  item *i;
  for (i = items[k]; i && strcmp(s, i->key); i = i->next)
    ;
  if (i == 0) { /* not found; create a new item */
    i = new item(s, default_value, items[k]);
    items[k] = i;
  }
  return i->val;
}

template<class T> inline const T & Table<T>::operator[](const char *s) const
{
  int k = hash(s) % n;
  for (item *i = items[k]; i; i = i->next)
    if (!strcmp(s, i->key))
      return i->val;
  return default_value;
}

template<class T> inline int Table<T>::exists(const char *s) const
{
  int k = hash(s) % n;
  for (item *i = items[k]; i; i = i->next)
    if (!strcmp(s, i->key))
      return 1;
  return 0;
}

template<class T> inline const char * Table<T>::key(const char *s) const
{
  int k = hash(s) % n;
  item *i;
  for (i = items[k]; i; i = i->next)
    if (!strcmp(s, i->key))
      return i->key;
  return 0;
}

template<class T> inline unsigned Table<T>::hash(const char *s) const
{
  unsigned k;
  for (k = 0; *s; s++)
    k = 31*k + (*s);
  return k;
}

template<class T> inline int Table<T>::entries() const
{
  int e = 0;
  for (int i = 0; i < n; i++)
    for (item *j = items[i]; j; j = j->next)
      e++;
  return e;
}

template<class T> inline void Table<T>::remove_all_entries()
{ 
  for (int i = 0; i < n; i++) {
    delete items[i];
    items[i] = 0;
  }
}

template<class T> class TableIterator
{
  Table<T> *t;
  int k;
  typename Table<T>::item *i;
public:
  TableIterator() : t(0) { }
  void init(Table<T> &tab)
  {
    t = &tab;
    k = 0;
    i = t->items[0];
    if (!i)
      next();
  }
  int ok() const
  { return k < t->n; }
  const char *key()
  { return i->key; }
  T & val()
  { return i->val; }
  
  T & value() { return val(); }
  
  void next()
  { 
    if (i == 0 || i->next == 0)
      for (k++ ; k < t->n && (i = t->items[k]) == 0; k++)
	;
    else
      i = i->next;
  }
};

template<class T> class ConstTableIterator
{
  const Table<T> *t;
  int k;
  const typename Table<T>::item *i;
public:
  ConstTableIterator() : t(0) { }
  void init(const Table<T> &tab)
  {
    t = &tab;
    k = 0;
    i = t->items[0];
    if (!i)
      next();
  }
  int ok() const
  { return k < t->n; }
  const char *key()
  { return i->key; }
  const T & val()
  { return i->val; }
  
  const T & value() { return val(); }
  
  void next()
  { 
    if (i == 0 || i->next == 0)
      for (k++ ; k < t->n && (i = t->items[k]) == 0; k++)
	;
    else
      i = i->next;
  }
};

template<class T> inline ostream & operator<<(ostream &s, const Table<T> &t)
{
  ConstTableIterator<T> ti;
  int i;
  int n = t.entries();
  const char **keys = new const char*[n];
  for (i = 0, ti.init(t); ti.ok(); i++, ti.next())
    keys[i] = ti.key();
  /* Bubble sort */
  int any = 1;
  while(any)
    for (any = 0, i = 1; i < n; i++)
      if (strcmp(keys[i-1], keys[i]) > 0) {
	const char *tmp = keys[i];
	keys[i] = keys[i-1];
	keys[i-1] = tmp;
	any = 1;
      }
  s << "{\n";
  IndentPush();
  for (i = 0; i < n; i++)
    s << Indent() << keys[i] << " " << t[keys[i]] << "\n";
  IndentPop();
  delete[] keys;
  keys = 0;
  return s << Indent() << "}\n";
}

template<class T> inline istream & operator>>(istream &s, Table<T> &t)
{
  char bracket;
  char key[256];
  s >> bracket;
  if (bracket != '{')
    QCrash("Table::operator>>: expecting '{', found something else");
  while(s) {
    s >> bracket;
    if (bracket == '}')
      return s;
    else
      s.putback(bracket);
    s >> key;
    s >> t[key];
  }
  QCrash("Table::operator>>: EOF reached before reading '}'");
  return s;
}

// By Shenglong Wang

template<class T> inline int Table<T>::is_empty() const
{
  return entries() > 0 ? 0 : 1;
}

template<class T> inline void Table<T>::combine(Table<T> b)
{
  TableIterator<T> bi;
  for(bi.init(b); bi.ok(); bi.next()) {
    const char *b_key = bi.key();
    (*this)[b_key] = bi.val();
  }
}

#endif /* TABLE_H */
