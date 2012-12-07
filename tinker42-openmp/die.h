 
/* $Id$ */

#ifndef DIE_H
#define DIE_H

void die_at(const char *s, const char *file, int line);

#define die(x) die_at(x, __FILE__, __LINE__)

#define insist(x) if (!(x)) die_at("insist failed: " #x, __FILE__, __LINE__)

#endif /* DIE_H */
