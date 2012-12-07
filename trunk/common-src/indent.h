

/* $Id: indent.h 2 2007-10-20 14:59:10Z wangsl $ */

#ifndef INDENT_H
#define INDENT_H

#ifdef __cplusplus
extern "C" {
#endif

  void IndentPush();
  void IndentPop();
  const char *Indent();

#ifdef __cplusplus
}
#endif

#endif /* INDENT_H */
