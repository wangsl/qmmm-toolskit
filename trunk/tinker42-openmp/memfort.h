 
/* $Id: memfort.h 2 2008-02-28 22:52:52Z wangsl $ */

#ifndef FORTDEF_H
#define FORTDEF_H

#define Crash(x) CrashLoc(x, __FILE__, __LINE__)
#define crash(x) CrashLoc(x, __FILE__, __LINE__)
#define CRASH(x) CrashLoc(x, __FILE__, __LINE__)

#define Warn(x) WarnLoc(x, __FILE__, __LINE__)
#define warn(x) WarnLoc(x, __FILE__, __LINE__)
#define WARN(x) WarnLoc(x, __FILE__, __LINE__)

#define FORTAllocReal8(x, y) FORTAllocReal8Loc(x, y, __FILE__, LinEnd(__FILE__), __LINE__)
#define fortallocreal8(x, y) FORTAllocReal8Loc(x, y, __FILE__, LinEnd(__FILE__), __LINE__)
#define FORTALLOCREAL8(x, y) FORTAllocReal8Loc(x, y, __FILE__, LinEnd(__FILE__), __LINE__)

#define FORTAllocInt(x, y) FORTAllocIntLoc(x, y, __FILE__, LinEnd(__FILE__), __LINE__)
#define fortallocint(x, y) FORTAllocIntLoc(x, y, __FILE__, LinEnd(__FILE__), __LINE__)
#define FORTALLOCINT(x, y) FORTAllocIntLoc(x, y, __FILE__, LinEnd(__FILE__), __LINE__)

#define FORTAllocChar(x, y) FORTAllocCharLoc(x, y, __FILE__, LinEnd(__FILE__), __LINE__)
#define fortallocchar(x, y) FORTAllocCharLoc(x, y, __FILE__, LinEnd(__FILE__), __LINE__)
#define FORTALLOCCHAR(x, y) FORTAllocCharLoc(x, y, __FILE__, LinEnd(__FILE__), __LINE__)

#define FORTFree(x) FORTFreeLoc(x, __FILE__, LinEnd(__FILE__), __LINE__)
#define fortfree(x) FORTFreeLoc(x, __FILE__, LinEnd(__FILE__), __LINE__)
#define FORTFREE(x) FORTFreeLoc(x, __FILE__, LinEnd(__FILE__), __LINE__)

#define MemoryInformation FortranMemInf(__FILE__, LinEnd(__FILE__), __LINE__)
#define MEMORYINFORMATION FortranMemInf(__FILE__, LinEnd(__FILE__), __LINE__)
#define memoryinformation FortranMemInf(__FILE__, LinEnd(__FILE__), __LINE__)

#define QCrash(x) Crash(x)
#define FORTQFree(x) FORTFree(x)
#define FORTQAllocReal8(x, y) FORTAllocReal8(x, y)

#endif /* FORTDEF_H */


