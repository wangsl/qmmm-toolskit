 
C     $Id: crashF.F 2 2008-02-28 22:52:52Z wangsl $


      SubRoutine CrashLoc(Message, Module, Line)
      Implicit None
      Include 'iounit.i'
      Character*(*) Message, Module
      Integer Line
      
 1000 Format(/, ' Fatal error occured in module: ', A,
     $     ', line ', I4, ':', /, 2X, A,/)
      
      Write(IOut, 1000) Module, Line, Message
      Stop
      Return
      End
      
      SubRoutine WarnLoc(Message, Module, Line)
      Implicit None
      Include 'iounit.i'
      Character*(*) Message, Module
      Integer Line
      
 1000 Format(/, ' Warning occured in module: ', A,
     $     ', line ', I4, ':', /, 2X, A,/)
      
      Write(IOut, 1000) Module, Line, Message
      Return
      End
