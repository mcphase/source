DECLARE SUB iterate (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "SPLIT1 SPLIT1 SPLIT1 SPLIT1 SPLIT1 SPLIT1 SPLIT1 SPLIT1 SPLIT1 SPLIT1"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program SPLIT1 - use it like 'SPLIT1 *.* RRRR'                       "
DATA "    (means SPLIT file *.* into files  named *01.* *02.* ..."
DATA "     by looking for RRRR (not case sensitive) at line begin)"
DATA "*************************************************************************"
DIM text$(300), xm#(30)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(a$, INSTR(a$, " "))): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))
code$ = UCASE$(LTRIM$(RTRIM$(a$)))

999 CALL iterate(filename$)

'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'oooooooooooooooooooooooooooooooo OPEN FILES oooooooooooooooooooooooooooo
'open input file and input file header
OPEN "i", 1, filename$

i% = 0
WHILE EOF(1) = 0: i% = i% + 1
pp% = INSTR(filename$, ".")
IF pp% = 0 THEN
pp% = LEN(filename$) + 3: ext$ = ""
trunc$ = LEFT$(filename$, pp% - 3)
ELSE
 IF pp% > 7 THEN
   trunc$ = LEFT$(filename$, pp% - 3)
   ext$ = MID$(filename$, pp%)
  ELSE
   trunc$ = LEFT$(filename$, pp% - 1)
   ext$ = MID$(filename$, pp%)
  END IF
END IF

no$ = LTRIM$(RTRIM$(STR$(i%))): IF i% < 10 THEN no$ = "0" + no$
outputfile$ = trunc$ + no$ + ext$



' open output file and write fileheader
OPEN "o", 2, outputfile$
'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

REM input data columns>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

DO
  PRINT #2, a$
  IF EOF(1) = 0 THEN LINE INPUT #1, a$
LOOP WHILE LEFT$(a$, LEN(code$)) <> UCASE$(code$) AND EOF(1) = 0

'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 CLOSE 2
WEND
 CLOSE 1
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END

333 FOR i = 1 TO 5: READ a$: PRINT a$: NEXT: END

SUB iterate (file1$)
STATIC washere%

  IF INSTR(COMMAND$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere% = 0 THEN
      washere% = 1
      'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /b > SPLIT1.dir"
      OPEN "i", 9, "SPLIT1.dir"
   END IF
   IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del SPLIT1.dir": END
   INPUT #9, file1$
   IF LCASE$(file1$) = "split1.dir" THEN
    IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del split1.dir": END
    INPUT #9, file1$
   END IF
   IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del SPLIT1.dir": END
  END IF


END SUB

