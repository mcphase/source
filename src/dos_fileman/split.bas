DECLARE SUB iterate (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "SPLIT SPLIT SPLIT SPLIT SPLIT SPLIT SPLIT SPLIT SPLIT SPLIT"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program SPLIT - use it like 'SPLIT *.* 32'                       "
DATA "    (means split file *.* into 32 files of equal size named *01.* *02.* ..."
DATA "*************************************************************************"
DIM text$(300), xm#(30)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(a$, INSTR(a$, " "))): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))
noofblocks% = VAL(a$)  'get

999 CALL iterate(filename$)

'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'oooooooooooooooooooooooooooooooo OPEN FILES oooooooooooooooooooooooooooo
'open input file and input file header
OPEN "b", 1, filename$

FOR i% = 1 TO noofblocks%
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
FOR j = 1 TO LOF(1) / noofblocks% * 1#
a$ = INPUT$(1, #1)
' GET 1, a
PRINT #2, a$;
NEXT j

'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 CLOSE 2: NEXT i%: CLOSE 1
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END

333 FOR i = 1 TO 4: READ a$: PRINT a$: NEXT: END

SUB iterate (file1$)
STATIC washere%

  IF INSTR(COMMAND$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere% = 0 THEN
      washere% = 1
      'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /b > split.dir"
      OPEN "i", 9, "split.dir"
   END IF
   IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del split.dir": END
   INPUT #9, file1$
   IF LCASE$(file1$) = "split.dir" THEN
    IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del split.dir": END
    INPUT #9, file1$
   END IF
   IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del split.dir": END
  END IF


END SUB

