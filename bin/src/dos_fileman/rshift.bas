DECLARE SUB inputline (n!, d1#(), col1%)
DECLARE SUB iteratefile (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "RSHIFT RSHIFT RSHIFT RSHIFT RSHIFT RSHIFT RSHIFT RSHIFT RSHIFT RSHIFT"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA " program to shift the a function relativ to another function             "
DATA " command: RSHIFT *1.* 23 *2.* 14 24.4                                    "
DATA " function y(x)=col3(col2)  in file *1.* is shifted in y by adding a constant"
DATA " the constant is chosen such that the minimum difference between           "
DATA " y=col3 (in file *1.*) and y=col4(col1) (in file *2.*) is 24.4            "
DATA " options: /@list instead of *1.* takes filenames out of list file                          "
DATA " format of file                                                             "
DATA " { header: this is the                                                      "
DATA "   file header delimited                                                    "
DATA "   by brackets after this header there follow 3 or more data columns        "
DATA " time[s] T[K] dl/l}                                                         "
DATA " 11 3.14235 65367                                                           "
DATA "  .    .     .                                                              "
DATA "  .    .     .                                                              "
DATA " 32 2412.34 324.2                                                           "
DATA "                                                                            "
DATA "programm needs to call: shiftc.bas fact.bas add.bas                         "
DATA "*************************************************************************   "
DIM text$(300), x#(30), d1$(30)

' aaaaaaaaaaaa ANALYZE COMMAND$ aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
Ala$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))

 col1% = 0: WHILE LEN(Ala$) > 0  'split$
    col1% = col1% + 1
    Ala$ = LTRIM$(Ala$) + " "
    d1$(col1%) = LEFT$(Ala$, INSTR(Ala$, " "))
    Ala$ = LTRIM$(RIGHT$(Ala$, LEN(Ala$) - INSTR(Ala$, " ")))
 WEND

file1$ = d1$(1)
file2$ = d1$(3)
ij1$ = d1$(2)
ij2$ = d1$(4)
mindiff = VAL(LTRIM$(d1$(5)))
nooffiles% = 0

999 CALL iteratefile(file1$)
i1% = ASC(LEFT$(ij1$, 1)) MOD 48: j1% = ASC(MID$(ij1$, 2, 1)) MOD 48'get columns
i2% = ASC(LEFT$(ij2$, 1)) MOD 48: j2% = ASC(MID$(ij2$, 2, 1)) MOD 48'get columns

maximum# = -1E+37
minimum# = 1E+37
nooffiles% = nooffiles% + 1
SHELL "copy " + file2$ + " rshift.d" + LTRIM$(STR$(nooffiles%))

FOR jj% = 1 TO nooffiles%
 SHELL "copy " + file1$ + " rshift.dum"
 PRINT
 SHELL "fact rshift.dum " + STR$(j1%) + " -1"
 PRINT
 SHELL "add rshift.dum " + ij1$ + " rshift.d" + LTRIM$(STR$(jj%)) + " " + ij2$ + " /noex0"
 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
 ' ffffffffffffffffffff FIND maximum of x to determine shift constant ffffffffffffffff
 OPEN "i", 1, "rshift.dum": CALL headerinput(text$(), j, 1)' open file+ input file header
 WHILE EOF(1) = 0
  CALL inputline(1, x#(), col%)
  IF -x#(j1%) > maximum# THEN maximum# = -x#(j1%)
  IF -x#(j1%) < minimum# THEN minimum# = -x#(j1%)
  WEND  'ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
5 CLOSE 1
NEXT jj%

IF ABS(maximum#) > 1E+36 OR ABS(minimum#) > 1E+36 THEN PRINT "ERROR PROGRAM RSHIFT I could not calculate the difference between "; file1$; " and "; file2$: END

' add constant to column j1% using program shift
IF mindiff > 0 THEN
 constant = mindiff - minimum#
ELSE
 constant = mindiff - maximum#
END IF

 SHELL "shiftc " + file1$ + " " + STR$(j1%) + " " + STR$(constant)
PRINT
PRINT "END RSHIFT"
IF INSTR(COMMAND$, "*") <> 0 OR INSTR(COMMAND$, "/@") <> 0 THEN
 file2$ = file1$: ij2$ = ij1$
 GOTO 999
END IF
 SHELL "del rshift.dum"
 SHELL "del rshift.*"
PRINT
END

333 FOR I = 1 TO 19: READ a$: PRINT a$: NEXT I: END

SUB headerinput (text$(), j, n)
'**********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'**********************************************************************

1 INPUT #n, a$
   I = INSTR(a$, "{"): IF I > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(a$, LEN(a$) - I)
  j = 1: I = INSTR(text$(j), "}"): IF I > 0 GOTO 3  'look for "}" in first line
                                   
   FOR j = 2 TO 300
   INPUT #n, text$(j)
   I = INSTR(text$(j), "}"): IF I > 0 GOTO 3      'look for "}"
   NEXT j: PRINT "text in data file too long": END
3 text$(j) = LEFT$(text$(j), I - 1)

END SUB

SUB inputline (n, d1#(), col1%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)

LINE INPUT #n, Ala$
col1% = 0
 WHILE LEN(Ala$) > 0
    col1% = col1% + 1
    Ala$ = LTRIM$(Ala$) + " "
    d1#(col1%) = VAL(LEFT$(Ala$, INSTR(Ala$, " ")))
    Ala$ = LTRIM$(RIGHT$(Ala$, LEN(Ala$) - INSTR(Ala$, " ")))
 WEND


END SUB

SUB iteratefile (file1$)
STATIC washere%

  IF INSTR(COMMAND$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere% = 0 THEN
      washere% = 1
      'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /b > rshift.dir"
      OPEN "i", 9, "rshift.dir"
   END IF
   IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del rshift.dir": END
   INPUT #9, file1$
   IF LCASE$(file1$) = "fact.dir" THEN
    IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del fact.dir": END
    INPUT #9, file1$
   END IF
   IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del rshift.dir": END
  END IF
 
  a% = INSTR(COMMAND$, "/@")
  IF a% <> 0 THEN
   IF washere% = 0 THEN
      washere% = 1
      'get filenames to be operated on as file1$
      b$ = LTRIM$(MID$(COMMAND$, a% + 2))
      b% = INSTR(b$, " "): IF b% > 0 THEN b$ = LEFT$(b$, b%)
      OPEN "i", 9, b$
   END IF
   IF EOF(9) <> 0 THEN CLOSE 9: END
   INPUT #9, file1$
   IF LTRIM$(file1$) = "" THEN CLOSE 9: END
  END IF


END SUB

