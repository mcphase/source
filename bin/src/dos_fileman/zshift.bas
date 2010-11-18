DECLARE SUB inputline (n!, d1#(), col1%)
DECLARE SUB analizecommand (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "ZSHIFT ZSHIFT ZSHIFT ZSHIFT ZSHIFT ZSHIFT ZSHIFT ZSHIFT ZSHIFT ZSHIFT"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA " program to shift the zero of a function                                 "
DATA " command: zshift *.* 23 24.4                                             "
DATA " function y(x)=col3(col2)  y values are shifted bz a constant such that
DATA " for x=24.4 y(x)=0"
DATA " format of file                                                             "
DATA "                                                                            "
DATA " { header: this is the                                                      "
DATA "   file header delimited                                                    "
DATA "   by brackets after this header there follow 3 or more data columns        "
DATA " time[s] T[K] dl/l}                                                         "
DATA " 11 3.14235 65367                                                           "
DATA "  .    .     .                                                              "
DATA "  .    .     .      .  .  .                                                 "
DATA "  .    .     .                                                              "
DATA "  .    .     .                                                              "
DATA " 32 2412.34 324.2                                                           "
DATA "                                                                            "
DATA "                                                                            "
DATA "*************************************************************************   "
DIM text$(300), x#(30), xold#(30)

' aaaaaaaaaaaa ANALYZE COMMAND$ aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))
I = INSTR(a$, " "): filename$ = LEFT$(a$, I): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - I))
ii% = ASC(LEFT$(a$, 1)) MOD 48: jj% = ASC(MID$(a$, 2, 1)) MOD 48'get columns
tnorm = VAL(RIGHT$(a$, LEN(a$) - 2))
999 CALL analizecommand(filename$)
'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
CONSTANT = 1.3E+38
' ffffffffffffffffffff FIND X CORRESPONDING TO ZERO Y fffffffffffffffffffff
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)' open file+ input file header

IF EOF(1) <> 0 GOTO 5
CALL inputline(1, x#(), col%)

WHILE EOF(1) = 0  ' REM input columns
 FOR coll% = 1 TO col%: xold#(coll%) = x#(coll%): NEXT

CALL inputline(1, x#(), col%)

 IF xold#(ii%) <= tnorm + .000001 AND x#(ii%) >= tnorm - .000001 THEN CONSTANT = -x#(jj%): GOTO 5
 IF xold#(ii%) >= tnorm - .000001 AND x#(ii%) <= tnorm + .000001 THEN CONSTANT = -x#(jj%): GOTO 5
WEND  'ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

5 CLOSE 1
IF CONSTANT = 1.3E+38 THEN PRINT filename$; ":"; tnorm; "out of range - nothing has been shifted": GOTO 6

OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)' open file+ input file header
OPEN "o", 2, "fact.fac"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " column"; jj%;
PRINT #2, "shiftec by "; CONSTANT; " using program zshift.bas}"

WHILE EOF(1) = 0  ' REM input columns

CALL inputline(1, x#(), col%)
 IF col% = -1 GOTO 7
 x#(jj%) = x#(jj%) + CONSTANT
 FOR coll% = 1 TO col%: PRINT #2, x#(coll%); : NEXT
PRINT #2,

WEND

7 PRINT
CLOSE 1, 2

SHELL "copy fact.fac " + filename$
SHELL "del fact.fac"

PRINT "END ZSHIFT file "; filename$; ": function y(x)=col"; jj%; "(col"; ii%; ") - zero-y has been shifted "
PRINT "to x="; tnorm; " by adding a const.(lin.interpol. of data points)"
6 IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END

333 FOR I = 1 TO 20: READ a$: PRINT a$: NEXT I: END

SUB analizecommand (file1$)
STATIC washere%

  IF INSTR(COMMAND$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere% = 0 THEN
      washere% = 1
      'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /b > fact.dir"
      OPEN "i", 9, "fACT.dir"
   END IF
   IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del fact.dir": END
   INPUT #9, file1$
   IF LCASE$(file1$) = "fact.dir" THEN
    IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del fact.dir": END
    INPUT #9, file1$
   END IF
   IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del fact.dir": END
  END IF


END SUB

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

LINE INPUT #n, ALA$
col1% = 0
 WHILE LEN(ALA$) > 0
    col1% = col1% + 1
    ALA$ = LTRIM$(ALA$) + " "
    d1#(col1%) = VAL(LEFT$(ALA$, INSTR(ALA$, " ")))
    ALA$ = LTRIM$(RIGHT$(ALA$, LEN(ALA$) - INSTR(ALA$, " ")))
 WEND


END SUB

