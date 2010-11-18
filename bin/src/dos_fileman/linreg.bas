DECLARE SUB analizecommand (file1$)
DECLARE SUB inputline (n!, d1#(), col1%)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "LINREG LINREG LINREG LINREG LINREG LINREG LINREG LINREG LINREG LINREG"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program LINREG - use it like LINREG *.* 23"            
DATA "    (23 means calc.LINREG of 3rd columns with respect to second column) - "
DATA " ---> the result is written in file *.*,                                 "
DATA " the LINear REGression is calculated by least squares method              "
DATA " the linear regression is written to the last column of file *.*         "
DATA " format of file:                                                          "
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "*************************************************************************"
DIM text$(300), x#(30, 100), d#(30)

'aaaaa analyse command string (get filename$ - ii% - jj% and interval) aaaaaa
A$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))

filename$ = RTRIM$(LEFT$(A$, LEN(A$) - 2))
jj% = ASC(RIGHT$(A$, 1)) MOD 48: ii% = ASC(MID$(A$, LEN(A$) - 1, 1)) MOD 48  'get column numbers

999 CALL analizecommand(filename$)

IF RIGHT$(filename$, 4) = ".mrc" OR RIGHT$(filename$, 4) = ".rcp" THEN
   PRINT "you should never ever change data in *.rcp files": PLAY "dgdgdg"
22 INPUT "do you really want to continue (Y/N)"; ALA$: IF LCASE$(ALA$) = "n" THEN END
   IF LCASE$(ALA$) <> "y" GOTO 22
END IF
'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'oooooooooooooooooo OPEN FILES oooooooooooooooooooooooooooooooooooooooooooooo
'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)
'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


REM input data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
xxm# = 0: xym# = 0: ym# = 0: xm# = 0: I = 0

'input data and calculate linear regression  iiiiiiiiiiiiiiiiiiiii
WHILE EOF(1) = 0
CALL inputline(1, d#(), col%)
 I = I + 1
xm# = xm# + d#(ii%): ym# = ym# + d#(jj%)
xym# = xym# + d#(ii%) * d#(jj%)
xxm# = xxm# + d#(ii%) * d#(ii%)
WEND: CLOSE
xm# = xm# / I: ym# = ym# / I: xym# = xym# / I: xxm# = xxm# / I
' here calculate the y=offset+slope*x  of the column jj% with respect to column ii%
slope# = (xym# - ym# * xm#) / (xxm# - xm# * xm#)
offset# = ym# - slope# * xm#
'iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

'  write to output file  ddddddddddddddddddddddddddd
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

' open output file and write fileheader
OPEN "o", 2, "LINREG.lin"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " column"; jj%; "(av="; ym#; "LINREG calc. with respect to column"; ii%; "(av="; xm#;
PRINT #2, " slope="; slope#; "offset="; offset#; " regression put into column"; col% + 1
PRINT #2, "by program LINREG.bas}"

WHILE EOF(1) = 0
CALL inputline(1, d#(), col%)
FOR coll% = 1 TO col%: PRINT #2, d#(coll%); : NEXT: PRINT #2, offset# + slope# * d#(ii%)
WEND
'ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd

'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

CLOSE 1, 2
SHELL "copy LINREG.lin " + filename$
SHELL "del LINREG.lin"
PRINT
PRINT "END LINREG in file "; filename$;
PRINT " column"; jj%; "(av="; ym#; "LINREG calc. with respect to column"; ii%; "(av="; xm#;
PRINT " slope="; slope#; "offset="; offset#; " regression put into column"; col% + 1
IF INSTR(COMMAND$, "*") <> 0 GOTO 999

END
333 FOR I = 1 TO 14: READ A$: PRINT A$: NEXT I: END

SUB analizecommand (file1$)
STATIC washere%
  'this is for cumulative action on more files
  IF INSTR(COMMAND$, "*") <> 0 THEN
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
'*********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'*********************************************************************

1 INPUT #n, A$
   I = INSTR(A$, "{"): IF I > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(A$, LEN(A$) - I)
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

INPUT #n, ALA$
col1% = 0
 WHILE LEN(ALA$) > 0
    col1% = col1% + 1
    ALA$ = LTRIM$(ALA$) + " "
    d1#(col1%) = VAL(LEFT$(ALA$, INSTR(ALA$, " ")))
    ALA$ = LTRIM$(RIGHT$(ALA$, LEN(ALA$) - INSTR(ALA$, " ")))
 WEND
END SUB

