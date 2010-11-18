DECLARE SUB Inputline (n!, D1#(), col1%, d2$(), col2%)
DECLARE SUB analizecommand (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "FACT FACT FACT FACT FACT FACT FACT FACT FACT FACT FACT FACT FACT"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program fact - use it like FACT *.* 2 3.47"           
DATA "    (means multiply column 2 with 3.47) -                                "
DATA " ---> the result is written in file *.*                                  "
DATA " format of file                                                          "
DATA ""                                                        
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA "  .    .     .                                                           "
DATA "  .    .     .    .  .   .                                               "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA ""                                                        
DATA "*************************************************************************"
DIM text$(400), xm#(30), d2$(30)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(a$, INSTR(a$, " "))): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))
jj% = ASC(LEFT$(a$, 1)) MOD 48: a$ = LTRIM$(RIGHT$(a$, LEN(a$) - 1))'get column
factor = VAL(a$)  'get factor

999 CALL analizecommand(filename$)

IF RIGHT$(filename$, 4) = ".rcp" OR RIGHT$(filename$, 4) = ".mrc" THEN
   SHELL "echo you should never ever change data in *.rcp/*.mrc files": PLAY "dgdgdg"
222 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 222
END IF
'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'open input file and input file header
OPEN "i", 1, filename$:


' open output file and write fileheader
OPEN "o", 2, "fact.fac"

22 CALL headerinput(text$(), j, 1)
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " column"; jj%;
PRINT #2, "multiplied by "; factor; " using program fact.bas}"


 REM input data columns
WHILE EOF(1) = 0
 CALL Inputline(1, xm#(), col%, d2$(), col2%)
 IF col% = -1 GOTO 22

'factorize column jj%
 xm#(jj%) = xm#(jj%) * factor

 'write result to file
 FOR coll% = 1 TO col%: PRINT #2, xm#(coll%); : NEXT
 FOR coll% = 1 TO col2%: PRINT #2, " {" + d2$(coll%) + "} "; : NEXT: PRINT #2,
WEND

 CLOSE 1, 2
SHELL "copy fact.fac " + filename$
SHELL "del fact.fac"
PRINT
SHELL "echo END FACT in file " + filename$ + " column " + STR$(jj%) + " has been multiplied by " + STR$(factor)

IF INSTR(COMMAND$, "*") <> 0 THEN GOTO 999
END
333 FOR i = 1 TO 17: READ a$: PRINT a$: NEXT i: END

SUB analizecommand (file1$)
STATIC washere%
'*****************************************************************
' this sub analizes the command$ and detects the filename
'*****************************************************************

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
'**********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'**********************************************************************

1 INPUT #n, a$
   i = INSTR(a$, "{"): IF i > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(a$, LEN(a$) - i)
  j = 1: i = INSTR(text$(j), "}"): IF i > 0 GOTO 3  'look for "}" in first line
                                           
   FOR j = 2 TO 400
   INPUT #n, text$(j)
   i = INSTR(text$(j), "}"): IF i > 0 GOTO 3      'look for "}"
   NEXT j: SHELL "echo text in data file too long": END
3 text$(j) = LEFT$(text$(j), i - 1)

END SUB

SUB Inputline (n, D1#(), col1%, d2$(), col2%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in D1#(1...col1%)
'comments in {} are stored in d2$
'if a comment is started somewhere in this line by "{" and not finished
'then col1% is set -1 and the filepointer is set to the beginning of the
'line (with seek)

a$ = INKEY$: IF a$ <> "" THEN IF ASC(a$) = 27 THEN END


aa = SEEK(n)
LINE INPUT #n, ala$
WHILE INSTR(ala$, CHR$(9)) > 0  'abandon tabs
 i% = INSTR(ala$, CHR$(9))
 ala$ = LEFT$(ala$, i% - 1) + " " + MID$(ala$, i% + 1)
WEND
'treat comments in input line
klauf% = INSTR(ala$, "{")
klzu% = INSTR(ala$, "}")
col2% = 0
WHILE klauf% < klzu% AND klauf% > 0   'take out closed bracketed expressions
 col2% = col2% + 1
 d2$(col2%) = MID$(ala$, klauf% + 1, klzu% - klauf% - 1)
 ala$ = LEFT$(ala$, klauf% - 1) + " " + MID$(ala$, klzu% + 1)
 klauf% = INSTR(ala$, "{")
 klzu% = INSTR(ala$, "}")
WEND


IF klauf% > 0 THEN
 col1% = -1: SEEK n, aa 'a comment bracket is not closed ... no data read
ELSE
 col1% = 0
 WHILE LEN(ala$) > 0
    col1% = col1% + 1
    ala$ = LTRIM$(ala$) + " "
    D1#(col1%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
    ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " ")))
 WEND
END IF

END SUB

