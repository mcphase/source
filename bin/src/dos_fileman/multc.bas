DECLARE SUB inputline (n!, d1#(), col1%, d2$(), col2%)
DECLARE SUB analizecommand (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "MULTC MULTC MULTC MULTC MULTC MULTC MULTC MULTC MULTC MULTC MULTC MULTC MULTC"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program MULTC - use it like MULTC *.* 23"                     
DATA "    (means multiply col 2 and 3 point by point) -                        "
DATA " ---> the result is written in file *.* column 2                         "
DATA " format of file                                                          "
DATA ""                                                                   
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA "  .    .     .     .  .  .                                               "
DATA "  .    .     .                                                           "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA ""                                                                   
DATA "*************************************************************************"
DIM text$(300), xm#(30), d2$(30)
' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(a$, LEN(a$) - 2))
jj% = ASC(RIGHT$(a$, 1)) MOD 48: ii% = ASC(MID$(a$, LEN(a$) - 1, 1)) MOD 48

999 CALL analizecommand(filename$)

IF RIGHT$(filename$, 4) = ".rcp" THEN
   PRINT "you should never ever change data in *.rcp files": PLAY "dgdgdg"
22 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 22
END IF
'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'oooooooooooooooooooo OPEN FILES ooooooooooooooooooooooooooooooooooooooooo
'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

' open output file and write fileheader
OPEN "o", 2, "MULTC.mlc"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " col "; ii%; "and col"; jj%; "have been multiplied";
PRINT #2, " using program MULTC.bas}"
'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

REM input data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
WHILE EOF(1) = 0

 CALL inputline(1, xm#(), col%, d2$(), col2%)

 'mul and write result to file
 xm#(ii%) = xm#(jj%) * xm#(ii%)
 FOR coll% = 1 TO col%: PRINT #2, xm#(coll%); : NEXT:
 FOR coll% = 1 TO col2%: PRINT #2, " {"; d2$(coll%); "} "; : NEXT:
 PRINT #2,
WEND'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 CLOSE 1, 2
SHELL "copy MULTC.mlc " + filename$
SHELL "del MULTC.mlc"
PRINT
PRINT "END MULTC in file "; filename$; " column "; jj%; "and "; ii%; " have been multiplied"
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END
333 FOR i = 1 TO 17: READ a$: PRINT a$: NEXT i: END

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
   i = INSTR(a$, "{"): IF i > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(a$, LEN(a$) - i)
  j = 1: i = INSTR(text$(j), "}"): IF i > 0 GOTO 3  'look for "}" in first line
                                      
   FOR j = 2 TO 300
   INPUT #n, text$(j)
   i = INSTR(text$(j), "}"): IF i > 0 GOTO 3      'look for "}"
   NEXT j: PRINT "text in data file too long": END
3 text$(j) = LEFT$(text$(j), i - 1)

END SUB

SUB inputline (n, d1#(), col1%, d2$(), col2%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
'comments in {} are stored in d2$
'if a comment is started somewhere in this line by "{" and not finished
'then col1% is set -1 and the filepointer is set to the beginning of the
'line (with seek)

a$ = INKEY$: IF a$ <> "" THEN IF ASC(a$) = 27 THEN END
IF SCREEN(24, 10) <> 45 THEN PRINT "-";  ELSE LOCATE 24, 1: PRINT SPACE$(15); : LOCATE 24, 1

aa = SEEK(n)
LINE INPUT #n, ala$
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
    d1#(col1%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
    ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " ")))
 WEND
END IF

END SUB

