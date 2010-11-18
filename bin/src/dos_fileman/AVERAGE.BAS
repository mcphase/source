DECLARE SUB analizecommand (file1$)
DECLARE SUB inputline (n!, d1#(), col1%)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE AVERAGE"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*********************************************************************"
DATA "this program is designed to AVERAGE data by deleting close data points"
DATA "use as:
DATA " AVERAGE *.* 4 [options] .... takes 4 points and averages data
DATA "			       option
DATA "			                middle point is taken
DATA "			       /av      points are averaged
DATA "			       /med     median of points is calculated and kept
DATA " format of file:                                                          "
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "*************************************************************************"
DIM text$(300), x#(30, 100), xm#(30), in#(30)

'aaaaa analyse command string (get filename$ - ii%) aaaaaa
A$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))

filename$ = RTRIM$(LEFT$(A$, INSTR(A$, " "))): A$ = LTRIM$(RIGHT$(A$, LEN(A$) - INSTR(A$, " ")))
ii% = ASC(LEFT$(A$, 1)) MOD 48               'get number of points

999 CALL analizecommand(filename$)

IF RIGHT$(filename$, 4) = ".rcp" OR RIGHT$(filename$, 4) = ".mrc" THEN
   PRINT "you should never ever change data in *.rcp/mrc files": PLAY "dgdgdg"
22 INPUT "do you really want to continue (Y/N)"; ALA$: IF LCASE$(ALA$) = "n" THEN END
   IF LCASE$(ALA$) <> "y" GOTO 22
END IF
'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'oooooooooooooooooo OPEN FILE oooooooooooooooooooooooooooooooooooooooooooooo
'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

' open output file and write fileheader
OPEN "o", 2, "AVERAGE.red"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " "; ii%; "points AVERAGED  "; in; col; "; ii%;"
PRINT #2, "by program AVERAGE.bas "; COMMAND$; "}"
'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


REM input data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
' x#(ii%, 0) = 1E+30: x#(ii%, 1) = 1E+30

'input column y until abs(x#(ii%,i)-x#(ii%,1))>dmin iiiiiiiiiiii
4 FOR i = 1 TO ii%
     IF EOF(1) <> 0 GOTO 5
 IF i > 99 THEN PRINT "ERROR - number of points too big - program AVERAGE terminated !!!!": END
 CALL inputline(1, in#(), col%): FOR coll% = 1 TO col%: x#(coll%, i) = in#(coll%): NEXT coll%
 NEXT i
'iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

' AVERAGE and write to output file  ddddddddddddddddddddddddddd
FOR coll% = 1 TO col%
 'if no options just take a point in about the middle of the interval
  xm#(coll%) = x#(coll%, INT(ii% / 2) + 1)
 
 'if option /av just take the average of the interval dmin
  IF INSTR(A$, "/av") > 0 THEN
   xm#(coll%) = 0: FOR k% = 1 TO ii%
   xm#(coll%) = xm#(coll%) + x#(coll%, k%): NEXT
   xm#(coll%) = xm#(coll%) / ii%
  END IF

 'if option /med just take the median of the interval dmin
  IF INSTR(A$, "/med") > 0 THEN
  
   i = ii%
   WHILE i > 2
    max = -1E+10: min = 1E+10
    FOR k% = 1 TO ii%
     IF x#(coll%, k%) < min AND x#(coll%, k%) > -1E+10 THEN min = x#(coll%, k%): x#(coll%, k%) = -1E+10
     IF x#(coll%, k%) > max AND x#(coll%, k%) < 1E+10 THEN max = x#(coll%, k%): x#(coll%, k%) = 1E+10
    NEXT
    i = i - 2
   WEND
   
    FOR k% = 1 TO ii%
     IF x#(coll%, k%) > -1E+10 AND x#(coll%, k%) < 1E+10 THEN xm#(coll%) = x#(coll%, k%)
    NEXT
  END IF

NEXT coll%

'write result to file
 FOR coll% = 1 TO col%
nn$ = STR$(xm#(coll%)): MID$(nn$, INSTR(nn$, "D"), 1) = "E"
PRINT #2, " " + nn$; : NEXT: PRINT #2,
'ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
                       
GOTO 4 '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

5 CLOSE 1, 2
SHELL "copy AVERAGE.red " + filename$
SHELL "del AVERAGE.red"
PRINT
PRINT "END AVERAGE in file "; filename$;
PRINT " column"; ii%; "AVERAGED to minimal step "; dmin; " in col "; ii%;
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END
333 FOR i = 1 TO 16: READ A$: PRINT A$: NEXT i: END

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
   IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del fact.dir": END
  END IF


END SUB

SUB headerinput (text$(), j, n)
'*********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'*********************************************************************

1 INPUT #n, A$
   i = INSTR(A$, "{"): IF i > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(A$, LEN(A$) - i)
  j = 1: i = INSTR(text$(j), "}"): IF i > 0 GOTO 3  'look for "}" in first line
                                            
   FOR j = 2 TO 300
   INPUT #n, text$(j)
   i = INSTR(text$(j), "}"): IF i > 0 GOTO 3      'look for "}"
   NEXT j: PRINT "text in data file too long": END
3 text$(j) = LEFT$(text$(j), i - 1)

END SUB

SUB inputline (n, d1#(), col1%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
A$ = INKEY$: IF A$ <> "" THEN IF ASC(A$) = 27 THEN END

INPUT #n, ALA$
col1% = 0
 WHILE LEN(ALA$) > 0
    col1% = col1% + 1
    ALA$ = LTRIM$(ALA$) + " "
    d1#(col1%) = VAL(LEFT$(ALA$, INSTR(ALA$, " ")))
    ALA$ = LTRIM$(RIGHT$(ALA$, LEN(ALA$) - INSTR(ALA$, " ")))
 WEND

END SUB

