DECLARE SUB inputline (n!, d1#(), col1%, d2$(), colt%)
DECLARE SUB gobackonestrg (n!, ala$)
DECLARE SUB analizecommand (file1$, file2$, x1%, y1%, X2%, y2%)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "ADD ADD ADD ADD ADD ADD ADD ADD ADD ADD ADD ADD ADD ADD ADD ADD"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA " this programm is used to add functions y1(x1) and y2(x2)                "
DATA " command format: add *1.* 13 *2.* 21                                     "
DATA "     means function y1(x1)=col3(col1) in file *1.* is added              "
DATA "      to function y2(x2)=col1(col2) in file *2.* ...                     "
DATA "     the result is written into the col3 of file *1.*                    "
DATA " linear interpolation is used to match the x1 and x2 axis !!!!!          "
DATA "options:  /noex0   ... do not extrapolate and delete datapoints "                                         
DATA "          /noex1   ... do not extrapolate and add nothing to datapoints"                                 
DATA " file formats                                                            "
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns     "
DATA "   time[s]  T(or R)  capacitance[pF] }                                   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA "  .    .     .      .  .  .                                              "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA ""
DATA "*************************************************************************"
DIM text1$(300), text2$(300), d1#(30), D2A#(30), d2b#(30), d2$(10)

'ooooooooooooooooooooooooo OPEN FILES oooooooooooooooooooooooooooooooooooo
999 CALL analizecommand(file1$, file2$, x1%, y1%, X2%, y2%)

OPEN "i", 1, file1$: CALL headerinput(text1$(), j1, 1)' open file1 and input header
OPEN "i", 2, file2$: CALL headerinput(text2$(), j2, 2)' open file2 and input header
OPEN "o", 3, "add.add"                         ' open output file

PRINT #3, "{";                    'print file header of output file
FOR i = 1 TO j1: PRINT #3, text1$(i): NEXT
PRINT #3, DATE$; " "; TIME$; " Column"; y1%; " has been added to column "; y2%; "of file "; file2$
PRINT #3, " the common axis was "; x1%; "in file "; file1$; " and "; X2%; " in "; file2$
PRINT #3, "   (((file "; file2$; " information: ";
FOR i = 1 TO j2: PRINT #3, LCASE$(text2$(i)): NEXT
PRINT #3, "   ))) addition done by add.bas} "
'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

'******************************************************************
REM input columns of #1 and #2 and write result to #3
'******************************************************************

    'input a data point on #1
IF EOF(1) <> 0 GOTO 5 ELSE CALL inputline(1, d1#(), col1%, d2$(), colt%)
                                 
   'input first 2 values on file 2
IF EOF(2) <> 0 GOTO 5 ELSE CALL inputline(2, D2A#(), col2%, d2$(), colt%)
IF EOF(2) <> 0 GOTO 5 ELSE CALL inputline(2, d2b#(), col2%, d2$(), colt%)

' decide wether values in file 2 increase or decrease
sign% = SGN(d2b#(X2%) - D2A#(X2%))

' input d1 values until d1#(x1%)*sign% is bigger than d2a#(x2%)*sign%
'WHILE d1#(x1%) * sign% < D2A#(X2%) * sign%
'      IF EOF(1) <> 0 GOTO 5 ELSE CALL inputline(1, d1#(), col1%, d2$(), colt%)
'WEND

' \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \
'> > > > > > > > > > > > > >ADDING LOOP> > > > > > > > > > > > > > > > > >
' / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
PRINT "adding": PRINT
4 LOCATE 23, 1: PRINT USING "###%"; 100 * SEEK(1) / LOF(1); : LOCATE 24, 1
'-----------------------------------------------------------------------
'input d2a- d2b values going backward in file until d1 is in the interval
      WHILE d1#(x1%) * sign% < D2A#(X2%) * sign%
                
                 ' go back three strings in file2
  FOR cr% = 1 TO 3: CALL gobackonestrg(2, ala$)
  NEXT cr%
  IF INSTR(ala$, "}") <> 0 THEN
   sign% = SGN(d2b#(X2%) - D2A#(X2%))
   INPUT #2, a$: INPUT #2, a$: INPUT #2, a$
   IF INSTR(LCASE$(COMMAND$), "/noex") = 0 THEN
     PRINT "out of range datapoint x="; d1#(x1%); "linear extrapolated": GOTO 15
   ELSE
     IF INSTR(LCASE$(COMMAND$), "/noex0") > 0 THEN
      PRINT "out of range datapoint x="; d1#(x1%); "skipped and removed": GOTO 16
     ELSE
      PRINT "out of range datapoint x="; d1#(x1%); "skipped and unchanged": GOTO 17
     END IF
   END IF
  END IF
                              'we are at the beginning of the file
                             
 CALL inputline(2, D2A#(), col2%, d2$(), colt%): CALL inputline(2, d2b#(), col2%, d2$(), colt%)
                        ' input d2a-d2b interval of before
        WEND
'------------------------------------------------------------------------

'........................................................................
'input d2b values until d1#(x1%) is between d2a#(x2%) and d2b#(x2%)
WHILE d2b#(X2%) * sign% < d1#(x1%) * sign%

  IF EOF(2) <> 0 THEN
   sign% = SGN(d2b#(X2%) - D2A#(X2%))
   IF INSTR(LCASE$(COMMAND$), "/noex") = 0 THEN
    PRINT "out of range datapoint x="; d1#(x1%); "linear extrapolated": GOTO 15
   ELSE
     IF INSTR(LCASE$(COMMAND$), "/noex0") > 0 THEN
      PRINT "out of range datapoint x="; d1#(x1%); "skipped and removed": GOTO 16
     ELSE
      PRINT "out of range datapoint x="; d1#(x1%); "skipped and unchanged": GOTO 17
     END IF
   END IF
  END IF
           'if d1 exceeds d2-range get next d1 value
  FOR coll% = 1 TO col2%: D2A#(coll%) = d2b#(coll%): NEXT
  CALL inputline(2, d2b#(), col2%, d2$(), colt%)
WEND
'........................................................................

'!!!!!!!!!!!!!!!!!!!!! DO ADDITION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
' linear interpolate d2a and d2b
15 interpoly2 = D2A#(y2%) + (d2b#(y2%) - D2A#(y2%)) / (d2b#(X2%) - D2A#(X2%)) * (d1#(x1%) - D2A#(X2%))

' DO addition
d1#(y1%) = d1#(y1%) + interpoly2

'save the datapoint
17 FOR coll% = 1 TO col1%:
nn$ = STR$(d1#(coll%)): IF INSTR(nn$, "D") > 0 THEN MID$(nn$, INSTR(nn$, "D"), 1) = "E"
PRINT #3, " " + nn$; : NEXT: PRINT #3,
'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

16 ' input next d1 value
IF EOF(1) <> 0 THEN
 GOTO 5
ELSE
21  CALL inputline(1, d1#(), col1%, d2$(), colt%)
 IF col1% < 0 AND EOF(1) = 0 THEN
  CALL headerinput(text2$(), j1, 1)' open file2 and input header
  PRINT #3, "{";                    'print file header of output file
    FOR i = 1 TO j1: PRINT #3, text1$(i): NEXT
    PRINT #3, "}"
  GOTO 21
  END IF
END IF

GOTO 4
' / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
'< < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
' \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \

5 CLOSE 1, 2, 3
SHELL "copy add.add " + file1$
SHELL "del add.add"
PRINT
PRINT "END ADD col"; y2%; " of file "; file2$; " has been successfully added to";
PRINT "col"; y1%; " of file "; file1$: PRINT

IF INSTR(COMMAND$, "*") <> 0 GOTO 999

END
333 FOR i = 1 TO 21: READ a$: PRINT a$: NEXT: END

SUB analizecommand (file1$, file2$, x1%, y1%, X2%, y2%)
STATIC washere%
'*****************************************************************
' this sub analizes the command$ and detects the two filenames and
' 2 columns for each file
'*****************************************************************

a$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))

'detect file 1:
i = INSTR(a$, " "): file1$ = RTRIM$(LEFT$(a$, i)): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - i))
x1% = ASC(LEFT$(a$, 1)) MOD 48: y1% = ASC(MID$(a$, 2, 1)) MOD 48: a$ = LTRIM$(RIGHT$(a$, LEN(a$) - 2))
'detect file 2:
i = INSTR(a$, " "): file2$ = RTRIM$(LEFT$(a$, i)): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - i))
X2% = ASC(LEFT$(a$, 1)) MOD 48: y2% = ASC(MID$(a$, 2, 1)) MOD 48: a$ = LTRIM$(RIGHT$(a$, LEN(a$) - 2))

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



' check if file is *.rcp measuring file and warn if it is
IF LCASE$(RIGHT$(file1$, 4)) = ".rcp" OR LCASE$(RIGHT$(file1$, 4)) = ".mrc" THEN
   PRINT "you should never ever change data in *.rcp or *.mrc files": PLAY "dgdgdg"
22 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 22
END IF

END SUB

SUB gobackonestrg (n, ala$)
'this sub sets the filepointer of file #n one string back and
' puts the corresponding string into ala$
filepointer = SEEK(n)
FOR x = 2 TO filepointer STEP 30
SEEK #n, filepointer - x
fpn = -1: INPUT #n, a$
122 WHILE SEEK(n) < filepointer
       fpn = SEEK(n):  INPUT #n, a$
    WEND
IF fpn > 0 GOTO 233
NEXT x: x = filepointer - 1: SEEK #n, 1: GOTO 122

233 SEEK #n, fpn: ala$ = a$
END SUB

SUB headerinput (text$(), j, n)
'*****************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'*****************************************************************

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
'
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
a$ = INKEY$: IF a$ <> "" THEN IF ASC(a$) = 27 THEN END
IF SCREEN(24, 10) <> 45 THEN PRINT "-";  ELSE LOCATE 24, 1: PRINT SPACE$(15); : LOCATE 24, 1

'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
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
    d1#(col1%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
    ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " ")))
 WEND
END IF

END SUB

