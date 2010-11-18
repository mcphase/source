DECLARE SUB gobackonestrg (n!, ALA$)
DECLARE SUB inputline (n!, d1#(), col1%)
DECLARE SUB iterate (file1$, file2$, x1%, y1%, X2%, y2%)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "MULT MULT MULT MULT MULT MULT MULT MULT MULT MULT MULT MULT MULT"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA " this programm is intended to multiply functions y1(x1) and y2(x2)       "
DATA " command format: mult *1.* 13 *2.* 21                                    "
DATA "     means function y1(x1)=col3(col1) in file *1.* is multiplied         "
DATA "      with function y2(x2)=col1(col2) in file *2.* ...                   "
DATA "     the result is written into the col3 of file *1.*                    "
DATA " linear interpolation is used to match the x1 and x2 axis !!!!!          "
DATA ""                                                         
DATA " file formats                                                           "
DATA ""
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns     "
DATA "   time[s]  T(or R)  capacitance[pF] }                                   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA "  .    .     .       .   .   .                                           "
DATA "  .    .     .                                                           "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "                                                                      "
DATA "*************************************************************************"
DIM text1$(300), text2$(300), d1#(30), D2A#(30), d2b#(30)

'ooooooooooooooooooooooooo OPEN FILES oooooooooooooooooooooooooooooooooooo
999 CALL iterate(file1$, file2$, x1%, y1%, X2%, y2%)

OPEN "i", 1, file1$: CALL headerinput(text1$(), j1, 1)' open file1 and input header
OPEN "i", 2, file2$: CALL headerinput(text2$(), j2, 2)' open file2 and input header

OPEN "o", 3, "mult.mul"                        ' open output file
PRINT #3, "{";                    'print file header of output file
FOR i = 1 TO j1: PRINT #3, text1$(i): NEXT
PRINT #3, DATE$; " "; TIME$; " Column"; y1%; " has been multiplied with column "; y2%; "of file "; file2$
PRINT #3, " the common axis was "; x1%; "in file "; file1$; " and "; X2%; " in "; file2$
PRINT #3, "   (((file "; file2$; " information: ";
FOR i = 1 TO j2: PRINT #3, LCASE$(text2$(i)): NEXT
PRINT #3, "   ))) multiplication done by mult.bas} "
'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

'******************************************************************
REM input columns of #1 and #2 and write result to #3
'******************************************************************

    'input a data point on #1
IF EOF(1) <> 0 GOTO 5 ELSE CALL inputline(1, d1#(), col1%)

   'input first 2 values on file 2
IF EOF(2) <> 0 GOTO 5 ELSE CALL inputline(2, D2A#(), col2%)
IF EOF(2) <> 0 GOTO 5 ELSE CALL inputline(2, d2b#(), col2%)

' decide wether values in file 2 increase or decrease
sign% = SGN(d2b#(X2%) - D2A#(X2%))

' input d1 values until d1#(x1%)*sign% is bigger than d2a#(x2%)*sign%
WHILE d1#(x1%) * sign% < D2A#(X2%) * sign%
      IF EOF(1) <> 0 GOTO 5 ELSE CALL inputline(1, d1#(), col1%)
WEND

' \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \
'> > > > > > > > > > > > > >ADDING LOOP> > > > > > > > > > > > > > > > > >
' / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
PRINT "multiplying": PRINT
4 LOCATE 23, 1: PRINT USING "###%"; 100 * SEEK(1) / LOF(1); : LOCATE 24, 1
'-----------------------------------------------------------------------
'input d2a- d2b values going backward in file until d1 is in the interval
      WHILE d1#(x1%) * sign% < D2A#(X2%) * sign%
                
                 ' go back three strings in file2
  FOR cr% = 1 TO 3: CALL gobackonestrg(2, ALA$)
  NEXT cr%
   IF INSTR(ALA$, "}") <> 0 THEN INPUT #2, A$: INPUT #2, A$: INPUT #2, A$: PRINT "out of range datapoint x="; d1#(x1%); "linear extrapolated": GOTO 15
                              'we are at the beginning of the file
                             
 CALL inputline(2, D2A#(), col2%): CALL inputline(2, d2b#(), col2%)
                        ' input d2a-d2b interval of before
        WEND
'------------------------------------------------------------------------

'........................................................................
'input d2b values until d1#(x1%) is between d2a#(x2%) and d2b#(x2%)
WHILE d2b#(X2%) * sign% < d1#(x1%) * sign%

  IF EOF(2) <> 0 THEN PRINT "out of range datapoint x="; d1#(x1%); "linear extrapolated": GOTO 15
  'if d1 exceeds d2-range get next d1 value
  FOR coll% = 1 TO col2%: D2A#(coll%) = d2b#(coll%): NEXT
  CALL inputline(2, d2b#(), col2%)
WEND
'........................................................................

'!!!!!!!!!!!!!!!!!!!!!! DO MULTIPLICATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
' linear interpolate d2a and d2b
15 interpoly2 = D2A#(y2%) + (d2b#(y2%) - D2A#(y2%)) / (d2b#(X2%) - D2A#(X2%)) * (d1#(x1%) - D2A#(X2%))

' DO MULTIPLICATION
d1#(y1%) = d1#(y1%) * interpoly2

'save the datapoint
FOR coll% = 1 TO col1%: PRINT #3, d1#(coll%); : NEXT: PRINT #3,
'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

' input next d1 value
IF EOF(1) <> 0 GOTO 5 ELSE CALL inputline(1, d1#(), col1%)

GOTO 4
' / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
'< < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
' \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \

5 CLOSE 1, 2, 3

SHELL "copy mult.mul " + file1$
SHELL "del mult.mul"

PRINT
PRINT "END MULT col"; y1%; " of file "; file1$; " has been successfully multiplied with"
PRINT "col"; y2%; " of file "; file2$
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END

333 FOR i = 1 TO 22: READ A$: PRINT A$: NEXT: END

SUB gobackonestrg (n, ALA$)
'this sub sets the filepointer of file #n one string back and
' puts the corresponding string into ala$
filepointer = SEEK(n)
FOR x = 2 TO filepointer STEP 30
SEEK #n, filepointer - x
fpn = -1: INPUT #n, A$
122 WHILE SEEK(n) < filepointer
       fpn = SEEK(n):  INPUT #n, A$
    WEND
IF fpn > 0 GOTO 233
NEXT x: x = filepointer - 1: SEEK #n, 1: GOTO 122

233 SEEK #n, fpn: ALA$ = A$
END SUB

SUB headerinput (text$(), j, n)
'**********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'**********************************************************************

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
'*************************************************************
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
'*************************************************************
A$ = INKEY$: IF A$ <> "" THEN IF ASC(A$) = 27 THEN END
IF SCREEN(24, 10) <> 45 THEN PRINT "-";  ELSE LOCATE 24, 1: PRINT SPACE$(15); : LOCATE 24, 1

INPUT #n, ALA$
col1% = 0
 WHILE LEN(ALA$) > 0
    col1% = col1% + 1
    ALA$ = LTRIM$(ALA$) + " "
    d1#(col1%) = VAL(LEFT$(ALA$, INSTR(ALA$, " ")))
    ALA$ = LTRIM$(RIGHT$(ALA$, LEN(ALA$) - INSTR(ALA$, " ")))
 WEND
END SUB

SUB iterate (file1$, file2$, x1%, y1%, X2%, y2%)
STATIC washere1%, washere2%
'*****************************************************************
' this sub analizes the command$ and detects the two filenames and
' 2 columns for each file
'*****************************************************************

A$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))

'detect file 1:
i = INSTR(A$, " "): file1$ = LEFT$(A$, i): A$ = LTRIM$(RIGHT$(A$, LEN(A$) - i))
x1% = ASC(LEFT$(A$, 1)) MOD 48: y1% = ASC(MID$(A$, 2, 1)) MOD 48: A$ = LTRIM$(RIGHT$(A$, LEN(A$) - 2))
'detect file 2:
i = INSTR(A$, " "): file2$ = LEFT$(A$, i): A$ = LTRIM$(RIGHT$(A$, LEN(A$) - i))
X2% = ASC(LEFT$(A$, 1)) MOD 48: y2% = ASC(MID$(A$, 2, 1)) MOD 48: A$ = LTRIM$(RIGHT$(A$, LEN(A$) - 2))

IF INSTR(file1$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere1% = 0 THEN
      washere1% = 1
      'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /ON /b > file1.dir"
      OPEN "i", 9, "file1.dir"
   END IF
   IF EOF(9) <> 0 THEN CLOSE : SHELL "del file1.dir": END
   INPUT #9, file1$
   IF LCASE$(file1$) = "file1.dir" THEN
    IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del file1.dir": END
    INPUT #9, file1$
   END IF
   IF LTRIM$(file1$) = "" THEN CLOSE : SHELL "del file1.dir": END
  END IF

IF INSTR(file2$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere2% = 0 THEN
      washere2% = 1
      'get filenames to be operated on as file2$
      SHELL "dir " + file2$ + " /ON /b > file2.dir"
      OPEN "i", 8, "file2.dir"
   END IF
   IF EOF(8) <> 0 THEN CLOSE : SHELL "del file2.dir": END
   INPUT #8, file2$
   IF LTRIM$(file2$) = "" THEN CLOSE : SHELL "del file2.dir": END
  END IF



' check if file is *.rcp measuring file and warn if it is
IF LCASE$(RIGHT$(file1$, 4)) = ".mrc" OR LCASE$(RIGHT$(file1$, 4)) = ".rcp" THEN
   PRINT "you should never ever change data in *.rcp/mrc files": PLAY "dgdgdg"
22 INPUT "do you really want to continue (Y/N)"; ALA$: IF LCASE$(ALA$) = "n" THEN END
   IF LCASE$(ALA$) <> "y" GOTO 22
END IF

END SUB

