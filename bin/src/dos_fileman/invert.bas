DECLARE SUB analizecommand (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
DECLARE SUB pointdraw3d (x!, xmin!, xmax!, y!, ymin!, ymax!, z!, zmin!, zmax!, size!)
PRINT "INVERT INVERT INVERT INVERT INVERT INVERT INVERT INVERT INVERT"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program invert - use it like invert *.* 2"                        
DATA "    (2 means invert second column) -                                     "
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
DIM text$(300), xm#(30)

' analyse command string (get filename$) aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
A$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(A$, LEN(A$) - 1))
jj% = ASC(RIGHT$(A$, 1)) MOD 48  'get column

999 CALL analizecommand(filename$)

IF RIGHT$(filename$, 4) = ".rcp" THEN
   PRINT "you should never ever change data in *.rcp files": PLAY "dgdgdg"
22 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 22
END IF
'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'ooooooooooooooooooooo OPEN FILES ooooooooooooooooooooooooooooooooooooooo
'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

' open output file and write fileheader
OPEN "o", 2, "invert.inv"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " column"; jj%;
PRINT #2, "inverted by program invert.bas}"
'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

 REM input data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
WHILE EOF(1) = 0
 INPUT #1, ala$        'input data point line as string and split into numbers
            col% = 0: WHILE LEN(ala$) > 0: ala$ = LTRIM$(ala$) + " ": col% = col% + 1
 xm#(col%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
            ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " "))): WEND

 'invert columnd jj%
 IF xm#(jj%) = 0 THEN PRINT "INVERT division by zero at data point ": FOR alaa% = 1 TO col%: PRINT xm#(alaa%); : NEXT: PRINT "data point skipped": GOTO 12
 xm#(jj%) = 1 / xm#(jj%)

 'write result to file
 FOR alaa% = 1 TO col%: PRINT #2, xm#(alaa%); : NEXT: PRINT #2,
12 WEND '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 CLOSE 1, 2
SHELL "copy invert.inv " + filename$
SHELL "del invert.inv"
PRINT
PRINT "END INVERT in file "; filename$; " column "; jj%; "has been inverted"
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END

333 FOR I = 1 TO 17: READ A$: PRINT A$: NEXT I: END

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

