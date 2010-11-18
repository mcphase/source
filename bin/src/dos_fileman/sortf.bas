DECLARE SUB analizecommand (file1$)
DECLARE SUB recinput (n%, reclen%, startofrec%, x#(), col%)
DECLARE SUB lineinput (n!, x#(), col%)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "SORTF SORTF SORTF SORTF SORTF SORTF SORTF SORTF SORTF SORTF"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA " program to sort a file according to ascending values                    "
DATA " command:  sortf *.* 2 (2 means second column)                           "
DATA " format of file                                                          "
DATA "                                                                         "
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA "  .    .     .                                                           "
DATA "  .    .     .     .   .   .                                             "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "                                                                         "
DATA "*************************************************************************"
DIM text$(300), xx#(30), x#(30)

'aaaaaaaaaaaaaaaa ANALYZE COMMAND aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = COMMAND$: a$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(a$, LEN(a$) - 1))
ii% = ASC(RIGHT$(a$, 1)) MOD 48

999 CALL analizecommand(filename$)

' check if file is *.rcp measuring file and warn if it is
IF LCASE$(RIGHT$(filename$, 4)) = ".rcp" THEN
   PRINT "you should nevger change data in *.rcp files": PLAY "dgdgdg"
22 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 22
END IF
'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'rrrrrrrrrrrrrrrrr REFORMAT FILE rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1) 'open file and input file header
 numberofrec% = 0
OPEN "o", 2, "sortf.sor"
PRINT #2, "{"; : FOR i = 1 TO j - 1: PRINT #2, text$(i): NEXT: PRINT #2, text$(j) + "}"
WHILE EOF(1) = 0
INPUT #1, ala$        'input data point line as string and split into numbers
          col% = 0: WHILE LEN(ala$) > 0: ala$ = LTRIM$(ala$) + " ": col% = col% + 1
x#(col%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
          ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " "))): WEND
numberofrec% = numberofrec% + 1
IF numberofrec% = 1 THEN reclen% = col% * 18 + 2: coll% = col%
FOR i = 1 TO coll%: PRINT #2, USING "##.##########^^^^ "; x#(i); : NEXT: PRINT #2,
WEND
CLOSE 1, 2
'rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr

OPEN "i", 2, "sortf.sor": CALL headerinput(text$(), j, 2) 'open file and input file header
startofrec% = SEEK(2)

' sort according to column ii% sssssssssssssssssssssssssssssssssssssssssssss
PRINT "sorting ... "

 FOR n% = 1 TO numberofrec% - 1: PRINT USING "###%"; SQR(n% / numberofrec%) * 100;
  CALL recinput(n%, reclen%, startofrec%, x#(), col%)
  FOR nn% = n% + 1 TO numberofrec%
  CALL recinput(nn%, reclen%, startofrec%, xx#(), col%)
  IF x#(ii%) > xx#(ii%) THEN
                                'swap nn% and n% record
   CLOSE 2: OPEN "a", 1, "sortf.sor"  'open formatted file to sort it
   SEEK #1, startofrec% + reclen% * 1# * (n% - 1)
   FOR i = 1 TO coll%: PRINT #1, USING "##.##########^^^^ "; xx#(i); : NEXT: PRINT #1,
   SEEK #1, startofrec% + reclen% * 1# * (nn% - 1)
   FOR i = 1 TO coll%: PRINT #1, USING "##.##########^^^^ "; x#(i); : x#(i) = xx#(i): NEXT: PRINT #1,
   CLOSE 1: OPEN "i", 2, "sortf.sor"
  END IF
  NEXT nn%
 NEXT n%
'sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
CLOSE 2

'  create sorted output file ccccccccccccccccccccccccccccccccccccccccccc
SHELL "copy sortf.sor " + filename$
SHELL "del sortf.sor"
PRINT : PRINT "END SORTF file "; filename$; " sorted according to ascending column"; ii%
'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END
333 FOR i = 1 TO 16: READ a$: PRINT a$: NEXT: END

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

SUB recinput (n%, reclen%, startofrec%, x#(), col%)
'************************************************************************
' this sub is to read a record number n% and split into col% columns and store
' in x#(1...col%)
'************************************************************************
SEEK #2, startofrec% + reclen% * 1# * (n% - 1)'goto record in file
INPUT #2, ala$
        'input data point line as string and split into numbers
         
          col% = 0: WHILE LEN(ala$) > 0: ala$ = LTRIM$(ala$) + " ": col% = col% + 1
x#(col%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
          ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " "))): WEND
END SUB

