DECLARE SUB headerinput (text$(), j!, n!)
PRINT "SWAPF SWAPF SWAPF SWAPF SWAPF SWAPF SWAPF SWAPF SWAPF SWAPF SWAPF "
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program swapf - use it like 'SWAPF *.* *.*'                            "
DATA "    (means swap col 1 of two files of the same length) -                 "
DATA " format of file                                                          "
DATA "                                                                         "
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA "  .    .     .         .  .  .                                           "
DATA "  .    .     .                                                           "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "                                                                         "
DATA "*************************************************************************"
DIM text1$(300), text2$(300), x#(30), y#(30)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
i = INSTR(a$, " "): file1$ = LTRIM$(LEFT$(a$, i))
                    file2$ = LTRIM$(RIGHT$(a$, LEN(a$) - i))

IF RIGHT$(file1$, 4) = ".rcp" THEN
   PRINT "you should never ever change data in *.rcp files": PLAY "dgdgdg"
23 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 23
END IF
IF LCASE$(RIGHT$(file2$, 4)) = ".rcp" THEN
   PRINT "you should never ever change data in *.rcp files": PLAY "dgdgdg"
22 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 22
END IF
'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'ooooooooooooooooooooo OPEN FILES oooooooooooooooooooooooooooooooooooooooo
'open input file1 and input file header
OPEN "i", 1, file1$: CALL headerinput(text1$(), j1, 1)

' open output file1 and write fileheader
OPEN "o", 2, "swapf.sw1"
PRINT #2, "{"; : FOR iii = 1 TO j1: PRINT #2, text1$(iii): NEXT
PRINT #2, DATE$; " "; TIME$; " col1 "; file1$; "and col1"; file2$; "have been swapped";
PRINT #2, " using program swapf.bas}"

'open input file2 and input file header
OPEN "i", 3, file2$: CALL headerinput(text2$(), j2, 3)

' open output file2 and write fileheader
OPEN "o", 4, "swapf.sw2"
PRINT #4, "{"; : FOR iii = 1 TO j2: PRINT #4, text2$(iii): NEXT
PRINT #4, DATE$; " "; TIME$; " col1 "; file1$; "and col1"; file2$; "have been swapped";
PRINT #4, " using program swapf.bas}"
'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


REM input data columns >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
WHILE EOF(1) = 0 AND EOF(3) = 0
 INPUT #1, ala$        'input data point line as string and split into numbers
          col1% = 0: WHILE LEN(ala$) > 0: ala$ = LTRIM$(ala$) + " ": col1% = col1% + 1
 x#(col1%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
          ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " "))): WEND

 INPUT #3, ala$        'input data point line as string and split into numbers
          col2% = 0: WHILE LEN(ala$) > 0: ala$ = LTRIM$(ala$) + " ": col2% = col2% + 1
 y#(col2%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
          ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " "))): WEND

 'write result to file
 PRINT #2, y#(1); : FOR coll% = 2 TO col1%: PRINT #2, x#(coll%); : NEXT: PRINT #2,
 PRINT #4, x#(1); : FOR coll% = 2 TO col2%: PRINT #4, y#(coll%); : NEXT: PRINT #4,
WEND '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 CLOSE
SHELL "copy swapf.sw1 " + file1$: SHELL "del swapf.sw1"
SHELL "copy swapf.sw2 " + file2$: SHELL "del swapf.sw2"
PRINT
PRINT "END SWAPF files "; file1$; " and "; file2$; "1st columns have been swapped"
END
333 FOR i = 1 TO 16: READ a$: PRINT a$: NEXT i: END

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

