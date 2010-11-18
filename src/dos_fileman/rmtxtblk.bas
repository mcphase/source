DECLARE SUB countchar (t$, a$, cnt!)
DECLARE SUB detendofblock (i!, op%, t$)
DECLARE SUB blockinput (text$(), j!, n!, a$, e$)
DECLARE SUB inputline (n!, d1#(), col1%, d2$(), col2%)
DECLARE SUB checkrepeat (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "RMTXTBLK RMTXTBLK RMTXTBLK RMTXTBLK RMTXTBLK RMTXTBLK RMTXTBLK RMTXTBLK RMTXTBLK"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program RMTXTBLK - use it like RMTXTBLK *.* 'dsdad'                    "
DATA "    (remove {text blocks} from file *.* unless they contain dsdad "
DATA "    no case match is done"
DATA "format of file"
DATA " { header: this is the                                                  "
DATA "   file header delimited                                                "
DATA "   by brackets after this header there follow 3 or more data columns }  "
DATA " 11 3.14235 65367                                                       "
DATA "  .    .     .                                                          "
DATA " {text blocks}                                                          "
DATA "  .    .     .    .  .   .                                              "
DATA "  .    .     .                                                          "
DATA " 32 2412.34 324.2                                                       "
DATA ""
DATA "*************************************************************************"
DIM text$(300), xm#(30), d2$(30)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
filename$ = RTRIM$(LEFT$(a$, INSTR(a$, " "))): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))
aa% = INSTR(a$, "'"): bb% = INSTR(aa% + 1, a$, "'")
jj$ = LCASE$(MID$(a$, aa% + 1, bb% - aa% - 1))


999 CALL checkrepeat(filename$)

IF RIGHT$(filename$, 4) = ".rcp" OR RIGHT$(filename$, 4) = ".mrc" THEN
   PRINT "you should never ever change data in *.rcp/mrc files": PLAY "dgdgdg"
222 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 222
END IF
'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'open input file and input file header
OPEN "i", 1, filename$
' open output file
OPEN "o", 2, "RMTXTBLK.cut"

22 CALL blockinput(text$(), j, 1, a$, e$)

  IF j > 0 THEN
   PRINT #2, a$; "{"; : FOR iii = 1 TO j - 1: PRINT #2, text$(iii): NEXT
   PRINT #2, text$(j); "}";
   IF INSTR(e$, "{") = 0 THEN PRINT #2, e$
  ELSE
  PRINT #2, a$;
  IF INSTR(e$, "{") = 0 THEN PRINT #2, e$
  END IF
IF EOF(1) = 0 GOTO 22


23 CLOSE 1, 2
SHELL "copy RMTXTBLK.cut " + filename$
SHELL "del RMTXTBLK.cut"

PRINT
PRINT "END RMTXTBLK in file "; filename$; "all {text} without "; jj$; "has been removed"

IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END

333 FOR i = 1 TO 16: READ a$: PRINT a$: NEXT i: END

SUB blockinput (text$(), j, n, a$, e$)
SHARED jj$
'**********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'**********************************************************************

1 IF EOF(n) = 0 THEN ss = SEEK(n): LINE INPUT #n, a$ ELSE e$ = "": GOTO 4
   i = INSTR(a$, "{"): IF i > 0 GOTO 2 ELSE PRINT #2, a$: GOTO 1   'look for "{"
2 text$(1) = RIGHT$(a$, LEN(a$) - i): a$ = LEFT$(a$, i - 1)
  j = 1
CALL countchar(text$(j), "{", cntopen)
CALL countchar(text$(j), "}", cntclose)
oldlevel% = 0
level% = cntopen - cntclose
oo% = oldlevel%: CALL detendofblock(i, oo%, text$(j))
IF level% < -1 THEN PRINT "error - more } then { at "; text$(j): END
IF level% = -1 OR i < LEN(text$(j)) + 1 THEN
CALL detendofblock(i, oldlevel%, text$(j)): GOTO 3
ELSE
   FOR j = 2 TO 300
   ss = SEEK(n)
   LINE INPUT #n, text$(j)
  CALL countchar(text$(j), "{", cntopen)
  CALL countchar(text$(j), "}", cntclose)
  oldlevel% = level%
  level% = level% + cntopen - cntclose
  oo% = oldlevel%: CALL detendofblock(i, oo%, text$(j))
  IF level% < -1 THEN PRINT "error - more } then { at "; text$(j - 1): PRINT text$(j): END
  IF level% = -1 OR i < LEN(text$(j)) + 1 THEN CALL detendofblock(i, oldlevel%, text$(j)): GOTO 3
  NEXT j: PRINT "text in data file too long": END
END IF
3  e$ = MID$(text$(j), i + 1): text$(j) = LEFT$(text$(j), i - 1)
  IF INSTR(e$, "{") > 0 THEN SEEK #n, ss + i

jj = j    'cutheader if it does not contain jj$
FOR i = 1 TO jj
IF INSTR(LCASE$(text$(i)), jj$) > 0 GOTO 4
NEXT i: j = 0


4 END SUB

SUB checkrepeat (file1$)
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

SUB countchar (t$, a$, cnt)
'counts the times of appearance of a$ in t$
cnt = 0: p% = 1
WHILE INSTR(p%, t$, a$) > 0
cnt = cnt + 1
p% = INSTR(p%, t$, a$) + 1
WEND

END SUB

SUB detendofblock (i, op%, t$)
'determines the position of } (i) corresponding to end of block
' int t$, op% denotes on input the number of brackets } which
'are opened in addition to the } denoting the end of the block
cl% = 0
i = 0
WHILE op% > -1
o% = INSTR(i + 1, t$, "{")
c% = INSTR(i + 1, t$, "}")
IF o% < i + 1 THEN
 IF c% < i + 1 THEN
  i = LEN(t$) + 1: GOTO 21
 ELSE
  i = c%
  op% = op% - 1
 END IF
ELSE
 IF c% < i + 1 THEN
  op% = op% + 1
  i = o%
 ELSE
  IF o% < c% THEN
   op% = op% + 1
   i = o%
  ELSE
   i = c%
   op% = op% - 1
  END IF
 END IF
END IF
WEND


21 END SUB

SUB inputline (n, d1#(), col1%, d2$(), col2%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
'comments in {} are stored in comment$
'if a comment is started somewhere in this line by "{" and not finished
'then col1% is set -1

a$ = INKEY$: IF a$ <> "" THEN IF ASC(a$) = 27 THEN END
' IF SCREEN(24, 10) <> 45 THEN PRINT "-";  ELSE LOCATE 24, 1: PRINT SPACE$(15); : LOCATE 24, 1
PRINT USING "-###%"; 100 * SEEK(n) / LOF(n); : LOCATE 24, 1
'PRINT USING "-########"; SEEK(n); : LOCATE 24, 1

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

