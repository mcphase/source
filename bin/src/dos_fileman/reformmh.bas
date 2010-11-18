DECLARE SUB iteratefiles (file1$)
DECLARE SUB inputheader (text$(), j!, n!)
DECLARE SUB swapcol (text$, col1%, col2%)
DECLARE SUB getcol (text$, col$, col%)
DECLARE SUB inputline (n!, d1#(), col1%, d2$(), col2%)
DECLARE SUB splitstr (a$, s$(), no%)
DECLARE SUB getpar (x!, x$, a$)
DECLARE SUB analizecommand (file1$)
PRINT "REFORMMH REFORMMH REFORMMH REFORMMH REFORMMH REFORMMH REFORMMH REFORMMH REFORMMH"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program REFORMMH - use it like REFORMMH *.* "
DATA "             REFORMMH Herbie Muellers *.* files to filman format     "
DATA " format of output file                                                          "
DATA ""
DATA "#{ header: this is the                                                  "
DATA "#  file header delimited                                                "
DATA "#  by brackets after this header there follow 3 or more data columns }  "
DATA " 11 3.14235 65367                                                       "
DATA "  .    .     .                                                          "
DATA "  .    .     .                                                          "
DATA "  .    .     .    .  .   .                                              "
DATA "  .    .     .                                                          "
DATA " 32 2412.34 324.2                                                       "
DATA ""
DATA "*************************************************************************"
DIM text$(300), xm#(30), x$(30), s$(80), xmm$(30)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

filename$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))

999 CALL iteratefiles(filename$)

'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'open input file and input file header
OPEN "i", 1, filename$: CALL inputheader(text$(), j, 1)
IF LOF(1) - SEEK(1) < 10 THEN CLOSE 1: GOTO 11
id$ = ""

' open output file and write fileheader
OPEN "o", 2, "REFORMMH.ffo"
PRINT #2, "#{"; : FOR iii = 1 TO j: PRINT #2, "#"; text$(iii): NEXT
PRINT #2, "#"; DATE$; " "; TIME$;
PRINT #2, " REFORMMHted using program REFORMMH.bas}"


 REM input data columns
WHILE EOF(1) = 0
 CALL inputline(1, xm#(), col1%, xmm$(), col2%)
         
 'write result to file
 FOR coll% = 1 TO col1%: PRINT #2, xm#(coll%); : NEXT: PRINT #2,
WEND

22 CLOSE 1, 2

SHELL "copy REFORMMH.ffo " + filename$
SHELL "del REFORMMH.ffo"

PRINT
PRINT "END REFORMMH file "; filename$; "has been REFORMMHatted"


11 IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END

333 FOR i = 1 TO 16: READ a$: PRINT a$: NEXT i: END

SUB getcol (text$, col$, col%)
'in: text$(...column header), col$ ... col ID,
'out: col% ... number of column
DIM s$(30)

   CALL splitstr(text$, s$(), no%)
   FOR ss% = 1 TO no%
    IF INSTR(s$(ss%), col$) > 0 THEN col% = ss%: GOTO 21
   NEXT ss%

21 END SUB

SUB getpar (x, x$, a$)
'gets parameter x (name: x$) out of a$ (if present) - case match)
x = 0
in% = INSTR(a$, x$)
IF in% > 0 THEN
 in% = in% + LEN(x$)
 b$ = LTRIM$(MID$(a$, in%))
 IF LEFT$(b$, 1) = "=" THEN b$ = MID$(b$, 2)
 x = VAL(b$)
END IF

END SUB

SUB inputheader (text$(), j, n)
'**********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
' id$ contains info if it is a hscan kscan lscan ...
'**********************************************************************
DIM s$(80)

j = 0
DO
j = j + 1: LINE INPUT #n, text$(j)
'PRINT text$(j); : INPUT ala

LOOP WHILE INSTR(text$(j), "@EOH") = 0 AND EOF(n) = 0


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

SUB iteratefiles (file1$)
STATIC washere%, path$

  IF INSTR(COMMAND$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere% = 0 THEN
      washere% = 1
      path$ = ""
      IF INSTR(file1$, "\") > 0 THEN
       i% = INSTR(file1$, "\")
       WHILE INSTR(i% + 1, file1$, "\") <> 0: i% = INSTR(i% + 1, file1$, "\"): WEND
       path$ = LEFT$(file1$, i%)
      END IF
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
   file1$ = path$ + file1$
  END IF



END SUB

SUB splitstr (a$, s$(), no%)
' splits a$ into substrings
b$ = a$: no% = 0
WHILE LTRIM$(b$) <> ""
 b$ = LTRIM$(b$): no% = no% + 1
 IF INSTR(b$, " ") > 0 THEN
  ss% = INSTR(b$, " ")
  s$(no%) = LEFT$(b$, ss% - 1)
  b$ = MID$(b$, ss%)
 ELSE
  s$(no%) = b$: b$ = ""
 END IF
WEND
END SUB

SUB swapcol (text$, col1%, col2%)
' swap columns col1% and col2% in header text$
 
DIM s$(30)

   CALL splitstr(text$, s$(), no%)

SWAP s$(col1%), s$(col2%)

text$ = ""
FOR ss% = 1 TO no%
 text$ = text$ + " " + s$(ss%)
NEXT ss%


END SUB

