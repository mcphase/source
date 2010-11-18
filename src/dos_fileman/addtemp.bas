DECLARE SUB calctemp (er%, tempfile$, mtime$, mdate$, T1!, T2!, T3!, T4!)
DECLARE SUB splitstr (a$, s$(), no%)
DECLARE SUB getpar (x!, x$, a$)
DECLARE SUB inputline (n!, d1#(), col1%, d2$(), col2%)
DECLARE SUB iteratefiles (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "ADDTEMP ADDTEMP ADDTEMP ADDTEMP ADDTEMP ADDTEMP ADDTEMP ADDTEMP ADDTEMP ADDTEMP"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program ADDTEMP - use it like ADDTEMP *.* temp.txt "
DATA "    (means calc.temperatures from HMI T sensor file temp.txt "
DATA " ---> the result is written in file *.*,     "
DATA " option:  synctim=65000=18:00:22"
DATA ""
DATA " format of file: "
DATA " { header: this is the                    "
DATA "   file header delimited                  "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "*************************************************************************"
DIM text$(300), x#(30, 100), d1#(30), d2$(10), par#(30), av#(30), s$(30)

'aaaaa analyse command string (get filename$ - ii% - jj% and interval) aaaaaa
a$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))

CALL splitstr(a$, s$(), no%)
filename$ = s$(1)
tempfile$ = s$(2)


999 CALL iteratefiles(filename$)

'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'oooooooooooooooooo OPEN FILES oooooooooooooooooooooooooooooooooooooooooooooo
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

'calculate temperatures  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
FOR ii% = 1 TO j
 IF INSTR(text$(ii%), "DATE") > 0 THEN datum$ = text$(ii%)
 IF INSTR(text$(ii%), "TIME") > 0 THEN mtime$ = text$(ii%)
 IF INSTR(text$(ii%), "NR") > 0 THEN
  text$(ii%) = text$(ii%) + " T1S  T2S  T3S  T4S"
 END IF
NEXT ii%
CALL calctemp(er%, tempfile$, mtime$, datum$, T1, T2, T3, T4)
IF er% = 1 THEN CLOSE 1: PRINT "error calculating temp for file "; filename$: GOTO 222
'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

' open output file and write fileheader
OPEN "o", 2, "ADDTEMP.lin"
PRINT #2, "{";
FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, USING "T1=####.### K T2=####.### K T3=####.### K T4=####.### K "; T1; T2; T3; T4
PRINT #2, DATE$; " "; TIME$; " column"; jj%; "ADDTEMP calc. with respect to column"; ii%;
PRINT #2, "temperatures  put into column"; col% + 1; "to"; col% + 4
PRINT #2, "calculation done by program ADDTEMP.bas";
PRINT #2, "}"

i = 0:  WHILE EOF(1) = 0
CALL inputline(1, d1#(), col%, d2$(), col2%)
 i = i + 1: FOR coll% = 1 TO col%: PRINT #2, d1#(coll%);
NEXT

PRINT #2, T1; T2; T3; T4;
FOR coll% = 1 TO col2%: PRINT #2, " {" + d2$(coll%) + "}"; : NEXT: PRINT #2,
WEND
'ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd

'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

CLOSE 1, 2

SHELL "copy ADDTEMP.lin " + filename$
SHELL "del ADDTEMP.lin"
PRINT
PRINT "END ADDTEMP in file "; filename$;
PRINT "TEMP calc. with "; tempfile$

222 IF INSTR(COMMAND$, "*") <> 0 GOTO 999

END
333 FOR i = 1 TO 14: READ a$: PRINT a$: NEXT i: END

SUB calctemp (er%, tempfile$, mtime$, mdate$, T1, T2, T3, T4)
'calculate temperature using input file T1 T2 T3 T4
DIM tt$(30), s$(30), d1#(100), d2$(10)
er% = 0


OPEN "i", 3, tempfile$
CALL headerinput(tt$(), jt, 3)
'calculate measuring time - time from tt$(6) in seconds
'elapsedtime= ?
'format mtime$ 'TIME 11:45:20'
'format mdate$ 'DATE  5-OCT-99'
'format tt$(5) ... zero time 'Mo, 4. Okt. 99 00:00:00-35234'
 CALL splitstr(tt$(6), s$(), no%)
 ttday% = VAL(s$(1))
 tthour% = VAL(s$(4))
 ttmin% = VAL(MID$(LTRIM$(s$(4)), 4))

 CALL splitstr(mdate$, s$(), no%)
  mday% = VAL(s$(2))'day calculation not intended to cross month or year boundary
 CALL splitstr(mtime$, s$(), no%)
 mhour% = VAL(s$(2))
 mmin% = VAL(MID$(LTRIM$(s$(2)), 4))

 elapsedtime# = (((mday% - ttday%) * 24# + mhour% - tthour%) * 60# + mmin% - ttmin%) * 60#
IF INSTR(COMMAND$, "SYNC") > 0 THEN
 CALL getpar(ttsyn, "SYNC", COMMAND$)
 CALL splitstr(COMMAND$, s$(), no%)
 FOR i% = 1 TO no%
  IF INSTR(s$(i%), "SYNC") > 0 THEN
   ii% = INSTR(s$(i%), "=")
   ii% = INSTR(ii% + 1, s$(i%), "=")
   msyntime$ = LTRIM$(MID$(s$(i%), ii% + 1))
  END IF
 NEXT i%
 msynhour% = VAL(msyntime$)
 msynmin% = VAL(MID$(LTRIM$(msyntime$), 4))
 diff = ttsyn - (msynhour% * 60# + msynmin%) * 60
 PRINT "difference in computer times is "; diff; "seconds"
 elapsedtime# = elapsedtime# + diff
END IF



 IF elapsedtime# < 0 THEN er% = 1

WHILE EOF(3) = 0 AND elapsedtime# > timein#
CALL inputline(3, d1#(), col1%, d2$(), col2%)
timein# = d1#(1)
WEND
IF elapsedtime# <= timein# THEN
 T1 = d1#(2)
 T2 = d1#(3)
 T3 = d1#(4)
 T4 = d1#(5)

ELSE
 er% = 1  'elapsedtime is bigger than maximum time in tempfile
END IF

CLOSE 3
END SUB

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

SUB headerinput (text$(), j, n)
'*********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'*********************************************************************

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

'
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

