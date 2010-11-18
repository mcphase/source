DECLARE SUB inputline (n!, d1#(), col1%, d2$(), col2%)
DECLARE SUB getpar (er%, x!, x$, a$)
DECLARE SUB headerinput (text$(), j!, n!)
DECLARE SUB iteratefiles (file1$)
PRINT "VIEWPR VIEWPR VIEWPR VIEWPR VIEWPR VIEWPR VIEWPR VIEWPR VIEWPR"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program VIEWPR - use it like VIEWPR *.*sc "
DATA "                  VIEW and PRint HMI *.*sc files to filman format     "
DATA "                  acceptable options: xmin=  xmax= ymin= ymax=   C...continue"
DATA " format of input file                                                           "
DATA " { header: this is the   'SGEN dh=  dk= ...  '                                             "
DATA "   file header delimited                                                "
DATA "   by brackets after this header there follow 3 or more data columns }  "
DATA " 11 3.14235 65367                                                       "
DATA "  .    .     .                                                          "
DATA "  .    .     .                                                          "
DATA "  .    .     .    .  .   .                                              "
DATA "  .    .     .                                                          "
DATA " 32 2412.34 324.2                                                       "
DATA ""
DATA "*************************************************************************"
DIM text$(300), xm#(30), x$(30), d1#(30), d2$(30)

' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

filename$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
IF INSTR(filename$, " ") > 0 THEN
 i% = INSTR(filename$, " ")
 filename$ = RTRIM$(LEFT$(filename$, i%))
END IF
CALL getpar(erymin%, ymin, "ymin", LCASE$(COMMAND$))
CALL getpar(erymax%, ymax, "ymax", LCASE$(COMMAND$))
CALL getpar(erxmin%, xmin, "xmin", LCASE$(COMMAND$))
CALL getpar(erxmax%, xmax, "xmax", LCASE$(COMMAND$))

'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

OPEN "o", 1, "viewpr.bat": CLOSE 1

999 CALL iteratefiles(filename$)

IF RIGHT$(filename$, 4) = ".rcp" OR RIGHT$(filename$, 4) = ".mrc" THEN
   PRINT "you should never ever change data in *.rcp/mrc files": PLAY "dgdgdg"
222 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 222
END IF
'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

title$ = ""
'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)
CALL inputline(1, d1#(), col1%, d2$(), col2%)

fit$ = CHR$(col1% + 96)   'fit is always last column
fitm$ = CHR$(col1% - 1 + 96) 'fit is always last column

'get some scan parameters needed for plotting ggggggggggggggggggggggggggggg
FOR i% = 1 TO j
 
 
              'from line SGEN
 IF INSTR(text$(i%), "SGEN") > 0 THEN
  'for hscans
  CALL getpar(er%, dh, "dh", LCASE$(text$(i%))): IF dh <> 0 THEN xtitle$ = "H": view$ = "2": sdet$ = "j"
  CALL getpar(er%, h, "h", LCASE$(text$(i%)))
 
  'for kscans
  CALL getpar(er%, dk, "dk", LCASE$(text$(i%))): IF dk <> 0 THEN xtitle$ = "K": view$ = "3": sdet$ = "j"
  CALL getpar(er%, k, "k", LCASE$(text$(i%)))
 
  'for lscans
  CALL getpar(er%, dl, "dl", LCASE$(text$(i%))): IF dl <> 0 THEN xtitle$ = "L": view$ = "4": sdet$ = "j"
  CALL getpar(er%, l, "l", LCASE$(text$(i%)))
 
  'for omegascans
  CALL getpar(er%, domgs, "domgs", LCASE$(text$(i%))): IF domgs <> 0 THEN xtitle$ = "omega": view$ = "2": sdet$ = "4"
  CALL getpar(er%, dtths, "dtths", LCASE$(text$(i%))): IF dtths <> 0 THEN xtitle$ = "2theta": view$ = "2": sdet$ = "4"
  CALL getpar(er%, dchis, "dchis", LCASE$(text$(i%))): IF dchis <> 0 THEN xtitle$ = "chi": view$ = "2": sdet$ = "4"
  CALL getpar(er%, dphis, "dphis", LCASE$(text$(i%))): IF dphis <> 0 THEN xtitle$ = "phi": view$ = "2": sdet$ = "4"
 
  
  'set titles
  IF h <> 0 OR k <> 0 OR l <> 0 THEN title$ = "(" + STR$(h) + "  " + STR$(k) + "  " + STR$(l) + ")" + "         " + xtitle$ + "-scan"
  ytitle$ = "CountS"
 
 END IF

              'from line MON=
 IF INSTR(text$(i%), "MON=") > 0 THEN
  CALL getpar(er%, mon, "mon", LCASE$(text$(i%)))
  title$ = title$ + "      MON=" + STR$(mon)
 END IF

 IF INSTR(text$(i%), "T3=") > 0 THEN
  CALL getpar(er%, T3, "t3", LCASE$(text$(i%)))
  title$ = title$ + " T3=" + STR$(T3) + "K"
 END IF

 IF INSTR(text$(i%), "T4=") > 0 THEN
  CALL getpar(er%, T4, "t4", LCASE$(text$(i%)))
  title$ = title$ + " T4=" + STR$(T4) + "K"
 END IF


NEXT i%
CLOSE 1



'ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg

'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
'view the file
cont$ = "": IF INSTR(COMMAND$, "/C") > 0 THEN cont$ = " /c"
SHELL "view " + filename$ + " " + view$ + sdet$ + " " + filename$ + " " + view$ + fit$ + cont$
'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


'modify output batch file mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
OPEN "a", 1, "viewprn.bat"
IF LOF(1) < 10 THEN PRINT filename$ + " error viewpr: viewprn.bat too short": INPUT ala: CLOSE 1: GOTO 9966
CLOSE 1
OPEN "i", 1, "viewprn.bat"
OPEN "o", 2, "viewpr.out"
WHILE EOF(1) = 0
 LINE INPUT #1, a$
 a$ = a$ + " "
 ii% = INSTR(a$, "title="): IF ii% > 0 THEN a$ = LEFT$(a$, ii% + 5) + title$ + " " + MID$(a$, ii% + 6)
 ii% = INSTR(a$, "xaxis="): IF ii% > 0 THEN a$ = LEFT$(a$, ii% + 5) + xtitle$ + " " + MID$(a$, ii% + 6)
 ii% = INSTR(a$, "yaxis="): IF ii% > 0 THEN a$ = LEFT$(a$, ii% + 5) + ytitle$ + " " + MID$(a$, ii% + 6)

 ii% = INSTR(a$, "MINY="): IF ii% > 0 AND erymin% = 0 THEN a$ = LEFT$(a$, ii% + 4) + STR$(ymin)
 ii% = INSTR(a$, "MAXY="): IF ii% > 0 AND erymax% = 0 THEN a$ = LEFT$(a$, ii% + 4) + STR$(ymax)
 ii% = INSTR(a$, "MINX="): IF ii% > 0 AND erxmin% = 0 THEN a$ = LEFT$(a$, ii% + 4) + STR$(xmin)
 ii% = INSTR(a$, "MAXX="): IF ii% > 0 AND erxmax% = 0 THEN a$ = LEFT$(a$, ii% + 4) + STR$(xmax)
 CALL getpar(er%, d, "MINX=", a$): IF er% = 0 THEN xmin = d
 CALL getpar(er%, d, "MAXX=", a$): IF er% = 0 THEN xmax = d
 CALL getpar(er%, d, "MINY=", a$): IF er% = 0 THEN ymin = d
 CALL getpar(er%, d, "MAXY=", a$): IF er% = 0 THEN ymax = d

 'spline and plot fit
 ii% = INSTR(a$, view$ + fit$ + " symbol=")
 IF ii% > 0 THEN
  PRINT #2, "copy " + filename$ + " viewpr.fit"
  PRINT #2, "delcol viewpr.fit " + LTRIM$(RTRIM$(STR$(VAL(view$) + 1))) + fitm$
  PRINT #2, "delcol viewpr.fit 1" + LTRIM$(RTRIM$(STR$(VAL(view$) - 1)))
  PRINT #2, "spline viewpr.fit 1 spl=2"
   stps = (xmax - xmin) / 200
  PRINT #2, "interpol viewpr.fit 1 stp="; stps; " spl=2"
  'a$ = LEFT$(a$, ii% + 9) + "- " + MID$(a$, ii% + 10)
  a$ = "plot viewpr.fit 12 symbol=-"
 END IF

 'plot data
 ii% = INSTR(a$, view$ + sdet$ + " symbol="): IF ii% > 0 THEN a$ = LEFT$(a$, ii% + 9) + "O " + MID$(a$, ii% + 10)

 a$ = LTRIM$(RTRIM$(a$))

 PRINT #2, a$
 
WEND
newfilename$ = filename$: bs% = INSTR(newfilename$, "\"): bs% = bs% + 1
dot% = INSTR(bs%, newfilename$, ".")
IF dot% > 0 THEN newfilename$ = LEFT$(newfilename$, dot% - 1)
newfilename$ = newfilename$ + ".ps"
PRINT #2, "copy plot.ps " + newfilename$

CLOSE 1, 2
SHELL "del viewprn.bat"

SHELL "copy viewpr.bat+viewpr.out viewpr.sav"
SHELL "del viewpr.out"

SHELL "copy viewpr.sav viewpr.bat"
SHELL "del viewpr.sav"

'now we have extended viewpr.bat ...

PRINT
PRINT "END VIEWPR file "; filename$; "has been VIEWPRatted"


9966 IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END

333 FOR i = 1 TO 16: READ a$: PRINT a$: NEXT i: END

SUB getpar (er%, x, x$, a$)
'gets parameter x (name: x$) out of a$ (if present) - case match)
x = 0
in% = INSTR(a$, x$)
IF in% > 0 THEN
 in% = in% + LEN(x$)
 b$ = LTRIM$(MID$(a$, in%))
 IF LEFT$(b$, 1) = "=" THEN b$ = MID$(b$, 2)
 x = VAL(LTRIM$(b$))
 er% = 0
ELSE
 er% = 1
END IF

END SUB

SUB headerinput (text$(), j, n)
'**********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'**********************************************************************

1 LINE INPUT #n, a$
   i = INSTR(a$, "{"): IF i > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(a$, LEN(a$) - i)
  j = 1: i = INSTR(text$(j), "}"): IF i > 0 GOTO 3  'look for "}" in first line
                                         
   FOR j = 2 TO 300
   LINE INPUT #n, text$(j)
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

