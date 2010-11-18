DECLARE SUB inputline (n!, d1(), col1%, d2$(), col2%)
DECLARE SUB maxmincint (minxii!, maxxii!, minxjj!, maxxjj!, zoom%)
DECLARE SUB initviewprn (minxii!, maxxii!, minxjj!, maxxjj!)
DECLARE SUB plot (filename$, maxxii!, minxii!, maxxjj!, minxjj!, z!, delpoint%, colorswitch!, ii%, jj%, pr$)
DECLARE SUB iteratefilename (file1$)
DECLARE SUB getminmax (filename$, zoom%, init%, minxii!, maxxii!, minxjj!, maxxjj!, ii%, jj%)
DECLARE SUB headerinput (text$(), j!, n!)
DECLARE SUB pointdraw3d (x!, xmin!, xmax!, y!, ymin!, ymax!, z!, zmin!, zmax!, size!)
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA " this programm is intended to VIEW data on some files                    "
DATA " it views the i column vs the j column                                   "
DATA " command format: VIEW *1.* 23 *2.* 14 ... /C                             "
DATA " means view col 2 vs col 3 of *1.* and col 1 vs col 4 of *1.* and ..."
DATA " and continue (/C))                                                  "
DATA " format of file                                                      "
DATA "                                                                     "
DATA " { header: this is the                                               "
DATA "   file header delimited                                             "
DATA "   by brackets after this header there follow 3 or more data columns }"
DATA " 11 3.14235 65367                                                     "
DATA "  .    .     .                                                        "
DATA "  .    .     .    .. . . . .                                          "
DATA "  .    .     .                                                        "
DATA "  .    .     .                                                        "
DATA " 32 2412.34 324.2                                                     "
DATA "                                                                      "
DATA "                                                                      "
DATA "*************************************************************************"
DIM text$(300), x(30)
SCREEN 12, 6: WIDTH 80, 60: CLS ' clear screen
' aaaaaaaaaaaaaa ANALYSE COMMAND$ aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
a$ = LTRIM$(RTRIM$(LCASE$(COMMAND$)))
' look if option /C (continue) has been entered
continue% = 0: IF LCASE$(RIGHT$(a$, 2)) = "/c" THEN continue% = 1: a$ = RTRIM$(LEFT$(a$, LEN(a$) - 2))
save$ = a$  ' save all the filenames as save$
'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa


4 '  L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L  >>>
'   L   L   L   L start look for extremum   L   L   L   L   L   L   L   L >>>
'  L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L  >>>
a$ = save$: init% = 0
WHILE a$ <> ""
 'get next filename+columns gggggggggggggggggggggggggggggggggggggg
 IF INSTR(a$, " ") = 0 THEN PRINT "illegal file col1col2"; a$: INPUT ala: END
 filename$ = LEFT$(a$, INSTR(a$, " "))          'look up filename
 a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))  'rmove filename from a$
 ii% = ASC(LEFT$(a$, 1)) MOD 48: jj% = ASC(MID$(a$, 2, 1)) MOD 48 'get columns
 a$ = LTRIM$(RIGHT$(a$, LEN(a$) - 2))               'remove columns from a$
 'ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
 IF INSTR(filename$, "*") > 0 THEN 'see if multiple files are to be viewed
  CALL iteratefilename(filename$)
  WHILE filename$ <> ""
   CALL getminmax(filename$, zoom%, init%, minxii, maxxii, minxjj, maxxjj, ii%, jj%)
   CALL iteratefilename(filename$)
  WEND
 ELSE
  CALL getminmax(filename$, zoom%, init%, minxii, maxxii, minxjj, maxxjj, ii%, jj%)
 END IF
WEND
 
CALL maxmincint(minxii, maxxii, minxjj, maxxjj, zoom%)

' L  L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L<<<
'L  L   L   L   L   L   L   Lend look for extremum  L   L   L   L   L   L <<<
' L  L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L   L<<<


'colors used for plotting
CLS 0
11 DATA 14,13,10,11,12,15,9,7,3,2,6,5,4,8,1
RESTORE 11:  a$ = save$: READ colorswitch
CALL initviewprn(minxii, maxxii, minxjj, maxxjj)
'P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P>>>>
' P   P   P   P   P   P  START PLOT   P   P   P   P   P   P   P   P   P   P>>>
'P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P>>>>
nof% = 0: WHILE a$ <> "": nof% = nof% + 1
 'get next filename+columns gggggggggggggggggggggggggggggggggggggg
 filename$ = LEFT$(a$, INSTR(a$, " "))          'look up filename
 a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))  'rmove filename from a$
 pr$ = LEFT$(a$, 2) 'set column id for file viewprn.bat
 ii% = ASC(LEFT$(a$, 1)) MOD 48: jj% = ASC(MID$(a$, 2, 1)) MOD 48 'get columns
 a$ = LTRIM$(RIGHT$(a$, LEN(a$) - 2))               'remove columns from a$
 'ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg

 IF INSTR(filename$, "*") > 0 THEN 'see if multiple files are to be viewed
  CALL iteratefilename(filename$)
  WHILE filename$ <> ""
   CALL plot(filename$, maxxii, minxii, maxxjj, minxjj, z, delpoint%, colorswitch, ii%, jj%, pr$)
   CALL iteratefilename(filename$)
  WEND
 ELSE
  CALL plot(filename$, maxxii, minxii, maxxjj, minxjj, z, delpoint%, colorswitch, ii%, jj%, pr$)
 END IF

WEND: IF delpoint% > 0 THEN delpoint% = 0: SHELL "copy VIEW.del VIEW.sav": save$ = "VIEW.sav " + CHR$(ii% + 96) + CHR$(jj% + 96)
     
      PRINT #8, "plot /s": CLOSE 8   'endof file viewprn.bat
' P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P <<<<<
'P   P   P   P   P   P   P  END PLOT P   P   P   P   P   P   P   P   P  <<<<<
' P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P   P <<<<<




'eeeeeeeeeeeeeeeeeeeeeeeeee END PROGRAM (OR REPLOT)eeeeeeeeeeeeeeeeeeeeeee
COLOR 15 'switch to white color again
LOCATE 1, 1: PRINT " "; : VIEW PRINT 1 TO 2
IF continue% = 0 THEN
   PRINT "HOME+-_*" + CHR$(17) + CHR$(16) + CHR$(24) + CHR$(25) + "zoom,spc/del..VIEW,PRESS ESC TO CONTINUE!"; save$;
10 a$ = "": WHILE a$ = "": a$ = INKEY$: WEND 'get char from keybord

IF a$ = " " THEN delpoint% = 1: LOCATE 1, 1: PRINT "0"; : VIEW PRINT: GOTO 11
' save range to dosvariables to make plotting easier (with program
' plot.bas and graph.bas) and end program VIEW
 IF ASC(a$) = 27 THEN SHELL "del VIEW.sav": END
 
                'zoom functions
 rangex = maxxii - minxii: rangey = maxxjj - minxjj
 IF a$ = "+" THEN maxxii = minxii + rangex / 3: zoom% = 1: GOTO 4
 IF a$ = "-" THEN maxxii = minxii + rangex * 3: zoom% = 1: GOTO 4
 IF a$ = "*" THEN maxxjj = minxjj + rangey / 3: zoom% = 2: GOTO 4
 IF a$ = "_" THEN maxxjj = minxjj + rangey * 3: zoom% = 2: GOTO 4
 IF a$ = CHR$(0) + "M" THEN minxii = minxii + 1.01 * rangex / 5: maxxii = maxxii + .99 * rangex / 5: zoom% = 1: GOTO 4
 IF a$ = CHR$(0) + "K" THEN minxii = minxii - .99 * rangex / 5: maxxii = maxxii - 1.01 * rangex / 5: zoom% = 1: GOTO 4
 'up-down arrows
 IF a$ = CHR$(0) + "H" THEN minxjj = minxjj + .99 * rangey / 5: maxxjj = maxxjj + 1.01 * rangey / 5: zoom% = 2: GOTO 4
 IF a$ = CHR$(0) + "P" THEN minxjj = minxjj - .99 * rangey / 5: maxxjj = maxxjj - 1.01 * rangey / 5: zoom% = 2: GOTO 4
 'home
 IF a$ = CHR$(0) + "G" THEN zoom% = 0: GOTO 4

 GOTO 10        'right and left arrow
END IF

END 'eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
333 FOR i = 1 TO 20: READ a$: PRINT a$: NEXT i: END

SUB getminmax (filename$, zoom%, init%, minxii, maxxii, minxjj, maxxjj, ii%, jj%)
 DIM x(30), text$(300), d2$(10)
  OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1) 'open file and get fileheader

  ' input columns and look for ex

   WHILE EOF(1) = 0 '-----------------------------------------
21    CALL inputline(1, x(), col%, d2$(), col2%)
 IF col% < 0 AND EOF(1) = 0 THEN
  CALL headerinput(text$(), j, 1)' open file2 and input header
  GOTO 21
  END IF
    IF col% > 1 THEN
     IF zoom% = 0 THEN   'change maxxii and minxii and maxxjj and minxjj
      IF init% = 0 THEN init% = 1: minxii = x(ii%): maxxii = x(ii%): minxjj = x(jj%): maxxjj = x(jj%)
                                  'initialize ii and jj extremums
      IF x(ii%) < minxii THEN minxii = x(ii%)'hor axis
      IF x(ii%) > maxxii THEN maxxii = x(ii%)
     
      IF x(jj%) < minxjj THEN minxjj = x(jj%)'vert axis
      IF x(jj%) > maxxjj THEN maxxjj = x(jj%)
     
     END IF
    
     IF zoom% = 1 THEN  'change only maxxjj and minxjj
      IF x(ii%) >= minxii AND x(ii%) <= maxxii THEN
       IF init% = 0 THEN init% = 1:  minxjj = x(jj%): maxxjj = x(jj%)
                                   'initialize jj extremums
       IF x(jj%) < minxjj THEN minxjj = x(jj%)'vert axis
       IF x(jj%) > maxxjj THEN maxxjj = x(jj%)
      END IF
     END IF
    END IF
   WEND: CLOSE 1 '--------------------------------------------

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

SUB initviewprn (minxii, maxxii, minxjj, maxxjj)
      OPEN "o", 8, "viewprn.bat"
      PRINT #8, "set MINX=" + STR$(minxii)
      PRINT #8, "set MAXX=" + STR$(maxxii)
      PRINT #8, "set MINY=" + STR$(minxjj)
      PRINT #8, "set MAXY=" + STR$(maxxjj)
      PRINT #8, "del plot.ps"
PRINT #8, "set locx = 120"
PRINT #8, "set locy = 160"
PRINT #8, "set width=150"
PRINT #8, "set height=220"

PRINT #8, "del graph.ps"
PRINT #8, "graph xaxis= yaxis= title="
PRINT #8, "rename graph.ps plot.ps"

PRINT #8, "rem plot file using standard symbols"

END SUB

SUB inputline (n, d1(), col1%, d2$(), col2%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1(1...col1%)
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
    d1(col1%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
    ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " ")))
 WEND
END IF

'

END SUB

SUB iteratefilename (file1$)
STATIC path$
'*****************************************************************
' this sub iterates the filename if there is a * in it ...
' it returns "" if no new filename is found any more
'*****************************************************************

  IF INSTR(file1$, "*") <> 0 THEN
    'look for path
    IF INSTR(file1$, "\") > 0 THEN
     p% = 1: WHILE p% > 0: pl% = p%: p% = INSTR(pl% + 1, file1$, "\"): WEND
     path$ = LCASE$(LEFT$(file1$, pl% - 1)) + "\"
    ELSE
     path$ = ""
    END IF
     
      'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /b > view.dir"
      OPEN "i", 9, "view.dir"
   END IF

   IF EOF(9) <> 0 THEN
    CLOSE 9: SHELL "del view.dir"
    file1$ = ""
   ELSE
    INPUT #9, file1$
    WHILE LCASE$(file1$) = "view.dir" OR LCASE$(file1$) = "viewprn.bat"
     IF EOF(9) <> 0 THEN
      CLOSE 9: SHELL "del view.dir": file1$ = ""
     ELSE
      INPUT #9, file1$
     END IF
    WEND
   
    file1$ = path$ + file1$: file1$ = LTRIM$(file1$)
    IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del view.dir"
   END IF


END SUB

SUB maxmincint (minxii, maxxii, minxjj, maxxjj, zoom%)
'dddddddddddddddddddddddd  make sensible range ddddddddddddddddddddddddd
 IF zoom% = 0 OR zoom% = 1 THEN
  IF minxii = maxxii THEN maxxii = maxxii + 10 ^ 20 'for x axis
   minxsav = minxii
   maxxsav = maxxii
  delta = (maxxii - minxii) / 10: exponent% = INT(LOG(delta) / LOG(10)) - 1
  delta = (INT(delta / 10 ^ exponent%)) * 10 ^ exponent%  'round delta to sensible number
  IF minxii <> 0 THEN exponentmin% = INT(LOG(ABS(minxii)) / LOG(10)) ELSE exponentmin% = -100
  IF maxxii <> 0 THEN exponentmax% = INT(LOG(ABS(maxxii)) / LOG(10)) ELSE exponentmax% = -100
  expon% = exponentmin%: IF exponentmin% < exponentmax% THEN expon% = exponentmax%
  IF expon% < eponent% + 5 THEN minxii = INT(minxii / delta) * delta
  maxxii = minxii + 10 * delta
  WHILE maxxii < maxxsav: maxxii = maxxii + delta: WEND
  WHILE minxii > minxsav: minxii = minxii - delta: WEND
 END IF
 IF zoom% = 0 OR zoom% = 1 OR zoom% = 2 THEN
  IF minxjj = maxxjj THEN maxxjj = maxxjj + 10 ^ 20  'for y axis
   minysav = minxjj
   maxysav = maxxjj
  delta = (maxxjj - minxjj) / 10: exponent% = INT(LOG(delta) / LOG(10)) - 1
  delta = (INT(delta / 10 ^ exponent%)) * 10 ^ exponent% 'round delta to sensible number
  IF minxjj <> 0 THEN exponentmin% = INT(LOG(ABS(minxjj)) / LOG(10)) ELSE exponentmin% = -100
  IF maxxjj <> 0 THEN exponentmax% = INT(LOG(ABS(maxxjj)) / LOG(10)) ELSE exponentmax% = -100
  expon% = exponentmin%: IF exponentmin% < exponentmax% THEN expon% = exponentmax%
  IF expon% < eponent% + 5 THEN minxjj = INT(minxjj / delta) * delta
  maxxjj = minxjj + 10 * delta
  WHILE maxxjj < maxysav: maxxjj = maxxjj + delta: WEND
  WHILE minxjj > minysav: minxjj = minxjj - delta: WEND
 END IF
'dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd

END SUB

SUB plot (filename$, maxxii, minxii, maxxjj, minxjj, z, delpoint%, colorswitch, ii%, jj%, pr$)
STATIC plotsymbol%
 DIM text$(300), x(30), d2$(30)
 OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1) 'open file and get fileheader

 IF delpoint% = 1 THEN OPEN "o", 2, "VIEW.del": PRINT #2, "{"; : FOR k% = 1 TO j - 1: PRINT #2, text$(k%): NEXT: PRINT #2, text$(j) + "}"

 ' input columns and plot data on screen ----------------
 WHILE EOF(1) = 0
22  CALL inputline(1, x(), col%, d2$(), col2%)
 IF col% < 0 AND EOF(1) = 0 THEN
  CALL headerinput(text$(), j, 1)' open file2 and input header
  GOTO 22
  END IF
            'input dataline and plot
  IF x(ii%) < maxxii AND x(ii%) > minxii AND x(jj%) < maxxjj AND x(jj%) > minxjj THEN
      CALL pointdraw3d(x(ii%), minxii, maxxii, x(jj%), minxjj, maxxjj, z, 1, 1, colorswitch + .02 * delpoint%)
    del% = 0
    IF delpoint% = 1 THEN
      IF LEN(d2$(1)) > 40 THEN d2$(1) = LEFT$(d2$(1), 40)
      COLOR colorswitch: LOCATE 6, 1: PRINT USING "(x/y)=(##.#####^^^^/##.#####^^^^) "; x(ii%); x(jj%); : PRINT d2$(1);
      i$ = "": WHILE i$ = "": i$ = INKEY$: WEND 'get char from keybord
      IF ASC(i$) = 27 THEN delpoint% = 2: del% = 0: i$ = " "
      WHILE i$ <> " "
       IF i$ = CHR$(0) + "R" THEN CALL pointdraw3d(x(ii%), minxii, maxxii, x(jj%), minxjj, maxxjj, z, 1, 1, colorswitch + .02 * delpoint%): del% = 0
                                 'insert
       IF i$ = CHR$(0) + "S" THEN CALL pointdraw3d(x(ii%), minxii, maxxii, x(jj%), minxjj, maxxjj, z, 1, 1, colorswitch - 1 + .02 * delpoint%): del% = 1
                                 'delpoint
       i$ = "": WHILE i$ = "": i$ = INKEY$: WEND'get char from keybord
       IF ASC(i$) = 27 THEN delpoint% = 2: del% = 0: i$ = " "
      WEND
    END IF
  END IF
    IF delpoint% > 0 AND del% = 0 THEN FOR k% = 1 TO col%: PRINT #2, x(k%); : NEXT: PRINT #2,
 WEND: CLOSE 1, 2
 IF delpoint% = 2 THEN delpoint% = 0: SHELL "copy VIEW.del VIEW.sav": save$ = "VIEW.sav " + CHR$(ii% + 96) + CHR$(jj% + 96)
 '-------------------------------------------------------------

  VIEW PRINT 53 TO 60:  COLOR colorswitch ' print text of file
  tt$ = "": FOR kkk% = 1 TO j: IF LEN(tt$) < 80 THEN tt$ = tt$ + RTRIM$(LTRIM$(text$(kkk%)))
  NEXT kkk%
  LOCATE 60, 1: PRINT RIGHT$(filename$, 12) + LEFT$(tt$, 68): VIEW PRINT
  IF colorswitch = 1 THEN RESTORE 11
  READ colorswitch
  plotsymbol% = plotsymbol% + 1
  IF plotsymbol% > 24 THEN plotsymbol% = 1
  PRINT #8, "plot " + filename$ + " " + pr$ + " symbol=" + CHR$(plotsymbol% + 96)
END SUB

SUB pointdraw3d (x, xmin, xmax, y, ymin, ymax, z, zmin, zmax, size)
STATIC mins(), maxs(), i, s()
REM*************************************************************************
REM*************************************************************************
REM diese routine zeichnet einen punkt der groesse size
REM in einer 3dgrafik und mahlt die skala, wenn sie das erste mal aufgerufen
REM wird bzw. wenn die skala sich aendert werden auch alle anderen punkte
REM neu skaliert und gemahlt
REM nach dem aufruf erscheint in der linke oberen ecke des bildschrims
REM eine 0, wird diese ueberschrieben, so wird bei einem neuerlichen
REM aufruf der routine ein neuer plot angefangen und die alten punkte
REM vergessen
REM
REM eingabe
REM x,y,z               punktkoordinaten
REM x,y,z,min,xyzmax       bereich (ist xyzmin>xyzmax, so erfolgt automatische
REM                                     bereichswahl,
REM                                     ist zmin=zmax so folgt 2d darstellung)
REM size                punktgroesse (>0), (<0..zeichnet linie von vorhergehendem punkt)
      'z.B 3.14 means color 3 and circle diameter 14

REM some fix parameters of the routine XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
alpha = 45 / 180 * 3.1415: REM winkel des schraegrisses
zsmax = 210: REM laenge der z-achse
x$ = " ": y$ = " ": z$ = " " 'achsenbeschriftung
height = 360: widt = 620: REM hoehe und breite des bildes
h0 = 40: w0 = 0: REM linke obere ecke des bildes
'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

REM*************************************************************************
REM*************************************************************************
DIM min(3), max(3): min(1) = xmin: min(2) = ymin: min(3) = zmin 'rename xyzmin
                    max(1) = xmax: max(2) = ymax: max(3) = zmax 'rename xyzmax

'ccccccccc calculate maximum length of x,y (z) axis on the screen cccccccccc
IF min(3) = max(3) THEN zsmax = 20 / COS(alpha)  'wenn 2dim plot dann zsmax fast 0
xsmax = widt - zsmax * SIN(alpha)
ysmax = height - 10 - zsmax * COS(alpha)
'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

'iiiii this is only done once on first call(initialize plotroutine) iiiiiiiiiiiiii
IF i = 0 THEN
 DIM s(3, 1000), mins(3), maxs(3)' dimension array(static) for saving data points
 SCREEN 12, 6:  WIDTH 80, 60: CLS
 VIEW SCREEN (w0, h0 + 15)-(w0 + widt + 11, h0 + height + 15), 0, 3
                       'setup output screen
END IF 'iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

'look if pointdraw should make new plot or add point to old one
IF SCREEN(1, 1, 0) <> 48 THEN makenewgraph% = 1 ELSE makenewgraph% = 0
IF makenewgraph% = 1 THEN i = 0

i = i + 1: 'a new point is to be drawn, give it number i
     'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
     i = 2:  ' this line is ONLY FOR PROGRAM VIEW.BAS it gives all points to number
        ' one - that means no points are kept in memory !!!!!!!
        ' for normal operation of this routine delete this paragraph
     'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
s(1, i) = x: s(2, i) = y: s(3, i) = z: s(0, i) = size   'store the new point

REM abfrage automatic scale for x y and z axis  sssssssssssssssssssssssssssss
FOR axis% = 1 TO 3
IF min(axis%) > max(axis%) THEN
 min(axis%) = s(axis%, 1): max(axis%) = s(axis%, 1)
 FOR j = 1 TO i 'calculate extremums
  IF s(axis%, j) > max(axis%) THEN max(axis%) = s(axis%, j)
  IF s(axis%, j) < min(axis%) THEN min(axis%) = s(axis%, j)
 NEXT: IF min(axis%) = max(axis%) THEN max(axis%) = max(axis%) + 10 ^ 20
 delta = (max(axis%) - min(axis%)) / 9: exponent% = INT(LOG(delta) / LOG(10))
 delta = (1 + INT(delta / 10 ^ exponent%)) * 10 ^ exponent%'round delta to sensible number
 IF min(axis%) <> 0 THEN exponentmin% = INT(LOG(ABS(min(axis%))) / LOG(10)) ELSE exponentmin% = -100
 IF max(axis%) <> 0 THEN exponentmax% = INT(LOG(ABS(max(axis%))) / LOG(10)) ELSE exponentmax% = -100
 expon% = exponentmin%: IF exponentmin% < exponentmax% THEN expon% = exponentmax%
 IF expon% < eponent% + 5 THEN min(axis%) = INT(min(axis%) / delta) * delta
 max(axis%) = min(axis%) + 10 * delta
END IF
NEXT axis% 'Sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

redrawscale% = makenewgraph% 'redrawscale% tells if the scales has to be redrawn-------
FOR axis% = 1 TO 3
IF min(axis%) <> mins(axis%) OR max(axis%) <> maxs(axis%) THEN redrawscale% = 1
  mins(axis%) = min(axis%): maxs(axis%) = max(axis%)
NEXT axis% '----------------------------------------------------------------

IF redrawscale% = 1 THEN
 CLS : COLOR 15: VIEW PRINT: LOCATE 1, 1: PRINT "0";  'skala machen
              
 ' labels and scale for x,y axis xxxxxyyyyyyyyyyyyyyyxxxxxxxxxxxxyyyyyyyyyy
 LOCATE CINT((h0 + height + 8) * .125), 1: PRINT USING "##.##^^^^"; min(1);
  FOR j = .2 TO 1.1 STEP .2
   LOCATE CINT((h0 + height + 8) * .125), CINT(.12 * w0 + .12 * j * xsmax): PRINT USING "##.##^^^^"; min(1) + j * (max(1) - min(1));
  NEXT: PRINT x$;
 LOCATE CINT((h0 + height - 8) * .125), 1 + CINT(.12 * w0): PRINT USING "##.####^^^^"; min(2);
 LOCATE CINT((h0 + height) * .125) - CINT(.125 * ysmax), 1 + CINT(.12 * w0): PRINT USING "##.####^^^^"; max(2); : PRINT y$;

 LINE (5 + w0, h0 + height - 5)-(w0 + xsmax, h0 + height - 5), 3
 LINE (5 + w0, h0 + height - 5 - ysmax)-(w0 + 5, h0 + height - 5), 3
 FOR j = 0 TO 1.01 STEP .1
 LINE (w0 + xsmax * j + 5, h0 + height - 5)-(w0 + xsmax * j + 5, h0 + height - 2), 3
 LINE (w0 + 2, h0 + height - 5 - j * ysmax)-(w0 + 5, h0 + height - 5 - j * ysmax), 3
 NEXT 'xxxxxxxxxyyyyyyyyyyyyyxxxxxxxxxxxxxxxxyyyyyyyyyyyyyyyxxxxxxxxxxyyyyy

 IF min(3) <> max(3) THEN  'scale z-axis zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
  LOCATE CINT((h0 + height) * .132) + 1 - CINT(.16 * zsmax * COS(alpha)), CINT(.12 * (w0 + zsmax * SIN(alpha)))
  PRINT USING "##.#^^^^"; max(3); : PRINT z$;
  LOCATE CINT((h0 + height) * .132) + 1 - CINT(.1 * zsmax * COS(alpha)), 1 + CINT(.12 * w0): 'CINT(.07 * zsmax * SIN(alpha)) - 7
  PRINT USING "##.#^^^^"; min(3) + .5 * (max(3) - min(3))
  LINE (w0 + 5, h0 + height - 5)-(w0 + 5 + SIN(alpha) * zsmax, h0 + height - 5 - COS(alpha) * zsmax), 3
  FOR j = 0 TO 1.01 STEP .1
  LINE (w0 + 2 + j * zsmax * SIN(alpha), h0 + height - 5 - COS(alpha) * j * zsmax)-(5 + j * zsmax * SIN(alpha), height - 5 - COS(alpha) * j * zsmax), 3
  NEXT
 END IF 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

  FOR j = 1 TO i: GOSUB 3199: NEXT j 'redraw old points if there are
END IF

j = i: GOSUB 3199  'draw new point i if there are already more than one
LOCATE CINT(height * .123) + 1, 1
GOTO 3200   'goto end sub
 
'$$$$$$$$$$$$$$$$$$$$$$$$$  GOSUB drawpoint $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
3199  '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
REM zeichnen des punktes nummer j drawpoint $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
     'calculate position of point xp,yp on screen
IF min(3) <> max(3) THEN    'make 3d-plot
 xp = 5 + (s(1, j) - mins(1)) / (maxs(1) - mins(1)) * xsmax + (s(3, j) - mins(3)) / (maxs(3) - mins(3)) * SIN(alpha) * zsmax
 yp = height - 5 - (s(2, j) - mins(2)) / (maxs(2) - mins(2)) * ysmax - (s(3, j) - mins(3)) / (maxs(3) - mins(3)) * COS(alpha) * zsmax
ELSE   'make 2d plot
 xp = 5 + (s(1, j) - mins(1)) / (maxs(1) - mins(1)) * xsmax
 yp = height - 5 - (s(2, j) - mins(2)) / (maxs(2) - mins(2)) * ysmax
END IF

IF s(0, j) >= 0 THEN      'draw circle
CIRCLE (w0 + xp, h0 + yp), (s(0, j) - INT(s(0, j))) * 100, INT(s(0, j))
ELSE
 IF j > 1 THEN            'draw line
  IF min(3) <> max(3) THEN
   xpm = 5 + (s(1, j - 1) - mins(1)) / (maxs(1) - mins(1)) * xsmax + (s(3, j - 1) - mins(3)) / (maxs(3) - mins(3)) * SIN(alpha) * zsmax
   ypm = height - 5 - (s(2, j - 1) - mins(2)) / (maxs(2) - mins(2)) * ysmax - (s(3, j - 1) - mins(3)) / (maxs(3) - mins(3)) * COS(alpha) * zsmax
  ELSE
   xpm = 5 + (s(1, j - 1) - mins(1)) / (maxs(1) - mins(1)) * xsmax
   ypm = height - 5 - (s(2, j - 1) - mins(2)) / (maxs(2) - mins(2)) * ysmax
  END IF
  LINE (w0 + xp, h0 + yp)-(w0 + xpm, h0 + ypm), INT(s(0, j))
 END IF
END IF
RETURN '$$$$$$$$$$$$$$$$$$ end sub drawpoint $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

3200 END SUB

