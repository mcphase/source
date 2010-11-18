DECLARE SUB inputline (n!, d1#(), col1%, d2$(), col2%)
DECLARE SUB checkrepeat (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
DECLARE SUB pointdraw3d (x!, xmin!, xmax!, y!, ymin!, ymax!, z!, zmin!, zmax!, size!)
PRINT "PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA " **********************************************************************"
DATA " program PLOT: use as"
DATA " plot *.* 13 [size=7][shade=7][xbars=4][xlow=2 xup=4][ybars=5][ylow=2 yup=3][symbol=a][/s/tx/ty/dx=6/dy=8]"
DATA " ->plots data stored in file *.*, xaxis is column1, yaxis is column 3"
DATA " postscript file plot.ps is created
DATA " symbol: any character, 0 for empty circle, . for full circel or - for line"
DATA " size: size of symbol(<1), shade: intensity of symbol(<1)"
DATA " xbars: x-error bar column   ybars: y-errorbars column"
DATA " xlow/xup: lower and upper measurement limits in x direction"
DATA " ylow/yup: lower and upper measurement limits in y direction"
DATA " option /dx=6 /dy=8 displays numbers on plot corresponding to col6 and 8"
DATA " option /tx, /ty displays comment text for each data point (horizontal,vertical)"
DATA " remember: only when the option /s is added the        "
DATA " postscript file is finished using the 'showpage' command                  "
DATA " you must have defined dos variables minx maxx miny ymaxy (range of plot)"
DATA " and locx locy heigth width (location and size of the plot (in mm))      "
DATA " use program graph to create the scales for such a plot                  "
DATA "                                                                         "
DATA " the filename of the text is just written, where the cursor is located   "
DATA "**************************************************************************"
IF INSTR(COMMAND$, "*") > 0 THEN SHELL "copy plot.ps plot.dus"
REM**********einlesen der dos umgebungsvariablen des graphs****************
DIM check%(8)
SHELL "set > plot.dum": OPEN "i", 2, "plot.dum"
WHILE EOF(2) = 0
INPUT #2, a$: a$ = LCASE$(a$)
 IF LEFT$(a$, 4) = "minx" THEN xmin = VAL(RIGHT$(a$, LEN(a$) - INSTR(a$, "="))): check%(1) = 1
 IF LEFT$(a$, 4) = "maxx" THEN xmax = VAL(RIGHT$(a$, LEN(a$) - INSTR(a$, "="))): check%(2) = 1
 IF LEFT$(a$, 4) = "miny" THEN ymin = VAL(RIGHT$(a$, LEN(a$) - INSTR(a$, "="))): check%(3) = 1
 IF LEFT$(a$, 4) = "maxy" THEN ymax = VAL(RIGHT$(a$, LEN(a$) - INSTR(a$, "="))): check%(4) = 1
 IF LEFT$(a$, 5) = "width" THEN gwidth = VAL(RIGHT$(a$, LEN(a$) - INSTR(a$, "="))): check%(5) = 1
 IF LEFT$(a$, 6) = "height" THEN gheight = VAL(RIGHT$(a$, LEN(a$) - INSTR(a$, "="))): check%(6) = 1
 IF LEFT$(a$, 4) = "locx" THEN locx = VAL(RIGHT$(a$, LEN(a$) - INSTR(a$, "="))): check%(7) = 1
 IF LEFT$(a$, 4) = "locy" THEN locy = VAL(RIGHT$(a$, LEN(a$) - INSTR(a$, "="))): check%(8) = 1
WEND
FOR i% = 1 TO 8: IF check%(i%) <> 1 THEN PRINT "ERROR PROGRAM PLOT - dos variable "; i%; " not set correctly": END
NEXT i%
CLOSE 2: SHELL "del plot.dum"
'**************************************************************************
' analyse command$
a$ = COMMAND$: a$ = LCASE$(LTRIM$(RTRIM$(a$)))
filename$ = LEFT$(a$, INSTR(a$, " ")): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))

IF a$ = "" OR a$ = "/s" THEN OPEN "a", 1, "plot.ps": GOTO 133
colx% = ASC(LEFT$(a$, 1)) MOD 48: coly% = ASC(MID$(a$, 2, 1)) MOD 48
a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))
size = 0: IF INSTR(a$, "size=") > 0 THEN size = ASC(MID$(a$, INSTR(a$, "size=") + 5, 1)) MOD 48
shade = 0: IF INSTR(a$, "shade=") > 0 THEN shade = ASC(MID$(a$, INSTR(a$, "shade=") + 6, 1)) MOD 48
xbars = 0: IF INSTR(a$, "xbars=") > 0 THEN xbars = ASC(MID$(a$, INSTR(a$, "xbars=") + 6, 1)) MOD 48
xlow = 0: IF INSTR(a$, "xlow=") > 0 THEN xlow = ASC(MID$(a$, INSTR(a$, "xlow=") + 5, 1)) MOD 48
xup = 0: IF INSTR(a$, "xup=") > 0 THEN xup = ASC(MID$(a$, INSTR(a$, "xup=") + 4, 1)) MOD 48
ybars = 0: IF INSTR(a$, "ybars=") > 0 THEN ybars = ASC(MID$(a$, INSTR(a$, "ybars=") + 6, 1)) MOD 48
ylow = 0: IF INSTR(a$, "ylow=") > 0 THEN ylow = ASC(MID$(a$, INSTR(a$, "ylow=") + 5, 1)) MOD 48
yup = 0: IF INSTR(a$, "yup=") > 0 THEN yup = ASC(MID$(a$, INSTR(a$, "yup=") + 4, 1)) MOD 48
dx = 0: IF INSTR(a$, "/dx=") > 0 THEN dx = ASC(MID$(a$, INSTR(a$, "/dx=") + 4, 1)) MOD 48
dy = 0: IF INSTR(a$, "/dy=") > 0 THEN dy = ASC(MID$(a$, INSTR(a$, "/dy=") + 4, 1)) MOD 48
symbol$ = "": IF INSTR(a$, "symbol=") > 0 THEN symbol$ = MID$(a$, INSTR(a$, "symbol=") + 7, 1)

999 CALL checkrepeat(filename$)

OPEN "a", 1, "plot.ps"

'**************************************************************************

DIM text$(1000), d#(30), d2$(30)

IF LOF(1) < 1 THEN newfile% = 1

REM header einlesen*****************************************************
OPEN "i", 2, filename$
22 CALL headerinput(text$(), j, 2)
'**************************************************************************

PRINT #1, "/mm {72 mul 25.4 div} bind def"
PRINT #1, "/Helvetica findfont"
PRINT #1, "6 scalefont setfont"
PRINT #1, "0 setgray"
IF newfile% = 1 THEN
PRINT #1, locx - gwidth / 2 - 7; " mm "; locy - gheight / 2 - 14; " mm moveto"
ELSE
PRINT #1, "0 -3 mm rmoveto"
PRINT #1, "gsave"   '-store xy position on stack
END IF
  tt$ = "": FOR kkk% = 1 TO j: IF LEN(tt$) < 200 THEN tt$ = tt$ + RTRIM$(LTRIM$(text$(kkk%)))
  NEXT kkk%

t$ = LEFT$(COMMAND$ + tt$, INT(gwidth * 1.1))
'remove brackets () from text
WHILE INSTR(t$, "(") <> 0: t$ = LEFT$(t$, INSTR(t$, "(") - 1) + RIGHT$(t$, LEN(t$) - INSTR(t$, "(")): WEND
WHILE INSTR(t$, ")") <> 0: t$ = LEFT$(t$, INSTR(t$, ")") - 1) + RIGHT$(t$, LEN(t$) - INSTR(t$, ")")): WEND
PRINT #1, "("; t$; ") show"
'**************************************************************************


REM reflexe einlesen*******************************************
111 IF EOF(2) <> 0 GOTO 110
   CALL inputline(2, d#(), col%, d2$(), col2%)
 IF col% = -1 THEN PRINT #1, "grestore": newfile% = 0: GOTO 22 'a comment line has started and has to be treated
 comment$ = "": IF col2% > 0 THEN FOR i% = 1 TO col2%: comment$ = comment$ + d2$(i%): NEXT i%

REM pixel ausdrucken
 x = d#(colx%) - xmin: y = d#(coly%) - ymin
 tx = x: ty = y
IF d#(colx%) <= xmax AND x >= 0 AND y >= 0 AND d#(coly%) < ymax THEN

 x0 = locx - gwidth / 2: y0 = locy - gheight / 2
 scalex = gwidth / (xmax - xmin): scaley = gheight / (ymax - ymin)

 'determine size
 IF size > 0 THEN ssize = d#(size) ELSE ssize = SQR(gwidth * gheight) / 100: IF ssize < .6 THEN ssize = 1
 IF ssize < .4 THEN ssize = .4  'lower limit for ssize


 'set shade
 IF shade > 0 THEN sshade = d#(shade) ELSE sshade = 1
 IF sshade < .1 THEN sshade = .1'lower limit for sshade

 ' ALL SHADED THINGS
 IF sshade < 1 THEN PRINT #1, 1 - sshade; "setgray"
  
  'draw symbol
  IF symbol$ <> "0" AND symbol$ <> "-" AND symbol$ <> "" THEN
   PRINT #1, "/Helvetica findfont"
   PRINT #1, STR$(1 + INT(ssize)) + " mm scalefont setfont"
   PRINT #1, x0 + x * scalex - ssize / 4; " mm "; y0 + y * scaley - ssize / 4; "mm moveto"
   PRINT #1, "("; symbol$; ") show"
  END IF
  
   
  ' draw empty circle when necessary
  IF symbol$ = "0" OR symbol$ = "" THEN
   PRINT #1, x0 + x * scalex + ssize / 2; " mm "; y0 + y * scaley; "mm moveto"
   PRINT #1, x0 + x * scalex; " mm "; y0 + y * scaley; " mm "; ssize / 2; " mm 0 360 arc"
   PRINT #1, "stroke"
  END IF
  
  ' draw full circle when necessary
  IF symbol$ = "." THEN
   PRINT #1, x0 + x * scalex + ssize / 2; " mm "; y0 + y * scaley; "mm moveto"
   PRINT #1, x0 + x * scalex; " mm "; y0 + y * scaley; " mm "; ssize / 2; " mm 0 360 arc"
   PRINT #1, "fill"
  END IF
  
  'draw lines if necessary
  IF symbol$ = "-" AND washere% = 1 AND xsave >= 0 AND ysave >= 0 AND xsave <= xmax - xmin AND ysave < ymax - ymin THEN
   PRINT #1, "0.03 setlinewidth"
   PRINT #1, x0 + xsave * scalex; " mm "; y0 + ysave * scaley; " mm moveto"
   PRINT #1, x0 + x * scalex; " mm "; y0 + y * scaley; " mm lineto"
   PRINT #1, "stroke"
  END IF
  
  
  REM  error bars if necessary
  IF xbars > 0 AND d#(xbars) > 0 THEN
   PRINT #1, "1.5 setlinewidth"
   PRINT #1, x0 + (x + d#(xbars) / 2) * scalex; " mm "; y0 + y * scaley; " mm moveto"
   PRINT #1, x0 + (x - d#(xbars) / 2) * scalex; " mm "; y0 + y * scaley; " mm lineto"
   PRINT #1, "stroke"
   PRINT #1, "0.5 setlinewidth"
   tx = x + d#(xbars) / 2
  END IF
  IF ybars > 0 AND d#(ybars) > 0 THEN
   PRINT #1, "1.5 setlinewidth"
   PRINT #1, x0 + x * scalex; " mm "; y0 + (y + d#(ybars) / 2) * scaley; " mm moveto"
   PRINT #1, x0 + x * scalex; " mm "; y0 + (y - d#(ybars) / 2) * scaley; " mm lineto"
   PRINT #1, "stroke"
   PRINT #1, "0.5 setlinewidth"
   ty = y + d#(ybars) / 2
  END IF
  
 'ALL UNSHADED THINGS
 IF sshade < 1 THEN PRINT #1, "0 setgray"
  

 'draw x - scanrange (draw a thin line according to scan range)
 IF xlow > 0 AND xup > 0 THEN
  IF d#(xup) > xmax THEN d#(xup) = xmax
  IF d#(xlow) < xmin THEN d#(xlow) = xmin
  IF d#(xlow) < d#(xup) THEN
   PRINT #1, "0.1 setlinewidth"
   PRINT #1, x0 + (d#(xlow) - xmin) * scalex - 1; " mm "; y0 + y * scaley; " mm "; 1; " mm 270 450 arc"
   PRINT #1, "stroke"
   PRINT #1, x0 + (d#(xlow) - xmin) * scalex; " mm "; y0 + y * scaley; " mm  moveto"
   PRINT #1, x0 + (d#(xup) - xmin) * scalex; " mm "; y0 + y * scaley; " mm lineto"
   PRINT #1, "stroke"
   PRINT #1, x0 + (d#(xup) - xmin) * scalex + 1; " mm "; y0 + y * scaley; " mm "; 1; " mm 90 270 arc"
   PRINT #1, "stroke"
   PRINT #1, "0.5 setlinewidth"
   IF tx < d#(xup) THEN tx = d#(xup) - xmin
  END IF
 END IF

 'draw y - scanrange
 IF ylow > 0 AND yup > 0 THEN
  IF d#(yup) > ymax THEN d#(yup) = ymax
  IF d#(ylow) < ymin THEN d#(ylow) = ymin
  IF d#(ylow) < d#(yup) THEN
   PRINT #1, x0 + x * scalex; " mm "; y0 + (d#(ylow) - ymin) * scaley - 1; " mm "; 1; " mm 0 180 arc"
   PRINT #1, "stroke"
   PRINT #1, x0 + x * scalex; " mm "; y0 + (d#(ylow) - ymin) * scaley; " mm  moveto"
   PRINT #1, x0 + x * scalex; " mm "; y0 + (d#(yup) - ymin) * scaley; " mm lineto"
   PRINT #1, "stroke"
   PRINT #1, x0 + x * scalex; " mm "; y0 + (d#(yup) - ymin) * scaley + 1; " mm "; 1; " mm 180 360 arc"
   PRINT #1, "stroke"
   PRINT #1, "0.5 setlinewidth"
   IF ty < d#(yup) THEN ty = d#(yup) - ymin
  END IF
 END IF


 x$ = "": y$ = ""

 'draw comment  horizontal if option /tx is entered
 IF INSTR(LCASE$(COMMAND$), "/tx") > 0 THEN
  x$ = comment$
 END IF

 'draw comment  vertical if option /ty is entered
 IF INSTR(LCASE$(COMMAND$), "/ty") > 0 THEN
  y$ = comment$
 END IF

 'draw data from column in x direction
 IF dx > 0 THEN
  x$ = x$ + STR$(d#(dx))
 END IF

 'draw data from column in y direction
 IF dy > 0 THEN
   y$ = y$ + STR$(d#(dy))
 END IF


 
  'draw texts
  IF x$ <> "" THEN
   PRINT #1, "/Helvetica findfont"
   PRINT #1, "7 scalefont setfont"
   PRINT #1, x0 + tx * scalex + 3; " mm "; y0 + y * scaley - 1; " mm moveto"
   PRINT #1, "("; x$; ") show"
  END IF

  IF y$ <> "" THEN
   PRINT #1, "/Helvetica findfont"
   PRINT #1, "7 scalefont setfont"
   PRINT #1, x0 + x * scalex + 1; " mm "; y0 + ty * scaley + 3; " mm  moveto"
   PRINT #1, "90 rotate"
   PRINT #1, "("; y$; ") show"
   PRINT #1, "-90 rotate"
  END IF

END IF
 washere% = 1: xsave = x: ysave = y

110 IF EOF(2) = 0 GOTO 111
CLOSE 2
PRINT #1, "grestore" ' get xyposition from stack

133 IF INSTR(LCASE$(COMMAND$), "/s") > 0 THEN PRINT #1, "showpage"
CLOSE 1: PRINT "END PLOT"
IF INSTR(COMMAND$, "*") <> 0 THEN
 GOTO 999
END IF
END
333 FOR i = 1 TO 15: READ a$: PRINT a$: NEXT: END

SUB checkrepeat (file1$)
STATIC washere%

  IF INSTR(COMMAND$, "*") <> 0 THEN  'this is for cumulative action on more files
   IF washere% = 0 THEN
    'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /b > fact.dir"
      OPEN "i", 9, "fACT.dir"
   END IF
   IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del fact.dir": SHELL "del plot.dus": END
   INPUT #9, file1$
   IF LCASE$(file1$) = "fact.dir" THEN
    IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del fact.dir": SHELL "del plot.dus": END
    INPUT #9, file1$
   END IF
   IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del fact.dir": SHELL "del plot.dus": END
   IF washere% = 1 THEN
    IF INSTR(LCASE$(COMMAND$), "/s") > 0 THEN
     SHELL "copy plot.ps+plot.dus plot.dut"
     SHELL "copy plot.dut plot.ps"
     SHELL "del plot.dut"
    END IF
   END IF
   washere% = 1
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
   NEXT j: SHELL "echo text in data file too long": END
3 text$(j) = LEFT$(text$(j), i - 1)

END SUB

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

