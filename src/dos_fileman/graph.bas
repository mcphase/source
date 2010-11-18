DECLARE SUB calcticklin (min!, max!, axislength!, n!, tick!(), m!, mintick!())
DECLARE SUB checktext (text$)
DECLARE SUB inputline (n!, d1#(), col1%)
DECLARE SUB headerinput (text$(), j!, n!)
DECLARE SUB pointdraw3d (x!, xmin!, xmax!, y!, ymin!, ymax!, z!, zmin!, zmax!, size!)
PRINT "GRAPH GRAPH GRAPH GRAPH GRAPH GRAPH GRAPH GRAPH GRAPH GRAPH GRAPH GRAPH GRAPH"
DATA " **********************************************************************"
DATA " program GRAPH:use as GRAPH [xaxis=blablabla] [yaxis=blanblA] [title=blabla]"
DATA " this program makes the scales and the axis TITLES FOR in the postscript"
DATA " file graph.ps                                                          "
DATA " file GRAPH.ps ... remember only when the option /s is added the        "
DATA " postscript file is finished using the 'showpage' command               "
DATA " you must have defined dos variables minx maxx miny maxy (range of graph)"
DATA " and locx locy heigth width (location and size of the GRAPH (in mm))     "
DATA "**************************************************************************"

REM**********einlesen der dos umgebungsvariablen des graphs****************
DIM check%(8)
SHELL "set > GRAPH.dum": OPEN "i", 1, "GRAPH.dum"
WHILE EOF(1) = 0
INPUT #1, A$: A$ = LCASE$(A$)
 IF LEFT$(A$, 4) = "minx" THEN xmin = VAL(RIGHT$(A$, LEN(A$) - INSTR(A$, "="))): check%(1) = 1
 IF LEFT$(A$, 4) = "maxx" THEN xmax = VAL(RIGHT$(A$, LEN(A$) - INSTR(A$, "="))): check%(2) = 1
 IF LEFT$(A$, 4) = "miny" THEN ymin = VAL(RIGHT$(A$, LEN(A$) - INSTR(A$, "="))): check%(3) = 1
 IF LEFT$(A$, 4) = "maxy" THEN ymax = VAL(RIGHT$(A$, LEN(A$) - INSTR(A$, "="))): check%(4) = 1
 IF LEFT$(A$, 5) = "width" THEN gwidth = VAL(RIGHT$(A$, LEN(A$) - INSTR(A$, "="))): check%(5) = 1
 IF LEFT$(A$, 6) = "height" THEN gheight = VAL(RIGHT$(A$, LEN(A$) - INSTR(A$, "="))): check%(6) = 1
 IF LEFT$(A$, 4) = "locx" THEN locx = VAL(RIGHT$(A$, LEN(A$) - INSTR(A$, "="))): check%(7) = 1
 IF LEFT$(A$, 4) = "locy" THEN locy = VAL(RIGHT$(A$, LEN(A$) - INSTR(A$, "="))): check%(8) = 1
WEND
FOR i% = 1 TO 8: IF check%(i%) <> 1 THEN PRINT "ERROR PROGRAM GRAPH - dos variable "; i%; " not set correctly": GOTO 333
NEXT i%
CLOSE 1: SHELL "del GRAPH.dum"
'**************************************************************************
' analyse command$

A$ = COMMAND$: A$ = LCASE$(LTRIM$(RTRIM$(A$))): Aa$ = LTRIM$(RTRIM$(COMMAND$))
xtext$ = "": IF INSTR(A$, "xaxis=") > 0 THEN xtext$ = MID$(Aa$, INSTR(A$, "xaxis=") + 6)
ytext$ = "": IF INSTR(A$, "yaxis=") > 0 THEN ytext$ = MID$(Aa$, INSTR(A$, "yaxis=") + 6)
title$ = "": IF INSTR(A$, "title=") > 0 THEN title$ = MID$(Aa$, INSTR(A$, "title=") + 6)
CALL checktext(xtext$)
CALL checktext(ytext$)
CALL checktext(title$)
'**************************************************************************

'**************************************************************************
REM postscriptfile kreieren***********************************************

REM text ausdrucken *****************************************************
OPEN "a", 1, "GRAPH.ps"
PRINT #1, "/mm {72 mul 25.4 div} bind def"

PRINT #1, "/RightAllign{"
PRINT #1, "dup"
PRINT #1, "stringwidth pop"
PRINT #1, "neg 0 rmoveto"
PRINT #1, "show"
PRINT #1, "}bind def"

PRINT #1, "/Helvetica findfont"
PRINT #1, "12 scalefont setfont"
PRINT #1, "0 setgray"

'rem locate to title
IF title$ <> "" THEN PRINT #1, locx - gwidth / 2; " mm "; locy + gheight / 2 + 10; " mm moveto": PRINT #1, "("; title$; ") show"

PRINT #1, "/Helvetica findfont"
PRINT #1, "12 scalefont setfont"
PRINT #1, "0 setgray"

'print date and time on right top
PRINT #1, 205; " mm "; 290; " mm moveto": PRINT #1, "("; DATE$; " "; TIME$; ") RightAllign":

'rem locate to xtext
IF xtext$ <> "" THEN PRINT #1, locx - LEN(xtext$) * 1; " mm "; locy - gheight / 2 - 10; " mm moveto": PRINT #1, "("; xtext$; ") show"
'rem locate to ytext vertical writing
IF ytext$ <> "" THEN
   PRINT #1, locx - gwidth / 2 - 25; " mm "; locy - LEN(ytext$) * 1; " mm moveto"
   PRINT #1, "90 rotate"
   PRINT #1, "("; ytext$; ") RightAllign"
   PRINT #1, "-90 rotate"
END IF


'**************************************************************************

REM koordinatensystem ------------------------------------------------
x0 = locx - gwidth / 2: y0 = locy - gheight / 2
scalex = gwidth / (xmax - xmin): scaley = gheight / (ymax - ymin)
PRINT #1, x0; " mm "; y0; " mm moveto": PRINT #1, x0; " mm "; y0 + gheight; " mm lineto"
PRINT #1, "stroke"
PRINT #1, x0; " mm "; y0 + gheight; " mm "; " moveto": PRINT #1, x0 + gwidth; " mm "; y0 + gheight; " mm "; " lineto"
PRINT #1, "stroke"
PRINT #1, x0 + gwidth; " mm "; y0 + gheight; " mm "; " moveto"
PRINT #1, x0 + gwidth; " mm "; y0; " mm "; " lineto"
PRINT #1, "stroke"
PRINT #1, x0 + gwidth; " mm "; y0; " mm "; " moveto"
PRINT #1, x0; " mm "; y0; " mm "; " lineto"
PRINT #1, "stroke"
'----------------------------------------------------------------------
'               SCALES SCALES
'----------------------------------------------------------------------
DIM xtick(100), ytick(100)
DIM xminortick(110), yminortick(110)
CALL calcticklin(xmin, xmax, gwidth, nx, xtick(), nxm, xminortick())
CALL calcticklin(ymin, ymax, gheight, ny, ytick(), nym, yminortick())
x0 = locx - gwidth / 2: y0 = locy - gheight / 2

REM xticks und scale
FOR i = 1 TO nx
x = x0 + (xtick(i) - xmin) * scalex
PRINT #1, x + 5; " mm "; y0 - 5; " mm "; " moveto"
PRINT #1, "(";
PRINT #1, USING "##.##^^^^"; xtick(i);
PRINT #1, ") RightAllign"
PRINT #1, x; " mm "; y0; " mm "; " moveto"
PRINT #1, x; " mm "; y0 + gheight / 50; " mm "; " lineto"
PRINT #1, "stroke"
PRINT #1, x; " mm "; y0 + gheight; " mm "; " moveto"
PRINT #1, x; " mm "; y0 + gheight - gheight / 50; " mm "; " lineto"
PRINT #1, "stroke"
   
NEXT i
   REM xminorticks und scale
        FOR Ii = 1 TO nxm
         IF xminortick(Ii) < xmax THEN
          x = x0 + (xminortick(Ii) - xmin) * scalex
          PRINT #1, x; " mm "; y0; " mm "; " moveto"
          PRINT #1, x; " mm "; y0 + gheight / 100; " mm "; " lineto"
          PRINT #1, "stroke"
          PRINT #1, x; " mm "; y0 + gheight; " mm "; " moveto"
          PRINT #1, x; " mm "; y0 + gheight - gheight / 100; " mm "; " lineto"
          PRINT #1, "stroke"
         END IF
        NEXT Ii


' yticks und scale
FOR i = 1 TO ny
y = y0 + (ytick(i) - ymin) * scaley
PRINT #1, x0; " mm "; y; " mm "; " moveto"
PRINT #1, x0 + gwidth / 50; " mm "; y; " mm "; " lineto"
PRINT #1, "stroke"
PRINT #1, x0 + gwidth; " mm "; y; " mm "; " moveto"
PRINT #1, x0 + gwidth - gwidth / 50; " mm "; y; " mm "; " lineto"
PRINT #1, "stroke"
PRINT #1, x0; " mm "; y - 1; " mm "; " moveto"
PRINT #1, "(";
PRINT #1, USING "##.####^^^^"; ytick(i);
PRINT #1, ") RightAllign"
NEXT i
    ' yminorticks und scale
        FOR Ii = 1 TO nym
         IF yminortick(Ii) < ymax THEN
          y = y0 + (yminortick(Ii) - ymin) * scaley
          PRINT #1, x0; " mm "; y; " mm "; " moveto"
          PRINT #1, x0 + gwidth / 100; " mm "; y; " mm "; " lineto"
          PRINT #1, "stroke"
          PRINT #1, x0 + gwidth; " mm "; y; " mm "; " moveto"
          PRINT #1, x0 + gwidth - gwidth / 100; " mm "; y; " mm "; " lineto"
          PRINT #1, "stroke"
         END IF
        NEXT Ii

'----------------------------------------------------------------------
PRINT #1, x0 - 20; " mm "; y0 - 10; " mm moveto"

IF INSTR(A$, "/s") > 0 THEN PRINT #1, "showpage"
CLOSE 1: PRINT "END GRAPH"
END
333 FOR i = 1 TO 9: READ A$: PRINT A$: NEXT i: END

SUB calcticklin (min, max, axislength, n, tick(), m, mintick())
' this sub takes the minimum and the maximum and the length[mm] of an axis
' and calculates a convenient number of ticks n and their
' position tick(1....n)

 IF min >= max THEN PRINT "ERROR GRAPH sub calcticklin min>=max": END
 IF axislength < 30 THEN PRINT "ERROR GRAPH axislength<30mm": END

 delta = (max - min) / axislength * 20: exponent% = INT(LOG(delta) / LOG(10))
 i% = (1 + INT(delta / 10 ^ exponent%))
 delta = i% * 10 ^ exponent%'round delta to sensible number
 

 IF min <> 0 THEN exponentmin% = INT(LOG(ABS(min)) / LOG(10)) ELSE exponentmin% = -100
 IF max <> 0 THEN exponentmax% = INT(LOG(ABS(max)) / LOG(10)) ELSE exponentmax% = -100
 expon% = exponentmin%: IF exponentmin% < exponentmax% THEN expon% = exponentmax%

 n = 1: tick(1) = min: IF expon% < eponent% + 5 THEN tick(1) = INT(min / delta) * delta
 IF tick(1) < min THEN tick(1) = tick(1) + delta

 WHILE tick(n) + delta < max
 n = n + 1: tick(n) = tick(n - 1) + delta
 WEND: tick(n + 1) = tick(n) + delta

mdelta = .1
 IF i% = 1 THEN mdelta = delta / 10
 IF i% = 2 THEN mdelta = delta / 10
 IF i% = 3 THEN mdelta = delta / 6
 IF i% = 4 THEN mdelta = delta / 8
 IF i% = 5 THEN mdelta = delta / 10
 IF i% = 6 THEN mdelta = delta / 6
 IF i% = 7 THEN mdelta = delta / 7
 IF i% = 8 THEN mdelta = delta / 8
 IF i% = 9 THEN mdelta = delta / 9
 IF i% = 10 THEN mdelta = delta / 10

 m = 1: mintick(1) = tick(1)
 WHILE mintick(m) + mdelta < max
 m = m + 1: mintick(m) = mintick(m - 1) + mdelta
 WEND: mintick(m + 1) = mintick(m) + mdelta






END SUB

SUB checktext (text$)

IF INSTR(LCASE$(text$), "title=") > 0 THEN text$ = LEFT$(text$, INSTR(LCASE$(text$), "title=") - 1)
IF INSTR(LCASE$(text$), "yaxis=") > 0 THEN text$ = LEFT$(text$, INSTR(LCASE$(text$), "yaxis=") - 1)
IF INSTR(LCASE$(text$), "xaxis=") > 0 THEN text$ = LEFT$(text$, INSTR(LCASE$(text$), "xaxis=") - 1)

END SUB

