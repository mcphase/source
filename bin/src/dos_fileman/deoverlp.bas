DECLARE SUB addpoint (er%, q!(), x!, y!, j%, n%)
DECLARE SUB checkiflinescross (true%, line1point1x!, line1point1y!, line1point2x!, line1point2y!, line2point1x!, line2point1y!, line2point2x!, line2point2y!, pintx!, pinty!)
DECLARE SUB detintersectionpts (square!(), squarechk!(), intersct!(), noofintersect%)
DECLARE SUB area (a1!, squarechk!())
DECLARE SUB checkifptisinside (yes%, x!, y!, squarechk!())
DECLARE SUB detoverlap (overlap!, square!(), squarechk!())
DECLARE SUB getsquare (d#(), col%, d2$(), col2%, squarepts!())
DECLARE SUB inputline (n!, d1#(), col1%, d2$(), col2%)
DECLARE SUB checkrepeat (file1$)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "DEOVERLP DEOVERLP DEOVERLP DEOVERLP DEOVERLP DEOVERLP DEOVERLP DEOVERLP DEOVERLP DEOVERLP DEOVERLP DEOVERLP DEOVERLP"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA " **********************************************************************"
DATA " program DEOVERLP: use as"
DATA " DEOVERLP *.* 13 [size=7][xbars=4][xlow=2 xup=4][ybars=5][ylow=2 yup=3][/tx/ty/dx=6/dy=8/b]"
DATA " ->DEOVERLaPs data stored in file *.*, xaxis is column1, yaxis is column 3"
DATA " files DEOVERLP.n are created such that the data points do not overlap in plot"
DATA " size: size of symbol(<1)"
DATA " xbars: x-error bar column   ybars: y-errorbars column"
DATA " xlow/xup: lower and upper measurement limits in x direction"
DATA " ylow/yup: lower and upper measurement limits in y direction"
DATA " option /dx=6 /dy=8 displays numbers on plot corresponding to col6 and 8"
DATA " option /tx, /ty displays comment text for each data point (horizontal,vertical)"
DATA " option /b allows overlap within data blocks separated by {multiline comments}"
DATA "**************************************************************************"
DIM tdum$(1000), text$(1000), d#(30), d2$(30), dchk#(30), dchk2$(30), square(10, 3), squarechk(10, 3)
DIM wfileheader%(999), headerpos(999)
maxoverlap = .1'maximal overlap of datapoints 50%
maxnooffiles% = 10 'maximal number of output files
noofoutputfiles% = 1: n% = 1: newfile% = 1
REM**********einlesen der dos umgebungsvariablen des graphs****************
DIM check%(8)
SHELL "set > deoverlp.dum": OPEN "i", 2, "deoverlp.dum"
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
FOR i% = 1 TO 8: IF check%(i%) <> 1 THEN PRINT "ERROR PROGRAM deoverlp - dos variable "; i%; " not set correctly": END
NEXT i%
CLOSE 2: SHELL "del deoverlp.dum"
'**************************************************************************
'**************************************************************************
' analyse command$
a$ = COMMAND$: a$ = LCASE$(LTRIM$(RTRIM$(a$)))
filename$ = LEFT$(a$, INSTR(a$, " ")): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))
999 CALL checkrepeat(filename$)

IF a$ = "" GOTO 133
colx% = ASC(LEFT$(a$, 1)) MOD 48: coly% = ASC(MID$(a$, 2, 1)) MOD 48
a$ = LTRIM$(RIGHT$(a$, LEN(a$) - INSTR(a$, " ")))
size = 0: IF INSTR(a$, "size=") > 0 THEN size = ASC(MID$(a$, INSTR(a$, "size=") + 5, 1)) MOD 48
xbars = 0: IF INSTR(a$, "xbars=") > 0 THEN xbars = ASC(MID$(a$, INSTR(a$, "xbars=") + 6, 1)) MOD 48
xlow = 0: IF INSTR(a$, "xlow=") > 0 THEN xlow = ASC(MID$(a$, INSTR(a$, "xlow=") + 5, 1)) MOD 48
xup = 0: IF INSTR(a$, "xup=") > 0 THEN xup = ASC(MID$(a$, INSTR(a$, "xup=") + 4, 1)) MOD 48
ybars = 0: IF INSTR(a$, "ybars=") > 0 THEN ybars = ASC(MID$(a$, INSTR(a$, "ybars=") + 6, 1)) MOD 48
ylow = 0: IF INSTR(a$, "ylow=") > 0 THEN ylow = ASC(MID$(a$, INSTR(a$, "ylow=") + 5, 1)) MOD 48
yup = 0: IF INSTR(a$, "yup=") > 0 THEN yup = ASC(MID$(a$, INSTR(a$, "yup=") + 4, 1)) MOD 48
dx = 0: IF INSTR(a$, "/dx=") > 0 THEN dx = ASC(MID$(a$, INSTR(a$, "/dx=") + 4, 1)) MOD 48
dy = 0: IF INSTR(a$, "/dy=") > 0 THEN dy = ASC(MID$(a$, INSTR(a$, "/dy=") + 4, 1)) MOD 48
'**************************************************************************

REM header einlesen*****************************************************
OPEN "i", 2, filename$
REM reflexe einlesen*******************************************
IF INSTR(LCASE$(COMMAND$), "/b") = 0 THEN  'no option b has been entered bbbbbbbbbbbb
111 CALL inputline(2, d#(), col%, d2$(), col2%)
    PRINT USING "-###%"; 100 * SEEK(2) / LOF(2); : LOCATE 24, 1
     IF col% = 0 THEN 'no data ...
      IF EOF(2) = 0 GOTO 111 ELSE GOTO 133
     END IF
    IF col% = -1 THEN  'new text section
     CALL headerinput(text$(), j, 2) 'a comment line has started and has to be treated
     FOR i% = 1 TO noofoutputfiles%: wfileheader%(i%) = 1: NEXT i%
     IF EOF(2) = 0 GOTO 111 ELSE GOTO 133
    END IF

 IF washere% > 0 THEN  'check datapoint
  CALL getsquare(d#(), col%, d2$(), col2%, square())
  FOR n% = 1 TO noofoutputfiles%
   OPEN "i", 3, "deoverlp." + LTRIM$(STR$(n%))
   CALL headerinput(tdum$(), jdum, 3)
   WHILE EOF(3) = 0
25  CALL inputline(3, dchk#(), colchk%, dchk2$(), colchk2%)
    IF colchk% = -1 THEN CALL headerinput(tdum$(), jdum, 3): IF EOF(3) = 0 GOTO 25 ELSE CLOSE 3: GOTO 33
    CALL getsquare(dchk#(), colchk%, dchk2$(), colchk2%, squarechk())
    CALL detoverlap(overlap, square(), squarechk())
    IF overlap > maxoverlap GOTO 32
   WEND: CLOSE 3: GOTO 33
32 CLOSE 3: NEXT n%
   'no place found for the datapoint --- create new file if possible
   IF noofoutputfiles% >= maxnooffiles% THEN
    n% = noofoutputfiles%
   ELSE
    noofoutputfiles% = noofoutputfiles% + 1: newfile% = 1
    n% = noofoutputfiles%
    wfileheader%(n%) = 1
    PRINT "creating file deoverlap."; n%
   END IF
 END IF

33 IF newfile% = 1 THEN
    newfile% = 0: OPEN "o", 1, "deoverlp." + LTRIM$(STR$(n%))'open correct output file
   ELSE
    OPEN "a", 1, "deoverlp." + LTRIM$(STR$(n%)) 'open correct output file
   END IF
 IF wfileheader%(n%) = 1 THEN  'output fileheader
  PRINT #1, "{"; : FOR i% = 1 TO j - 1: PRINT #1, text$(i%): NEXT i%: PRINT #1, text$(j); "}"
  wfileheader%(n%) = 0: headerpos(n%) = SEEK(1)
 END IF
 washere% = 1              'printout datapoint on correct file
 IF col% > 0 OR col2% > 0 THEN
  FOR i% = 1 TO col%:
nn$ = STR$(d#(i%)): MID$(nn$, INSTR(nn$, "D"), 1) = "E"
PRINT #2, " " + nn$; : NEXT i%:
  FOR i% = 1 TO col2%: PRINT #1, "{"; d2$(i%); "}"; : NEXT i%: PRINT #1,
 END IF
 CLOSE 1
 IF EOF(2) = 0 GOTO 111 'pick next data point
ELSE  'option /b is treated here bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
  'when /b is entered only whole datasets should be moved
  'move next unempty data block to deoverlp.dum
212 IF EOF(2) <> 0 GOTO 133
    PRINT USING "-###%"; 100 * SEEK(2) / LOF(2); : LOCATE 24, 1
    CALL inputline(2, d#(), col%, d2$(), col2%)
    IF col% = -1 THEN CALL headerinput(text$(), j, 2): GOTO 212
    OPEN "o", 4, "deoverlp.dum"
    PRINT #4, "{"; : FOR i% = 1 TO j - 1: PRINT #4, text$(i%): NEXT i%: PRINT #4, text$(j); "}"
    FOR i% = 1 TO col%: PRINT #4, d#(i%); : NEXT i%:
    FOR i% = 1 TO col2%: PRINT #4, "{"; d2$(i%); "}"; : NEXT i%: PRINT #4,
    WHILE EOF(2) = 0
     CALL inputline(2, d#(), col%, d2$(), col2%)
     IF col% = -1 GOTO 213
     FOR i% = 1 TO col%: PRINT #4, d#(i%); : NEXT i%:
     FOR i% = 1 TO col2%: PRINT #4, "{"; d2$(i%); "}"; : NEXT i%: PRINT #4,
    WEND
213 CLOSE 4 'now check with saved blocks
    'treat first data block
    IF washere% = 0 THEN washere% = 1: SHELL "del deoverlp.1": SHELL "copy deoverlp.dum deoverlp.1": GOTO 212
    FOR n% = 1 TO noofoutputfiles%
     OPEN "i", 3, "deoverlp." + LTRIM$(STR$(n%))
     CALL headerinput(tdum$(), jdum, 3)
     WHILE EOF(3) = 0
225   CALL inputline(3, dchk#(), colchk%, dchk2$(), colchk2%)
      IF colchk% = -1 THEN CALL headerinput(tdum$(), jdum, 3): IF EOF(3) = 0 GOTO 225 ELSE CLOSE 3: GOTO 233
      OPEN "i", 4, "deoverlp.dum"
      CALL headerinput(tdum$(), jdum, 4)
      WHILE EOF(4) = 0
       CALL inputline(4, d#(), col%, d2$(), col2%)
       CALL getsquare(dchk#(), colchk%, dchk2$(), colchk2%, squarechk())
       CALL getsquare(d#(), col%, d2$(), col2%, square())
       CALL detoverlap(overlap, square(), squarechk())
       IF overlap > maxoverlap THEN CLOSE 3, 4: GOTO 232
      WEND: CLOSE 4
     WEND: CLOSE 3: GOTO 233 '->ok
232 NEXT n%
   'no place found for the datapoint --- create new file if possible
   IF noofoutputfiles% >= maxnooffiles% THEN
    n% = noofoutputfiles%
   ELSE
    noofoutputfiles% = noofoutputfiles% + 1: newfile% = 1
    n% = noofoutputfiles%
    PRINT "creating file deoverlap."; n%
   END IF
233
  IF newfile% = 1 THEN
   'create new file
   SHELL "del deoverlp." + LTRIM$(STR$(n%))
   SHELL "copy deoverlp.dum deoverlp." + LTRIM$(STR$(n%))
   newfile% = 0
  ELSE
  'ok add data block to file deoverlp.n%
   SHELL "copy deoverlp." + LTRIM$(STR$(n%)) + "+deoverlp.dum deoverlp.du1"
   SHELL "del deoverlp." + LTRIM$(STR$(n%))
   SHELL "copy deoverlp.du1 deoverlp." + LTRIM$(STR$(n%))
   SHELL "del deoverlp.du1"
  END IF
  GOTO 212
END IF 'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
133 CLOSE 1: PRINT "END DEOVERLP": SHELL "del deoverlp.dum"
IF INSTR(COMMAND$, "*") <> 0 GOTO 999
END
333 FOR i = 1 TO 13: READ a$: PRINT a$: NEXT: END

SUB addpoint (er%, q(), x, y, j%, n%)
'add a point to the polygon q at position j%
' n%... noofpoints in q()
'check if new point is already in polygon
er% = 0
qi% = q(1, 0): DO
IF ABS(x - q(qi%, 1)) < 1E-10 AND ABS(y - q(qi%, 2)) < 1E-10 THEN
 er% = 1: GOTO 41
END IF
qi% = q(qi%, 0): LOOP WHILE qi% <> q(1, 0)


js% = j%
n% = n% + 1
q(n%, 0) = q(js%, 0) 'successor of point j% follows new point
q(js%, 0) = n% 'new point follows point j%
q(n%, 1) = x
q(n%, 2) = y
41 END SUB

SUB area (a1, squarechk())
dx = ABS(squarechk(1, 1) - squarechk(2, 1))
IF dx < .00001 THEN dx = ABS(squarechk(1, 1) - squarechk(3, 1))
dy = ABS(squarechk(1, 2) - squarechk(2, 2))
IF dy < .00001 THEN dy = ABS(squarechk(1, 2) - squarechk(3, 2))
a1 = dx * dy
END SUB

SUB checkiflinescross (true%, line1point1x, line1point1y, line1point2x, line1point2y, line2point1x, line2point1y, line2point2x, line2point2y, pintx, pinty)
'this sub checks, if two lines given by their endpoints intersect
'output: true%=1 ... lines cross
'        true%=0 ... lines do not cross
'        true%=-1... lines overlap one 'crossingpoint'


true% = 0

'setup vectors
vector1x = line1point2x - line1point1x
vector1y = line1point2y - line1point1y
vector2x = line2point2x - line2point1x
vector2y = line2point2y - line2point1y

'check if lines are really lines or rather close points
IF ABS(vector1x) + ABS(vector1y) < 1E-10 THEN true% = 0: GOTO 4354
IF ABS(vector2x) + ABS(vector2y) < 1E-10 THEN true% = 0: GOTO 4354

parallel% = 0
'check parallelity
IF ABS(vector1x * vector2y - vector1y * vector2x) < 1E-10 THEN
 vectorx = line2point2x - line1point1x
 vectory = line2point2y - line1point1y
 'check if parallel lines are offset
 IF ABS(vectorx * vector2y - vectory * vector2x) < 1E-10 THEN
  'check overlap of lines by checking if one of the two points
  'of line 2 lies on line one
  IF ABS(vector1x) > 1E-10 THEN
   lambdaa = (-line1point1x + line2point1x) / vector1x
   lambdab = (-line1point1x + line2point2x) / vector1x
  ELSE
   lambdaa = (-line1point1y + line2point1y) / vector1y
   lambdab = (-line1point1y + line2point2y) / vector1y
  END IF
  IF lambdaa >= -1E-10 AND lambdaa <= 1.00000001# THEN true% = -1: pintx = line2point1x: pinty = line2point1y: GOTO 4354
  IF lambdab >= -1E-10 AND lambdab <= 1.00000001# THEN true% = -1: pintx = line2point2x: pinty = line2point2y: GOTO 4354
  'check overlap of lines by checking if one of the two points
  'of line 1 lies on line 2
  IF ABS(vector2x) > 1E-10 THEN
   lambdaa = (-line2point1x + line1point1x) / vector2x
   lambdab = (-line2point1x + line1point2x) / vector2x
  ELSE
   lambdaa = (-line2point1y + line1point1y) / vector2y
   lambdab = (-line2point1y + line1point2y) / vector2y
  END IF
  IF lambdaa >= -1E-10 AND lambdaa <= 1.00000001# THEN true% = -1: pintx = line1point1x: pinty = line1point1y: GOTO 4354
  IF lambdab >= -1E-10 AND lambdab <= 1.00000001# THEN true% = -1: pintx = line1point2x: pinty = line1point2y: GOTO 4354
 END IF
 true% = 0: GOTO 4354
END IF

' calculate intersection point
det = -vector2x * vector1y + vector2y * vector1x
mu = (-vector1y * (line1point1x - line2point1x) + vector1x * (line1point1y - line2point1y)) / det
lambda = (vector2y * (line2point1x - line1point1x) - vector2x * (line2point1y - line1point1y)) / det
pintx = line2point1x + mu * vector2x
pinty = line2point1y + mu * vector2y
IF 0 < mu AND mu < 1 AND 0 < lambda AND lambda < 1 THEN true% = 1

4354 END SUB

SUB checkifptisinside (yes%, x, y, squarechk())
 'check if xypoint is in plygon squarechk()
'take a line through the xypoint and a distant point

 'check how often the polygon is crossed by this line
 crossings% = 0
j% = squarechk(1, 0): DO'check intersections with qi%
  i% = squarechk(j%, 0)
  CALL checkiflinescross(true%, x, y, 400, 200, squarechk(j%, 1), squarechk(j%, 2), squarechk(i%, 1), squarechk(i%, 2), pintx, pinty)
  IF true% = 1 THEN crossings% = crossings% + 1

j% = squarechk(j%, 0): LOOP WHILE j% <> squarechk(1, 0)


 IF INT(crossings% / 2) <> crossings% / 2 THEN yes% = 1 ELSE yes% = 0
 'if the number of crossings is odd then
 'YES, the point lies inside !!!

END SUB

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

SUB detintersectionpts (p(), q(), intersct(), n%)
' this sub determines the intersectionpoints of the polygons
' p() and q()
n% = 0: intersct(1, 0) = 1: 'initialize zero points polygon

pi% = p(1, 0): DO 'go around p()
line1point1x = p(pi%, 1)
line1point1y = p(pi%, 2)
line1point2x = p(p(pi%, 0), 1)
line1point2y = p(p(pi%, 0), 2)

'check if a point on q is equal to p - if yes move it a bit
qi% = q(1, 0): DO'check intersections with qi%
dx = q(qi%, 1) - p(pi%, 1)
dy = q(qi%, 2) - p(pi%, 2)
distance = SQR(dx * dx + dy * dy)
IF distance < .01 THEN q(qi%, 1) = q(qi%, 1) + .01: q(qi%, 2) = q(qi%, 2) + .01
dx = q(qi%, 1) - p(p(pi%, 0), 1)
dy = q(qi%, 2) - p(p(pi%, 0), 2)
distance = SQR(dx * dx + dy * dy)
IF distance < .01 THEN q(qi%, 1) = q(qi%, 1) + .01: q(qi%, 2) = q(qi%, 2) + .01
qi% = q(qi%, 0): LOOP WHILE qi% <> q(1, 0)

qi% = q(1, 0): DO'check intersections with qi%
line2point1x = q(qi%, 1)
line2point1y = q(qi%, 2)
line2point2x = q(q(qi%, 0), 1)
line2point2y = q(q(qi%, 0), 2)
 
 CALL checkiflinescross(true%, line1point1x, line1point1y, line1point2x, line1point2y, line2point1x, line2point1y, line2point2x, line2point2y, pintx, pinty)
 IF true% = 1 OR true% = -1 THEN  'if lines cross or overlap a point has been found
  CALL addpoint(er%, intersct(), pintx, pinty, n%, n%)
 END IF

qi% = q(qi%, 0): LOOP WHILE qi% <> q(1, 0)
pi% = p(pi%, 0): LOOP WHILE pi% <> p(1, 0)



END SUB

 SUB detoverlap (overlap, square(), squarechk())
'this sub determines the overlap of the two
'polygons square and squarechk
' overlap= max[area(square and squarechk)/(area(square),
'              area(square and squarechk)/(area(squarechk)]
DIM intersct(4, 3)

CALL detintersectionpts(square(), squarechk(), intersct(), nop%)


IF nop% < 4 THEN
 FOR i% = 1 TO 4
  CALL checkifptisinside(yes%, square(i%, 1), square(i%, 2), squarechk())
  IF yes% = 1 THEN CALL addpoint(er%, intersct(), square(i%, 1), square(i%, 2), 1, nop%)
 NEXT i%
 IF nop% < 4 THEN
  FOR i% = 1 TO 4
   CALL checkifptisinside(yes%, squarechk(i%, 1), squarechk(i%, 2), square())
   IF yes% = 1 THEN CALL addpoint(er%, intersct(), squarechk(i%, 1), squarechk(i%, 2), 1, nop%)
  NEXT i%
 END IF
END IF

'calculate areas

CALL area(a1, squarechk())
CALL area(a2, square())
CALL area(a3, intersct())

IF a1 = 0 OR a2 = 0 THEN
 overlap = 0
ELSE
 IF a3 / a2 > a3 / a1 THEN
  overlap = a3 / a2
 ELSE
  overlap = a3 / a1
 END IF
END IF
'PRINT "square"
'PRINT square(1, 1); square(2, 1); square(3, 1); square(4, 1)
'PRINT square(1, 2); square(2, 2); square(3, 2); square(4, 2)
'PRINT "squarechk"
'PRINT squarechk(1, 1); squarechk(2, 1); squarechk(3, 1); squarechk(4, 1)
'PRINT squarechk(1, 2); squarechk(2, 2); squarechk(3, 2); squarechk(4, 2)
'PRINT "intersections"
'PRINT intersct(1, 1); intersct(2, 1); intersct(3, 1); intersct(4, 1)
'PRINT intersct(1, 2); intersct(2, 2); intersct(3, 2); intersct(4, 2)


 END SUB

SUB getsquare (d#(), col%, d2$(), col2%, squarepts())
SHARED xmin, ymin, xmax, ymax, size, gheight, gwidth, xbars, ybars, xlow, xup, ylow, yup, colx%, coly%
'this sub determines the square of a data point used to determine
'the covered area
'squarepts(1-4,1-2)  1-4.... ptnumber, 1-2 x-y coordinates of the point
'squarepts(1-4,0) index of the following point, i.e. squarepts(1,0)=2
' (datenstruktur eines rings)
' (technically the square will be determined as sqxmin sqxmax sqymin sqymax
' and then inserted into squarepts()

'initialize:
sqxmin = 0: sqxmax = 0: sqymax = 0: sqymin = 0

comment$ = "": IF col2% > 0 THEN FOR i% = 1 TO col2%: comment$ = comment$ + d2$(i%): NEXT i%

REM pixel positions
 x = d#(colx%) - xmin: y = d#(coly%) - ymin
 tx = x: ty = y
 x0 = locx - gwidth / 2: y0 = locy - gheight / 2
 scalex = gwidth / (xmax - xmin): scaley = gheight / (ymax - ymin)

IF d#(colx%) <= xmax AND x >= 0 AND y >= 0 AND d#(coly%) < ymax THEN

   'determine size
   IF size > 0 THEN ssize = d#(size) ELSE ssize = SQR(gheight * gwidth) / 100: IF ssize < .6 THEN ssize = 1
   IF ssize < .4 THEN ssize = .4'lower limit


  'draw symbol
   sqxmin = x0 + x * scalex - ssize / 2
   sqxmax = x0 + x * scalex + ssize / 2
   sqymin = y0 + y * scaley - ssize / 2
   sqymax = y0 + y * scaley + ssize / 2
  
  REM  error bars if necessary
  IF xbars > 0 AND d#(xbars) > 0 THEN
   PRINT #1, "1.5 setlinewidth"
   sqxmax = x0 + (x + d#(xbars) / 2) * scalex
   sqxmin = x0 + (x - d#(xbars) / 2) * scalex
   tx = x + d#(xbars) / 2
  END IF
  IF ybars > 0 AND d#(ybars) > 0 THEN
   sqymax = y0 + (y + d#(ybars) / 2) * scaley
   sqymin = y0 + (y - d#(ybars) / 2) * scaley
   ty = y + d#(ybars) / 2
  END IF

 'draw x - scanrange (draw a thin line according to scan range)
 IF xlow > 0 AND xup > 0 THEN
  IF d#(xup) > xmax THEN d#(xup) = xmax
  IF d#(xlow) < xmin THEN d#(xlow) = xmin
  IF d#(xlow) < d#(xup) THEN
   sqxmin = x0 + (d#(xlow) - xmin) * scalex - 3
   sqxmax = x0 + (d#(xup) - xmin) * scalex
   IF tx < d#(xup) - xmin THEN tx = d#(xup) - xmin
  END IF
 END IF

 'draw y - scanrange
 IF ylow > 0 AND yup > 0 THEN
  IF d#(yup) > ymax THEN d#(yup) = ymax
  IF d#(ylow) < ymin THEN d#(ylow) = ymin
  IF d#(ylow) < d#(yup) THEN
   sqymin = y0 + (d#(ylow) - ymin) * scaley - 3
   sqymax = y0 + (d#(yup) - ymin) * scaley
   IF ty < d#(yup) - ymin THEN ty = d#(yup) - ymin
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
    sqxmax = sqxmax + LEN(x$)
  END IF

  IF y$ <> "" THEN
    sqymax = sqymax + LEN(y$)
  END IF

END IF
110  'set square
squarepts(1, 0) = 2
squarepts(1, 1) = sqxmin
squarepts(1, 2) = sqymin
squarepts(2, 0) = 3
squarepts(2, 1) = sqxmax
squarepts(2, 2) = sqymin
squarepts(3, 0) = 4
squarepts(3, 1) = sqxmax
squarepts(3, 2) = sqymax
squarepts(4, 0) = 1
squarepts(4, 1) = sqxmin
squarepts(4, 2) = sqymax

'IF ABS(d#(colx%) - .5) < .001 THEN STOP

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

