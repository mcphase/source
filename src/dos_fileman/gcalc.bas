DECLARE SUB inputline (n!, D1#(), col1%, D2$(), colII%)
DECLARE SUB analizecommand (c0#, b#, ra#, ri#, file1$, file2$)
DECLARE SUB gobackonestrg (n!, ala$)
DECLARE SUB headerinput (text$(), j!, n!)
version$ = "<GCALC v. 09.06.10 "
PRINT version$; COMMAND$;
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA " programm to calibrate gap d in full Compensated microdilatometer "
DATA " command format: "
DATA "  GCALC 3.2 9.794 2 6.25  *.ttc silver*.lit  [-r]"
DATA "   parameters mean: c0=3.2pF, b=9.794mm, ri=2, ra=6.25"
DATA "   option -r:  calculate (gap-k0) instead of gap"
DATA "  means function from file *.ttc (time vs temperature vs capacity)    "
DATA "  the capacity c is used to calculate the gap g using                 "
DATA "  the formula:c(T)=ca(T)-ci(T)
DATA ""
DATA "  ca(T)= Aa(T)*epsilon_0 / g(T)* 2*(1-sqr(1-gamma^2))/gamma^2"
DATA "  ci(T)= Ai(T)*epsilon_0 / g(T)* 2*(1-sqr(1-gammi^2))/gammi^2"
DATA "  Aa(T)=Aa0*(1+dl/l_Ag(T)) ... capacitance area"
DATA "  Ai(T)=Ai0*(1+dl/l_Ag(T)) ... capacitance area"
DATA "           note: we do not use ^2 in these formula, because the gap "
DATA "                 expands due to thermal expansion of silver bearings, too"
DATA "  gammi= ri/b*(k(T)/g(T)-1)"
DATA "  gamma= ra/b*(k(T)/g(T)-1)"
DATA "  k(T)=k0 ... k0:pivot distance at 300K "
DATA "      ...c0:=epsilon_0 (Aa(300K)-Ai(300K))/k0 ... corresponding capacitance"
DATA "  b#  'distance between center of capacitor and pivot in mm    "
DATA "  ra#  'outer plate radius in mm "
DATA "  ri# 'inner plate radius in mm"
DATA " "
DATA "Some parameters for different dilatometers:"
DATA "Dilatometer              b/ra/ri (mm)"
DATA "Dresden22-96            9.794/6.25/2.0"
DATA "Prague22-02             9.794/6.25/2.667"
DATA "Vienna20-02             9.794/6.25/2.333"
DATA ""
DATA "file format"
DATA " { header: this is the file header delimited         "
DATA "   by brackets after this header there follow 3 or more data columns     "
DATA "   time[s]  T[K]  capacitance[pF] }                                   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "*************************************************************************"
DIM text1$(300), text2$(300), text3$(300), D1#(30), D2$(30), D2A#(30), d2b#(30), D3A#(30), d3b#(30)

'ooooooooooooooooooooooooo OPEN FILES oooooooooooooooooooooooooooooooooooo
999 CALL analizecommand(c0#, b#, r#, ri#, file1$, file2$)
x1% = 2  '.... Temperature column  in *.ttc
y1% = 3  '.... capacity column in *.ttc
x2% = 1  '.... temperature column in silver*.lit
y2% = 2  '....dl/l column in silver*.lit
x3% = 1  '.... temperature column in saphir*.lit
y3% = 2  '....dl/l column in saphir*.lit

REM SHELL "copy " + file2$ + " GCALCAg.lit"

ddd = 1E+10
OPEN "i", 1, file2$
CALL headerinput(text1$(), j1, 1)' open file1 and input header
 WHILE EOF(1) = 0
 CALL inputline(1, D1#(), col1%, D2$(), colII%)
 IF col1% < 2 THEN INPUT "error gcalc:reading less than two columns in lit file"; dd$: END
 IF ABS(D1#(1) - 300) < ddd THEN ddd = ABS(D1#(1) - 300): shiftc = D1#(2)
 WEND
CLOSE 1

OPEN "i", 1, file2$
CALL headerinput(text1$(), j1, 1)' open file1 and input header
OPEN "o", 2, "GCALCAg.lit"
 PRINT #2, "{silver file}"
 WHILE EOF(1) = 0
  CALL inputline(1, D1#(), col1%, D2$(), colII%)
  D1#(2) = D1#(2) - shiftc
  PRINT #2, D1#(1); " "; D1#(2)
 WEND
CLOSE 1
CLOSE 2
file2$ = "GCALCAg.lit"

OPEN "i", 1, file1$: CALL headerinput(text1$(), j1, 1)' open file1 and input header
OPEN "i", 2, file2$: CALL headerinput(text2$(), j2, 2)' open file2 and input header

OPEN "o", 4, "GCALC.gca"                         ' open output file
'oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
Aa0# = 3.14159265# * r# * r#
Ai0# = 3.14159265# * ri# * ri#
eps0# = 8.854188000000001D-03 '[pF/mm]
k0# = eps0# * (Aa0# - Ai0#) / c0#


'hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh 'print file header of output file hhhhhhhhhhhh
PRINT #4, "{";
FOR i = 1 TO j1: PRINT #4, text1$(i): NEXT
PRINT #4, version$; COMMAND$; " "; DATE$; " "; TIME$;
IF INSTR(LCASE$(COMMAND$), "-r") > 0 THEN PRINT #4, "using option -r";
PRINT #4, "(pivot pos k0="; k0#; "mm c0#="; c0#; "pF)GCALC> "
PRINT #4, "   (((file "; file2$; " information: ";
FOR i = 1 TO j2: PRINT #4, LCASE$(text2$(i)): NEXT
PRINT #4, "   )))}"
'hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh

'******************************************************************
REM input columns of #1 and #2 and #3 write result to #4
'******************************************************************

    'input a data point on #1
IF EOF(1) <> 0 GOTO 5 ELSE CALL inputline(1, D1#(), col1%, D2$(), colII%)

   'input first 2 values on file 2
IF EOF(2) <> 0 GOTO 5 ELSE CALL inputline(2, D2A#(), col2%, D2$(), colII%)
IF EOF(2) <> 0 GOTO 5 ELSE CALL inputline(2, d2b#(), col2%, D2$(), colII%)
  

' decide wether values in file 2 increase or decrease
sign2% = SGN(d2b#(x2%) - D2A#(x2%))

' input d1 values until d1#(x1%)*sign2% is bigger than d2a#(x2%)*sign2%
WHILE D1#(x1%) * sign2% < D2A#(x2%) * sign2%
      IF EOF(1) <> 0 GOTO 5 ELSE CALL inputline(1, D1#(), col1%, D2$(), colII%)
WEND

' \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \
'> > > > > > > > > > > > > >CALIBRATION LOOP > > > > > > > > > > > > > > > >
' / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /

4 LOCATE 23, 1: PRINT USING "###%"; 100 * SEEK(1) / LOF(1); : LOCATE 24, 1
'-----------------------------------------------------------------------
'input d2a- d2b values going backward in file until d1 is in the interval
      WHILE D1#(x1%) * sign2% < D2A#(x2%) * sign2%
               
                 ' go back three strings in file2
  FOR cr% = 1 TO 3: CALL gobackonestrg(2, ala$)
  NEXT cr%
   IF INSTR(ala$, "}") <> 0 THEN INPUT #2, a$: INPUT #2, a$: INPUT #2, a$: PRINT "out of range datapoint x="; D1#(x1%); "skipped": GOTO 15
                              'we are at the beginning of the file
                            
 CALL inputline(2, D2A#(), col2%, D2$(), colII%): CALL inputline(2, d2b#(), col2%, D2$(), colII%)
                        ' input d2a-d2b interval of before
        WEND
'------------------------------------------------------------------------
'........................................................................
'input d2b values until d1#(x1%) is between d2a#(x2%) and d2b#(x2%)
WHILE d2b#(x2%) * sign2% < D1#(x1%) * sign2%

  IF EOF(2) <> 0 THEN PRINT "out of range datapoint x="; D1#(x1%); "skipped": GOTO 15
                      'if d1 exceeds d2-range get next d1 value
  FOR coll% = 1 TO col2%: D2A#(coll%) = d2b#(coll%): NEXT
  CALL inputline(2, d2b#(), col2%, D2$(), colII%)
WEND
'........................................................................
'!!!!!!!!!!!!!!!!!!!!! DO GAPCALIBRATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
' linear interpolate d2a and d2b
dllAg# = D2A#(y2%) + (d2b#(y2%) - D2A#(y2%)) / (d2b#(x2%) - D2A#(x2%)) * (D1#(x1%) - D2A#(x2%))

kT# = k0# 'removed * (1 + dllAg#) 'also the pivot changes with temperature - however
                           ' this compensates with the fact that we look
                           ' at a virtual gap = real distance /(1+dllAg#)
AaT# = Aa0# * (1 + dllAg#) ' square removed (1+dllag)nicht quadrieren, da capacitaetssensor
AiT# = Ai0# * (1 + dllAg#) ' einerseits durch vergroessern der platten mit
                           'steigender T um (1+dllag)^2 zu grosse cap werte
                           'liefert, andererseits aber wegen der thermischen
                           'verkuerzung des plattenabstands um (1-dllag) zu
                           ' kleine cap werte liefert - beide effekte zusammen
                           ' sollten der korrekturfaktor (1+dllag) am
                           ' besten beschreiben
' DO GapCalibration
c# = D1#(y1%): ccalc# = 1000!
IF c# <= 0 THEN PRINT "datapoint T="; D1#(x1%); "K skipped because C="; c#; "pF": GOTO 15
'initialize d#
d# = (AaT# - AiT#) * eps0# / c#: dstep# = .001
IF d# < kT# * r# / (b# + r#) THEN d# = 1.1 * kT# * r# / (b# + r#)
'iterate d#
WHILE ABS(c# - ccalc#) > 1E-09 AND dstep# / d# > 1E-15
gamma# = r# / b# * (kT# - d#) / d#
gamma2# = gamma# * gamma#
ccalc# = AaT# * eps0# / d# * 2 / gamma2# * (1 - SQR(1 - gamma2#))
gamma# = ri# / b# * (kT# - d#) / d#
gamma2# = gamma# * gamma#
ccalc# = ccalc# - AiT# * eps0# / d# * 2 / gamma2# * (1 - SQR(1 - gamma2#))

d# = d# + SGN(ccalc# - c#) * dstep#
IF oldsign <> SGN(ccalc# - c#) THEN oldsign = SGN(ccalc# - c#): dstep# = dstep# / 10
'PRINT d#, c#, ccalc#
WEND

IF ABS(c# - ccalc#) > .0001 THEN PRINT "error gcalc: iteration fails at datapoint T="; D1#(x1%); "K C="; c#; "pF": END

'option -r calculates gap relatively to c0
IF INSTR(LCASE$(COMMAND$), "-r") > 0 THEN d# = d# - k0#

D1#(y1%) = d#


'save the datapoint
FOR coll% = 1 TO col1%:
nn$ = STR$(D1#(coll%)): IF INSTR(nn$, "D") > 0 THEN MID$(nn$, INSTR(nn$, "D"), 1) = "E"
PRINT #4, " " + nn$; : NEXT: PRINT #4,
'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

' input next d1 value

15 IF EOF(1) <> 0 GOTO 5 ELSE CALL inputline(1, D1#(), col1%, D2$(), colII%)

GOTO 4
' / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /
'< < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
' \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \

5 CLOSE 1, 2, 3, 4
SHELL "copy GCALC.gca " + file1$
SHELL "del GCALC.gca"
SHELL "del GCALCAg.lit"
PRINT
PRINT "GCALC "; COMMAND$; " pivot pos k0="; k0#; "mm c0="; c0#; "pF";
IF INSTR(LCASE$(COMMAND$), "-r") > 0 THEN PRINT "using option -r";
PRINT "b="; b#; "ra="; r#; " ri="; ri#; "mm)"
PRINT "GCALC>"

IF INSTR(COMMAND$, "*") <> 0 GOTO 999 'do the same with more files
                                                             
END
333 FOR i = 1 TO 19: READ a$: PRINT a$: NEXT: INPUT ala
    FOR i = 1 TO 19: READ a$: PRINT a$: NEXT: END

SUB analizecommand (c0#, b#, ra#, ri#, file1$, file2$)
STATIC washere%
'*****************************************************************
' this sub analizes the command$ and detects the two filenames and
' 2 columns for each file
'*****************************************************************
a$ = LCASE$(RTRIM$(LTRIM$(COMMAND$)))




 'detect c0:
 i = INSTR(a$, " "): c0# = VAL(MID$(a$, 1, i)): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - i))
 i = INSTR(a$, " "): b# = VAL(MID$(a$, 1, i)): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - i))
 i = INSTR(a$, " "): ri# = VAL(MID$(a$, 1, i)): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - i))
 i = INSTR(a$, " "): ra# = VAL(MID$(a$, 1, i)): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - i))

 'detect file *.ttc:
 i = INSTR(a$, " "): file1$ = LEFT$(a$, i): a$ = LTRIM$(RIGHT$(a$, LEN(a$) - i))
  IF INSTR(COMMAND$, "*") <> 0 THEN
   IF washere% = 0 THEN
      washere% = 1
      'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /b > GCALC.dir"
      OPEN "i", 9, "GCALC.dir"
   END IF
   IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del GCALC.dir": END
   INPUT #9, file1$
   IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del GCALC.dir": END
  END IF

 'detect file *.lit:
 
 file2$ = a$
 i = INSTR(a$, " "): IF i > 0 THEN file2$ = LEFT$(a$, i)

check$ = LCASE$(RIGHT$(RTRIM$(file1$), 4))
IF check$ = ".rcp" OR check$ = ".mrc" THEN
   SHELL "echo you should never ever change data in *.rcp files": PLAY "dgdgdg"
222 INPUT "do you really want to continue (Y/N)"; ala$: IF LCASE$(ala$) = "n" THEN END
   IF LCASE$(ala$) <> "y" GOTO 222
END IF

END SUB

SUB gobackonestrg (n, ala$)
'this sub sets the filepointer of file #n one string back and
' puts the corresponding string into ala$
filepointer = SEEK(n)
FOR x = 2 TO filepointer STEP 30
SEEK #n, filepointer - x
fpn = -1: INPUT #n, a$
122 WHILE SEEK(n) < filepointer
       fpn = SEEK(n):  INPUT #n, a$
    WEND
IF fpn > 0 GOTO 233
NEXT x: x = filepointer - 1: SEEK #n, 1: GOTO 122

233 SEEK #n, fpn: ala$ = a$
END SUB

SUB headerinput (text$(), j, n)
'*****************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'*****************************************************************

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

SUB inputline (n!, D1#(), col1%, D2$(), colII%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in D1#(1...col1%)
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
colII% = 0
WHILE klauf% < klzu% AND klauf% > 0   'take out closed bracketed expressions
 colII% = colII% + 1
 D2$(colII%) = MID$(ala$, klauf% + 1, klzu% - klauf% - 1)
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
    D1#(col1%) = VAL(LEFT$(ala$, INSTR(ala$, " ")))
    ala$ = LTRIM$(RIGHT$(ala$, LEN(ala$) - INSTR(ala$, " ")))
 WEND
END IF
END SUB

