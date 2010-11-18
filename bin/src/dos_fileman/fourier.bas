DECLARE SUB analizecommand (file1$, colx%, xmax!, deltax!)
DECLARE SUB transform (freqx!, x!, f!, sumreal#, sumimag#, xmax!, deltax!)
DECLARE SUB gobackonestrg (n!, ALA$)
DECLARE SUB inputline (n!, d1#(), col1%)
DECLARE SUB headerinput (text$(), j!, n!)
PRINT "FOURIER FOURIER FOURIER FOURIER FOURIER FOURIER FOURIER FOURIER FOURIER"
IF LTRIM$(COMMAND$) = "" GOTO 333
DATA "*************************************************************************"
DATA "  program FOURIER - use it like"
DATA "         FOURIER *.* 12 xmax=3.47 deltax=0.1"
DATA "(means fouriertransform functions coln(x)=(other columns) x=col1 "
DATA " ---> in the range x=[-3.47;3.47] y=[-10,10] datapoint sampling is done  "
DATA " intervals deltax and linear interpolation"
DATA " the complex fouriertransform is written into files"
DATA " FOURIER.re and FOURIER.im []..[]..xcolumn[1/[x]]..[]..."
DATA "
DATA " format of file:                                                         "
DATA " { header: this is the                                                   "
DATA "   file header delimited                                                 "
DATA "   by brackets after this header there follow 3 or more data columns }   "
DATA " 11 3.14235 65367                                                        "
DATA "  .    .     .                                                           "
DATA "  .    .     .                                                           "
DATA "  .    .     .    .  .   .                                               "
DATA "  .    .     .                                                           "
DATA " 32 2412.34 324.2                                                        "
DATA "*************************************************************************"
DIM text$(300), d1a#(30), d1b#(30), sumreal#(30), sumimag#(30)


' analyse command string aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
999 CALL analizecommand(filename$, colx%, xmax, deltax)
'Aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)

' open output file and write fileheader
OPEN "o", 2, "FOURIER.re"
PRINT #2, "{"; : FOR iii = 1 TO j: PRINT #2, text$(iii): NEXT
PRINT #2, DATE$; " "; TIME$;
PRINT #2, "fouriertransform real part by command: "; COMMAND$; "}"
OPEN "o", 3, "FOURIER.im"
PRINT #3, "{"; : FOR iii = 1 TO j: PRINT #3, text$(iii): NEXT
PRINT #3, DATE$; " "; TIME$;
PRINT #3, "fouriertransform imaginary part by command: "; COMMAND$; "}"

 FOR freqx = -.5 / deltax TO .5 / deltax STEP .5 / xmax
 x = -xmax: PRINT USING "###%"; (freqx + .5 / deltax) * deltax * 100;
 'input first 2 values on file          
 IF EOF(1) <> 0 GOTO 5 ELSE CALL inputline(1, d1a#(), col1%)
 IF EOF(1) <> 0 GOTO 5 ELSE CALL inputline(1, d1b#(), col1%)
 ' decide wether values in file increase or decrease
 sign% = SGN(d1b#(colx%) - d1a#(colx%))

4 'input d1a- d1b values going backward in file until d1 is in the interval
      WHILE x * sign% < d1a#(colx%) * sign%
              
                 ' go back three strings in file1
  FOR cr% = 1 TO 3: CALL gobackonestrg(1, ALA$)
  NEXT cr%
   IF INSTR(ALA$, "}") <> 0 THEN INPUT #1, A$: INPUT #1, A$: INPUT #1, A$: PRINT "out of range datapoint x="; x; "linear extrapolated": GOTO 15
                               'we are at the beginning of the file
                            
  CALL inputline(1, d1a#(), col1%): CALL inputline(1, d1b#(), col1%)
                         ' input d1a-d1b interval of before
         WEND
 '------------------------------------------------------------------------
 'input d1b values until d1#(colx%) is between d1a#(colx%) and d1b#(colx%)
 WHILE d1b#(colx%) * sign% < x * sign%
 
   IF EOF(1) <> 0 THEN PRINT "out of range datapoint x="; x; "linear extrapolated": GOTO 15
            'if d1 exceeds d1-range get next d1 value
   FOR coll% = 1 TO col1%: d1a#(coll%) = d1b#(coll%): NEXT
   CALL inputline(1, d1b#(), col1%)
 WEND
 '........................................................................

 '!!!!!!!!!!!!!!!!!!!!! DO ADDITION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ' linear interpolate d1a and d1b
15 FOR coll% = 1 TO col1%
  IF coll% <> colx% THEN
   interpoly = d1a#(coll%) + (d1b#(coll%) - d1a#(coll%)) / (d1b#(colx%) - d1a#(colx%)) * (x - d1a#(colx%))
   ' DO transform
   CALL transform(freqx, x, interpoly, sumreal#(coll%), sumimag#(coll%), xmax, deltax)
  END IF
  NEXT coll%
 '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ' input next d1 value
 x = x + deltax: IF x < xmax + deltax * .1 GOTO 4
5 CLOSE 1
FOR coll% = 1 TO col1%  'write output file
  IF coll% <> colx% THEN
   PRINT #2, sumreal#(coll%);
   PRINT #3, sumimag#(coll%);
   sumreal#(coll%) = 0: sumimag#(coll%) = 0
  ELSE
   PRINT #2, freqx; : PRINT #3, freqx;
  END IF
  NEXT coll%: PRINT #2, : PRINT #3,

'open input file and input file header
OPEN "i", 1, filename$: CALL headerinput(text$(), j, 1)
NEXT freqx: CLOSE 1

 CLOSE 2, 3
PRINT
SHELL "echo END FOURIER file " + filename$ + " has been fouriertransformed"

IF INSTR(COMMAND$, "*") <> 0 THEN GOTO 999
END
333 FOR i = 1 TO 23: READ A$: PRINT A$: NEXT i: END

SUB analizecommand (file1$, colx%, xmax, deltax)
STATIC washere%
'*****************************************************************
' this sub analizes the command$ and detects the filename
'*****************************************************************
A$ = LCASE$(LTRIM$(RTRIM$(COMMAND$)))
file1$ = RTRIM$(LEFT$(A$, INSTR(A$, " "))): A$ = LTRIM$(RIGHT$(A$, LEN(A$) - INSTR(A$, " ")))
colx% = ASC(LEFT$(A$, 1)) MOD 48: A$ = LTRIM$(RIGHT$(A$, LEN(A$) - 1))'get xcolumn

n = INSTR(A$, "xmax=")
IF n > 0 THEN xmax = VAL(MID$(A$, n + 5)) ELSE PRINT "error program fourier - parameter xmax not found": END
n = INSTR(A$, "deltax=")
IF n > 0 THEN deltax = VAL(MID$(A$, n + 7)) ELSE PRINT "error program fourier - parameter deltax not found": END


 
  IF INSTR(COMMAND$, "*") <> 0 THEN
   IF washere% = 0 THEN
      washere% = 1
      'get filenames to be operated on as file1$
      SHELL "dir " + file1$ + " /b > FOURIER.dir"
      OPEN "i", 9, "FOURIER.dir"
   END IF
   IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del FOURIER.dir": END
   INPUT #9, file1$
   IF LCASE$(file1$) = "fact.dir" THEN
    IF EOF(9) <> 0 THEN CLOSE 9: SHELL "del fact.dir": END
    INPUT #9, file1$
   END IF
   IF LTRIM$(file1$) = "" THEN CLOSE 9: SHELL "del FOURIER.dir": END
  END IF

END SUB

SUB gobackonestrg (n, ALA$)
'this sub sets the filepointer of file #n one string back and
' puts the corresponding string into ala$
filepointer = SEEK(n)
FOR x = 2 TO filepointer STEP 30
SEEK #n, filepointer - x
fpn = -1: INPUT #n, A$
122 WHILE SEEK(n) < filepointer
       fpn = SEEK(n):  INPUT #n, A$
    WEND
IF fpn > 0 GOTO 233
NEXT x: x = filepointer - 1: SEEK #n, 1: GOTO 122

233 SEEK #n, fpn: ALA$ = A$

END SUB

SUB headerinput (text$(), j, n)
'**********************************************************************
' input file header of measurement file opened with #n
' the file header inputted has then j lines and is stored in text$(1-j)
'**********************************************************************

1 INPUT #n, A$
   i = INSTR(A$, "{"): IF i > 0 GOTO 2 ELSE GOTO 1     'look for "{"
2 text$(1) = RIGHT$(A$, LEN(A$) - i)
  j = 1: i = INSTR(text$(j), "}"): IF i > 0 GOTO 3  'look for "}" in first line
                                           
   FOR j = 2 TO 300
   INPUT #n, text$(j)
   i = INSTR(text$(j), "}"): IF i > 0 GOTO 3      'look for "}"
   NEXT j: SHELL "echo text in data file too long": END
3 text$(j) = LEFT$(text$(j), i - 1)

END SUB

SUB inputline (n, d1#(), col1%)
'input data point line on #n as string and split into numbers
' determine col1% and save data columns in d1#(1...col1%)
A$ = INKEY$: IF A$ <> "" THEN IF ASC(A$) = 27 THEN END

INPUT #n, ALA$
col1% = 0
 WHILE LEN(ALA$) > 0
    col1% = col1% + 1
    ALA$ = LTRIM$(ALA$) + " "
    d1#(col1%) = VAL(LEFT$(ALA$, INSTR(ALA$, " ")))
    ALA$ = LTRIM$(RIGHT$(ALA$, LEN(ALA$) - INSTR(ALA$, " ")))
 WEND

END SUB

SUB transform (freqx, x, f, sumreal#, sumimag#, xmax, deltax)

'sum up transform sum(freqx)=summe[x] f(x)*exp(i*freqx*2*pi*x)
IF ABS(ABS(x) - xmax) > deltax / 2 THEN
 sumreal# = sumreal# + f * COS(2 * 3.141573 * freqx * x) / SQR(2 * xmax / deltax)
 sumimag# = sumimag# + f * SIN(2 * 3.141573 * freqx * x) / SQR(2 * xmax / deltax)
ELSE
 sumreal# = sumreal# + f * COS(2 * 3.141573 * freqx * x) * .5 / SQR(2 * xmax / deltax)
 sumimag# = sumimag# + f * SIN(2 * 3.141573 * freqx * x) * .5 / SQR(2 * xmax / deltax)
END IF

'REM nql 2dimensional
'nql(l%) = 0
'FOR m = 0 TO nanz%
'IF na%(m) = 0 AND nb%(m) = 0 THEN nql(l%) = nql(l%) + n(m)
'IF na%(m) <> 0 AND nb%(m) = 0 THEN nql(l%) = nql(l%) + 2 * n(m) * COS(na%(m) * 2 * 3.141573 * ql(l%, 1))
'IF na%(m) = 0 AND nb%(m) <> 0 THEN nql(l%) = nql(l%) + 2 * n(m) * COS(nb%(m) * 2 * 3.141573 * ql(l%, 2))
'IF na%(m) <> 0 AND nb%(m) <> 0 THEN
'nql(l%) = nql(l%) + 2 * n(m) * COS(2 * 3.141573 * (na%(m) * ql(l%, 1) + nb%(m) * ql(l%, 2)))
'nql(l%) = nql(l%) + 2 * n(m) * COS(2 * 3.141573 * (na%(m) * ql(l%, 1) - nb%(m) * ql(l%, 2)))
'END IF
'NEXT m

END SUB

