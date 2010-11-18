
SUB getpath (path$)
'get current path
SHELL "dir > getpath.dum": OPEN "i", 1, "getpath.dum"
WHILE INSTR(a$, ":\") = 0: INPUT #1, a$: WEND
path$ = LCASE$(RIGHT$(a$, LEN(a$) - INSTR(a$, ":\") + 2)) + "\": CLOSE 1: SHELL "del getpath.dum"

END SUB

SUB getsetvar (a$)
REM**********einlesen einer dos umgebungsvariablen a$ ****************
'oon return a$ contains value
' on error program is ended

SHELL "set > getsetv.dum": OPEN "i", 2, "getsetv.dum"
WHILE EOF(2) = 0
INPUT #2, a$: a$ = LCASE$(a$)
 IF LEFT$(a$, 4) = a$ THEN a$ = RIGHT$(a$, LEN(a$) - INSTR(a$, "=")): check% = 1
WEND
IF check% <> 1 THEN PRINT "ERROR in sub getsetvar - dos variable "; a$; " not found": END
CLOSE 2: SHELL "del getsetv.dum"
'**************************************************************************

END SUB

