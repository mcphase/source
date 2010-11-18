OPEN "i", 1, "l:l990825.bib"
OPEN "i", 2, "l:l990823.bib"

WHILE EOF(1) = 0:
INPUT #1, a$
INPUT #2, b$
IF RTRIM$(LTRIM$(a$)) <> RTRIM$(LTRIM$(b$)) THEN PRINT a$: PRINT b$
WEND
CLOSE



