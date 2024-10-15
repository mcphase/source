/*-----------------------------------------------------------------------------
 
                                 M  O  D  U  L
 
                                    TYPES.C
 
-------------------------------------------------------------------------------
 
Aufgabe               :  Definition der benutzten Datentypen
 
-------------------------------------------------------------------------------
 
Definierte Datentypen :  INT      := Datentyp der ganzen Zahlen
                         CHAR     := Datentyp der Bytes
                         LONG     := Datentyp der grossen ganzen Zahlen
                         FLOAT    := Datentyp der reellen Zahlen
                         REAL     := FLOAT
                         DOUBLE   := Datentyp der reellen Zahlen mit
                                     doppelter Genauigkeit
 
                         KOMPLEX  := Datentyp der komplexen Zahlen
 
-----------------------------------------------------------------------------*/
 
 
 
 
/*-----------------------------------------------------------------------------
benutzte Datentypen definieren
-----------------------------------------------------------------------------*/
 
typedef int    INT;
typedef char   CHAR;
typedef long   LONG;
typedef float  FLOAT;
typedef double DOUBLE;
typedef float  REAL;
 
#ifdef _CRAY
    #define clearscreen  if(1)  /* clearscreen durch dummy */
   typedef void   SPEICHER;
   typedef size_t TYPUS;
   typedef void   FREE;
#else
    #define clearscreen  if(1){}  /* clearscreen durch dummy */
   typedef void   SPEICHER;
   typedef INT    TYPUS;
   typedef void   FREE;
#endif                                /* ersetzen                */
 
/*-----------------------------------------------------------------------------
extern definierte Funktionen
-----------------------------------------------------------------------------*/
extern CHAR *c_alloc(INT,INT);
/* extern FREE free(); */
#define free_(cp) free((SPEICHER*)(cp))
 
/*****************************************************************************/
/*   Version          ********************************************************/
/*****************************************************************************/
#define VERSION 5.6  /********************************************************/
/*****************************************************************************/
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _mnbrak{
    DOUBLE ax;
    DOUBLE bx;
    DOUBLE cx;
    DOUBLE fa;
    DOUBLE fb;
    DOUBLE fc;
}MNBRAK;
 
 
#define MNBRAK_ALLOC(nr)  ((MNBRAK*)c_alloc( (nr),sizeof(MNBRAK)))
#define AX(s) ((s)->ax)
#define BX(s) ((s)->bx)
#define CX(s) ((s)->cx)
#define FA(s) ((s)->fa)
#define FB(s) ((s)->fb)
#define FC(s) ((s)->fc)
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _read{
    CHAR   is_feld;
    CHAR   modus;
    CHAR   *einheit_in;
    CHAR   *einheit_out;
    DOUBLE zwei_j;
    DOUBLE temperatur;
    DOUBLE b1;
    DOUBLE b2;
    DOUBLE b3;
    INT    ionennr;
}READ;
 
 
#define READ_ALLOC(nr)  ((READ*)c_alloc( (nr),sizeof(READ)))
 
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _exponent{
    INT    B2i;
    INT    B4i;
    INT    B6i;
    INT    Bmol;
    INT    chi2;
    INT    x;
    INT    W;
    DOUBLE facB2i;
    DOUBLE facB4i;
    DOUBLE facB6i;
    DOUBLE facBmol;
    DOUBLE facchi2;
    DOUBLE facx;
    DOUBLE facW;
}EXPONENT;
 
#define EXP_ALLOC(nr)  ((EXPONENT*)c_alloc( (nr),sizeof(EXPONENT)))
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _quadrant{
    DOUBLE winkel;
    DOUBLE radius;
}QUADRANT;
 
#define WINKEL(q)    ( (q)->winkel  )
#define RADIUS(q)    ( (q)->radius  )
 
#define Q_ALLOC(nr)  (  (QUADRANT*)c_alloc( (nr),sizeof(QUADRANT) )  )
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _bruch{
    LONG zaehler;  /* Zaehler  */
    LONG nenner;   /* Nenner   */
}BRUCH;
 
#define ZAEHLER(b)  ( (b)->zaehler )
#define NENNER( b)  ( (b)->nenner  )
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _komplex{
    DOUBLE rt;     /* Realteil */
    DOUBLE it;     /* Imaginaerteil */
}KOMPLEX;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _tabelle{
    CHAR *t20n;
    CHAR *t20r;
    CHAR *t20c;
    CHAR *t21n;
    CHAR *t21r;
    CHAR *t21c;
    CHAR *t22n;
    CHAR *t22r;
    CHAR *t22c;
 
    CHAR *t40n;
    CHAR *t40r;
    CHAR *t40c;
    CHAR *t41n;
    CHAR *t41r;
    CHAR *t41c;
    CHAR *t42n;
    CHAR *t42r;
    CHAR *t42c;
    CHAR *t43n;
    CHAR *t43r;
    CHAR *t43c;
    CHAR *t44n;
    CHAR *t44r;
    CHAR *t44c;
 
    CHAR *t60n;
    CHAR *t60r;
    CHAR *t60c;
    CHAR *t61n;
    CHAR *t61r;
    CHAR *t61c;
    CHAR *t62n;
    CHAR *t62r;
    CHAR *t62c;
    CHAR *t63n;
    CHAR *t63r;
    CHAR *t63c;
    CHAR *t64n;
    CHAR *t64r;
    CHAR *t64c;
    CHAR *t65n;
    CHAR *t65r;
    CHAR *t65c;
    CHAR *t66n;
    CHAR *t66r;
    CHAR *t66c;
 
    CHAR *t20;
    CHAR *t21;
    CHAR *t22;
 
    CHAR *t40;
    CHAR *t41;
    CHAR *t42;
    CHAR *t43;
    CHAR *t44;
 
    CHAR *t60;
    CHAR *t61;
    CHAR *t62;
    CHAR *t63;
    CHAR *t64;
    CHAR *t65;
    CHAR *t66;
 
    CHAR *t11;
    CHAR *t12;
    CHAR *t13;
    CHAR *t14;
    CHAR *t15;
 
    CHAR *str;
    CHAR *tss;
}TABELLE;
 
#define TABELLE_ALLOC(nr)       (  (TABELLE* )c_alloc( (nr),sizeof(TABELLE) )  )
 
#define T20n  ( tabelle->t20n )
#define T20r  ( tabelle->t20r )
#define T20c  ( tabelle->t20c )
#define T21n  ( tabelle->t21n )
#define T21r  ( tabelle->t21r )
#define T21c  ( tabelle->t21c )
#define T22n  ( tabelle->t22n )
#define T22r  ( tabelle->t22r )
#define T22c  ( tabelle->t22c )
 
#define T40n  ( tabelle->t40n )
#define T40r  ( tabelle->t40r )
#define T40c  ( tabelle->t40c )
#define T41n  ( tabelle->t41n )
#define T41r  ( tabelle->t41r )
#define T41c  ( tabelle->t41c )
#define T42n  ( tabelle->t42n )
#define T42r  ( tabelle->t42r )
#define T42c  ( tabelle->t42c )
#define T43n  ( tabelle->t43n )
#define T43r  ( tabelle->t43r )
#define T43c  ( tabelle->t43c )
#define T44n  ( tabelle->t44n )
#define T44r  ( tabelle->t44r )
#define T44c  ( tabelle->t44c )
 
#define T60n  ( tabelle->t60n )
#define T60r  ( tabelle->t60r )
#define T60c  ( tabelle->t60c )
#define T61n  ( tabelle->t61n )
#define T61r  ( tabelle->t61r )
#define T61c  ( tabelle->t61c )
#define T62n  ( tabelle->t62n )
#define T62r  ( tabelle->t62r )
#define T62c  ( tabelle->t62c )
#define T63n  ( tabelle->t63n )
#define T63r  ( tabelle->t63r )
#define T63c  ( tabelle->t63c )
#define T64n  ( tabelle->t64n )
#define T64r  ( tabelle->t64r )
#define T64c  ( tabelle->t64c )
#define T65n  ( tabelle->t65n )
#define T65r  ( tabelle->t65r )
#define T65c  ( tabelle->t65c )
#define T66n  ( tabelle->t66n )
#define T66r  ( tabelle->t66r )
#define T66c  ( tabelle->t66c )
 
#define T20   ( tabelle->t20  )
#define T21   ( tabelle->t21  )
#define T22   ( tabelle->t22  )
 
#define T40   ( tabelle->t40  )
#define T41   ( tabelle->t41  )
#define T42   ( tabelle->t42  )
#define T43   ( tabelle->t43  )
#define T44   ( tabelle->t44  )
 
#define T60   ( tabelle->t60  )
#define T61   ( tabelle->t61  )
#define T62   ( tabelle->t62  )
#define T63   ( tabelle->t63  )
#define T64   ( tabelle->t64  )
#define T65   ( tabelle->t65  )
#define T66   ( tabelle->t66  )
 
#define T11   ( tabelle->t11  )
#define T12   ( tabelle->t12  )
#define T13   ( tabelle->t13  )
#define T14   ( tabelle->t14  )
#define T15   ( tabelle->t15  )
#define STR   ( tabelle->str  )
#define TSS   ( tabelle->tss  )
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _vektor{
    INT     _dimension;   /* Anzahl der Zeilen des komplexen Spaltenvektors */
    KOMPLEX *_zeile;      /* Komponenten       des komplexen Spaltenvektors */
}VEKTOR;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _matrix{
    INT    _dimension; /* Anzahl  der Spaltenvektoren der hermiteschen Matrix */
    VEKTOR *_spalte;   /* Spalten der hermiteschen Matrix */
}MATRIX;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _comhes{
    INT overwrite;/* entscheidet ,ob urspruengliche Matrix ueberschrieben wird*/
    MATRIX *matrix; /* urspruengliche Matrix */
    MATRIX *hesse;  /* Hessenbergmatrix ,merken wegen der Ruecktransf. */
    INT *swap;      /* Diese Matrix enthaelt die Zeilen- und Spalten- */
    INT swap_dim;   /* vertauschungen, welche bei der Reduktion der   */
}COMHES;            /* urspruenglichen komplexen    Matrix auf die    */
                    /* Hessenbergmatrix gemacht wurden                */
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _comrl2{
    COMHES *comhes;
    INT    *energie_nummer;  /* hier stehen die Nummern der Energieniveaus */
                             /* mit E( energnr.(1) )<= ... E( energnr.(n) )*/
    MATRIX *entartung;       /* siehe Fkt entartung() in DIAHERMX.c */
    MATRIX *eigenvektoren;   /* Eigenvektoren als Spalten gespeichert*/
    VEKTOR *eigenwerte;      /* Spaltenvektor mit den Eigenwerten    */
    INT    *gi;              /* entartungen der versch. Niveaus */
    INT    iterationsteps;
    CHAR   *iteration;      /* ist FAILED  oder SUCCESSFUL  */
    DOUBLE eps_machine;      /* MACHEPSFACT*Rechnergenauigkeit */
    DOUBLE eps_setup;        /* gesetzte Genauigkeit fuer Energieentartung*/
    DOUBLE shift;            /* kleinster Ew zu 0 geshiftet*/
}COMLR2,EWPROBLEM;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
typedef struct _stevens{
    INT    dimj;             /* Gesamtdrehimpuls J ,dimj := 2J+1*/
 
    MATRIX *p0p0;            /* Matrix (JM'� P0+0(J) �JM) */
 
    MATRIX *p2p2;            /* Matrix (JM'� P2+2(J) �JM) */
    MATRIX *p2p1;            /* Matrix (JM'� P2+1(J) �JM) */
    MATRIX *p2p0;            /* Matrix (JM'� P2+0(J) �JM) */
    MATRIX *p2m1;            /* Matrix (JM'� P2-1(J) �JM) */
    MATRIX *p2m2;            /* Matrix (JM'� P2-2(J) �JM) */
 
    MATRIX *p3p3;            /* Matrix (JM'| P3+3(J) |JM) */
    MATRIX *p3p2;            /* Matrix (JM'| P3+2(J) |JM) */
    MATRIX *p3p1;            /* Matrix (JM'| P3+1(J) |JM) */
    MATRIX *p3p0;            /* Matrix (JM'| P3+0(J) |JM) */
    MATRIX *p3m1;            /* Matrix (JM'| P3-1(J) |JM) */
    MATRIX *p3m2;            /* Matrix (JM'| P3-2(J) |JM) */
    MATRIX *p3m3;            /* Matrix (JM'| P3-3(J) |JM) */

 
    MATRIX *p4p4;            /* Matrix (JM'� P4+4(J) �JM) */
    MATRIX *p4p3;            /* Matrix (JM'� P4+3(J) �JM) */
    MATRIX *p4p2;            /* Matrix (JM'� P4+2(J) �JM) */
    MATRIX *p4p1;            /* Matrix (JM'� P4+1(J) �JM) */
    MATRIX *p4p0;            /* Matrix (JM'� P4+0(J) �JM) */
    MATRIX *p4m1;            /* Matrix (JM'� P4-1(J) �JM) */
    MATRIX *p4m2;            /* Matrix (JM'� P4-2(J) �JM) */
    MATRIX *p4m3;            /* Matrix (JM'� P4-3(J) �JM) */
    MATRIX *p4m4;            /* Matrix (JM'� P4-4(J) �JM) */
 
    MATRIX *p5p5;            /* Matrix (JM'| P5+5(J) |JM) */
    MATRIX *p5p4;            /* Matrix (JM'| P5+4(J) |JM) */
    MATRIX *p5p3;            /* Matrix (JM'| P5+3(J) |JM) */
    MATRIX *p5p2;            /* Matrix (JM'| P5+2(J) |JM) */
    MATRIX *p5p1;            /* Matrix (JM'| P5+1(J) |JM) */
    MATRIX *p5p0;            /* Matrix (JM'| P5+0(J) |JM) */
    MATRIX *p5m1;            /* Matrix (JM'| P5-1(J) |JM) */
    MATRIX *p5m2;            /* Matrix (JM'| P5-2(J) |JM) */
    MATRIX *p5m3;            /* Matrix (JM'| P5-3(J) |JM) */
    MATRIX *p5m4;            /* Matrix (JM'| P5-4(J) |JM) */
    MATRIX *p5m5;            /* Matrix (JM'| P5-5(J) |JM) */

 
    MATRIX *p6p6;            /* Matrix (JM'� P6+6(J) �JM) */
    MATRIX *p6p5;            /* Matrix (JM'� P6+5(J) �JM) */
    MATRIX *p6p4;            /* Matrix (JM'� P6+4(J) �JM) */
    MATRIX *p6p3;            /* Matrix (JM'� P6+3(J) �JM) */
    MATRIX *p6p2;            /* Matrix (JM'� P6+2(J) �JM) */
    MATRIX *p6p1;            /* Matrix (JM'� P6+1(J) �JM) */
    MATRIX *p6p0;            /* Matrix (JM'� P6+0(J) �JM) */
    MATRIX *p6m1;            /* Matrix (JM'� P6-1(J) �JM) */
    MATRIX *p6m2;            /* Matrix (JM'� P6-2(J) �JM) */
    MATRIX *p6m3;            /* Matrix (JM'� P6-3(J) �JM) */
    MATRIX *p6m4;            /* Matrix (JM'� P6-4(J) �JM) */
    MATRIX *p6m5;            /* Matrix (JM'� P6-5(J) �JM) */
    MATRIX *p6m6;            /* Matrix (JM'� P6-6(J) �JM) */
 
    MATRIX *o2p0;            /* Matrix (JM'� O20(J)  �JM) */
    MATRIX *o2p2;            /* Matrix (JM'� O22(J)  �JM) */
    MATRIX *o4p2;            /* Matrix (JM'� O42(J)  �JM) */
    MATRIX *o6p2;            /* Matrix (JM'� O62(J)  �JM) */
    MATRIX *o6p6;            /* Matrix (JM'� O66(J)  �JM) */
 
    MATRIX *o4p05;           /* Matrix (JM'� O40(J)+ 5O44(J)�JM) */
    MATRIX *o4m05;           /* Matrix (JM'� O40(J)- 5O44(J)�JM) */
 
    MATRIX *o6p21;           /* Matrix (JM'� O60(J)+21O64(J)�JM) */
    MATRIX *o6m21;           /* Matrix (JM'� O60(J)-21O64(J)�JM) */
 
    DOUBLE n_o2p0;           /*  || O20(J) || */
    DOUBLE n_o2p2;
    DOUBLE n_o4p2;
    DOUBLE n_o6p2;
    DOUBLE n_o6p6;
 
    DOUBLE n_o4p05;          /* �� O40(J) + 05*O44(J) �� */
    DOUBLE n_o4m05;          /* �� O40(J) - 05*O44(J) �� */
 
    DOUBLE n_o6p21;
    DOUBLE n_o6m21;
 
}STEVENS;
 
#define P0P0(s)     ( (s)->p0p0 )  /* Matrixelemente (JM'|Pkq(j)|JM) */
 
#define P2P2(s)     ( (s)->p2p2 )
#define P2P1(s)     ( (s)->p2p1 )
#define P2P0(s)     ( (s)->p2p0 )
#define P2M1(s)     ( (s)->p2m1 )
#define P2M2(s)     ( (s)->p2m2 )
 
#define P3P3(s)     ( (s)->p3p3 )
#define P3P2(s)     ( (s)->p3p2 )
#define P3P1(s)     ( (s)->p3p1 )
#define P3P0(s)     ( (s)->p3p0 )
#define P3M1(s)     ( (s)->p3m1 )
#define P3M2(s)     ( (s)->p3m2 )
#define P3M3(s)     ( (s)->p3m3 )


#define P4P4(s)     ( (s)->p4p4 )
#define P4P3(s)     ( (s)->p4p3 )
#define P4P2(s)     ( (s)->p4p2 )
#define P4P1(s)     ( (s)->p4p1 )
#define P4P0(s)     ( (s)->p4p0 )
#define P4M1(s)     ( (s)->p4m1 )
#define P4M2(s)     ( (s)->p4m2 )
#define P4M3(s)     ( (s)->p4m3 )
#define P4M4(s)     ( (s)->p4m4 )


#define P5P5(s)     ( (s)->p5p5 )
#define P5P4(s)     ( (s)->p5p4 )
#define P5P3(s)     ( (s)->p5p3 )
#define P5P2(s)     ( (s)->p5p2 )
#define P5P1(s)     ( (s)->p5p1 )
#define P5P0(s)     ( (s)->p5p0 )
#define P5M1(s)     ( (s)->p5m1 )
#define P5M2(s)     ( (s)->p5m2 )
#define P5M3(s)     ( (s)->p5m3 )
#define P5M4(s)     ( (s)->p5m4 )
#define P5M5(s)     ( (s)->p5m5 )

 
#define P6P6(s)     ( (s)->p6p6 )
#define P6P5(s)     ( (s)->p6p5 )
#define P6P4(s)     ( (s)->p6p4 )
#define P6P3(s)     ( (s)->p6p3 )
#define P6P2(s)     ( (s)->p6p2 )
#define P6P1(s)     ( (s)->p6p1 )
#define P6P0(s)     ( (s)->p6p0 )
#define P6M1(s)     ( (s)->p6m1 )
#define P6M2(s)     ( (s)->p6m2 )
#define P6M3(s)     ( (s)->p6m3 )
#define P6M4(s)     ( (s)->p6m4 )
#define P6M5(s)     ( (s)->p6m5 )
#define P6M6(s)     ( (s)->p6m6 )
 
#define O2P0(s)     ( (s)->o2p0 )
#define O2P2(s)     ( (s)->o2p2 )
#define O4P2(s)     ( (s)->o4p2 )
#define O6P2(s)     ( (s)->o6p2 )
#define O6P6(s)     ( (s)->o6p6 )
 
#define O4P05(s)     ( (s)->o4p05 )
#define O4M05(s)     ( (s)->o4m05 )
 
#define O6P21(s)     ( (s)->o6p21 )
#define O6M21(s)     ( (s)->o6m21 )
 
#define N_O2P0(s)     ( (s)->n_o2p0 )
#define N_O2P2(s)     ( (s)->n_o2p2 )
#define N_O4P2(s)     ( (s)->n_o4p2 )
#define N_O6P2(s)     ( (s)->n_o6p2 )
#define N_O6P6(s)     ( (s)->n_o6p6 )
 
#define N_O4P05(s)     ( (s)->n_o4p05 )
#define N_O4M05(s)     ( (s)->n_o4m05 )
 
#define N_O6P21(s)     ( (s)->n_o6p21 )
#define N_O6M21(s)     ( (s)->n_o6m21 )
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _info{
    CHAR *info;             /* Syntax  des Infos */
}INFO;
 
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _befehl{
    CHAR *befehlkommentar;
    CHAR *befehl;             /* Syntax  des Befehls */
}BEFEHL;
 
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _ludcmp{
    INT    invers;
    MATRIX *matrix;
    VEKTOR *index;
    DOUBLE d;
}LUDCMP;
 
#define LUDCMP_ALLOC(nr)    ( (LUDCMP*)c_alloc((nr),sizeof(LUDCMP)) )
 
/*-----------------------------------------------------------------------------
                        Implementierte Ionen
    weitere Implementationen durch Anfuegen einer Zeile in IONENIMP[]
-----------------------------------------------------------------------------*/
typedef struct ionen{              /* Ionen implementieren */
    CHAR    *ionname;              /* Name  des implementierten Ions     */
    INT     elektronen_in_vier_f;  /* Anzahl der Elektronen in 4f Schale */
    DOUBLE  gj;                    /* Landefaktor */
    INT     dimj;                  /* dimj = 2*J+1 */
    INT     f4;                    /* F(4) fuer x,W Parameter */
    INT     f6;                    /* F(6) fuer x,W Parameter */
    DOUBLE  r2;                    /* <r2> */
    DOUBLE  r4;                    /* <r4> */
    DOUBLE  r6;                    /* <r6> */
}IONEN;
 
#define E4f(implement_ionnr)   IONENIMP[implement_ionnr].elektronen_in_vier_f
#define F4( implement_ionnr)   IONENIMP[implement_ionnr].f4
#define F6( implement_ionnr)   IONENIMP[implement_ionnr].f6
#define r2( implement_ionnr)   IONENIMP[implement_ionnr].r2
#define r4( implement_ionnr)   IONENIMP[implement_ionnr].r4
#define r6( implement_ionnr)   IONENIMP[implement_ionnr].r6
#define e4f(umgebung)          E4f( IONENNR(umgebung) ) /* Elektronen in 4f  */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 typedef struct tafel{
    INT     anz_inaequivalenter_klassen;
    INT     *anz_aequivalenter_klassen;
    CHAR    **klassen;
    CHAR    **irreduzible_darstellungen;
    MATRIX  *charaktertafel;
}CHARAKTERTAFEL;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _einheit{
    CHAR *einheit;      /* implementierte Einheiten als string */
    CHAR c;             /* Erkennungszeichen einer implementierten Einheit */
    DOUBLE myB;         /* Bohrsches Magneton in der implement. Einheit */
    DOUBLE fke;         /* Umrechnungsfaktor von kelvin auf einheit */
    DOUBLE fek;         /* Umrechnungsfaktor von einheit auf Kelvin  */
}EINHEIT;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _setup{
    INT nummerierung   ;  /* Bestimmt Nummerierung der Datensaetze*/
    INT nachkommastelle;  /* Bestimmt die Nachkommastelle ab */
    INT anz_datenpunkte;  /* Anzahl der berechneten Datenpkte bei plots*/
    DOUBLE macheps;       /* gesetzte Genauigk. fuer Vergleich zweier  */
}SETUP;                   /* Energiewerte                              */
 
#define NUMMERIERUNG(i)    ( (i)->nummerierung    )
#define KOMMASTELLE(i)     ( (i)->nachkommastelle )
#define ANZ_DATENPUNKTE(i) ( (i)->anz_datenpunkte )
#define MACHEPS(i)         ( (i)->macheps )
#define SETUP_ALLOC(nr)    ( (SETUP*)c_alloc((nr),sizeof(SETUP)) )
 
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _umgebung{
    DOUBLE temperatur;             /* Temperatur der Probe    */
    INT    ionennr;                /* Nummer des implementierten Ions */
    INT    anz_nn;                 /* Anzahl der naechsten Nachbarn */
    CHAR   modus;                  /* ist 'r' oder 's' oder 'p'     */
    DOUBLE *q;                     /* Ladung der Umgebungsionen(Umi)*/
    DOUBLE *x1;                    /* x = (x1,x2,x2)  : Ort der Umi.*/
    DOUBLE *x2;
    DOUBLE *x3;
 
    DOUBLE b1s,b2s,b3s; /* Anisotropieparameter Sx^2, Sy^2, Sz^2 */
    DOUBLE b1ss,b2ss,b3ss; /* Anisotropieparameter Sx^4, Sy^4, Sz^4 */
 
    DOUBLE b1;     /* b = ( b1,b2,b3 ) angelegtes Magnetfeld */
    DOUBLE b2;     /*                                        */
    DOUBLE b3;     /*  H = H  + H    mit  H  = -g my J*b     */
}UMGEBUNG;         /*       0    mag       mag   J  B        */
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _nebenbedingung{
    INT    anz_par;            /* anzahl der Parameter im Hamiltonian*/
    INT    anz_var;            /* davon sind anz_var variabel        */
    INT    anzahl;             /* anzahl der Nebenbedingungen */
    VEKTOR *fix;               /* 9-dim, merken welche Parameter fest*/
    VEKTOR *eigenwerte;        /* gemessene Energieeigenwerte */
    VEKTOR *d_eigenwerte;      /* und deren angenommene Fehler */
    MATRIX *intensitaeten;     /* gemessene Uebergangsintensitaeten */
    MATRIX *d_intensitaeten;   /* und deren angenommene Fehler */
}NEBENBEDINGUNG;
 
#define NEBEN_ALLOC(nr) ((NEBENBEDINGUNG*)c_alloc( (nr),sizeof(NEBENBEDINGUNG)))
#define NEBENBEDINGUNG(i)  ((i)->nebenbedingung)
#define   ANZAHL(n)          ((n)->  anzahl    )
#define   ANZ_PAR(n)         ((n)->  anz_par   )
#define   ANZ_VAR(n)         ((n)->  anz_var   )
#define   FIX(    n)         ((n)->  fix       )
#define   EIGENWERTE(n)      ((n)->  eigenwerte)
#define D_EIGENWERTE(n)      ((n)->d_eigenwerte)
#define   INTENSITAETEN(n)   ((n)->  intensitaeten)
#define D_INTENSITAETEN(n)   ((n)->d_intensitaeten)
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _iteration{
    INT     einheitnr_in;       /* Nr der Implementierten Einheit */
    INT     einheitnr_out;      /* Nr der Implementierten Einheit */
    INT     symmetrienr;     /* Symmetrienr des Kristallfeldes */
    DOUBLE  temperatur;      /* Temperatur der Probe    */
 
    CHAR    *ionname;   /* Name des Ions dessen Umgebung betrachtet wird */
    INT     ionennummer;  /* number of the ion */
    INT     anz_nn;     /* Anzahl der naechsten Nachbarn */
    INT     dimj;       /* Gesamtdrehimpuls J : dimj = 2J+1 */
    DOUBLE  gj;         /* Landefaktor von ionname */
 
    CHAR   modus;       /* ist 'r' oder 's' oder 'p'     */
    CHAR   if_feld;     /* ist 'n' oder 'a' */
    DOUBLE efversion;   /* Versionsnummer des Eingabefiles  */
 
    DOUBLE  *q;         /* Ladungen der Umgebungsionen */
    DOUBLE  *r2_3;      /*    k           k+1  */
    DOUBLE  *r4_5;      /*  <r >  /  | R |     */
    DOUBLE  *r6_7;      /*                     */
 
    DOUBLE b1s,b2s,b3s; /* Anisotropieparameter fuer Sx^2, Sy^2, Sz^2 operatoren */
    DOUBLE b1ss,b2ss,b3ss; /* Anisotropieparameter fuer Sx^4, Sy^4, Sz^4 operatoren */

    DOUBLE b1;     /* b = ( b1,b2,b3 ) externes   Magnetfeld */
    DOUBLE b2;     /*                                        */
    DOUBLE b3;     /*  H = H  + H    mit  H  =  g my J*b     */
                   /*   1   0    mag1      mag1  J  B        */
 
    DOUBLE b1mol;  /* b = ( b1,b2,b3 ) Molekularfeld         */
    DOUBLE b2mol;  /*                                        */
    DOUBLE b3mol;  /*  H = H  + H    mit  H  =   2(g -1) my J*b     */
                   /*       1    mag2      mag2     J      B        */
 
    DOUBLE bmol;   /*  (b1mol,b2mol,b3mol) in sphaerischen */
    DOUBLE phi;    /*  Koordinaten  geschrieben            */
    DOUBLE theta;  /*                                      */
 
    KOMPLEX **c20;  /*                           */
    KOMPLEX **c21;  /*                           */
    KOMPLEX **c22;  /*   C   der Umgebungsionen  */
    KOMPLEX **c40;  /*    kq                     */
    KOMPLEX **c41;
    KOMPLEX **c42;
    KOMPLEX **c43;
    KOMPLEX **c44;
    KOMPLEX **c60;
    KOMPLEX **c61;
    KOMPLEX **c62;
    KOMPLEX **c63;
    KOMPLEX **c64;
    KOMPLEX **c65;
    KOMPLEX **c66;
 
    MATRIX  *h_magnetfeld;     /* = (JM| g my J H | M'J) */
 
    STEVENS *Pkq;  /* Verall. Operatoren P (J)  */
                   /*                     kq    */
    /*--------------------------------------------------------
 
           --               --                  *    +
    H   =  >  V  * P (J)  + >   V  * P (J)  +  V  * P (J)
     CF    --  k0   k0      --   kq   kq        kq   kq
            k               k>=1
                            q>0
 
                             2S+1
    V   =  beta   *  alpha (     L  )  *  B
     kq        kq         k       J        kq
 
              k
    B   =  < r  >  *  Q
     kq                kq
 
 
                 /   Q(R)       *            3
    Q   =  -|e|  |  ------  *  C ( omega )  d R
     kq          /     k+1      kq
                      R
    ---------------------------------------------------------*/
 
    KOMPLEX *V20;
    KOMPLEX *V21;
    KOMPLEX *V22;
    KOMPLEX *V40;
    KOMPLEX *V41;
    KOMPLEX *V42;
    KOMPLEX *V43;
    KOMPLEX *V44;
    KOMPLEX *V60;
    KOMPLEX *V61;
    KOMPLEX *V62;
    KOMPLEX *V63;
    KOMPLEX *V64;
    KOMPLEX *V65;
    KOMPLEX *V66;
 
    MATRIX  *hamiltonian; /*  =   H              +   H           */
                          /*       kristallfeld       Magnetfeld */
 
    NEBENBEDINGUNG *nebenbedingung;    /* CHI-Quadrat Nebenbedingungen */
    INT            fitroutinennr;      /* Nummer der Fitroutine        */
    INT            fitmax;             /* Abbruch der Iteration bei fitmax  */
    VEKTOR         *v0;                /* v0 merken  */
 
 
    INT            zeile;              /* nummeriert Parametersaetze durch */
    INT            max_zeile;          /* nummeriert Parametersaetze durch */
    INT            is_intensitaet;
    INT            is_matrixelement;
    DOUBLE         inttemp;            /* Intensitaeten zur Temp. inttemp  */
    INT            is_eigenwert;
 
    INT            is_menue;
    INT            menue;
 
    INT            is_susfit;
    INT            is_susa;
    INT            is_susb;
    INT            is_susc;
    INT            is_susp;
 
    INT            anzsa;
    INT            anzsb;
    INT            anzsc;
    INT            anzsp;
 
    INT            nrsusa;
    INT            nrsusb;
    INT            nrsusc;
    INT            nrsusp;
    CHAR           *susname;
    INT            formatpa;
    INT            formatpb;
    INT            formatna;
    INT            formatnb;
 
    DOUBLE         *susat;
    DOUBLE         *susasus;
    DOUBLE         *susbt;
    DOUBLE         *susbsus;
    DOUBLE         *susct;
    DOUBLE         *suscsus;
    DOUBLE         *suspt;
    DOUBLE         *suspsus;
 
    DOUBLE         lambda;
 
    INT            is_magfit;
    INT            is_maga;
    INT            is_magb;
    INT            is_magc;
    INT            is_magp;
    DOUBLE         magatemp;           /* m_a           zur Temp. inttemp  */
    DOUBLE         magbtemp;           /* m_b           zur Temp. inttemp  */
    DOUBLE         magctemp;           /* m_c           zur Temp. inttemp  */
    DOUBLE         magptemp;           /* m_p           zur Temp. inttemp  */
 
    INT            anzma;
    INT            anzmb;
    INT            anzmc;
    INT            anzmp;
 
    INT            nrmaga;
    INT            nrmagb;
    INT            nrmagc;
    INT            nrmagp;
    CHAR           *magname;
 
    DOUBLE         *magaB;
    DOUBLE         *magamag;
    DOUBLE         *magbB;
    DOUBLE         *magbmag;
    DOUBLE         *magcB;
    DOUBLE         *magcmag;
    DOUBLE         *magpB;
    DOUBLE         *magpmag;
 
    INT            is_posfit;
    CHAR           *posname;
    INT            posdatanz;
    INT            nrposa;
    INT            nrpose;
    DOUBLE         *post;
    DOUBLE         *posicin;
    DOUBLE         *posicqe;
    INT            *posanz;
    DOUBLE         **pose;
    DOUBLE         **posi;
 
    INT            r1;
    INT            r2;
    INT            r3;
 
    DOUBLE         wew;
    DOUBLE         wpos;
    DOUBLE         wint;
    DOUBLE         wmat;
    DOUBLE         wsus;
    DOUBLE         wmag;
 
    DOUBLE         vor20;
    DOUBLE         vor21;
    DOUBLE         vor22;
    DOUBLE         vor40;
    DOUBLE         vor41;
    DOUBLE         vor42;
    DOUBLE         vor43;
    DOUBLE         vor44;
    DOUBLE         vor60;
    DOUBLE         vor61;
    DOUBLE         vor62;
    DOUBLE         vor63;
    DOUBLE         vor64;
    DOUBLE         vor65;
    DOUBLE         vor66;
 
    DOUBLE         thetaa;
    INT            lesethetafile;
    CHAR          *namethetafile;
    DOUBLE         *thetax;
    DOUBLE         *thetaf;
    INT            thetaanz;
 
}ITERATION;
#define THETAA(i)   ((i)->thetaa)
#define THETAX(i)   ((i)->thetax)
#define THETAF(i)   ((i)->thetaf)
#define THETAANZ(i) ((i)->thetaanz)
 
#define BMOL(i)   ((i)->bmol)
#define THETA(i)  ((i)->theta)
#define PHI(i)    ((i)->phi  )
 
#define VOR20(i)  ((i)->vor20 )
#define VOR21(i)  ((i)->vor21 )
#define VOR22(i)  ((i)->vor22 )
#define VOR40(i)  ((i)->vor40 )
#define VOR41(i)  ((i)->vor41 )
#define VOR42(i)  ((i)->vor42 )
#define VOR43(i)  ((i)->vor43 )
#define VOR44(i)  ((i)->vor44 )
#define VOR60(i)  ((i)->vor60 )
#define VOR61(i)  ((i)->vor61 )
#define VOR62(i)  ((i)->vor62 )
#define VOR63(i)  ((i)->vor63 )
#define VOR64(i)  ((i)->vor64 )
#define VOR65(i)  ((i)->vor65 )
#define VOR66(i)  ((i)->vor66 )
 
#define INTTEMP(i)  ((i)->inttemp)
#define MAGATEMP(i) ((i)->magatemp)
#define MAGBTEMP(i) ((i)->magbtemp)
#define MAGCTEMP(i) ((i)->magctemp)
#define MAGPTEMP(i) ((i)->magptemp)
 
#define WEW(i)  ((i)->wew)
#define WPOS(i) ((i)->wpos)
#define WINT(i) ((i)->wint)
#define WMAT(i) ((i)->wmat)
#define WSUS(i) ((i)->wsus)
#define WMAG(i) ((i)->wmag)
 
#define IS_FELD(i) ((i)->is_feld)
#define IS_EIGENWERT(i) ((i)->is_eigenwert)
#define R1(i) ((i)->r1 )
#define R2(i) ((i)->r2 )
#define R3(i) ((i)->r3 )
 
#define FITROUTINENNR(i) ((i)->fitroutinennr)
#define FITMAX(i) ((i)->fitmax)
#define V0(       i) ((i)->v0 )
#define ZEILE(    i) ((i)->zeile )
#define MAX_ZEILE(i) ((i)->max_zeile)
#define IS_INTENSITAET(i) ((i)->is_intensitaet)
#define IS_MATRIXELEMENT(i) ((i)->is_matrixelement)
 
#define IS_MENUE(i) ((i)->is_menue)
#define MENUE(i) ((i)->menue)
 
#define IS_SUSFIT(i) ((i)->is_susfit)
#define IS_SUSA(i) ((i)->is_susa)
#define IS_SUSB(i) ((i)->is_susb)
#define IS_SUSC(i) ((i)->is_susc)
#define IS_SUSP(i) ((i)->is_susp)
#define NRSUSA(i) ((i)->nrsusa)
#define NRSUSB(i) ((i)->nrsusb)
#define NRSUSC(i) ((i)->nrsusc)
#define NRSUSP(i) ((i)->nrsusp)
#define SUSNAME(i) ((i)->susname)
#define FORMATPA(i) ((i)->formatpa)
#define FORMATPB(i) ((i)->formatpb)
#define FORMATNA(i) ((i)->formatna)
#define FORMATNB(i) ((i)->formatnb)
 
#define ANZSA(i) ((i)->anzsa)
#define ANZSB(i) ((i)->anzsb)
#define ANZSC(i) ((i)->anzsc)
#define ANZSP(i) ((i)->anzsp)
 
#define SUSAT(i) ((i)->susat)
#define SUSBT(i) ((i)->susbt)
#define SUSCT(i) ((i)->susct)
#define SUSPT(i) ((i)->suspt)
#define SUSASUS(i) ((i)->susasus)
#define SUSBSUS(i) ((i)->susbsus)
#define SUSCSUS(i) ((i)->suscsus)
#define SUSPSUS(i) ((i)->suspsus)
 
 
#define IS_MAGFIT(i) ((i)->is_magfit)
#define IS_MAGA(i) ((i)->is_maga)
#define IS_MAGB(i) ((i)->is_magb)
#define IS_MAGC(i) ((i)->is_magc)
#define IS_MAGP(i) ((i)->is_magp)
 
#define NRMAGA(i) ((i)->nrmaga)
#define NRMAGB(i) ((i)->nrmagb)
#define NRMAGC(i) ((i)->nrmagc)
#define NRMAGP(i) ((i)->nrmagp)
 
#define MAGNAME(i) ((i)->magname)
 
#define ANZMA(i) ((i)->anzma)
#define ANZMB(i) ((i)->anzmb)
#define ANZMC(i) ((i)->anzmc)
#define ANZMP(i) ((i)->anzmp)
 
#define MAGAB(i) ((i)->magaB)
#define MAGBB(i) ((i)->magbB)
#define MAGCB(i) ((i)->magcB)
#define MAGPB(i) ((i)->magpB)
#define MAGAMAG(i) ((i)->magamag)
#define MAGBMAG(i) ((i)->magbmag)
#define MAGCMAG(i) ((i)->magcmag)
#define MAGPMAG(i) ((i)->magpmag)
 
#define IS_POSFIT(i) ((i)->is_posfit)
#define POSNAME(i) ((i)->posname)
#define NRPOSA(i) ((i)->nrposa)
#define NRPOSE(i) ((i)->nrpose)
#define POST(i) ((i)->post)
#define POSICIN(i) ((i)->posicin)
#define POSICQE(i) ((i)->posicqe)
#define POSANZ(i) ((i)->posanz)
#define POSDATANZ(i) ((i)->posdatanz)
#define POSE(i) ((i)->pose)
#define POSI(i) ((i)->posi)
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _amoeba{
    CHAR   *text1;           /* Ist Iteration gelungen ? */
    CHAR   *text2;           /* Iterationsverhalten: stabil? */
    INT    iter_steps;      /* benoetigte zahl der Iterationsschritte */
    MATRIX *matrix;         /* Matrix p */
    VEKTOR *p;              /* Vektor y */
    VEKTOR *xi;             /* Vektor y */
    DOUBLE r;               /* r,phi1,...,phi8 */
    DOUBLE phi1;
    DOUBLE phi2;
    DOUBLE phi3;
    DOUBLE phi4;
    DOUBLE phi5;
    DOUBLE phi6;
    DOUBLE phi7;
    DOUBLE phi8;
 
    DOUBLE chi2;
    INT    anz_wiederholung;
    INT    max_wiederholung;
    ITERATION *iteration;
 
    DOUBLE brent;
    DOUBLE fret;
    DOUBLE xmin;
 
}MINIMUM;
 
#define MINIMUM_ALLOC(nr) ((MINIMUM*)c_alloc( (nr),sizeof(MINIMUM)))
#define TEXT1(s)          ( (s)->text1 )
#define TEXT2(s)          ( (s)->text2 )
#define MAX_WIEDER(s)     ( (s)->max_wiederholung )
#define ITER_STEPS(s)     ( (s)->iter_steps )
#define MATRIX(s)         ( (s)->matrix )
#define VEKTOR(s)         ( (s)->p      )
#define BRENT( s)         ( (s)->brent )
#define FRET( s)          ( (s)->fret  )
#define XMIN(  s)         ( (s)->xmin  )
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _linmin{
    VEKTOR *p;
    VEKTOR *xi;
    DOUBLE fret;
}LINMIN;
 
#define LINMIN_ALLOC(nr) ((LINMIN*)c_alloc( (nr),sizeof(LINMIN)))
#define FRET(s)        ( (s)->fret )
#define P_VEKTOR( s)         ( (s)->p )
#define XI_VEKTOR(s)         ( (s)->xi)
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _sphaere{
    DOUBLE r;        /* Radius         */
    DOUBLE theta;    /* Polarwinkel    */
    DOUBLE phi;      /* Azimutwinkel   */
}SPHAERE;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _symmetrie{
    INT    symnummer; /* Symmetrienummer,welche zum Symmetrienamen gehoe rt */
    CHAR   *symname;  /* Symmetriename */
}SYMMETRIE;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _fit{
    INT    fitnr;     /* Fitroutinennummer */
    CHAR   *fitname;  /* Fitroutinenname   */
}FIT;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _kristallfeld{
    ITERATION *iteration;   /* Umgebungsionen vorbereitet fuer Iteration */
    INT       symmetrienr;  /* Symmetrienr des Kristallfeldes */
    INT       ionennummer;  /* number of the ion */
    INT       is_suszept;   /* Soll Suszeptibilitaet berechnet werden ?*/
    INT       is_magnetm;   /* Magnetisches RE3+ Moment berechnen?*/
    INT       is_kpoly;     /* Magnetisches RE3+ Moment mitteln (kubisch)*/
    INT       is_ortho;     /* (Ortho)rhombischer Kristallfeldfit ?*/
    DOUBLE anfang_temperatur;
    DOUBLE end_temperatur;
    DOUBLE lambda;
    DOUBLE theta;
    DOUBLE anfang_feld;
    DOUBLE end_feld;
    DOUBLE temp;
    CHAR      parameterart; /* eingabeparameterart merken */
    CHAR     *filename;     /* filename des outputs */
    CHAR     *namethetafile;/* filename von theta(T)*/
    CHAR     *infile;       /* filename des intputfiles */
    INT      lesethetafile;
}KRISTALLFELD;              /* nach Tabelle  -i s             */
 
#define INFILE(           s)  ((s)-> infile  )
#define FILENAME(           s)  ((s)-> filename  )
#define EFVERSION(          s)  ((s)-> efversion)
#define IS_SUSZEPT(         s)  ((s)-> is_suszept)
#define ANFANG_TEMPERATUR(  s)  ((s)-> anfang_temperatur)
#define END_TEMPERATUR(     s)  ((s)-> end_temperatur)
#define NAMETHETAFILE(      s)  ((s)-> namethetafile )
#define LESETHETAFILE(      s)  ((s)-> lesethetafile )
#define LAMBDA(             s)  ((s)-> lambda)
#define IS_MAGNETM(         s)  ((s)-> is_magnetm)
#define IS_KPOLY(           s)  ((s)-> is_kpoly  )
#define IS_ORTHO(           s)  ((s)-> is_ortho  )
#define ANFANG_FELD(  s)  ((s)-> anfang_feld)
#define END_FELD(     s)  ((s)-> end_feld)
#define TEMP(         s)  ((s)-> temp)
#define EINGABEPARAMETERART(s)  ((s)-> parameterart)
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* commons fuer c - fortran interlanguages  */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
typedef struct _calfun{ /* struct _CALFUN wird von CALFUN gebraucht */
   DOUBLE   (*funk)();  /* CALFUN wird in MINIMA.C  definiert       */
   SETUP     *setup;
   EWPROBLEM *ewproblem;
   ITERATION *iteration;
   VEKTOR    *vektor;
   VEKTOR    *start;
}_CALFUN;
 
/*-----------------------------------------------------------------------------
benutzte Defines definieren
-----------------------------------------------------------------------------*/
 
#define MXSP(m,sp)     (  (m)->_spalte +(sp) )  /*Spaltenvektor der Matrix m*/
#define ANZ_SP(m)      (  (m)->_dimension    )  /*Spaltenzahl   der Matrix m*/
#define ANZ_ZN(m,sp)   (  MXSP(m,sp) -> _dimension ) /*Zeilenzahl des sp-ten*/
                                                     /*Spaltenvektors von m */
#define ANZ_ZE(m)         ANZ_ZN(m,0) /*Zeilenzahl 0.Spalte von m */
#define MXDIM(m)          ANZ_SP(m)   /*Dimension einer quadratischen Matrix*/
#define MXZE0(m,sp)    (  MXSP(m,sp) -> _zeile  ) /*0.Zeile der sp-ten Spalte*/
#define MX(m,ze,sp)    (  MXZE0(m,sp) + (ze)  )/*ze-te Zeile der sp-ten Spalte*/
#define RT(z)          (  (z) -> rt  )      /* Realteil der komplexen Zahl z*/
#define IT(z)          (  (z) -> it  )      /* Imagteil der komplexen Zahl z*/
#define R(m,ze,sp)        RT(MX(m,ze,sp))/* Realtl des ze-sp-ten Elements  */
#define I(m,ze,sp)        IT(MX(m,ze,sp))/* Imagtl des ze-sp-ten Elements  */
 
#define VRDIM(v)       (  (v)->_dimension     ) /* Spaltenvektorlaenge */
#define VR(v,ze)       (  (v)->_zeile + (ze)  ) /* ze-te Zeile von v   */
#define RV(v,ze)          RT(VR(v,ze))  /* Realtl der ze-ten Zeile von v */
#define IV(v,ze)          IT(VR(v,ze))  /* Imagtl der ze-ten Zeile von v */
 
/* Speicherplatz fuer die verschiedenen Strukturen */
#define KX_P_ALLOC(nr)        (  (KOMPLEX**)c_alloc( (nr),sizeof(KOMPLEX*) )  )
#define BRUCH_ALLOC(nr)       (  (BRUCH*   )c_alloc( (nr),sizeof(BRUCH   ) )  )
#define KX_ALLOC(nr)          (  (KOMPLEX* )c_alloc( (nr),sizeof(KOMPLEX ) )  )
#define VR_ALLOC(nr)          (  (VEKTOR* )c_alloc( (nr),sizeof(VEKTOR ) )  )
#define MX_ALLOC(nr)          (  (MATRIX* )c_alloc( (nr),sizeof(MATRIX ) )  )
#define COMHES_ALLOC(nr)      (  (COMHES* )c_alloc( (nr),sizeof(COMHES ) )  )
#define COMLR2_ALLOC(nr)      (  (COMLR2* )c_alloc( (nr),sizeof(COMLR2 ) )  )
#define INT_ALLOC(nr)         (  (INT*    )c_alloc( (nr+1),sizeof(INT  ) )  )
#define INTP_ALLOC(nr)        (  (INT**   )c_alloc( (nr+1),sizeof(INT* ) )  )
#define EWPROBLEM_ALLOC(nr)   COMLR2_ALLOC(nr)
#define STEVENS_ALLOC(nr)     (  (STEVENS*)c_alloc( (nr),sizeof(STEVENS) )  )
#define VALUE(i,n)            ( *((i) + (n)) )
 
#define STRINGP_ALLOC(nr)     ( (CHAR**  )c_alloc( (nr+1),sizeof(CHAR*) )  )
#define STRING_ALLOC(nr)      ( (CHAR*   )c_alloc( (nr+1),sizeof(CHAR ) )  )
 
#define INFO(nr)              INFOLISTE[nr].info
 
#define BEFEHL(nr)            BEFEHLLISTE[nr].befehl
#define BEFEHLKOM(nr)         BEFEHLLISTE[nr].befehlkommentar
 
#define UMGEBUNG_ALLOC(nr)    (  (UMGEBUNG*)c_alloc( (nr),sizeof(UMGEBUNG) )  )
#define DOUBLE_ALLOC(nr)      (  (DOUBLE* )c_alloc((nr+1),sizeof(DOUBLE  ) )  )
#define REAL_ALLOC(nr)        (  (REAL*   )c_alloc((nr+1),sizeof(REAL    ) )  )
#define DOUBLEP_ALLOC(nr)     (  (DOUBLE**)c_alloc((nr+1),sizeof(DOUBLE* ) )  )
#define IONENNR(i)            (   (i)-> ionennr  )
#define EINHEITNRIN(i)        (   (i)-> einheitnr_in)
#define EINHEITNROUT(i)       (   (i)-> einheitnr_out)
#define ANZ_NN(i)             (   (i)-> anz_nn   )
#define MODUS(i)              (   (i)-> modus    )
#define Q_P(i)                (   (i)-> q        )
#define X1_P(i)               (   (i)-> x1       )
#define X2_P(i)               (   (i)-> x2       )
#define X3_P(i)               (   (i)-> x3       )
#define Q(i,nr)               VALUE(Q_P(i)  ,nr)
#define X1(i,nr)              VALUE(X1_P(i) ,nr)
#define X2(i,nr)              VALUE(X2_P(i) ,nr)
#define X3(i,nr)              VALUE(X3_P(i) ,nr)
#define B1SS(i)                 (  (i) -> b1ss  )
#define B2SS(i)                 (  (i) -> b2ss  )
#define B3SS(i)                 (  (i) -> b3ss  )
#define B1S(i)                 (  (i) -> b1s  )
#define B2S(i)                 (  (i) -> b2s  )
#define B3S(i)                 (  (i) -> b3s  )
#define B1(i)                 (  (i) -> b1  )
#define B2(i)                 (  (i) -> b2  )
#define B3(i)                 (  (i) -> b3  )
#define B1MOL(i)                 (  (i) -> b1mol  )
#define B2MOL(i)                 (  (i) -> b2mol  )
#define B3MOL(i)                 (  (i) -> b3mol  )
 
#define SPHAERE_ALLOC(nr)     ( (SPHAERE* )c_alloc( (nr),sizeof(SPHAERE  ) )  )
#define SPH_R(i)              (  (i) -> r      )
#define SPH_THETA(i)          (  (i) -> theta  )
#define SPH_PHI(i)            (  (i) -> phi    )
 
#define ITERATION_ALLOC(nr)   ( (ITERATION*)c_alloc( (nr),sizeof(ITERATION) )  )
#define IONNAME(i)            (   (i)-> ionname)
#define IONENNUMMER(i)        (   (i)-> ionennummer)
#define DIMJ(i)               (   (i)-> dimj   )
#define TEMPERATUR(i)         (   (i)-> temperatur )
#define GJ(i)                 (   (i)-> gj     )
#define EINHEIT(i)            (   (i)-> energieeinheit )
#define Q_R_P(i)              (   (i)-> q      )
#define Q_R(i,nr)             VALUE(Q_R_P(i),nr)
#define R2_3_P(i)              (   (i)-> r2_3   )
#define R4_5_P(i)              (   (i)-> r4_5   )
#define R6_7_P(i)              (   (i)-> r4_5   )
#define R2_3(i,nr)            VALUE(R2_3_P(i),nr)
#define R4_5(i,nr)            VALUE(R4_5_P(i),nr)
#define R6_7(i,nr)            VALUE(R6_7_P(i),nr)
#define C20_P(i)              (   (i)-> c20    )
#define C21_P(i)              (   (i)-> c21    )
#define C22_P(i)              (   (i)-> c22    )
#define C40_P(i)              (   (i)-> c40    )
#define C41_P(i)              (   (i)-> c41    )
#define C42_P(i)              (   (i)-> c42    )
#define C43_P(i)              (   (i)-> c43    )
#define C44_P(i)              (   (i)-> c44    )
#define C60_P(i)              (   (i)-> c60    )
#define C61_P(i)              (   (i)-> c61    )
#define C62_P(i)              (   (i)-> c62    )
#define C63_P(i)              (   (i)-> c63    )
#define C64_P(i)              (   (i)-> c64    )
#define C65_P(i)              (   (i)-> c65    )
#define C66_P(i)              (   (i)-> c66    )
 
#define C20FS(i,nr)           VALUE(  C20_P(i) , nr  )
#define C21FS(i,nr)           VALUE(  C21_P(i) , nr  )
#define C22FS(i,nr)           VALUE(  C22_P(i) , nr  )
#define C40FS(i,nr)           VALUE(  C40_P(i) , nr  )
#define C41FS(i,nr)           VALUE(  C41_P(i) , nr  )
#define C42FS(i,nr)           VALUE(  C42_P(i) , nr  )
#define C43FS(i,nr)           VALUE(  C43_P(i) , nr  )
#define C44FS(i,nr)           VALUE(  C44_P(i) , nr  )
#define C60FS(i,nr)           VALUE(  C60_P(i) , nr  )
#define C61FS(i,nr)           VALUE(  C61_P(i) , nr  )
#define C62FS(i,nr)           VALUE(  C62_P(i) , nr  )
#define C63FS(i,nr)           VALUE(  C63_P(i) , nr  )
#define C64FS(i,nr)           VALUE(  C64_P(i) , nr  )
#define C65FS(i,nr)           VALUE(  C65_P(i) , nr  )
#define C66FS(i,nr)           VALUE(  C66_P(i) , nr  )
 
#define V20(i)               (  (i) -> V20  )
#define V21(i)               (  (i) -> V21  )
#define V22(i)               (  (i) -> V22  )
#define V40(i)               (  (i) -> V40  )
#define V41(i)               (  (i) -> V41  )
#define V42(i)               (  (i) -> V42  )
#define V43(i)               (  (i) -> V43  )
#define V44(i)               (  (i) -> V44  )
#define V60(i)               (  (i) -> V60  )
#define V61(i)               (  (i) -> V61  )
#define V62(i)               (  (i) -> V62  )
#define V63(i)               (  (i) -> V63  )
#define V64(i)               (  (i) -> V64  )
#define V65(i)               (  (i) -> V65  )
#define V66(i)               (  (i) -> V66  )
 
#define HMAG(i)               (   (i)-> h_magnetfeld   )
#define HAMILTONIAN(i)        (   (i)-> hamiltonian    )
#define PKQ(i)                (   (i)-> Pkq            )
 
 
#define KRISTALLFELD(nr) ((KRISTALLFELD*)c_alloc((nr),sizeof(KRISTALLFELD)))
#define ITERATION(u)          (   (u)-> iteration    )
#define SYMMETRIENR(u)        (   (u)-> symmetrienr  )
 
#define ABSD(x)         ( ((x)>=0.0) ? (x) : (-(x))  )
#define ABS(x)          ( ((x)>=0  ) ? (x) : (-(x))  )
 
#define ABSDR(m,ze,sp)  ABSD( R(m,ze,sp) )
#define ABSDI(m,ze,sp)  ABSD( I(m,ze,sp) )
 
#define ABSDRV(v,ze)    ABSD(RV(v,ze))
#define ABSDIV(v,ze)    ABSD(IV(v,ze))
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 
 
/*---------------------------------------------------------------------------*/
#define JA     1
#define NEIN   0
#define NOSPACE 99
 
#define NICHTIMP -1
 
#define TEST   1
#define NOTEST 0
 
#define AKQNAME     "akq.parmeter"
#define BKQNAME     "bkq.parmeter"
#define DKQNAME     "dkq.parmeter"
#define LKQNAME     "lkq.parmeter"
#define VKQNAME     "vkq.parmeter"
#define WKQNAME     "wkq.parmeter"
#define XWNAME      "xw.parmeter"
#define CFIELDSETUP "results/digits.setup"
#define ORTHO       "results/cfield.ortho"
#define CHI2        "results/cfield.chi2"
#define SUSZEPT     "results/suszept.rtplot"
#define MAGNETM     "results/moment.rtplot"
#define KPOLY       "results/kpoly.rtplot"
 
#define LOOPMAX      100  /*  Maximale Anzahl der moeglichen */
                          /*  Datensaetzen in Datenfiles     */
 
#define ANZ_NIVEAUS  17   /*  wird fuer CHI2  in Ortho.c */
#define ANZ_TRANS    17   /*  gebraucht.                 */
#define FITNR         1   /*  Defaultnummer der Fitroutine */
#define MACHEPSFACT 10000 /*  d.h macheps = (MACHEPSFACT*accuraccy())*/
#define MAXDOWNHILL 300
#define MAXPOWELL   100
#define MAXGRADIENT 100
#define MAXVA05A    100
#define MAXDUMMY    100
 
#define MDOWNHILL  10   /* zur menuesteuerung */
#define MPOWELL     2
#define MGRADIENT   2
#define MVA05A      2
#define MDUMMY      2
 
#define SINGLEION 's'
#define SIN       's'
#define AKQ       'a'
#define BKQ       'b'
#define DKQ       'd'
#define LKQ       'l'
#define VKQ       'v'
#define WKQ       'w'
#define XW        'x'
 
#define GGT_  'G'
#define GGT   'g'
#define NORM_ 'N'
#define NORM  'n'
 
#define NOMAG  ' '  /* kein Magnetfeld erwuenscht*/
 
#define NICHTANGELEGT  "not applied   "
#define ANGELEGT       "applied       "
#define FAILED         "failed        "
#define SUCCESSFUL     "successful    "
#define STABLE         "stable        "
#define UNSTABLE       "unstable      "
 
 
#define EPS1 0.001   /* es muss im kubischen Falle dann |V44- 5V40|<= EPS1 */
                     /*                            und  |V64+21V60|<= EPS1 */
/*---------------------------------------------------------------------------*/
 
/*=============================================================================
                             Benutzte Konstanten
=============================================================================*/
#define _k  1.38062      /* x10**(-23) J/K      Boltzmannkonstante k_B       */
#define _h  6.626075540  /* x10**(-34) Jsec     Planksches Wirkungsquantum h */
#define _hq (_h/2/pi)
#define _e  1.6021773349 /* x10**(-19) Coulomb  Elektronenladung e           */
#define _c  2.99792458   /* x10**(  8) m/sec    Lichtgeschwindigkeit c       */
#define _m  9.109389754  /* x10**(-31) Kg       Elektronenmasse m            */
#define _NA 6.022045     /* x10**( 23) Teilchen/mol Avogardokonstante*/
 
/*----------------------------------------------------------------------------
                  von den Konstanten abgeleitete Groessen
-----------------------------------------------------------------------------*/
                                          /*            2    */
                                          /*          hq     */
#define _a0 (100*_hq*_hq/_m/_e/_e/_c/_c ) /*  a0 =  ------   */ /* Bohrsche  */
                                          /*             2   */ /*  Radius   */
                                          /*          m*e    */ /* in Angstr.*/
 
 
#define _E0in_eV        ( _m*_e*_e*_e*_c*_c*_c*_c/_hq/_hq/100)    /* = e*e/a0 */
#define _E0programm     (- _E0in_eV)           /* Energieeinheit des Programms*/
 
#define _myBplus        (_hq/_m/2/1000)                      /* in eV / Tesla */
#define _myBprogramm    (_myBplus/_E0programm)               /* in  1 / Tesla */
 
/*
         h*c                -24
  E   = ----- = _h * _c * 10   J
   cm    1cm                                      _h *_c      -5
                                   => E  / E  =  -------- * 10
                       -19             cm   eV      _e
  E   =  1eV  = _e * 10    J
   eV
*/
 
#define FEcmeV  ( _h*_c/_e * 1e-5 )     /* 1Ecm = FEcmeV * 1eV   */
#define FEcmJ   ( _h * _c  * 1e-24  )   /* 1Ecm = FecmJ  * 1J    */
#define FeVJ    ( _e *       1e-19  )   /* 1eV  = FeVJ   * 1J    */
#define FeVEcm  ( 1/FEcmeV       )      /* 1eV  = FeVEcm * 1Ecm  */
 
 
 
#define FeVKelvin  ( 11604.9     )
#define FKelvineV  ( 1/FeVKelvin )
 
#define FcmKelvin  ( FEcmeV*FeVKelvin )
#define FKelvincm  ( 1/FcmKelvin )
 
#define FmeVKelvin  ( FeVKelvin/1000 )
#define FKelvinmeV  ( 1/FmeVKelvin   )
#define FKelvinKelvin  1
 
#define FThzmeV   ( 4.13570   )
#define FmeVThz   (1/FThzmeV  )
 
#define FThzKelvin  ( FThzmeV*FmeVKelvin )
#define FKelvinThz  ( 1/FThzKelvin )
 
#define FeVThz    (1/FThzmeV*1000  )
#define FThzeV    (FThzmeV/1000)
 
#define _my0   ( 1 ) /* x10**(-7) Vs/(Am)   = my0/4pi*/
 
#define myBeV      _myBplus            /* myB in eV/Tesla */
#define myBmeV    ( myBeV * 1000     ) /* myB in meV    /Tesla */
#define myBcm     ( myBeV * FeVEcm   ) /* myB in 1/Ecm  /Tesla */
#define myBKelvin ( myBeV * FeVKelvin) /* myB in Kelvin /Tesla */
#define myBThz    ( myBeV * FeVThz   ) /* myB in Thz    /Tesla */
#define myBJoule  ( myBeV * FeVJ     ) /* myB in Joule  /Tesla */
 
/* Spezifische Suszeptibilitaet ist Vielfaches  von XHI_0 */
/* Einheit von XHI_0 ist   emu  /mol                            */
/*       2                                                      */
/*  X = g * XHI_0 * Matrixelement , [X] = [XHI_0] = cm**3 / mol */
/*       j                                                      */
/*                                                              */
/*                      2     6                                 */
/*              my0  myB  * 10                                  */
/*  XHI_0   =   ---  -----------------  * N                     */
/*              4pi      k                 A                    */
/*                        B                                     */
/*                                                              */
/*                                                              */
/*          =  0.375150...  emu / mol                           */
/*                                                              */
 
#define XHI_0  ( myBeV*myBeV*_e*_e/_k*_NA*_my0*1e+7)
 
#define MYB        _myBprogramm
#define A0_BOHR    _a0
#define E0_EINHEIT _E0programm
 
 
#define delta(nj,mj) ( ((nj)==(mj))?  1 : 0  )
 
#define DDSIGN(a,b) ( ((b) > (0.0))?  ABSD(a) : (-ABSD(a))  )
 
#define DSIGN(n) ( ((n)>=(0.0))?  1.0 : -1.0  )
#define ISIGN(n) ( ((n)>=(0  ))?  1   : -1  )
#define MAX(a,b) ( ((a)>=(b  ))?  (a) : (b) )
#define MIN(a,b) ( ((a)<=(b  ))?  (a) : (b) )

/* For unused parameters, cast to void - to suppress warnings in GCC>4.5 */ 
#define UNUSED_PARAMETER(a) (void)a
 
/*------------------------------------------------------------------------------
ENDEMODUL    T Y P E S   C
------------------------------------------------------------------------------*/
