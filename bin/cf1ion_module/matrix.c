/*-----------------------------------------------------------------------------
 
                                M   O  D  U  L
 
                                   MATRIX C
 
-------------------------------------------------------------------------------
 
Aufgabe               :  Funktionen fuer den Umgang mit komplexen
                         Matrizen definieren
 
-------------------------------------------------------------------------------
 
Definierte Funktionen :
 
-----------------------
 
c_alloc()     : Kontroliert Speicherplatz holen
error()       : Fehlerbehandlung
mx_alloc()    : Reserviert Speicherplatz fuer eine komplexe Matrix
mx_copy()     : Dubliziert eine komplexe Matrix
testij()      : Testet Spalten- und Zeilenindex
show_mij()    : Zeigt kontrolliert ij-tes Element der Matrix m
                unkontrolliert :    MX(...)  bzw    MX(...)
set_mij()     : Setzt kontrolliert ij-tes Element der Matrix m
                unkontrolliert : RT(MX(...)) bzw IT(MX(...))
vr_alloc()    : Reserviert Speicher fuer komplexen Vektor
warning()     : Warnung ausgeben
 
 
 
-----------------------------------------------------------------------------*/
 
 
 
 
/*----------------------------------------------------------------------------
Includedateien holen
-----------------------------------------------------------------------------*/
#include <stdio.h>          /* damit FILE definiert wird               */
#include <stdlib.h>
#include <math.h>           /* damit sqrt in define_j.c definiert wird */
#include <limits.h>
#define pi (4.0*atan(1.0))  /* atan() braucht <math.h>                 */
#include "types.c"          /* benutze Datentypen laden                */
 
/*----------------------------------------------------------------------------
Extern definierte Funktionen
-----------------------------------------------------------------------------*/
extern DOUBLE is_equal();
extern DOUBLE accuracy();
extern SPEICHER *calloc();
/*----------------------------------------------------------------------------
                                   c_alloc()
------------------------------------------------------------------------------*/
CHAR* c_alloc( anzahl,laenge )
    INT anzahl;
    INT laenge;
{
    CHAR* cp;
    INT error();
 
    if( (cp=(CHAR*)calloc((TYPUS)anzahl,(TYPUS)laenge)) == (CHAR*)0 )
         error("No more space avaliable on disk - Programm aborted !");
 
    return(cp);
}
/*----------------------------------------------------------------------------
                                   warning()
------------------------------------------------------------------------------*/
void warning(s)  /* Warnung ausgeben */
    CHAR *s;
{
    printf("\nWarning : ");
    printf("%s,",s);
    printf("\n");
}
/*----------------------------------------------------------------------------
                                    error()
------------------------------------------------------------------------------*/
INT error(s)  /* Ausstieg beim Auftreten eines Fehlers */
    CHAR *s;
{
    printf("\nError : ");
    printf("%s",s);
    printf("\n");
    exit(1);
}
/*----------------------------------------------------------------------------
                                   mx_copy()
------------------------------------------------------------------------------*/
MATRIX *mx_copy(a)      /* kopiert Matrix a in neue Matrix m */
    MATRIX *a;
{
    MATRIX *m;          /* neue Matrix m */
    MATRIX *mx_alloc(); /* Funktion,welche Speicherplatz fuer eine Matrix holt*/
    INT    ze,sp;
    INT    anz_ze;
    INT    anz_sp;
 
    anz_sp = ANZ_SP(a);              /* Anzahl der Spalten der alten Matrix a */
    anz_ze = ANZ_ZE(a);              /* Anzahl der Zeilen  der alten Matrix a */
    m      = mx_alloc(anz_ze,anz_sp);/* Platz fuer neue Matrix m holen */
 
    for(sp=0;sp<=anz_sp;++sp)         /* Matrix von a nach m kopieren */
         for(ze=0;ze<=anz_ze;++ze){
              R(m,ze,sp) = R(a,ze,sp);
              I(m,ze,sp) = I(a,ze,sp);
         }
 
    return( m );        /* neue Matrix m zurueckgeben */
}
/*----------------------------------------------------------------------------
                                   mx_add()
------------------------------------------------------------------------------*/
MATRIX *mx_add(fac,a,b)      /* c :=  fac*(a + b)  */
    DOUBLE fac;
    MATRIX *a;
    MATRIX *b;
{
    MATRIX *c;          /* neue Matrix c */
    MATRIX *mx_alloc(); /* Funktion,welche Speicherplatz fuer eine Matrix holt*/
    INT    ze,sp;
    INT    anz_ze ;
    INT    anz_sp ;
 
    anz_sp = ANZ_SP(a);              /* Anzahl der Spalten der alten Matrix a */
    anz_ze = ANZ_ZE(a);              /* Anzahl der Zeilen  der alten Matrix a */
 
 
 
    c      = mx_alloc(anz_ze,anz_sp);/* Platz fuer neue Matrix c holen */
 
    for(sp=0;sp<=anz_sp;++sp)         /* Matrix von a nach m kopieren */
         for(ze=0;ze<=anz_ze;++ze){
              R(c,ze,sp) = fac*( R(a,ze,sp) + R(b,ze,sp) );
              I(c,ze,sp) = fac*( I(a,ze,sp) + I(b,ze,sp) );
         }
 
    return( c );        /* neue Matrix c zurueckgeben */
}
/*----------------------------------------------------------------------------
                                   mx_sub()
------------------------------------------------------------------------------*/
MATRIX *mx_sub(fac,a,b)      /* c :=  fac*(a - b)  */
    DOUBLE fac;
    MATRIX *a;
    MATRIX *b;
{
    MATRIX *c;          /* neue Matrix c */
    MATRIX *mx_alloc(); /* Funktion,welche Speicherplatz fuer eine Matrix holt*/
    INT    ze,sp;
    INT    anz_ze ;
    INT    anz_sp ;
 
    anz_sp = ANZ_SP(a);              /* Anzahl der Spalten der alten Matrix a */
    anz_ze = ANZ_ZE(a);              /* Anzahl der Zeilen  der alten Matrix a */
 
 
 
    c      = mx_alloc(anz_ze,anz_sp);/* Platz fuer neue Matrix c holen */
 
    for(sp=0;sp<=anz_sp;++sp)         /* Matrix von a nach m kopieren */
         for(ze=0;ze<=anz_ze;++ze){
              R(c,ze,sp) = fac*( R(a,ze,sp) - R(b,ze,sp) );
              I(c,ze,sp) = fac*( I(a,ze,sp) - I(b,ze,sp) );
         }
 
    return( c );        /* neue Matrix c zurueckgeben */
}
/*----------------------------------------------------------------------------
                                   mx_addf()
------------------------------------------------------------------------------*/
MATRIX *mx_addf(f1,f2,a,b)      /* c := f1*a + f2*b  */
    DOUBLE f1,f2;
    MATRIX *a;
    MATRIX *b;
{
    MATRIX *c;          /* neue Matrix c */
    MATRIX *mx_alloc(); /* Funktion,welche Speicherplatz fuer eine Matrix holt*/
    INT    ze,sp;
    INT    anz_ze ;
    INT    anz_sp ;
 
    anz_sp = ANZ_SP(a);              /* Anzahl der Spalten der alten Matrix a */
    anz_ze = ANZ_ZE(a);              /* Anzahl der Zeilen  der alten Matrix a */
 
 
 
    c      = mx_alloc(anz_ze,anz_sp);/* Platz fuer neue Matrix c holen */
 
    for(sp=0;sp<=anz_sp;++sp)         /* Matrix von a nach m kopieren */
         for(ze=0;ze<=anz_ze;++ze){
              R(c,ze,sp) = f1*R(a,ze,sp) + f2*R(b,ze,sp) ;
              I(c,ze,sp) = f1*I(a,ze,sp) + f2*I(b,ze,sp) ;
         }
 
    return( c );        /* neue Matrix c zurueckgeben */
}
/*----------------------------------------------------------------------------
                                   mx_mult()
------------------------------------------------------------------------------*/
MATRIX *mx_mult(a,b)      /* c := a * b  */
    MATRIX *a;
    MATRIX *b;
{
    MATRIX *c;          /* neue Matrix c */
    MATRIX *mx_alloc(); /* Funktion,welche Speicherplatz fuer eine Matrix holt*/
    INT    ze,sp,k;
    INT    dim;
 
    dim    = MXDIM(a);              /* Anzahl der Spalten der alten Matrix a */
 
 
 
    c      = mx_alloc(dim,dim);/* Platz fuer neue Matrix c holen */
 
    for(sp=0;sp<=dim;++sp)
         for(ze=0;ze<=dim;++ze){
              R(c,ze,sp) = 0.0;
              I(c,ze,sp) = 0.0;
              for(k=0;k<=dim;++k){
                 R(c,ze,sp) += R(a,ze,k)*R(b,k,sp) - I(a,ze,k)*I(b,k,sp);
                 I(c,ze,sp) += I(a,ze,k)*R(b,k,sp) + R(a,ze,k)*I(b,k,sp);
              }
 
         }
 
    return( c );        /* neue Matrix c zurueckgeben */
}
/*----------------------------------------------------------------------------
                                   mx_ewev()
------------------------------------------------------------------------------*/
MATRIX *mx_ewev(ew,ev,shift)/*       |         º                º   | */
    VEKTOR *ew;             /*  c := | ew(1)*ev(1)  ... ew(n)*ev(n) | */
    MATRIX *ev;             /*       |         º                º   | */
    DOUBLE shift;
{
    MATRIX *c;          /* neue Matrix c */
    MATRIX *mx_alloc(); /* Funktion,welche Speicherplatz fuer eine Matrix holt*/
    INT    ze,sp;
    INT    dim;
 
    dim    = MXDIM(ev);             /* Anzahl der Spalten der alten Matrix a */
 
 
 
    c      = mx_alloc(dim,dim);/* Platz fuer neue Matrix c holen */
 
    for(sp=1;sp<=dim;++sp)
         for(ze=1;ze<=dim;++ze){
            R(c,ze,sp) = (RV(ew,sp)+shift)*R(ev,ze,sp) - IV(ew,sp)*I(ev,ze,sp);
            I(c,ze,sp) = IV(ew,sp)*R(ev,ze,sp) + (RV(ew,sp)+shift)*I(ev,ze,sp);
 
         }
 
    return( c );        /* neue Matrix c zurueckgeben */
}
/*----------------------------------------------------------------------------
                                  mx_alloc()
------------------------------------------------------------------------------*/
MATRIX *mx_alloc(anz_ze,anz_sp)   /* Holt Speicherplatz fuer eine komplexe */
                                  /* Matrix der Dimension anz_ze x anz_sp  */
    INT anz_ze;         /* Anzahl der Zeilen  der benoetigten Matrix */
    INT anz_sp;         /* Anzahl der Spalten der benoetigten Matrix */
{
    MATRIX *m;          /* neue Matrix m */
    INT    ze,sp;
/*  INT    warning(); */
 
    if( anz_ze==0 || anz_sp==0 ){
         warning("Es wurde eine Matrix mit Dimension 0 definiert");
         return( (MATRIX*)0 );
    }
 
/* Platz fuer Matrix m holen */
 
    m          = MX_ALLOC(1);       /* Platz fuer MATRIX-Struktur */
    ANZ_SP(m)  = anz_sp;            /* Anzahl der Spalten von m in m notieren*/
    m->_spalte = VR_ALLOC(anz_sp+1);/* Platz fuer Pointer auf die Spalten von*/
    for(sp=0;sp<=anz_sp;++sp){       /* m holen */
          ANZ_ZN(m,sp) = anz_ze;
          MXZE0(m,sp)  = KX_ALLOC(anz_ze+1);/* Spaltenweise Platz fuer die*/
    }                                       /* einzelen Zeilen von m holen*/
 
/* Alle Komponenten der Matrix m auf 0 setzen */
 
    for(sp=0;sp<=ANZ_SP(m);++sp)
         for(ze=0;ze<=ANZ_ZN(m,sp);++ze)
              R(m,ze,sp) = I(m,ze,sp) = 0.0;
 
    return( m );   /* die komplexe Nullmatrix m zurueckgeben */
}
/*----------------------------------------------------------------------------
                                  free_mx()
------------------------------------------------------------------------------*/
void free_mx(mx)  /* Speicherplatz der Matrix mx freigeben */
  MATRIX *mx;
{
    INT    sp;
 
    for(sp=0;sp<=ANZ_SP(mx);++sp)
          free_( MXZE0(mx,sp) );
 
    free_(mx->_spalte);
    free_(mx);
}
/*----------------------------------------------------------------------------
                                  vr_alloc()
------------------------------------------------------------------------------*/
VEKTOR *vr_alloc(n) /* holt Platz fuer einen n-dimensionalen komplexen */
    INT n;          /* Spaltenvektor */
{
    INT    zeile;
    VEKTOR *v;
 
/* Platz fuer den komplexen Vektor holen */
    v        = VR_ALLOC(1);  /* Platz fuer die VEKTOR-Struktur holen */
    v->_zeile= KX_ALLOC(n+1);/* Platz fuer die n komplexen Komponenten holen*/
    VRDIM(v) = n;            /* merken ,dass es n Komponenten sind */
 
/* Vektorkomponenten auf Null setzen */
    for(zeile=0;zeile<=n;++zeile)
         RV(v,zeile) = IV(v,zeile) = 0.0;
 
    return( v );
}
/*----------------------------------------------------------------------------
                                 _vr_copy();
------------------------------------------------------------------------------*/
VEKTOR *_vr_copy(a,b)  /*  b nach a kopieren    */
       VEKTOR *a;
       VEKTOR *b;
{
   INT    zeile/*,i*/;
 
 
    for( zeile=0 ; zeile<=VRDIM(a) ; ++zeile){
         RV(a,zeile) = RV(b,zeile);
         IV(a,zeile) = IV(b,zeile);
    }
 
    return( a );
}
/*----------------------------------------------------------------------------
                                  vr_copy();
------------------------------------------------------------------------------*/
VEKTOR *vr_copy(a)  /*  a duplizieren */
       VEKTOR *a;
{
   VEKTOR *vr_alloc();
   VEKTOR *c;
   INT    zeile;
 
   c = vr_alloc( VRDIM(a) );
 
    for( zeile=0 ; zeile<=VRDIM(a) ; ++zeile){
         RV(c,zeile) = RV(a,zeile);
         IV(c,zeile) = IV(a,zeile);
    }
    return( c );
}
/*----------------------------------------------------------------------------
                                  vr_sub();
------------------------------------------------------------------------------*/
VEKTOR *vr_sub(a,b)  /*  c = a-b */
       VEKTOR *a,*b;
{
   VEKTOR *vr_alloc();
   VEKTOR *c;
   INT    zeile;
 
   c = vr_alloc( VRDIM(a) );
 
    for( zeile=0 ; zeile<=VRDIM(a) ; ++zeile){
         RV(c,zeile)  = RV(a,zeile) - RV(b,zeile);
         IV(c,zeile)  = IV(a,zeile) - IV(b,zeile);
    }
    return( c );
}
/*----------------------------------------------------------------------------
                                  free_vr()
------------------------------------------------------------------------------*/
void free_vr(v)        /* Speicherplatz von v freigeben */
 VEKTOR *v;
{
   free_(v->_zeile );
   free_( v );
}
/*----------------------------------------------------------------------------
                                   ludcmp()
aus numerical recipes, seite 35-36
------------------------------------------------------------------------------*/
LUDCMP *ludcmp(mx)
    MATRIX *mx;
{
   LUDCMP *ludcmp;
 
   MATRIX *a, *mx_copy();
   VEKTOR *vr_alloc(), *vr_copy(), *index, *vv;
   INT    n,i,j,k,invers,imax=-INT_MAX;
   DOUBLE accuracy(), tiny;
   DOUBLE aamax;
   DOUBLE d,sum,dum,is_equal();
 
/* INT    nmax = 100; / * maximale  Iterationsschritte */
 
 
   invers = JA;  /* inverse Matrix exestiert */
   tiny   = (MACHEPSFACT*accuracy());
   n      = MXDIM(mx);
   a      = mx_copy(mx);
 
   vv     = vr_alloc(n);
   index  = vr_alloc(n);
   d      = 1.0;
 
   for( i=1; i<=n; ++i ){
        aamax = 0.0;
        for( j=1; j<=n; ++j )
             if( ABSD(R(a,i,j)) > aamax ) aamax = ABSD(R(a,i,j));
        if( is_equal(aamax,0.0,tiny) ){
             invers = NEIN;
             goto ende;
        }
        RV(vv,i) = 1.0/aamax;
   }
 
   for( j=1; j<=n; ++j ){
        for( i=1; i<=j-1; ++i ){
             sum = R(a,i,j);
             for( k=1; k<=i-1; ++k )
                  sum -= R(a,i,k)*R(a,k,j);
             R(a,i,j) = sum;
        }
        aamax = 0.0;
        for( i=j; i<=n; ++i ){
             sum = R(a,i,j);
             for( k=1; k<=j-1; ++k )
                  sum -= R(a,i,k)*R(a,k,j);
             R(a,i,j) = sum;
             dum = RV(vv,i)*ABSD(sum);
             if( dum >= aamax ) {imax=i;aamax=dum;}
        }
        if( j != imax ){
            for( k=1; k<=n; ++k ){
                 dum         = R(a,imax,k);
                 R(a,imax,k) = R(a,j   ,k);
                 R(a,j   ,k) = dum;
            }
            d *= -1.0;
            RV(vv,imax) = RV(vv,j);
        }
        RV(index,j) = (DOUBLE)imax;
        if( is_equal(R(a,j,j),0.0,tiny) )  R(a,j,j) = tiny;
        if( j != n ){
            dum = 1.0/R(a,j,j);
            for( i=j+1; i<=n; ++i )
                 R(a,i,j) *= dum;
        }
   }
 
 
ende:
   free_vr(vv);
 
   ludcmp         = LUDCMP_ALLOC(1);
   ludcmp->invers = invers;
   ludcmp->matrix = a;
   ludcmp->index  = index;
   ludcmp->d      = d;
 
   return(ludcmp);
}
/*----------------------------------------------------------------------------
                                   lubksb()
aus numerical recipes, seite 37
------------------------------------------------------------------------------*/
VEKTOR *lubksb(ludcmp,v)
   LUDCMP *ludcmp;
   VEKTOR *v;
{
   VEKTOR *b, *vr_alloc(), *index;
   MATRIX *a;
   INT    n,i,j,ii,ll;
   DOUBLE /*d,*/sum,macheps,accuracy(),is_equal();
 
 
   a      = ludcmp->matrix;
   index  = ludcmp->index;
/* d      = ludcmp->d; */
   macheps= (MACHEPSFACT*accuracy());
 
   n = MXDIM(a);
   b = vr_copy(v);
 
   ii = 0;
   for( i=1; i<=n; ++i ){
        ll       = (INT)RV(index,i);
        sum      = RV(b,ll);
        RV(b,ll) = RV(b,i);
        if( ii!=0 ){
            for( j=ii; j<=i-1; ++j )
                 sum -= R(a,i,j)*RV(b,j);
        }
        else if( !is_equal(sum,0.0,macheps) )  ii = i;
        RV(b,i) = sum;
   }
 
   for( i=n; i>=1; --i ){
        sum = RV(b,i);
        if( i<= n ){
            for( j=i+1; j<=n; ++j )
                 sum -= R(a,i,j)*RV(b,j);
        }
        RV(b,i) = sum/R(a,i,i);
   }
 
   return(b);
}
/*----------------------------------------------------------------------------
                                   lubksb()
aus numerical recipes, seite 37
------------------------------------------------------------------------------*/
VEKTOR *loese(a,b)      /* loese  Ax = b */
   MATRIX a;
   VEKTOR b;
{
   LUDCMP *_ludcmp, *ludcmp();
   VEKTOR *x,       *lubksb();
 
   _ludcmp = ludcmp(a);
   x       = (VEKTOR*)0;
   if( _ludcmp->invers == JA ){  /* inverse matrix exestiert */
       x = lubksb(_ludcmp,b);
   }
 
   free_vr( _ludcmp->index );
   free_mx( _ludcmp->matrix);
   free_( _ludcmp );
 
   return( x );
}
/*------------------------------------------------------------------------------
ENDEMODUL    M A T R I X    C
------------------------------------------------------------------------------*/
