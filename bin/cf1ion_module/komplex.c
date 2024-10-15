/*-----------------------------------------------------------------------------
 
                                M   O  D  U  L
 
                                   KOMPLEX.C
 
-------------------------------------------------------------------------------
 
Aufgabe               :  Funktionen zum Umgang mit komplexen Zahlen
                         definieren
 
-------------------------------------------------------------------------------
 
Definierte Funktionen :
 
-----------------------
cdiv()
cabs()
csqroot()
-----------------------------------------------------------------------------*/
 
 
 
 
/*----------------------------------------------------------------------------
Includedateien holen
-----------------------------------------------------------------------------*/
#include <stdio.h>          /* damit FILE definiert wird               */
#include <stdlib.h>
#include <math.h>           /* damit sqrt in define_j.c definiert wird */
#define pi (4.0*atan(1.0))  /* atan() braucht <math.h>                 */
#include "types.c"          /* benutze Datentypen laden                */
/*----------------------------------------------------------------------------
Extern definierte Funktionen 
-----------------------------------------------------------------------------*/
extern INT    free_vr(VEKTOR *v);       /* definiert in MATRIX.C   */
extern VEKTOR *vr_alloc(INT n);     /* definiert in MATRIX.C   */
extern VEKTOR *_vr_copy(VEKTOR *a,VEKTOR *b);     /* definiert in MATRIX.C   */
 
/*----------------------------------------------------------------------------
                                     ckon()
------------------------------------------------------------------------------*/
                  /*                 *                   */
KOMPLEX *ckon(KOMPLEX *a)  /* Konjugation :  a   gespeichert in c */
{
   KOMPLEX *c;
 
   c = KX_ALLOC(1);
 
   RT(c) =  RT(a);
   IT(c) = -IT(a);
 
   return(c);
}
/*----------------------------------------------------------------------------
                                     cadd()
------------------------------------------------------------------------------*/
KOMPLEX *cadd(KOMPLEX *a,KOMPLEX *b)  /* komplexe Addition : a+b gespeichert in c */
{
   KOMPLEX *c;
 
   c = KX_ALLOC(1);
 
   RT(c) =  RT(a)+RT(b);
   IT(c) =  IT(a)+IT(b);
 
   return(c);
}
/*----------------------------------------------------------------------------
                                     csub()
------------------------------------------------------------------------------*/
KOMPLEX *csub(KOMPLEX *a,KOMPLEX *b)  /* komplexe Addition : a-b gespeichert in c */
{
   KOMPLEX *c;
 
   c = KX_ALLOC(1);
 
   RT(c) =  RT(a)-RT(b);
   IT(c) =  IT(a)-IT(b);
 
   return(c);
}
/*----------------------------------------------------------------------------
                                     cmult()
------------------------------------------------------------------------------*/
KOMPLEX *cmult(KOMPLEX *a,KOMPLEX *b)  /* komplexe Muliplikation: a*b gespeichert in c */
{
   KOMPLEX *c;
 
   c = KX_ALLOC(1);
 
   RT(c) =  RT(a)*RT(b) - IT(a)*IT(b);
   IT(c) =  RT(a)*IT(b) + IT(a)*RT(b);
 
   return(c);
}
/*----------------------------------------------------------------------------
                                     cdiv()
------------------------------------------------------------------------------*/
KOMPLEX *cdiv(DOUBLE xr,DOUBLE xi,DOUBLE yr,DOUBLE yi)    /* komplexe Division : x/y gespeichert in z */
{
    KOMPLEX *z;
    DOUBLE   h,nenner;
 
    z =   KX_ALLOC(1);  /* Platz fuer die komplexe Zahl z holen */
 
    if( ABSD(yr) > ABSD(yi) ){
         h     = yi / yr;
         nenner= yr * ( 1 + h*h );
         RT(z) = ( xr + h * xi )/nenner;
         IT(z) = ( xi - h * xr )/nenner;
    }
    else{
              h     = yr / yi;
              nenner= yi * ( 1 + h*h );
              RT(z) = ( xi + h * xr )/nenner;
              IT(z) = (-xr + h * xi )/nenner;
        }
 
    return(z);
}
/*----------------------------------------------------------------------------
                                     cabs()
------------------------------------------------------------------------------*/
DOUBLE cabs1(DOUBLE xr,DOUBLE xi)   /* cabs = +sqrt( xi*xi + xr*xr ) */
{
    DOUBLE h;
 
    xr = ABSD(xr);
    xi = ABSD(xi);
 
    if( xi > xr ){  /* xr = max( |xr| , |xi| )  */
         h  = xr;   /* xi = min( |xr| , |xi| )  */
         xr = xi;
         xi = h;
    }
 
    return(   (xi==0.0) ? xr : xr * sqrt( 1.0+(xi/xr)*(xi/xr) )   );
}
/*----------------------------------------------------------------------------
                                    csqroot()
------------------------------------------------------------------------------*/
KOMPLEX *csqroot(DOUBLE xr,DOUBLE xi)
{
    DOUBLE h,cabs1(DOUBLE xr,DOUBLE xi);
    KOMPLEX *y;
 
    y =  KX_ALLOC(1);  /* Platz fuer die komplexe Zahl y holen */
    h =  sqrt(  ( ABSD(xr)+ cabs1(xr,xi) )/2.0  );
 
    if( xi != 0.0 )
         xi /= 2.0 * h;
 
    if( xr >= 0.0 )
         xr = h;
    else if( xi >= 0.0 ){
              xr = xi;
              xi = h;
         }
         else {    xr = -xi;
                   xi = -h;
              }
 
    RT(y) = xr;
    IT(y) = xi;
 
    return( y );
}
/*----------------------------------------------------------------------------
                              vskalar()
-----------------------------------------------------------------------------*/
KOMPLEX *vskalar(VEKTOR *a,VEKTOR *b)  /* Multiplikation zweier komplexer */
     /* Vektoren <a| und |b> : <a|b>    */
{                      /*              +                  */
    KOMPLEX *c;        /* mit <a| = |a>                   */
    INT dim,n;         /* uebergeben werden 2 Spaltenvektoren */
                       /* a,b  = �a>,�b>                      */
                       /*                                     */
                       /*            2                        */
    dim = VRDIM(b);    /* mit || a ||  := <a|a>;              */
    c   = KX_ALLOC(1); /*                                     */
 
    RT(c) = 0.0;
    IT(c) = 0.0;
    for( n=dim  ; n>=1 ; --n  ){                   /*      +    */
       RT(c) += RV(a,n)*RV(b,n) + IV(a,n)*IV(b,n); /* (|a>) |b> */
       IT(c) += RV(a,n)*IV(b,n) - IV(a,n)*RV(b,n); /*           */
    }
 
    return(c);
}
/*----------------------------------------------------------------------------
                              cnorm()
-----------------------------------------------------------------------------*/
DOUBLE cnorm(VEKTOR *a)        /* Euklidnorm eines komplexen      */
        /* Vektors |a>                     */
{
    KOMPLEX *c;
    KOMPLEX *vskalar(VEKTOR *a,VEKTOR *b);
    DOUBLE  norm,sqrt(DOUBLE f);
 
    c     = vskalar(a,a);   /*  <a|a> */
    norm = RT(c);
    free_(c);
 
    return( sqrt(norm) );
 
}
/*----------------------------------------------------------------------------
                              cv_mult()
-----------------------------------------------------------------------------*/
VEKTOR *cv_mult( KOMPLEX *c,VEKTOR  *v)   /*  |w> := c * |v> */
{
       INT     zeile;
       KOMPLEX *d;
       KOMPLEX *cmult(KOMPLEX *a,KOMPLEX *b);
       VEKTOR  *vr_alloc(INT n),*w;
 
       w = vr_alloc( VRDIM(v) );
 
       for( zeile=1 ; zeile<=VRDIM(v) ; ++zeile ){
           d     = cmult( c,VR(v,zeile) );
           RV(w,zeile) = RT(d);
           IV(w,zeile) = IT(d);
           free_(d);
       }
 
       return( w );
}
/*----------------------------------------------------------------------------
                              vr_normalisieren()
-----------------------------------------------------------------------------*/
VEKTOR *vr_normalisieren(VEKTOR *v)
{
    DOUBLE  cnorm(VEKTOR *a),norm;
    KOMPLEX *c;
    VEKTOR  *cv_mult(KOMPLEX *c,VEKTOR  *v);
    VEKTOR  *_vr_copy(VEKTOR *a,VEKTOR *b);
    VEKTOR  *w;
/*  INT     i; */
 
    c     = KX_ALLOC(1);
    norm  = cnorm(v);
 
    if( norm==0.0 ){
 
       printf("\nError in vr_normalisieren in KOMPLEX.C .\n");
       printf("Norm of the Vectors is Null.\n\n");
       exit(1);
    }
 
    RT(c) = 1/norm;
    IT(c) = 0.0;
 
    w = cv_mult(c,v);
    v = _vr_copy(v,w);
 
    free_vr(w);
    free_(c);
 
    return(v);
}
/*------------------------------------------------------------------------------
ENDEMODUL    K O M P L E X    C
------------------------------------------------------------------------------*/
