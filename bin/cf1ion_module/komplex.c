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
extern INT    free_vr();       /* definiert in MATRIX.C   */
extern VEKTOR *vr_alloc();     /* definiert in MATRIX.C   */
extern VEKTOR *_vr_copy();     /* definiert in MATRIX.C   */
 
/*----------------------------------------------------------------------------
                                     ckon()
------------------------------------------------------------------------------*/
                  /*                 *                   */
KOMPLEX *ckon(a)  /* Konjugation :  a   gespeichert in c */
   KOMPLEX *a;
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
KOMPLEX *cadd(a,b)  /* komplexe Addition : a+b gespeichert in c */
   KOMPLEX *a;
   KOMPLEX *b;
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
KOMPLEX *csub(a,b)  /* komplexe Addition : a-b gespeichert in c */
   KOMPLEX *a;
   KOMPLEX *b;
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
KOMPLEX *cmult(a,b)  /* komplexe Muliplikation: a*b gespeichert in c */
   KOMPLEX *a;
   KOMPLEX *b;
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
KOMPLEX *cdiv(xr,xi,yr,yi)    /* komplexe Division : x/y gespeichert in z */
    DOUBLE   xr;
    DOUBLE   xi;
    DOUBLE   yr;
    DOUBLE   yi;
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
DOUBLE cabs1(xr,xi)   /* cabs = +sqrt( xi*xi + xr*xr ) */
    DOUBLE   xr;
    DOUBLE   xi;
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
KOMPLEX *csqroot(xr,xi)
    DOUBLE   xr;
    DOUBLE   xi;
{
    DOUBLE h,cabs1();
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
KOMPLEX *vskalar(a,b)  /* Multiplikation zweier komplexer */
    VEKTOR *a,*b;      /* Vektoren <a| und |b> : <a|b>    */
{                      /*              +                  */
    KOMPLEX *c;        /* mit <a| = |a>                   */
    INT dim,n;         /* uebergeben werden 2 Spaltenvektoren */
                       /* a,b  = ºa>,ºb>                      */
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
DOUBLE cnorm(a)        /* Euklidnorm eines komplexen      */
    VEKTOR *a;         /* Vektors |a>                     */
{
    KOMPLEX *c;
    KOMPLEX *vskalar();
    DOUBLE  norm,sqrt();
 
    c     = vskalar(a,a);   /*  <a|a> */
    norm = RT(c);
    free_(c);
 
    return( sqrt(norm) );
 
}
/*----------------------------------------------------------------------------
                              cv_mult()
-----------------------------------------------------------------------------*/
VEKTOR *cv_mult(c,v)   /*  |w> := c * |v> */
       KOMPLEX *c;
       VEKTOR  *v;
{
       INT     zeile;
       KOMPLEX *d;
       KOMPLEX *cmult();
       VEKTOR  *vr_alloc(),*w;
 
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
VEKTOR *vr_normalisieren(v)
       VEKTOR *v;
{
    DOUBLE  cnorm(),norm;
    KOMPLEX *c;
    VEKTOR  *cv_mult();
    VEKTOR  *_vr_copy();
    VEKTOR  *w;
    INT     i;
 
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
