/*-----------------------------------------------------------------------------
 
                                M  O  D  U  L
 
                                  STEVENS C
 
-------------------------------------------------------------------------------
 
Aufgabe               : verallgemeinerte Stevensoperatoren definieren
 
-------------------------------------------------------------------------------
Literatur             : [1] D.Smith and J.H.M. Thornley
                            "The use of operator equivalents",
                            Proc.Phys.1966,Vol.89
 
                        [2] Diplomarbeit Peter Hoffmann 1989
                            Kfa Juelich ,Iff
-------------------------------------------------------------------------------
 
Definierte Funktionen :
 
-----------------------
 
calc_stevens_op() : Matrixelemente (JM'|Pkq(J)|MJ) berechnen und speichern
info_stevens_op() : Matrixelemente (JM'|Pkq(J)|MJ) ausgeben
 
Clm()                : Tensor Clm berechnen
info_tensor_Clm()    : Tensor Clm ansehen
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
extern MATRIX *mx_alloc(INT anz_ze,INT anz_sp);
extern MATRIX *mx_add(DOUBLE z, MATRIX*a,MATRIX*b);
extern MATRIX *mx_sub(DOUBLE z, MATRIX*a,MATRIX*b);
extern MATRIX *mx_addf(DOUBLE f1,DOUBLE f2,MATRIX*m1,MATRIX*m2);
extern INT    free_mx(MATRIX * m);
extern INT    is_equal(DOUBLE a,DOUBLE b,DOUBLE macheps);
 
extern DOUBLE omegan0n(DOUBLE f);   /*           */
extern DOUBLE omegan1n(DOUBLE f);   /*  omega    */
extern DOUBLE omegan2n(DOUBLE f);   /*        kq */
extern DOUBLE omegan3n(DOUBLE f);
extern DOUBLE omegan4n(DOUBLE f);
extern DOUBLE omegan5n(DOUBLE f);
extern DOUBLE omegan6n(DOUBLE f);
 
extern IONEN  IONENIMP[];
 
extern FILE *fopen_errchk(const char * filename,const char * mode);         /* definiert in EINGABE.C*/ 

/*----------------------------------------------------------------------------
   Internal function declarations
-----------------------------------------------------------------------------*/
void showPkq(FILE *fp, MATRIX *mx, INT k, INT q, INT dimj, CHAR modus);
void showSTEVkq(FILE *fp, MATRIX *mx, INT k, INT q, INT dimj, CHAR modus);

/*----------------------------------------------------------------------------
                               norm_Pkq()
------------------------------------------------------------------------------*/
/*
 
                       J                              2
                    ------   || < JM | P (J) | JN > ||
                2    \                  kq
     || P (J) ||  :=  >      -----------------------
         kq          /
                    ------          2 * J + 1
                    N,M = -J
 
 
                                 2                         +
   mit || < JM | P (J) | JN > ||  = ( < JM | P (J) | JN > )  < JM | P (J) | JN >
                  kq                          kq                     kq
 
                                              +
                                  =   < JN | P (J) | JM >  < JM | P (J) | JN >
                                              kq                   kq
 
allgemein fuer einen nicht hermiteschen operator A:
 
            J                               J
           ---            +       +        ---
      2    \    0.5*<JN| A A + A A |JN>    \    <JN|Re(A)^2+Im(A)^2|JN>
 ��A��   := >   -----------------------  =  >   -----------------------
           /           2J+1                /          2J+1
           ---                             ---
           N=-J                            N=-J
 
--------------------------------------------------------------------------------
*/
DOUBLE norm_Pkq(INT k,INT q,INT dimj)   /*    || Pkq ||   definieren  */
{
   MATRIX *pkq_plus;
   MATRIX *pkq_minus;
   MATRIX *Pkq(INT k,INT q,INT dimj);
   DOUBLE sum=0.0,sqrt(DOUBLE f);
   DOUBLE exp,exp1,exp2,test_exp,j,pow_(DOUBLE x,INT n) ;
   INT n,m;
   DOUBLE dummy;
                               /*                       */
   pkq_plus  = Pkq(k, q,dimj); /*  < JN | P (J) | JM >  */
                               /*          kq           */
 
                               /*          +            */
   pkq_minus = Pkq(k,-q,dimj); /*  < JN | P (J) | JM >  */
                               /*          kq           */
 
   /* 10-er Potenz der Matrixelemente finden */
   exp1 = 0.0;
   for( n=dimj ; n>=1 ; --n )
      for( m=dimj ; m>=1 ; --m ){
           dummy = ABSDR( pkq_plus,n,m );
           if( dummy >= 10.0 ){
               test_exp = log( dummy )/log(10.0);
               if( test_exp > exp1 )  exp1 = test_exp;
           }
      }
   /* 10-er Potenz der Matrixelemente finden */
   exp2 = 0.0;
   for( n=dimj ; n>=1 ; --n )
      for( m=dimj ; m>=1 ; --m ){
           dummy = ABSDR( pkq_minus,n,m );
           if( dummy >= 10.0 ){
               test_exp = log( dummy )/log(10.0);
               if( test_exp > exp2 )  exp2 = test_exp;
           }
      }
 
    exp = exp1;
    if( exp2 > exp1 ) exp = exp2;
 
      if( exp>=10.0 ){
           j = ((DOUBLE)dimj-1)/2;
           printf("\n");
           printf("cannot calculate the Norm of P%1d%1d(%3.1f) ",k,q,j);
           printf(", \n");
           printf("because the matrix elements of the operators are too big.\n");
           printf("\n");
           exit(1);
      }
 
   for( n=dimj ; n>=1 ; --n )
      for( m=dimj ; m>=1 ; --m )
           sum +=  (R( pkq_minus,m,n )/pow_(10.0,(INT)exp))
                  *(R( pkq_plus ,n,m )/pow_(10.0,(INT)exp));
 
   sum /= dimj;
 
 
   free_mx( pkq_plus  );
   free_mx( pkq_minus );
 
 
   return( pow_(10.0,(INT)exp)*sqrt(sum) );
}
/*----------------------------------------------------------------------------
                               norm_stevkq()
------------------------------------------------------------------------------*/
/*
                                            +
    sei hier O (J) := STEV (J) = ( P (J) + P (J) )/2/theta
              kq          kq        kq      kq            kq
                       J
                                                      2
                    ------   || < JM | O (J) | JN > ||
                2    \                  kq
     || O (J) ||  :=  >      -----------------------
         kq          /
                    ------          2 * J + 1
                    N,M = -J
 
 
                                 2                         +
   mit || < JM | O (J) | JN > ||  = ( < JM | O (J) | JN > )  < JM | O (J) | JN >
                  kq                          kq                     kq
 
 
                                  =   < JN | O (J) | JM >  < JM | O (J) | JN >
                                              kq                   kq
 
 also :
 
 
 
                      J
                    ------   < JN | O (J) | JM >*< JM | O (J) | JN >
                2    \               kq                  kq
     || O (J) ||  :=  >      ---------------------------------------
         kq          /
                    ------                  2 * J + 1
                    M = -J
 
 
 
--------------------------------------------------------------------------------
*/
DOUBLE norm_Stevkq(INT k,INT q,INT dimj)   /*    || STEVkq ||   definieren  */
{
   MATRIX *stev;
   MATRIX *stevkq(INT k,INT q,INT dimj);
   DOUBLE norm_operator(MATRIX *mx),norm;
                               /*                       */
   stev = stevkq(k,q,dimj);    /*  < JN | O (J) | JM >  */
                               /*          kq           */
 
 
   norm = norm_operator(stev);
   free_mx( stev );
 
   return( norm );
}
/*----------------------------------------------------------------------------
                               norm_operator(MATRIX *mx)
------------------------------------------------------------------------------*/
DOUBLE norm_operator(MATRIX *mx)
  /* mx ist die Matrixdarstellung des operators dessen */
                /* Norm berechnet werden soll                        */
{  DOUBLE  macheps=0.0001;
   DOUBLE sum=0.0,sqrt(DOUBLE f);
   DOUBLE exp,test_exp/*,j*/,pow_(DOUBLE x,INT n) ;
   DOUBLE expr,expi;
   INT n,m,dimj,is_equal(DOUBLE a,DOUBLE b,DOUBLE macheps);
   DOUBLE dummy;
 
   dimj = MXDIM(mx);
 
   /* 10-er Potenz der Matrixelemente finden */
   expr = 0.0;
   for( n=dimj ; n>=1 ; --n )
      for( m=dimj ; m>=1 ; --m ){
           dummy = ABSDR( mx  ,n,m );
           if( dummy > 1.0 ) {
               test_exp = log( dummy )/log(10.0);
           if( test_exp > expr )  expr = test_exp; }
      }
   expi= 0.0;
   for( n=dimj ; n>=1 ; --n )
      for( m=dimj ; m>=1 ; --m ){
           dummy = ABSDI( mx  ,n,m );
           if( dummy > 1.0 ) {
               test_exp = log( dummy )/log(10.0);
           if( test_exp > expi )  expi = test_exp; }
      }
 
   exp = expr;
   if( expi > expr ) exp = expi;
 
   if( exp>=10.0 ){
           printf("\n");
           printf("Error in: norm_operator() in stevens.c !\n");
           printf("cannot calculate the norm, \n");
           printf(".\n");
           printf("because the matrix elements of the operators are too big\n");
           exit(1);
   }
 
   for( n=dimj ; n>=1 ; --n )
      for( m=dimj ; m>=1 ; --m ){
           sum += (R( mx ,n,m )/pow_(10.0,(INT)exp))
                 *(R( mx ,n,m )/pow_(10.0,(INT)exp));
           sum += (I( mx ,n,m )/pow_(10.0,(INT)exp))
                 *(I( mx ,n,m )/pow_(10.0,(INT)exp));
      }
 
/*
   if( sum<= 0.0 ){
           printf("\n");
           printf("Error in: norm_operator() in stevens.c !\n");
           printf("cannot calculate the norm, \n");
           printf("because sum < 0 !\n");                                     ;
           printf("\n");
           exit(1);
   }
*/
   sum /= dimj;
 
 
   sum =  pow_(10.0,(INT)exp)*sqrt(sum);
 
   if( is_equal( sum, 0.0, macheps) ) sum=1.0;
   return( sum );
}
/*----------------------------------------------------------------------------
                               fkq_tabelle(INT anz_ionen)
------------------------------------------------------------------------------*/
void fkq_tabelle(INT anz_ionen)/* Fkq Faktoren fuer die implementierten Ionen ausgeben */
{
 INT  ionennr;
 INT  dimj;
 CHAR *ionname;
 FILE    *fopen(const char * filename,const char * mode),*out;
 CHAR *filename = "Fkq.tabelle";
 CHAR *t01,*t02,*t03,*t04,*t05;
 LONG   fkq(INT k,INT q,INT dimj);
 LONG   f20,f21,f22;
 LONG   f40,f41,f42,f43,f44;
 LONG   f60,f61,f62,f63,f64,f65,f66;
 
 out = fopen_errchk(filename,"w");
 
 clearscreen;
printf("Information is contained in the file %s .\n",filename);
 
fprintf(out,"========================================================\n");
fprintf(out,"| Information on the factors  F                        |\n");
fprintf(out,"|                              kq                      |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"| F (J)  :=  GGT(  { <JM'|STEV (J) |MJ> }           )  |\n");
fprintf(out,"|  kq                         kq         M,M'=-J..J    |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"| und F (J)  ganz .                                    |\n");
fprintf(out,"|      kq                                              |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"========================================================\n");
fprintf(out,"\n");
 printf("calculating Fkq , please wait ... \n");
 printf("calculating the F2q ... \n");
 
t01="================================\n";
t02="| ion  |  J  | F20 | F21 | F22 |\n";
t03="|------+-----+-----+-----+-----|\n";
t04="| %4s | %3.1f |  %1d  |  %1d  |  %1d  |\n";
t05="|------+-----+-----+-----+-----|\n";
 
fprintf(out,"%s",t01);fprintf(out,"%s",t02);fprintf(out,"%s",t03);
 
 for( ionennr=0;ionennr< anz_ionen; ++ionennr ){
    dimj    = IONENIMP[ ionennr ].dimj;
    ionname = IONENIMP[ ionennr ].ionname;
 
 
    f20 = fkq(2,0,dimj);
    f21 = fkq(2,1,dimj);
    f22 = fkq(2,2,dimj);
 
    fprintf(out,t04,ionname,(DOUBLE)(dimj-1)/2,f20,f21,f22  );
    if( ionennr < anz_ionen-1 ) fprintf(out,"%s",t05);
    else { fprintf(out,"%s",t01);fprintf(out,"\n"); }
 }
 
 printf("calculating F4q ... \n");
t01="============================================\n";
t02="| ion  |  J  | F40 | F41 | F42 | F43 | F44 |\n";
t03="|------+-----+-----+-----+-----+-----+-----|\n";
t04="| %4s | %3.1f | %3d | %3d | %3d | %3d | %3d |\n";
t05="|------+-----+-----+-----+-----+-----+-----|\n";
 
fprintf(out,"%s",t01);fprintf(out,"%s",t02);fprintf(out,"%s",t03);
 
 for( ionennr=0;ionennr< anz_ionen; ++ionennr ){
    dimj    = IONENIMP[ ionennr ].dimj;
    ionname = IONENIMP[ ionennr ].ionname;
 
    f40 = fkq(4,0,dimj);
    f41 = fkq(4,1,dimj);
    f42 = fkq(4,2,dimj);
    f43 = fkq(4,3,dimj);
    f44 = fkq(4,4,dimj);
 
    fprintf(out,t04,ionname,(DOUBLE)(dimj-1)/2,f40,f41,f42,f43,f44 );
 
 
    if( ionennr < anz_ionen-1 ) fprintf(out,"%s",t05);
    else { fprintf(out,"%s",t01);fprintf(out,"\n"); }
 }
 
 printf("calculating F6q ... \n");
t01="======================================================================\n";
t02="| ion  |  J  |  F60  |  F61  |  F62  |  F63  |  F64  |  F65  |  F66  |\n";
t03="|------+-----+-------+-------+-------+-------+-------+-------+-------|\n";
t04="| %4s | %3.1f |%7d|%7d|%7d|%7d|%7d|%7d|%7d|\n";
t05="|------+-----+-------+-------+-------+-------+-------+-------+-------|\n";
 
fprintf(out,"%s",t01);fprintf(out,"%s",t02);fprintf(out,"%s",t03);
 
 for( ionennr=0;ionennr< anz_ionen; ++ionennr ){
    dimj    = IONENIMP[ ionennr ].dimj;
    ionname = IONENIMP[ ionennr ].ionname;
 
 
    f60 = fkq(6,0,dimj);
    f61 = fkq(6,1,dimj);
    f62 = fkq(6,2,dimj);
    f63 = fkq(6,3,dimj);
    f64 = fkq(6,4,dimj);
    f65 = fkq(6,5,dimj);
    f66 = fkq(6,6,dimj);
 
    fprintf(out,t04,ionname,(DOUBLE)(dimj-1)/2,f60,f61,f62,f63,f64,f65,f66);
    if( ionennr < anz_ionen-1 ) fprintf(out,"%s",t05);
    else { fprintf(out,"%s",t01);fprintf(out,"\n"); }
 }
 
 fclose(out);
}
/*----------------------------------------------------------------------------
                               gkq_tabelle(INT anz_ionen)
------------------------------------------------------------------------------*/
void gkq_tabelle(INT anz_ionen)/* Gkq Faktoren fuer die implementierten Ionen ausgeben */
{
 INT  ionennr;
 INT  dimj;
 CHAR *ionname;
 FILE    *fopen(const char * filename,const char * mode),*out;
 CHAR *filename = "Gkq.tabelle";
 CHAR *t01,*t02,*t03,*t04,*t05;
 LONG   gkq(INT k,INT q,INT dimj);
 LONG   g20,g21,g22;
 LONG   g40,g41,g42,g43,g44;
 LONG   g60,g61,g62,g63,g64,g65,g66;
 LONG   g2m1,g2m2;
 LONG   g4m1,g4m2,g4m3,g4m4;
 LONG   g6m1,g6m2,g6m3,g6m4,g6m5,g6m6;
 
 out = fopen_errchk(filename,"w");
 
 clearscreen;
printf("Information contained in the file %s.\n",filename);
 
fprintf(out,"========================================================\n");
fprintf(out,"| Information about the Factor   G                     |\n");
fprintf(out,"|                                 kq                   |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"| G (J)  :=  GGT(  { <JM'|P (J) |MJ> }           )     |\n");
fprintf(out,"|  kq                      kq         M,M'=-J..J       |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"| and G (J) entirely .                                    |\n");
fprintf(out,"|      kq                                              |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"========================================================\n");
fprintf(out,"\n");
 printf("calculating Gkq , please wait ... \n");
 printf("calculating G2q ... \n");
 
t01="============================================\n";
t02="| ion  |  J  | G2-2| G2-1| G20 | G21 | G22 |\n";
t03="|------+-----+-----+-----+-----+-----+-----|\n";
t04="| %4s | %3.1f |  %1d  |  %1d  |  %1d  |  %1d  |  %1d  |\n";
t05="|------+-----+-----+-----+-----+-----+-----|\n";
 
fprintf(out,"%s",t01);fprintf(out,"%s",t02);fprintf(out,"%s",t03);
 
 for( ionennr=0;ionennr< anz_ionen; ++ionennr ){
    dimj    = IONENIMP[ ionennr ].dimj;
    ionname = IONENIMP[ ionennr ].ionname;
 
 
    g2m2= gkq(2,-2,dimj);
    g2m1= gkq(2,-1,dimj);
    g20 = gkq(2, 0,dimj);
    g21 = gkq(2, 1,dimj);
    g22 = gkq(2, 2,dimj);
 
    fprintf(out,t04,ionname,(DOUBLE)(dimj-1)/2,g2m2,g2m1,g20,g21,g22  );
    if( ionennr < anz_ionen-1 ) fprintf(out,"%s",t05);
    else { fprintf(out,"%s",t01);fprintf(out,"\n"); }
 }
 
 printf("calculating G4q ... \n");
t01="====================================================================\n";
t02="| ion  |  J  | G4-4| G4-3| G4-2| G4-1| G40 | G41 | G42 | G43 | G44 |\n";
t03="|------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----|\n";
t04="| %4s | %3.1f | %3d | %3d | %3d | %3d | %3d | %3d | %3d | %3d | %3d |\n";
t05="|------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----|\n";
 
fprintf(out,"%s",t01);fprintf(out,"%s",t02);fprintf(out,"%s",t03);
 
 for( ionennr=0;ionennr< anz_ionen; ++ionennr ){
    dimj    = IONENIMP[ ionennr ].dimj;
    ionname = IONENIMP[ ionennr ].ionname;
 
    g4m4= gkq(4,-4,dimj);
    g4m3= gkq(4,-3,dimj);
    g4m2= gkq(4,-2,dimj);
    g4m1= gkq(4,-1,dimj);
    g40 = gkq(4, 0,dimj);
    g41 = gkq(4, 1,dimj);
    g42 = gkq(4, 2,dimj);
    g43 = gkq(4, 3,dimj);
    g44 = gkq(4, 4,dimj);
 
    fprintf(out,t04,ionname,(DOUBLE)(dimj-1)/2,g4m4,g4m3,g4m2,g4m1,
                                               g40,g41,g42,g43,g44 );
 
 
    if( ionennr < anz_ionen-1 ) fprintf(out,"%s",t05);
    else { fprintf(out,"%s",t01);fprintf(out,"\n"); }
 }
 
 printf("calculating G6q ... \n");
t01="======================================================================\n";
t02="| ion  |  J  |  G6-6 |  G6-5 |  G6-4 |  G6-3 |  G6-2 |  G6-1 |  G60  |\n";
t03="|------+-----+-------+-------+-------+-------+-------+-------+-------|\n";
t04="| %4s | %3.1f |%7d|%7d|%7d|%7d|%7d|%7d|%7d|\n";
t05="|------+-----+-------+-------+-------+-------+-------+-------+-------|\n";
 
fprintf(out,"%s",t01);fprintf(out,"%s",t02);fprintf(out,"%s",t03);
 
 for( ionennr=0;ionennr< anz_ionen; ++ionennr ){
    dimj    = IONENIMP[ ionennr ].dimj;
    ionname = IONENIMP[ ionennr ].ionname;
 
 
    g6m6= gkq(6,-6,dimj);
    g6m5= gkq(6,-5,dimj);
    g6m4= gkq(6,-4,dimj);
    g6m3= gkq(6,-3,dimj);
    g6m2= gkq(6,-2,dimj);
    g6m1= gkq(6,-1,dimj);
    g60 = gkq(6, 0,dimj);
 
    fprintf(out,t04,ionname,(DOUBLE)(dimj-1)/2,g6m6,g6m5,g6m4,g6m3,g6m2,g6m1,
                                               g60);
    if( ionennr < anz_ionen-1 ) fprintf(out,"%s",t05);
    else { fprintf(out,"%s",t01);fprintf(out,"\n"); }
 }
 
t01="==============================================================\n";
t02="| ion  |  J  |  G61  |  G62  |  G63  |  G64  |  G65  |  G66  |\n";
t03="|------+-----+-------+-------+-------+-------+-------+-------|\n";
t04="| %4s | %3.1f |%7d|%7d|%7d|%7d|%7d|%7d|\n";
t05="|------+-----+-------+-------+-------+-------+-------+-------|\n";
 
fprintf(out,"%s",t01);fprintf(out,"%s",t02);fprintf(out,"%s",t03);
 
 for( ionennr=0;ionennr< anz_ionen; ++ionennr ){
    dimj    = IONENIMP[ ionennr ].dimj;
    ionname = IONENIMP[ ionennr ].ionname;
 
 
    g61 = gkq(6, 1,dimj);
    g62 = gkq(6, 2,dimj);
    g63 = gkq(6, 3,dimj);
    g64 = gkq(6, 4,dimj);
    g65 = gkq(6, 5,dimj);
    g66 = gkq(6, 6,dimj);
 
    fprintf(out,t04,ionname,(DOUBLE)(dimj-1)/2,g61,g62,g63,g64,g65,g66);
    if( ionennr < anz_ionen-1 ) fprintf(out,"%s",t05);
    else { fprintf(out,"%s",t01);fprintf(out,"\n"); }
 }
 fclose(out);
}
/*----------------------------------------------------------------------------
                                  Pkq(INT k,INT q,INT dimj)
------------------------------------------------------------------------------*/
MATRIX *Pkq(INT k,INT q,INT dimj)  /* Matrix <JM'|Pkq(J)|MJ> zurueckgeben */
{
   MATRIX  *pn0n(INT k,INT dimj),*pn1n(INT k,INT dimj),*pn2n(INT k,INT dimj),*pn3n(INT k,INT dimj),*a;
   MATRIX  *pn4n(INT k,INT dimj),*pn5n(INT k,INT dimj),*pn6n(INT k,INT dimj);
 
   switch( ABS( k-ABS(q) ) ){
          case 0 : a=pn0n( q,dimj );
                   break;
          case 1 : a=pn1n( q,dimj );
                   break;
          case 2 : a=pn2n( q,dimj );
                   break;
          case 3 : a=pn3n( q,dimj );
                   break;
          case 4 : a=pn4n( q,dimj );
                   break;
          case 5 : a=pn5n( q,dimj );
                   break;
          case 6 : a=pn6n( q,dimj );
                   break;
          default: fprintf(stderr,"stevkq(): Invalid value of k=%i.\n",k); exit(-1);
   }
   return( a );
}
/*----------------------------------------------------------------------------
                               calc_Pkq(INT dimj)
------------------------------------------------------------------------------*/
STEVENS *calc_Pkq(INT dimj)        /* Stevensoperatoren initialisieren */             /* Gesamtdrehimpuls J, dimj := 2J+1 */
{
    MATRIX  *pn0n(INT k,INT dimj),*pn1n(INT k,INT dimj),*pn2n(INT k,INT dimj);
    MATRIX  *pn3n(INT k,INT dimj),*pn4n(INT k,INT dimj),*pn5n(INT k,INT dimj),*pn6n(INT k,INT dimj);
    STEVENS *stevens;
    MATRIX *Pkq(INT k,INT q,INT dimj);
 
    stevens       = STEVENS_ALLOC(1); /* Speicher fuer die Struktur STEVENS */
    DIMJ(stevens) = dimj;
    P0P0(stevens) = Pkq(0, 0,dimj);      /* Matrix (JN�P0+0(j)�JM) */
 
    P2P2(stevens) = Pkq(2, 2,dimj);      /* Matrix (JN�P2+2(j)�JM) */
    P2P1(stevens) = Pkq(2, 1,dimj);      /* Matrix (JN�P2+1(j)�JM) */
    P2P0(stevens) = Pkq(2, 0,dimj);      /* Matrix (JN�P2+0(j)�JM) */
    P2M1(stevens) = Pkq(2,-1,dimj);      /* Matrix (JN�P2-1(j)�JM) */
    P2M2(stevens) = Pkq(2,-2,dimj);      /* Matrix (JN�P2-2(j)�JM) */
 
    P3P3(stevens) = Pkq(3, 3,dimj);      /* Matrix (JN|P3+3(j)|JM) */
    P3P2(stevens) = Pkq(3, 2,dimj);      /* Matrix (JN|P3+2(j)|JM) */
    P3P1(stevens) = Pkq(3, 1,dimj);      /* Matrix (JN|P3+1(j)|JM) */
    P3P0(stevens) = Pkq(3, 0,dimj);      /* Matrix (JN|P3+0(j)|JM) */
    P3M1(stevens) = Pkq(3,-1,dimj);      /* Matrix (JN|P3-1(j)|JM) */
    P3M2(stevens) = Pkq(3,-2,dimj);      /* Matrix (JN|P3-2(j)|JM) */
    P3M3(stevens) = Pkq(3,-3,dimj);      /* Matrix (JN|P3-3(j)|JM) */

    P4P4(stevens) = Pkq(4, 4,dimj);      /* Matrix (JN�P4+4(j)�JM) */
    P4P3(stevens) = Pkq(4, 3,dimj);      /* Matrix (JN�P4+3(j)�JM) */
    P4P2(stevens) = Pkq(4, 2,dimj);      /* Matrix (JN�P4+2(j)�JM) */
    P4P1(stevens) = Pkq(4, 1,dimj);      /* Matrix (JN�P4+1(j)�JM) */
    P4P0(stevens) = Pkq(4, 0,dimj);      /* Matrix (JN�P4+0(j)�JM) */
    P4M1(stevens) = Pkq(4,-1,dimj);      /* Matrix (JN�P4-1(j)�JM) */
    P4M2(stevens) = Pkq(4,-2,dimj);      /* Matrix (JN�P4-2(j)�JM) */
    P4M3(stevens) = Pkq(4,-3,dimj);      /* Matrix (JN�P4-3(j)�JM) */
    P4M4(stevens) = Pkq(4,-4,dimj);      /* Matrix (JN�P4-4(j)�JM) */

    P5P5(stevens) = Pkq(5, 5,dimj);      /* Matrix (JN|P5+5(j)|JM) */
    P5P4(stevens) = Pkq(5, 4,dimj);      /* Matrix (JN|P5+4(j)|JM) */
    P5P3(stevens) = Pkq(5, 3,dimj);      /* Matrix (JN|P5+3(j)|JM) */
    P5P2(stevens) = Pkq(5, 2,dimj);      /* Matrix (JN|P5+2(j)|JM) */
    P5P1(stevens) = Pkq(5, 1,dimj);      /* Matrix (JN|P5+1(j)|JM) */
    P5P0(stevens) = Pkq(5, 0,dimj);      /* Matrix (JN|P5+0(j)|JM) */
    P5M1(stevens) = Pkq(5,-1,dimj);      /* Matrix (JN|P5-1(j)|JM) */
    P5M2(stevens) = Pkq(5,-2,dimj);      /* Matrix (JN|P5-2(j)|JM) */
    P5M3(stevens) = Pkq(5,-3,dimj);      /* Matrix (JN|P5-3(j)|JM) */
    P5M4(stevens) = Pkq(5,-4,dimj);      /* Matrix (JN|P5-4(j)|JM) */
    P5M5(stevens) = Pkq(5,-5,dimj);      /* Matrix (JN|P5-5(j)|JM) */
 
    P6P6(stevens) = Pkq(6, 6,dimj);      /* Matrix (JN�P6+6(j)�JM) */
    P6P5(stevens) = Pkq(6, 5,dimj);      /* Matrix (JN�P6+5(j)�JM) */
    P6P4(stevens) = Pkq(6, 4,dimj);      /* Matrix (JN�P6+4(j)�JM) */
    P6P3(stevens) = Pkq(6, 3,dimj);      /* Matrix (JN�P6+3(j)�JM) */
    P6P2(stevens) = Pkq(6, 2,dimj);      /* Matrix (JN�P6+2(j)�JM) */
    P6P1(stevens) = Pkq(6, 1,dimj);      /* Matrix (JN�P6+1(j)�JM) */
    P6P0(stevens) = Pkq(6, 0,dimj);      /* Matrix (JN�P6+0(j)�JM) */
    P6M1(stevens) = Pkq(6,-1,dimj);      /* Matrix (JN�P6-1(j)�JM) */
    P6M2(stevens) = Pkq(6,-2,dimj);      /* Matrix (JN�P6-2(j)�JM) */
    P6M3(stevens) = Pkq(6,-3,dimj);      /* Matrix (JN�P6-3(j)�JM) */
    P6M4(stevens) = Pkq(6,-4,dimj);      /* Matrix (JN�P6-4(j)�JM) */
    P6M5(stevens) = Pkq(6,-5,dimj);      /* Matrix (JN�P6-5(j)�JM) */
    P6M6(stevens) = Pkq(6,-6,dimj);      /* Matrix (JN�P6-6(j)�JM) */
 
 
    return(stevens);
}
/*----------------------------------------------------------------------------
                               free_Pkq(STEVENS * stevens )
------------------------------------------------------------------------------*/
void free_Pkq( STEVENS * stevens  )
{
 
    free_mx(P0P0(stevens));
 
    free_mx(P2P2(stevens));
    free_mx(P2P1(stevens));
    free_mx(P2P0(stevens));
    free_mx(P2M1(stevens));
    free_mx(P2M2(stevens));
 
    free_mx(P4P4(stevens));
    free_mx(P4P3(stevens));
    free_mx(P4P2(stevens));
    free_mx(P4P1(stevens));
    free_mx(P4P0(stevens));
    free_mx(P4M1(stevens));
    free_mx(P4M2(stevens));
    free_mx(P4M3(stevens));
    free_mx(P4M4(stevens));
 
    free_mx(P6P6(stevens));
    free_mx(P6P5(stevens));
    free_mx(P6P4(stevens));
    free_mx(P6P3(stevens));
    free_mx(P6P2(stevens));
    free_mx(P6P1(stevens));
    free_mx(P6P0(stevens));
    free_mx(P6M1(stevens));
    free_mx(P6M2(stevens));
    free_mx(P6M3(stevens));
    free_mx(P6M4(stevens));
    free_mx(P6M5(stevens));
    free_mx(P6M6(stevens));
 
 
 
  free_(stevens);
}
/*----------------------------------------------------------------------------
                               free_Okq(STEVENS * stevens )
------------------------------------------------------------------------------*/
void free_Okq(STEVENS * stevens )
{
 
    free_mx(O2P0(stevens));
    free_mx(O2P2(stevens));
    free_mx(O4P2(stevens));
    free_mx(O6P2(stevens));
    free_mx(O6P6(stevens));
    free_mx(O4P05(stevens));
    free_mx(O4M05(stevens));
    free_mx(O6P21(stevens));
    free_mx(O6M21(stevens));
 
    free_(stevens);
}
/*----------------------------------------------------------------------------
                               STEVkq_tabelle()
------------------------------------------------------------------------------*/
/* STEVkq_tabelle()
{
    DOUBLE  j;
    FILE    *fopen(const char * filename,const char * mode),*out;
    INT dimj,k,q;
    CHAR *name="STEVkq.tabelle";
 
    out = fopen_errchk(name,"w");
    clearscreen;
    printf("Inofmation is contained in the file %s .\n",name);
 
    fprintf(out,"========================================================\n");
    fprintf(out,"| Inofrmation about the Stevens operators              |\n");
    fprintf(out,"|                                                      |\n");
    fprintf(out,"| STEV    :=  {  P  ( J )   +  P  ( J )  }/2/theta     |\n");
    fprintf(out,"|     kq          kq  -         k,-q              kq   |\n");
    fprintf(out,"|                                                      |\n");
    fprintf(out,"========================================================\n");
    fprintf(out,"\n");
 
  for( dimj=1 ; dimj<=17 ; ++dimj ){
 
    fprintf(out,"=========================================================\n");
    fprintf(out,"|Matrix element (JM'|STEVkq(J)|MJ) for k=%2d q=%2d J=%5.1f|\n",
                k,q,j);
    fprintf(out,"=========================================================\n");
 
 
  }
    fclose(out);
}*/
/*----------------------------------------------------------------------------
                               info_Pkq(INT k,INT q,INT dimj,CHAR modus)
------------------------------------------------------------------------------*/
/*Matrixelemente (JM'|Pkq(J)|JM) auf file name ausgeben */
void info_Pkq(INT k,INT q,INT dimj,CHAR modus)
{
    CHAR    *name  = "results/Pkq.info";
    MATRIX  *pn0n(INT k,INT dimj),*pn1n(INT k,INT dimj),*pn2n(INT k,INT dimj),*pn3n(INT k,INT dimj),*mx;
    MATRIX  *pn4n(INT k,INT dimj),*pn5n(INT k,INT dimj),*pn6n(INT k,INT dimj),*Pkq(INT k,INT q,INT dimj);
 
    DOUBLE  j;
    FILE    *fopen(const char * filename,const char * mode),*out;
 
    out = fopen_errchk(name,"w");
    clearscreen;
    printf("Inofmation is contained in the file %s.\n",name);
 
    j = ((DOUBLE)dimj-1)/2;
fprintf(out,"========================================================\n");
fprintf(out,"| Information about the generalised                    |\n");
fprintf(out,"| Stevens operator P  ( J )  .                          |\n");
fprintf(out,"|                  kq  -                               |\n");
fprintf(out,"========================================================\n");
fprintf(out,"\n");
 
if(modus==GGT ){
fprintf(out,"\n");
fprintf(out,"========================================================\n");
fprintf(out,"| Matrix elements (JM'|Pkq(J)|MJ) in units of G (J)     |\n");
fprintf(out,"|                                                  kq  |\n");
fprintf(out,"| with                                                 |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"| G (J) := GGT(  (JM'|Pkq(J)|MJ)        ) , G (J)      |\n");
fprintf(out,"|  kq                           M,M'=-J..J   kq        |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"========================================================\n");
fprintf(out,"\n");
}
 
if(modus==NORM ){
fprintf(out,"\n");
fprintf(out,"========================================================\n");
fprintf(out,"| Matrix elements (JM'|Pkq(J)|MJ) in units of  Norm    |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"|                      J                         2     |\n");
fprintf(out,"|               2   -------  || (JM'|Pkq(J)|MJ)||      |\n");
fprintf(out,"| with ||Pkq(J)|| := >        ---------------------     |\n");
fprintf(out,"|                   -------           2*J+1            |\n");
fprintf(out,"|                   M,M'=-J                            |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"|                     2               +                |\n");
fprintf(out,"| and ||(JM'|Pkq|MJ)||  = (JM'|Pkq| MJ) * (JM'|Pkq| MJ)|\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"|                               +                      |\n");
fprintf(out,"|                       = (JM |Pkq|'MJ) * (JM'|Pkq| MJ)|\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"========================================================\n");
fprintf(out,"\n");
}
 
 
    fprintf(out,"========================================================\n");
fprintf(out,"|Matrix element   (JM'|Pkq(J)|MJ) for k=%2d q=%2d J=%5.1f|\n",k
             ,q,j);
    fprintf(out,"========================================================\n");
    fprintf(out,"Only the matrix elements other than zero are ");
    fprintf(out,"given .\n ");
    fprintf(out,"\n");
    fprintf(out,"\n");
 
   mx = Pkq(k,q,dimj);
   switch(modus){
   case GGT : showPkq(out,mx,k,q,dimj,GGT);
              break;
   case NORM: showPkq(out,mx,k,q,dimj,NORM);
              break;
   }
   fprintf(out,"\n");
 
   fclose(out);
}
/*----------------------------------------------------------------------------
                                  showPkq(FILE   *fp,MATRIX *mx,INT k,INT q,INT dimj,CHAR modus)
------------------------------------------------------------------------------*/
void showPkq(FILE   *fp,MATRIX *mx,INT k,INT q,INT dimj,CHAR modus) /* Matrix  (JM'|Pkq(J)|JM)  zeigen     */
              /* in Einheiten des GGT's oder der Norm*/
    /* dimj = 2J+1        als INTeger*/

{
CHAR *text1="(%3.1f %4.1f |P%2d%2d| %4.1f %3.1f) = %1d/%1d * %13.8f\n";
CHAR *text2="(%3.1f %4.1f |P%2d%2d| %4.1f %3.1f) = %1d * %13.8f\n";
CHAR *text3=
"(%3.1f %4.1f |P%2d%2d| %4.1f %3.1f) = %1d/%1d * %1csqrt(%1.0f)\n";
CHAR *text4=
"(%3.1f %4.1f |P%2d%2d| %4.1f %3.1f) = %1d * %1csqrt(%1.0f)\n";
 
CHAR
*text5="(%3.1f %4.1f|P%2d%2d|%4.1f %3.1f) = %13.8f * %13.8f\n";
CHAR *text6=
"(%3.1f %4.1f|P%2d%2d|%4.1f %3.1f)= sqrt(%1d/%1d)*%1csqrt(%1d/%1d)\n";
CHAR *text7=
"(%3.1f %4.1f|P%2d%2d|%4.1f %3.1f)= sqrt(%1d)*%1csqrt(%1d/%1d)\n";
CHAR *text8=
"(%3.1f %4.1f|P%2d%2d|%4.1f %3.1f)= sqrt(%1d/%1d)*%1csqrt(%1d)\n";
CHAR *text9=
"(%3.1f %4.1f|P%2d%2d|%4.1f %3.1f)= sqrt(%1d)*%1csqrt(%1d)\n";
 
CHAR *texts="Alle (%3.1f M'|P%2d%2d| M %3.1f) sind null.\n";
 
CHAR *tnorm="Die Norm || P%1d%1d(%3.1f) || ist null.\n";
 
    INT   m,n,flag=0;
 
    DOUBLE rt,j,nj,mj;
    DOUBLE norm,norm_Stevkq(INT k,INT q,INT dimj);
    BRUCH *fkq,*rt_b,*ggt_matrix( MATRIX *matrix,INT q),*is_rational(DOUBLE zahl);
    LONG  zr,nr,zb,nb;
    LONG  haupt_n;
    CHAR  c;
 
    j = ((DOUBLE)dimj-1)/2; nr = 1; zr = 1; norm = 0.0;
 
    switch(modus){
 
       case GGT :  fkq   = ggt_matrix(mx,q);
                   zr    = ZAEHLER( fkq );
                   nr    = NENNER(  fkq );
                   break;
 
       case NORM:
                   norm  = norm_Pkq(k,q,dimj);
                   if(norm==0.0){
                       fprintf(fp,texts,j,k,q,j);
                       fprintf(fp,"\n");
                       fprintf(fp,tnorm,k,q,j);
                       fclose(fp);
                       return;
                   }
                   if( k<6 ){
                      fkq   = is_rational( norm*norm );
                      zr    = ZAEHLER( fkq );
                      nr    = NENNER(  fkq ); /* norm = sqrt( zr/nr ) */
                      break;
                   }
                   break;
    }
 
    flag=0;
    for( n=dimj ; n>=1 ; --n )
       for( m=dimj ; m>=1 ; --m ){
 
         switch(modus){
            case GGT :
                        rt = R(mx,n,m)*(DOUBLE)nr/(DOUBLE)zr;
                        mj = ((DOUBLE)m-j-1);
                        nj = ((DOUBLE)n-j-1);
                        if( rt!= 0.0 ){
                           ++flag;
                           if( nr!=1 )
                               fprintf(fp,text1,j,nj,k,q,mj,j,zr,nr,rt);
                           else
                               fprintf(fp,text2,j,nj,k,q,mj,j,zr,rt);
                        }
                        break;
 
            case NORM:
                        rt = R(mx,n,m)/norm;
                        if( rt!= 0.0 ){
                            ++flag;
                            mj = ((DOUBLE)m-j-1);
                            nj = ((DOUBLE)n-j-1);
                            fprintf(fp,text5,j,nj,k,q,mj,j,norm,rt);
                        }
                        break;
         }
       }
    if( flag==0 )  fprintf(fp,texts,j,k,q,j);
    fprintf(fp,"\n");
 
 
 haupt_n = 1;
 if( modus==NORM && k<6 ){  /* Hauptnenner eines Matrixelementes */
    for( n=dimj ; n>=1 ; --n )
       for( m=dimj ; m>=1 ; --m ){
             rt = R(mx,n,m)*sqrt((DOUBLE)nr)/sqrt((DOUBLE)zr);
                  if( rt!= 0.0 ){
                        rt_b = is_rational( rt*rt );
                        nb = NENNER(  rt_b );
                        if( nb>haupt_n ) haupt_n = nb;
 
         }
       }
 }
 
 
 
 if( (modus==GGT && q!=0) || (modus==NORM && k<6 ) ){
    flag=0;
    for( n=dimj ; n>=1 ; --n )
       for( m=dimj ; m>=1 ; --m ){
 
         switch(modus){
            case GGT :
                        rt = R(mx,n,m)*(DOUBLE)nr/(DOUBLE)zr;
                        c=' ';if(rt<0) c='-';
                        mj = ((DOUBLE)m-j-1);
                        nj = ((DOUBLE)n-j-1);
                        if( rt!= 0.0 ){
                         ++flag;
                         if( nr!=1 )
                           fprintf(fp,text3,j,nj,k,q,mj,j,zr,nr,c,rt*rt);
                          else
                           fprintf(fp,text4,j,nj,k,q,mj,j,zr,c,rt*rt);
                        }
                        break;
 
            case NORM:
                        rt = R(mx,n,m)*sqrt((DOUBLE)nr)/sqrt((DOUBLE)zr);
                        c=' ';if(rt<0) c='-';
                        if( rt!= 0.0 ){
                            ++flag;
                            mj = ((DOUBLE)m-j-1);
                            nj = ((DOUBLE)n-j-1);
                            rt_b = is_rational( rt*rt );
                            zb = ZAEHLER( rt_b );
                            nb = NENNER(  rt_b );
                            zb *= (haupt_n/nb);
                            nb = haupt_n;
 
                          if( nr==1 && nb==1 )
                           fprintf(fp,text9,j,nj,k,q,mj,j,zr,c,zb);
                          if( nr==1 && nb!=1 )
                           fprintf(fp,text7,j,nj,k,q,mj,j,zr,c,zb,nb);
                          if( nr!=1 && nb==1 )
                           fprintf(fp,text8,j,nj,k,q,mj,j,zr,nr,c,zb);
                          if( nr!=1 && nb!=1 )
                           fprintf(fp,text6,j,nj,k,q,mj,j,zr,nr,c,zb,nb);
                        }
                        break;
         }
       }
    if( flag==0 )  fprintf(fp,texts,j,k,q,j);
    fprintf(fp,"\n");
 
 }
 
}
/*----------------------------------------------------------------------------
                               info_STEVkq(INT k,INT q,INT dimj,CHAR modus)
------------------------------------------------------------------------------*/
/*Matrixelemente (JM'|STEVkq(J)|JM) auf file name ausgeben */
void info_STEVkq(INT k,INT q,INT dimj,CHAR modus)
{
    CHAR    *name  = "results/STEVkq.info";
 
    DOUBLE  j;
    FILE    *fopen(const char * filename,const char * mode),*out;
    MATRIX  *stevkq(INT k,INT q,INT dimj),*mx;
 
    out = fopen_errchk(name,"w");
    clearscreen;
    printf("Information is contained in the file %s.\n",name);
 
    j = ((DOUBLE)dimj-1)/2;
 
fprintf(out,"==============================================================\n");
fprintf(out,"| Information about the  Stevens operators                   |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"| STEV     := {  P    (J)  +  P      (J)  }/2/omega     :q>=0|\n");
fprintf(out,"|     kq          k|q| -       k,-|q| -            k|q|      |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"| STEV  *i := {  P    (J)  -  P      (J)  }/2/omega     :q<0 |\n");
fprintf(out,"|     kq          k|q| -       k,-|q| -            k|q|      |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"==============================================================\n");
 
if(modus==GGT ){
fprintf(out,"\n");
fprintf(out,"==============================================================\n");
fprintf(out,"| Matrix elements (JM'|STEVkq(J)|MJ) bzw (JM'|i*STEVkq(J)|MJ) |\n");
fprintf(out,"| in units of F  (J)                                         |\n");
fprintf(out,"|                   kq                                       |\n");
fprintf(out,"| q >= 0 :                                                   |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"| F (J) := GGT(  (JM'|STEVkq(J)|MJ)        ) , F (J) .       |\n");
fprintf(out,"|  kq                              M,M'=-J..J   kq           |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"| q <  0 :                                                   |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"| F (J) := GGT((JM'|i*STEVkq(J)|MJ)        ) , F (J) .       |\n");
fprintf(out,"|  kq                              M,M'=-J..J   kq           |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"==============================================================\n");
fprintf(out,"\n");
}
if(modus==NORM ){
fprintf(out,"\n");
fprintf(out,"==============================================================\n");
fprintf(out,"| Matrix elements (JM'| A |MJ) in units of the Norm          |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"|                    J                   +   +               |\n");
fprintf(out,"|            2    -------  (JM | 0.5*(A*A + A*A ) |JM)       |\n");
fprintf(out,"| with || A ||  := >        ---------------------------       |\n");
fprintf(out,"|                 -------             2*J+1                  |\n");
fprintf(out,"|                  M = -J                                    |\n");
fprintf(out,"|                                                            |\n");
fprintf(out,"==============================================================\n");
fprintf(out,"\n");
}
fprintf(out,"\n");
fprintf(out,"=========================================================\n");
if( q>=0 )
fprintf(out,"|Matrix element (JM'|STEVkq(J)|MJ) for k=%2d q=%2d J=%5.1f|\n",k
             ,q,j);
else
fprintf(out,"|Matrix element (JM'|i*STEVkq |MJ) for k=%2d q=%2d J=%5.1f|\n",k
             ,q,j);
    fprintf(out,"=========================================================\n");
    fprintf(out,"Only the matrix elements other than zero are ");
    fprintf(out,"given .\n ");
    fprintf(out,"\n");
    fprintf(out,"\n");
 
    mx = stevkq(k,q,dimj);
    switch(modus){
    case GGT : showSTEVkq(out,mx,k,q,dimj,GGT);
               break;
    case NORM: showSTEVkq(out,mx,k,q,dimj,NORM);
               break;
    }
    fprintf(out,"\n");
 
    fclose(out);
}
/*----------------------------------------------------------------------------
                                  stevkq(INT k,INT q,INT dimj)
------------------------------------------------------------------------------*/
MATRIX *stevkq(INT k,INT q,INT dimj)  /* Matrix <JM'|STEVkq(J)|MJ> zurueckgeben */
{
   MATRIX  *pn0n(INT k,INT dimj),*pn1n(INT k,INT dimj),*pn2n(INT k,INT dimj),*pn3n(INT k,INT dimj),*a,*b,*c;
   MATRIX  *pn4n(INT k,INT dimj),*pn5n(INT k,INT dimj),*pn6n(INT k,INT dimj);
   MATRIX  *mx_add(DOUBLE z, MATRIX*a,MATRIX*b),*mx_sub(DOUBLE z, MATRIX*a,MATRIX*b);
   DOUBLE  z;
   switch( ABS( k-ABS(q) ) ){
          case 0 : a=pn0n( ABS(q),dimj);
                   b=pn0n(-ABS(q),dimj);
                   z=0.5/omegan0n(ABS(q));      /*         +   */
                   if( q >=0 ) c=mx_add(z,a,b); /* O  =(P+P )/2*/
                   else        c=mx_sub(z,a,b); /*         +   */
                   break;                       /* O*i=(P-P )/2*/
          case 1 : a=pn1n( ABS(q),dimj);
                   b=pn1n(-ABS(q),dimj);
                   z=0.5/omegan1n(ABS(q));
                   if( q >=0 ) c=mx_add(z,a,b);
                   else        c=mx_sub(z,a,b);
                   break;
          case 2 : a=pn2n( ABS(q),dimj);
                   b=pn2n(-ABS(q),dimj);
                   z=0.5/omegan2n(ABS(q));
                   if( q >=0 ) c=mx_add(z,a,b);
                   else        c=mx_sub(z,a,b);
                   break;
          case 3 : a=pn3n( ABS(q),dimj);
                   b=pn3n(-ABS(q),dimj);
                   z=0.5/omegan3n(ABS(q));
                   if( q >=0 ) c=mx_add(z,a,b);
                   else        c=mx_sub(z,a,b);
                   break;
          case 4 : a=pn4n( ABS(q),dimj);
                   b=pn4n(-ABS(q),dimj);
                   z=0.5/omegan4n(ABS(q));
                   if( q >=0 ) c=mx_add(z,a,b);
                   else        c=mx_sub(z,a,b);
                   break;
          case 5 : a=pn5n( ABS(q),dimj);
                   b=pn5n(-ABS(q),dimj);
                   z=0.5/omegan5n(ABS(q));
                   if( q >=0 ) c=mx_add(z,a,b);
                   else        c=mx_sub(z,a,b);
                   break;
          case 6 : a=pn6n( ABS(q),dimj);
                   b=pn6n(-ABS(q),dimj);
                   z=0.5/omegan6n(ABS(q));
                   if( q >=0 ) c=mx_add(z,a,b);
                   else        c=mx_sub(z,a,b);
                   break;
          default: fprintf(stderr,"stevkq(): Invalid value of k=%i.\n",k); exit(-1);
   }
 
   free_mx(a);
   free_mx(b);
   return( c );
}
/*----------------------------------------------------------------------------
                                  stevks(INT k,INT s,INT dimj,ITERATION *iter)
------------------------------------------------------------------------------*/
/* Matrix <JM'|STEVk0(J)+s*STEVk4(J)|MJ>zurueckgeben */
MATRIX *stevks(INT k,INT s,INT dimj,ITERATION *iter)
{
   MATRIX  *mx_addf(DOUBLE f1,DOUBLE f2,MATRIX*m1,MATRIX*m2),*stevkq(INT k,INT q,INT dimj);
   MATRIX  *o40,*o44,*o60,*o64,*o;
   DOUBLE  ss;
 
   ss = (DOUBLE)s;
   switch( k ){
          case 4 :
                   o40  = stevkq(4,0,dimj);
                   if( VOR44(iter)==1.0 ) o44  = stevkq(4, 4,dimj);
                   else                   o44  = stevkq(4,-4,dimj);
                   o    = mx_addf(1.0,ss,o40,o44);
                   free_mx( o40);
                   free_mx( o44);
                   break;
 
          case 6 :
                   o60  = stevkq(6,0,dimj);
                   if( VOR64(iter)==1.0 ) o64  = stevkq(6, 4,dimj);
                   else                   o64  = stevkq(6,-4,dimj);
                   o    = mx_addf(1.0,ss,o60,o64);
                   free_mx( o60);
                   free_mx( o64);
                   break;
          default: fprintf(stderr,"stevkq(): Invalid value of k=%i.\n",k); exit(-1);
   }
 
   return( o );
}
/*----------------------------------------------------------------------------
                                  showSTEVkq(FILE *fp, MATRIX *mx, INT k, INT q, INT dimj, CHAR modus)
------------------------------------------------------------------------------*/
void showSTEVkq(FILE *fp, MATRIX *mx, INT k, INT q, INT dimj, CHAR modus) /* Matrix  (JM'|STEVkq(J)|JM)  zeigen   */
                  /* in Einheiten des GGT's oder der Norm */
  /* dimj = 2J+1        als INTeger*/
{
CHAR *text1="(%3.1f %4.1f |STEV%2d%2d| %4.1f %3.1f) = %1d/%1d * %13.8f\n";
CHAR *text2="(%3.1f %4.1f |STEV%2d%2d| %4.1f %3.1f) = %1d * %13.8f\n";
CHAR *text3=
"(%3.1f %4.1f |STEV%2d%2d| %4.1f %3.1f) = %1d/%1d * %1csqrt(%1.0f)\n";
CHAR *text4=
"(%3.1f %4.1f |STEV%2d%2d| %4.1f %3.1f) = %1d * %1csqrt(%1.0f)\n";
 
CHAR
*text5="(%3.1f %4.1f|STEV%2d%2d|%4.1f %3.1f) = %13.8f * %13.8f\n";
CHAR *text6=
"(%3.1f %4.1f|STEV%2d%2d|%4.1f %3.1f)= sqrt(%1d/%1d)*%1csqrt(%1d/%1d)\n";
CHAR *text7=
"(%3.1f %4.1f|STEV%2d%2d|%4.1f %3.1f)= sqrt(%1d)*%1csqrt(%1d/%1d)\n";
CHAR *text8=
"(%3.1f %4.1f|STEV%2d%2d|%4.1f %3.1f)= sqrt(%1d/%1d)*%1csqrt(%1d)\n";
CHAR *text9=
"(%3.1f %4.1f|STEV%2d%2d|%4.1f %3.1f)= sqrt(%1d)*%1csqrt(%1d)\n";
 
CHAR *texts="Alle (%3.1f M'|STEV%2d%2d| M %3.1f) sind null.\n";
 
CHAR *tnorm="Die Norm || STEV%1d%1d(%3.1f) || ist null.\n";
 
    INT   m,n,flag=0;
    DOUBLE rt,j,nj,mj;
    DOUBLE norm,norm_Stevkq(INT k,INT q,INT dimj);
    BRUCH *fkq,*rt_b,*ggt_matrix( MATRIX *matrix,INT q),*is_rational(DOUBLE zahl);
    LONG  zr,nr,zb,nb;
    LONG  haupt_n;
    CHAR  c;
 
    j = ((DOUBLE)dimj-1)/2; nr=1; zr=1; norm=0.0;
 
    switch(modus){
 
       case GGT :  fkq   = ggt_matrix(mx,q);
                   zr    = ZAEHLER( fkq );
                   nr    = NENNER(  fkq );
                   break;
 
       case NORM:  norm  = norm_Stevkq(k,q,dimj);
                   if(norm==0.0){
                        fprintf(fp,texts,j,k,q,j);
                        fprintf(fp,"\n");
                        fprintf(fp,tnorm,k,q,j);
                        fclose(fp);
                        return;
                   }
                   if( k<6 ){
                      fkq   = is_rational( norm*norm );
                      zr    = ZAEHLER( fkq );
                      nr    = NENNER(  fkq ); /* norm = sqrt( zr/nr ) */
                      break;
                   }
    }
 
    flag=0;
    for( n=dimj ; n>=1 ; --n )
       for( m=dimj ; m>=1 ; --m ){
 
         switch(modus){
            case GGT :
                        rt = R(mx,n,m)*(DOUBLE)nr/(DOUBLE)zr;
                        mj = ((DOUBLE)m-j-1);
                        nj = ((DOUBLE)n-j-1);
                        if( rt!= 0.0 ){
                           ++flag;
                           if( nr!=1 )
                               fprintf(fp,text1,j,nj,k,q,mj,j,zr,nr,rt);
                           else
                               fprintf(fp,text2,j,nj,k,q,mj,j,zr,rt);
                        }
                        break;
 
            case NORM:
                        rt = R(mx,n,m)/norm;
                        if( rt!= 0.0 ){
                            ++flag;
                            mj = ((DOUBLE)m-j-1);
                            nj = ((DOUBLE)n-j-1);
                            fprintf(fp,text5,j,nj,k,q,mj,j,norm,rt);
                        }
                        break;
         }
       }
    if( flag==0 )  fprintf(fp,texts,j,k,q,j);
    fprintf(fp,"\n");
 
 
 
 haupt_n = 1;
 if( modus==NORM && k<6 ){  /* Hauptnenner eines Matrixelementes */
    for( n=dimj ; n>=1 ; --n )
       for( m=dimj ; m>=1 ; --m ){
             rt = R(mx,n,m)*sqrt((DOUBLE)nr)/sqrt((DOUBLE)zr);
                  if( rt!= 0.0 ){
                        rt_b = is_rational( rt*rt );
                        nb = NENNER(  rt_b );
                        if( nb>haupt_n ) haupt_n = nb;
 
         }
       }
 }
 
 
 
 
 if( (modus==GGT && q!=0) || (modus==NORM && k<6 ) ){
    flag=0;
    for( n=dimj ; n>=1 ; --n )
       for( m=dimj ; m>=1 ; --m ){
 
         switch(modus){
            case GGT :
                        rt = R(mx,n,m)*(DOUBLE)nr/(DOUBLE)zr;
                        c=' ';if(rt<0) c='-';
                        mj = ((DOUBLE)m-j-1);
                        nj = ((DOUBLE)n-j-1);
                        if( rt!= 0.0 ){
                         ++flag;
                         if( nr!=1 )
                           fprintf(fp,text3,j,nj,k,q,mj,j,zr,nr,c,rt*rt);
                          else
                           fprintf(fp,text4,j,nj,k,q,mj,j,zr,c,rt*rt);
                        }
                        break;
 
            case NORM:
                        rt = R(mx,n,m)*sqrt((DOUBLE)nr)/sqrt((DOUBLE)zr);
                        c=' ';if(rt<0) c='-';
                        if( rt!= 0.0 ){
                            ++flag;
                            mj = ((DOUBLE)m-j-1);
                            nj = ((DOUBLE)n-j-1);
                            rt_b = is_rational( rt*rt );
                            zb = ZAEHLER( rt_b );
                            nb = NENNER(  rt_b );
                            zb *= (haupt_n/nb);
                            nb = haupt_n;
                          if( nr==1 && nb==1 )
                           fprintf(fp,text9,j,nj,k,q,mj,j,zr,c,zb);
                          if( nr==1 && nb!=1 )
                           fprintf(fp,text7,j,nj,k,q,mj,j,zr,c,zb,nb);
                          if( nr!=1 && nb==1 )
                           fprintf(fp,text8,j,nj,k,q,mj,j,zr,nr,c,zb);
                          if( nr!=1 && nb!=1 )
                           fprintf(fp,text6,j,nj,k,q,mj,j,zr,nr,c,zb,nb);
                        }
                        break;
         }
       }
    if( flag==0 )  fprintf(fp,texts,j,k,q,j);
    fprintf(fp,"\n");
 
 }
 
}
/*----------------------------------------------------------------------------
                               info_Fkq(INT k,INT q,INT dimj)
------------------------------------------------------------------------------*/
/*
 
   F (J) :=  GGT(  { (JM'| STEV (J) |MJ) }              )
    kq                         kq         M,M'=-J...J
 
   mit GGT ganz ;
 
 
 
   beispiel :   ggt( 2/3*sqrt(2) ,4/3*sqrt(2) ) = 2/3
 
 
 
*/
void info_Fkq(INT k,INT q,INT dimj)
{
    CHAR    *name  = "results/Fkq.info";
    LONG    fkq(INT k,INT q,INT dimj),_fkq;
/*  DOUBLE  z; */
    DOUBLE  j;
/*  INT     n,m; */
    FILE    *fopen(const char * filename,const char * mode),*out;
 
    out = fopen_errchk(name,"w");
    clearscreen;
    printf("Information is contained in the file %s .\n",name);
 
 
fprintf(out,"========================================================\n");
fprintf(out,"| Information about the  Faktors F                     |\n");
fprintf(out,"|                                 kq                   |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"| F (J)  :=  GGT(  { <JM'|STEV (J) |MJ> }           )  |\n");
fprintf(out,"|  kq                         kq         M,M'=-J..J    |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"| and F (J)  entirely .                                |\n");
fprintf(out,"|      kq                                              |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"========================================================\n");
fprintf(out,"\n");
 
_fkq = fkq(k,q,dimj);
 
 
j = ((DOUBLE)dimj-1.0)/2;
fprintf(out,"========================================================\n");
    fprintf(out,"| F%2d%2d(%4.1f) = %9d                              |\n",
                  k,q,j, (int)_fkq  );
fprintf(out,"========================================================\n");
fprintf(out,"\n");
 
fclose(out);
 
 
}
/*----------------------------------------------------------------------------
                               info_Gkq(INT k,INT q,INT dimj)
------------------------------------------------------------------------------*/
/*
 
   G (J) :=  GGT(  { (JM'| P (J) |MJ) }              )
    kq                      kq         M,M'=-J...J
 
   mit GGT ganz ;
 
 
 
   beispiel :   ggt( 2/3*sqrt(2) ,4/3*sqrt(2) ) = 2/3
 
 
 
*/
void info_Gkq(INT k,INT q,INT dimj)
{
    CHAR    *name  = "results/Gkq.info";
    LONG    gkq(INT k,INT q,INT dimj),_gkq;
/*  DOUBLE  z; */
    DOUBLE  j;
/*  INT     n,m; */
    FILE    *fopen(const char * filename,const char * mode),*out;
 
    out = fopen_errchk(name,"w");
    clearscreen;
    printf("Information is contained in the file %s.\n",name);
 
 
fprintf(out,"========================================================\n");
fprintf(out,"| Information about the Faktors  G                     |\n");
fprintf(out,"|                                 kq                   |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"| G (J)  :=  GGT(  { <JM'|P (J) |MJ> }           )     |\n");
fprintf(out,"|  kq                      kq         M,M'=-J..J       |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"| and G (J)  entirely .                                    |\n");
fprintf(out,"|      kq                                              |\n");
fprintf(out,"|                                                      |\n");
fprintf(out,"========================================================\n");
fprintf(out,"\n");
 
_gkq = gkq(k,q,dimj);
 
 
j = ((DOUBLE)dimj-1.0)/2;
fprintf(out,"========================================================\n");
    fprintf(out,"| F%2d%2d(%4.1f) = %9d                              |\n",
                  k,q,j, (int)_gkq  );
fprintf(out,"========================================================\n");
fprintf(out,"\n");
 
fclose(out);
 
 
}
/*----------------------------------------------------------------------------
                             fkq(INT k,INT q,INT dimj)
------------------------------------------------------------------------------*/
LONG fkq(INT k,INT q,INT dimj)  /* Fkq fuer stevensoperatoren berechnen */
{
    MATRIX *stevkq(INT k,INT q,INT dimj), *mx;
    BRUCH  *ggt_matrix( MATRIX *matrix,INT q),*fkq;
 
    mx  = stevkq(k,q,dimj);
    fkq = ggt_matrix(mx,q);
    free_mx(mx);
 
    if( NENNER(fkq) != 1) {printf("denom!=1 in STEVENS.C in fkq()\n");exit(0);}
 
    return( ZAEHLER(fkq) );
}
/*----------------------------------------------------------------------------
                             gkq(INT k,INT q,INT dimj)
------------------------------------------------------------------------------*/
LONG gkq(INT k,INT q,INT dimj)  /* Gkq fuer stevensoperatoren berechnen */
{
    MATRIX *Pkq(INT k,INT q,INT dimj), *mx;
    BRUCH  *ggt_matrix( MATRIX *matrix,INT q),*gkq;
 
    mx  = Pkq(k,q,dimj);
    gkq = ggt_matrix(mx,q);
    free_mx(mx);
 
    if( NENNER(gkq) != 1) {printf("denom!=1 in STEVENS.C in gkq(INT k,INT q,INT dimj)\n");exit(0);}
 
    return( ZAEHLER(gkq) );
}
/*----------------------------------------------------------------------------
                             ggt_matrix( MATRIX *matrix,INT q)
------------------------------------------------------------------------------*/
BRUCH *ggt_matrix( MATRIX *matrix,INT q)
{
    BRUCH   *ggt_r(DOUBLE a,LONG z,LONG n),*bruch,*dummy,*suche_quadratzahl(BRUCH   *bruch);
    LONG    ggt_l(LONG a,LONG b);
    LONG    ggtl;
    INT     n,m,dimj;
    DOUBLE  rt2;
 
 
    bruch = BRUCH_ALLOC(1);
    ZAEHLER(bruch) = 0;
    NENNER( bruch) = 0;
 
    dimj = MXDIM(matrix);
    for( n=dimj ; n>=1 ; --n )
         for( m=dimj ; m>=1 ; --m ){
              if(q!=0 ) rt2 = R(matrix,n,m)*R(matrix,n,m);
              else      rt2 = R(matrix,n,m);
              if( rt2 != 0.0 ){
                 dummy = ggt_r( rt2 , ZAEHLER(bruch),NENNER(bruch) );
                 free_(bruch);
                 bruch = dummy;
              }
         }
 
    ZAEHLER(bruch)  = ABS( ZAEHLER(bruch) );
    NENNER( bruch)  = ABS( NENNER( bruch) );
    ggtl            = ggt_l( ZAEHLER(bruch) ,NENNER(bruch) );
    ZAEHLER(bruch) /= ggtl;
    NENNER( bruch) /= ggtl;
 
 
   if(q!=0)     bruch = suche_quadratzahl( bruch );
 
    if( ZAEHLER(bruch)==0 && NENNER(bruch)==0 ){
        ZAEHLER(bruch) = 1;
        NENNER(bruch)  = 1;
    }
 
    return( bruch );
}
/*----------------------------------------------------------------------------
                             pn0n(INT k,INT dimj)
------------------------------------------------------------------------------*/
MATRIX *pn0n(INT k,INT dimj)       /* MATRIX  (JN|Pk,k(J)|JM) */
         /* Gesamtdrehimpuls J , dimj := 2J+1 */
{
    INT    m,n;
    DOUBLE jzn,jzm,jpm;
    DOUBLE fac(INT n);
    MATRIX *okq,*mx_alloc(INT anz_ze,INT anz_sp);
 
    #include "define_j.c"      /* mj,J2,J4,... definieren */
    okq = mx_alloc(dimj,dimj); /* Speicher fuer Matrix (JN|Pk,k(J)|JM) */
 
    for( n=dimj ; n>=1 ; --n )
         for( m=dimj ; m>=1 ; --m ){
 
         jzn = 1;
         jzm = 1;
         jpm = JPM(k,mj)*D(nj,mj+k);
 
         R(okq,n,m) = jpm*(jzn + jzm )/2.0;
    }
 
    return( okq );
}
/*----------------------------------------------------------------------------
                             pn1n(INT k,INT dimj)
------------------------------------------------------------------------------*/
MATRIX *pn1n(INT k,INT dimj)       /* MATRIX  (JN|Pk+1,k(J)|JM) */
           /* Gesamtdrehimpuls J , dimj := 2J+1 */
{
    INT    m,n;
    DOUBLE jzn,jzm,jpm;
    DOUBLE fac(INT n);
    MATRIX *okq,*mx_alloc(INT anz_ze,INT anz_sp);
 
    #include "define_j.c"      /* mj,J2,J4,... definieren */
    okq = mx_alloc(dimj,dimj); /* Speicher fuer Matrix(JN|Pk+1,k(J)|JM)*/
 
    for( n=dimj ; n>=1 ; --n )
         for( m=dimj ; m>=1 ; --m ){
 
         jzn = nj;
         jzm = mj;
         jpm = JPM(k,mj)*D(nj,mj+k);
 
         R(okq,n,m) = jpm*(jzn + jzm )/2.0;
    }
 
    return( okq );
}
/*----------------------------------------------------------------------------
                             pn2n(INT k,INT dimj)
------------------------------------------------------------------------------*/
MATRIX *pn2n(INT k,INT dimj)       /* MATRIX  (JN|Pk+2,k(J)|JM) */
           /* Gesamtdrehimpuls J , dimj := 2J+1 */
{
    INT    m,n;
    INT    s;
    DOUBLE jzn,jzm/*,jpn*/,jpm;
    DOUBLE fac(INT n);
    MATRIX *okq,*mx_alloc(INT anz_ze,INT anz_sp);
 
    #include "define_j.c"      /* mj,J2,J4,... definieren */
    okq = mx_alloc(dimj,dimj); /* Speicher fuer Matrix(JN|Pk+2,k(J)|JM)*/
 
    s = ABS(k);
    for( n=dimj ; n>=1 ; --n )
         for( m=dimj ; m>=1 ; --m ){
 
         jzn = (2*s+3)*nj2 - J2 + s*P2(s+1)/2.0 -(2*s+3)*P2(s)/2.0;
         jzm = (2*s+3)*mj2 - J2 + s*P2(s+1)/2.0 -(2*s+3)*P2(s)/2.0;
         jpm = JPM(k,mj)*D(nj,mj+k);
 
         R(okq,n,m) = jpm*(jzn + jzm )/2.0;
    }
 
    return( okq );
}
/*----------------------------------------------------------------------------
                             pn3n(INT k,INT dimj)
------------------------------------------------------------------------------*/
MATRIX *pn3n(INT k,INT dimj)       /* MATRIX  (JN|Pk+3,k(J)|JM) */
          /* Gesamtdrehimpuls J , dimj := 2J+1 */
{
    INT    m,n;
    INT    s;
    DOUBLE jzn,jzm/*,jpn*/,jpm;
    DOUBLE fac(INT n);
    MATRIX *okq,*mx_alloc(INT anz_ze,INT anz_sp);
 
    #include "define_j.c"      /* mj,J2,J4,... definieren */
    okq = mx_alloc(dimj,dimj); /* Speicher fuer Matrix(JN|Pk+3,k(J)|JM)*/
 
    s = ABS(k);
    for( n=dimj ; n>=1 ; --n )
         for( m=dimj ; m>=1 ; --m ){
 
         jzn =( (2*s+5)*nj2
                -3*J2
                +P2(s+1)*(s+2)/2.0
                -(2*s+5)*P2(s)     )*nj;
 
         jzm =( (2*s+5)*mj2
                -3*J2
                +P2(s+1)*(s+2)/2.0
                -(2*s+5)*P2(s)     )*mj;
 
         jpm = JPM(k,mj)*D(nj,mj+k);
 
         R(okq,n,m) = jpm*(jzn + jzm )/2.0;
    }
 
    return( okq );
}
/*----------------------------------------------------------------------------
                             pn4n(INT k,INT dimj)
------------------------------------------------------------------------------*/
MATRIX *pn4n(INT k,INT dimj)       /* MATRIX  (JN|Pk+4,k(J)|JM) */
          /* Gesamtdrehimpuls J , dimj := 2J+1 */
{
    INT    m,n;
    INT    s;
    DOUBLE jzn,jzm/*,jpn*/,jpm;
    DOUBLE fac(INT n);
    MATRIX *okq,*mx_alloc(INT anz_ze,INT anz_sp);
 
    #include "define_j.c"      /* mj,J2,J4,... definieren */
    okq = mx_alloc(dimj,dimj); /* Speicher fuer Matrix(JN|Pk+4,k(J)|JM)*/
 
    s = ABS(k);
    for( n=dimj ; n>=1 ; --n )
         for( m=dimj ; m>=1 ; --m ){
 
         jzn = (2*s+5)*(2*s+7)*nj4
              -(2*s+5)*(  6*J2
                          +(s-1)*(3*P2(s)+12*s+5)
                       )*nj2
              -3*J2*( -P3(s) -2*P2(s) +2*s + 2 )
              +3*J4
              +s/4.0*( 5*P5(s) +27*P4(s) +23*P3(s) -39*P2(s)
                       -10*s +12 ) ;
 
         jzm = (2*s+5)*(2*s+7)*mj4
              -(2*s+5)*(  6*J2
                          +(s-1)*(3*P2(s)+12*s+5)
                       )*mj2
              -3*J2*( -P3(s) -2*P2(s) +2*s + 2 )
              +3*J4
              +s/4.0*( 5*P5(s) +27*P4(s) +23*P3(s) -39*P2(s)
                       -10*s +12 ) ;
 
         jpm = JPM(k,mj)*D(nj,mj+k);
 
         R(okq,n,m) = jpm*(jzn + jzm )/2.0;
    }
 
    return( okq );
}
/*----------------------------------------------------------------------------
                             pn5n(INT k,INT dimj)
------------------------------------------------------------------------------*/
MATRIX *pn5n(INT k,INT dimj)       /* MATRIX  (JN|Pk+5,k(J)|JM) */
          /* Gesamtdrehimpuls J , dimj := 2J+1 */
{
    INT    m,n;
    INT    s;
    DOUBLE jzn,jzm/*,jpn*/,jpm;
    DOUBLE fac(INT n);
    MATRIX *okq,*mx_alloc(INT anz_ze,INT anz_sp);
 
    #include "define_j.c"      /* mj,J2,J4,... definieren */
    okq = mx_alloc(dimj,dimj); /* Speicher fuer Matrix(JN|Pk+5,k(J)|JM)*/
 
    s = ABS(k);
    for( n=dimj ; n>=1 ; --n )
         for( m=dimj ; m>=1 ; --m ){
 
         jzn =(  (2*s+7)*(2*s+9)*nj4
                -(2*s+7)*(  10*J2
                            +5*P3(s)+20*P2(s)-15*s-15
                         )*nj2
                -5*J2*( -3*P3(s) -9*P2(s) +8*s + 10)
                +15*J4
                +1/4.0*( 25*P6(s)+185*P5(s)+255*P4(s)-365*P3(s)
                        -176*P2(s) +172*s +48 )
              )*nj;
 
         jzm =(  (2*s+7)*(2*s+9)*mj4
                -(2*s+7)*(  10*J2
                            +5*P3(s)+20*P2(s)-15*s-15
                         )*mj2
                -5*J2*( -3*P3(s) -9*P2(s) +8*s + 10)
                +15*J4
                +1/4.0*( 25*P6(s)+185*P5(s)+255*P4(s)-365*P3(s)
                        -176*P2(s) +172*s +48 )
              )*mj;
 
         jpm = JPM(k,mj)*D(nj,mj+k);
 
         R(okq,n,m) = jpm*(jzn + jzm )/2.0;
    }
 
    return( okq );
}
/*----------------------------------------------------------------------------
                             pn6n(INT k,INT dimj)
------------------------------------------------------------------------------*/
MATRIX *pn6n(INT k,INT dimj)       /* MATRIX  (JN|Pk+6,k(J)|JM) */           /* Gesamtdrehimpuls J , dimj := 2J+1 */
{
    INT    m,n;
/*  INT    s; */
    DOUBLE jzn,jzm/*,jpn*/,jpm;
    DOUBLE fac(INT n);
    MATRIX *okq,*mx_alloc(INT anz_ze,INT anz_sp);
 
    #include "define_j.c"      /* mj,J2,J4,... definieren */
    okq = mx_alloc(dimj,dimj); /* Speicher fuer Matrix(JN|Pk+6,k(J)|JM)*/
 
/*  s = ABS(k); */
    for( n=dimj ; n>=1 ; --n )
         for( m=dimj ; m>=1 ; --m ){
 
         jzn =(             231       *nj6
                -( 315*J2 -735)       *nj4
                +( 105*J4 -525*J2+294)*nj2
                -5*J6 + 40*J4 -60*J2
              );
 
         jzm =(             231       *mj6
                -( 315*J2 -735)       *mj4
                +( 105*J4 -525*J2+294)*mj2
                -5*J6 + 40*J4 -60*J2
              );
 
         jpm = JPM(k,mj)*D(nj,mj+k);
 
         R(okq,n,m) = jpm*(jzn + jzm )/2.0;
    }
 
    return( okq );
}
/*----------------------------------------------------------------------------
                               Clm(INT l,INT m,DOUBLE theta,DOUBLE phi)
------------------------------------------------------------------------------*/
KOMPLEX *Clm(INT l,INT m,DOUBLE theta,DOUBLE phi) /* Tensor Clm berechnen */
  /* in Grad */
{
    DOUBLE sin(DOUBLE f);
    DOUBLE cos(DOUBLE f);
    DOUBLE pow_(DOUBLE x,INT n) ;
    DOUBLE sqrt(DOUBLE f);
 /* DOUBLE pi;                      extern in <math.h> definiert */
    DOUBLE alpha_m,beta_lm,sin_m,g_lm,glm(INT l,INT m,DOUBLE x);
    DOUBLE fac(INT n);
 
    KOMPLEX *clm;
 
    clm     = KX_ALLOC(1);
 
    theta   *= pi/180.0;  /* Umrechnung auf Bogenmass */
    phi     *= pi/180.0;
 
    alpha_m = pow_( -1.0 , (m + ABS(m))/2 );
 
    beta_lm = fac(l)/pow_(2.0,l)/ sqrt( fac(l-ABS(m))*fac(l+ABS(m)) );
 
    sin_m   = pow_( sin(theta) , ABS(m) );
    g_lm    = glm( l , ABS(m) , cos(theta) );
 
    RT(clm)  = cos(m*phi) * alpha_m * beta_lm * sin_m * g_lm;
    IT(clm)  = sin(m*phi) * alpha_m * beta_lm * sin_m * g_lm;
 
    return( clm );
}
/*----------------------------------------------------------------------------
                               glm(INT l,INT m,DOUBLE x)
------------------------------------------------------------------------------*/
DOUBLE glm(INT l,INT m,DOUBLE x) /* darf nur fuer positive m kleiner als l aufgerufen werden*/
{
    DOUBLE binom1,binom2,glm=0.0;
    DOUBLE pow_(DOUBLE x,INT n) ;
    DOUBLE fac(INT n);
    INT i;
 
    for(i=0 ; i<=l-m ; ++i ){
       binom1 = fac(l+m)/fac(l-i)  /fac(i+m);
       binom2 = fac(l-m)/fac(l-m-i)/fac(i  );
       glm += pow_(x+1,i)*pow_(x-1,l-m-i)*binom1*binom2                  ;
    }
 
    return( glm );
}
/*----------------------------------------------------------------------------
                               fac(INT n)
------------------------------------------------------------------------------*/
DOUBLE fac(INT n)    /* Fakultaet von n berechnen */
{
    return( (n==0) ? 1.0 : (DOUBLE)n*fac(n-1) );
}
/*----------------------------------------------------------------------------
                          info_tensor_Clm(INT l,INT m,DOUBLE theta,DOUBLE phi)
------------------------------------------------------------------------------*/
void info_tensor_Clm(INT l,INT m,DOUBLE theta,DOUBLE phi) /* Tensor Clm  auf file name ausgeben */
{
    CHAR *name = "results/Clm.info";
 
    FILE    *fopen(const char * filename,const char * mode),*out;
    KOMPLEX *Clm(INT l,INT m,DOUBLE theta,DOUBLE phi),*z;
    DOUBLE  round(DOUBLE digit);
    CHAR *textr="Realteil      von C%1d%1d(%6.2f, %6.2f) : %13.6e\n";
    CHAR *texti="Imaginaerteil von C%1d%1d(%6.2f, %6.2f) : %13.6e\n";
 
    out = fopen_errchk(name,"w");
    clearscreen;
    printf("Informationen werden auf File %s ausgegeben.\n",name);
 
 
  fprintf(out,"==========================================================\n");
  fprintf(out,"|The tensor Clm(theta,phi) is calculated for             |\n");
  fprintf(out,"|l = %3d                                                 |\n",l);
  fprintf(out,"|m = %3d                                                 |\n",m);
  fprintf(out,"==========================================================\n");
  fprintf(out,"|theta := Polangle   phi := Azimutal angle   so that:   |\n");
  fprintf(out,"|x := sin(theta) * cos(phi)                              |\n");
  fprintf(out,"|y := sin(theta) * sin(phi)                              |\n");
  fprintf(out,"|z := cos(theta)                                         |\n");
  fprintf(out,"----------------------------------------------------------\n");
  fprintf(out,"\n");
 
    z =  Clm( l , m , theta , phi );
    RT(z) = round( RT(z) );
    IT(z) = round( IT(z) );
    fprintf(out,textr,l,m,theta,phi,RT(z) );
    fprintf(out,texti,l,m,theta,phi,IT(z) );
    printf(textr,l,m,theta,phi,RT(z) );
    printf(texti,l,m,theta,phi,IT(z) );
    free_(z);
 
    fclose(out);
}
/*----------------------------------------------------------------------------
                                  round(DOUBLE digit)
------------------------------------------------------------------------------*/
DOUBLE round(DOUBLE digit)  /* signifikante Stellen holen */
{
 return( 1000*((digit/1000 +1.0) - 1.0)  );
}
/*----------------------------------------------------------------------------
                                  pow_(DOUBLE x,INT n) 
------------------------------------------------------------------------------*/
DOUBLE pow_(DOUBLE x,INT n)   /*    x  hoch  n  , mit n>=0  */
{
    return(   (n==0)? 1.0 : (x==0.0)? 0.0 : x*pow_(x,n-1)   );
}
/*----------------------------------------------------------------------------
                                  ggt_r(DOUBLE a,LONG z,LONG n)
------------------------------------------------------------------------------*/
BRUCH   *ggt_r(DOUBLE a,LONG z,LONG n)  /* groessten gemeinsamen Teiler von a und b bestimmen */
              /* a,b=z/n  rational */
{
   BRUCH   *ar/*,*br*/,*is_rational(DOUBLE zahl),*bruch;
   LONG    a_nenner,a_zaehler;
   LONG    b_nenner,b_zaehler;
   LONG      nenner,  zaehler;
   LONG    ggt_l(LONG a,LONG b);
 
   ar = is_rational( a );
   a_zaehler = ZAEHLER( ar );
   a_nenner  = NENNER( ar );
 
   b_zaehler = z;
   b_nenner  = n;
 
   zaehler   = ggt_l( a_zaehler , b_zaehler );
   nenner    = ggt_l( a_nenner  , b_nenner  );
 
   bruch     = BRUCH_ALLOC(1);
   ZAEHLER(bruch) = zaehler;
   NENNER( bruch) = nenner ;
 
   return( bruch );
}
/*----------------------------------------------------------------------------
                              suche_quadratzahl(BRUCH   *bruch)
------------------------------------------------------------------------------*/
BRUCH   *suche_quadratzahl(BRUCH   *bruch)
{
   LONG zaehler,nenner,sucheq(LONG zahl);
/* LONG z =1;
   LONG n =1; */
 
    zaehler = ZAEHLER(bruch);
    nenner  = NENNER(bruch);
 
    zaehler = sucheq(zaehler);
    nenner  = sucheq(nenner );
 
    ZAEHLER(bruch) = zaehler;
    NENNER( bruch) = nenner;
 
    return(bruch);
}
/*----------------------------------------------------------------------------
                           is_quadrat(LONG zahl)
------------------------------------------------------------------------------*/
LONG is_quadrat(LONG zahl)    /* ist zahl = q * q ? */
{
    LONG q;
 
    for( q=0 ; q*q <= zahl ; ++q )
           if( q*q == zahl )
                  return( q );
    return( -(--q) );
}
/*----------------------------------------------------------------------------
                              sucheq(LONG zahl)
------------------------------------------------------------------------------*/
LONG sucheq(LONG zahl)    /* suche maximales q mit :   q*q teilt zahl */
{
    LONG i,q/*,s*/;
 
    if( (q=is_quadrat(zahl)) >= 0)
         return( q );
    for(i=-q;i>=1;--i)
        if( zahl % (i*i) == 0 )
               return(i);
 
    return( 1 );
}
/*----------------------------------------------------------------------------
                              is_rational(DOUBLE zahl)
------------------------------------------------------------------------------*/
BRUCH   *is_rational(DOUBLE zahl)
{
       LONG    bereich = 1000000;
       LONG    i       = 0 ;
       DOUBLE  rest;
       DOUBLE  macheps=0.0001;
       DOUBLE  test;
       DOUBLE  sign=1.0;
       BRUCH   *bruch;
 
 
       if( zahl< 0 ) { sign = -1.0; zahl *= -1.0; }
       do{
            test  = zahl * ++i;
            rest  = test - (LONG)(test+macheps);
            rest *= 1/macheps;
            rest  = (LONG)(rest);
            if( rest <  1.0 ){
                bruch = BRUCH_ALLOC(1);
                ZAEHLER(bruch)  = test+macheps;
                ZAEHLER(bruch) *= sign;
                NENNER( bruch)  = i;
                return(bruch);
            }
 
       }while(i<=bereich );
 
       printf(" \nError in STEVENS.C in is_rational() :\n");
       printf(" the number %f ist not rational .\n\n",zahl);
       exit(1);
}
/*----------------------------------------------------------------------------
                                  ggt_l(LONG a,LONG b)
------------------------------------------------------------------------------*/
LONG ggt_l(LONG a,LONG b) /* groessten gemeinsamen Teiler von a und b bestimmen */
{
   LONG   c,rest;
/* LONG   zaehler,nenner; */
 
   if( a==0 && b==0 )   return( 1 );
   if( a==0 )           return( b );
   if( b==0 )           return( a );
   if( a==b )           return( a );
 
 
 
   if( b > a ) { c=a; a=b; b=c; }
 
   do{
         rest = a%b;
         if( rest!=0 ) { a=b; b=rest; }
 
   }while( rest!=0 );
 
   return( b );
}
/*------------------------------------------------------------------------------
ENDEMODUL    S T E V E N S    C
------------------------------------------------------------------------------*/
