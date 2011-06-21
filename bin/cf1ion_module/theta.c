/*-----------------------------------------------------------------------------
 
                                M  O  D  U  L
 
                                   THETA C
 
-------------------------------------------------------------------------------
 
Aufgabe               : alpha_J , beta_J , gamma_J  und
 
                          n
                        <r > definieren  ( n = 2 , 4 , 6 )
 
-------------------------------------------------------------------------------
 
Literatur zur Berechnung der alpha_J , beta_J , gamma_J :
---------------------------------------------------------
[0] Diplomarbeit Peter Hoffmann : Dort wird gezeigt werden ,dass
                                  die alpha_J,beta_J und gamma_J
                                  fuer die allgemeinen Operatoren
                                  die Gleichen sind ,wie fuer die
                                  Stevensoperatoren ,welche in [1]
                                  tabelliert sind.
 
[1] Abragam und Bleaney "Electronic Paramagnetic Resonance of Transition Ions",
                        1970 ,Table 20
 
                               2      4      6
Literatur zur Berechnung der <r > , <r > , <r > :
-------------------------------------------------
[2] A.J. Freeman and J.P. Desclaux: Journal of Magnetism and Magnetic Materials
                                    12 (1979) 11
[3] W. B. Lewis et al.: The Journal of Chemical Physics 53 (1970) p 809
 
-------------------------------------------------------------------------------
 
Bemerkungen           :
 
-------------------------------------------------------------------------------
 
Definierte Funktionen :
 
-----------------------
 
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
                                define
-----------------------------------------------------------------------------*/
 
#define NR_RE 15
 
/*----------------------------------------------------------------------------
Extern definierte Funktionen
-----------------------------------------------------------------------------*/
extern FILE *fopen_errchk();         /* definiert in EINGABE.C*/ 

/*ionenanzahl muss mit ANZ_IONEN in cfield.c und cfieldrout.c uebereinstimmen !!!*/

IONEN IONENIMP [] = {
/*           (for stevensfactors)                                                  */
/*                       V                                                         */
/*            ion    e_in_4f  g     2J+1  F(4)   F(6)   <r^2>    <r^4>    <r^6>    */
/*                        5f   J          for x-W pars  (a0^2)   (a0^4)   (a0^6)   */
/*                                                                                 */
/*  0 */   { "Ce3+"  ,  1 , 6.0/7  ,  6  , 60 ,     1 ,1.309    ,3.964   ,23.31},
/*  1 */   { "Pr3+"  ,  2 , 4.0/5  ,  9  , 60 ,  1260 ,1.1963   ,3.3335  ,18.353},/* r^n von U.Walter  Diss.   */
/*  2 */   { "Nd3+"  ,  3 , 8.0/11 , 10  , 60 ,  2520 ,1.114    ,2.910   ,15.03},
/*  3 */   { "Pm3+"  ,  4 , 3.0/5  ,  9  , 60 ,  1260 ,1.0353   ,2.5390  ,12.546},/* r^n von U.Walter  Diss.   */
/*  4 */   { "Sm3+"  ,  5 , 2.0/7  ,  6  , 60 ,     1 ,0.9743   ,2.260   ,10.55},
/*  5 */   { "Eu3+"  ,  6 , 0.0    ,  1  ,  1 ,     1 ,0.9175   ,2.020   ,9.039},
/*  6 */   { "Gd3+"  ,  7 , 2.0    ,  8  , 60 ,  1260 ,0.8671   ,1.820   ,7.831},
/*  7 */   { "Tb3+"  ,  8 , 3.0/2  , 13  , 60 ,  7560 ,0.8220   ,1.651   ,6.852},
/*  8 */   { "Dy3+"  ,  9 , 4.0/3  , 16  , 60 , 13860 ,0.7814   ,1.505   ,6.048},
/*  9 */   { "Ho3+"  , 10 , 5.0/4  , 17  , 60 , 13860 ,0.7446   ,1.379   ,5.379},
/* 10 */   { "Er3+"  , 11 , 6.0/5  , 16  , 60 , 13860 ,0.7111   ,1.270   ,4.816},
/* 11 */   { "Tm3+"  , 12 , 7.0/6  , 13  , 60 ,  7560 ,0.6804   ,1.174   ,4.340},
/* 12 */   { "Yb3+"  , 13 , 8.0/7  ,  8  , 60 ,  1260 ,0.6522   ,1.089   ,3.932},
/* 13 */   { "Nd2+"  ,  4 , 3.0/5  ,  9  , 60 ,  1260 ,1.392    ,5.344   ,45.450},
/* 14 */   { "Sm2+"  ,  6 , 0.0    ,  1  ,  1 ,     1 ,1.197    ,3.861   ,28.560},
/* 15 */   { "Eu2+"  ,  7 , 2.0    ,  8  , 60 ,  1260 ,1.098    ,3.368   ,23.580},
/* 16 */   { "Gd2+"  ,  8 , 3.0/2  , 13  , 60 ,  7560 ,1.028    ,2.975   ,19.850},
/* 17 */   { "Tb2+"  ,  9 , 4.0/3  , 16  , 60 , 13860 ,0.968    ,2.655   ,16.980},
/* 18 */   { "Dy2+"  , 10 , 5.0/4  , 17  , 60 , 13860 ,0.913    ,2.391   ,14.730},
/* 19 */   { "Ho2+"  , 11 , 6.0/5  , 16  , 60 , 13860 ,0.866    ,2.169   ,12.920},
/* 20 */   { "Er2+"  , 12 , 7.0/6  , 13  , 60 ,  7560 ,0.824    ,1.979   ,11.450},
/* 21 */   { "Tm2+"  , 13 , 8.0/7  ,  8  , 60 ,  1260 ,0.785    ,1.819   ,10.240},
/* 22 */   { "U4+"   ,  2 , 4.0/5  ,  9  , 60 ,  1260 ,2.042    ,7.632   ,47.774}, /* r^n Freeman et al. PRB 13 (1976) 1168 */
/* 23 */   { "U3+"   ,  3 , 8.0/11 , 10  , 60 ,  2520 ,2.346    ,10.906  ,90.544}, /* r^n Freeman et al. PRB 13 (1976) 1168 */
/* 24 */   { "U2+"   ,  4 , 3.0/5  ,  9  , 60 ,  1260 ,3.257    ,26.82   ,462.85},  /* r^n Lewis et al. J. Chem Phys. 53 (1970) 809 average j=5/2 and j=7/2 in DS method*/
/* 25 */   { "Np4+"  ,  3 , 8.0/11 , 10  , 60 ,  2520 ,1.884    ,6.504   ,37.80},  /* r^n Lewis et al. J. Chem Phys. 53 (1970) 809 average j=5/2 and j=7/2 in DF method*/
/* 26 */   { "Np3+"  ,  4 , 3.0/5  ,  9  , 60 ,  1260 ,2.297    ,11.00   ,98.63},  /* r^n Lewis et al. J. Chem Phys. 53 (1970) 809 average j=5/2 and j=7/2 in DS method*/
/* 27 */   { "Pu4+"  ,  4 , 3.0/5  ,  9  , 60 ,  1260 ,1.838    ,6.401   ,38.77},  /* r^n Lewis et al. J. Chem Phys. 53 (1970) 809 average j=5/2 and j=7/2 in DS method*/
/* 28 */   { "Pu3+"  ,  5 , 2.0/7  ,  6  , 60 ,     1 ,2.1025   ,9.1775  ,75.30},  /* r^n Lewis et al. J. Chem Phys. 53 (1970) 809 average j=5/2 and j=7/2 in DS method*/
/* a pure spin S=int or half int problem */
/* 29 */   { "S="    ,  0 , 2.0    ,  0  ,  0 ,     0 ,1e300    ,1e300   ,1e300},  /* r^n Lewis et al. J. Chem Phys. 53 (1970) 809 average j=5/2 and j=7/2 in DS method*/

}; 
#define ANZ_IONEN              (  sizeof(IONENIMP)/sizeof(IONEN)  )

/*----------------------------------------------------------------------------
                                alpha_J
-----------------------------------------------------------------------------*/
 
DOUBLE alpha_J[NR_RE]={
/* 4f_0        */    1.0 * 0                ,
/* 4f_1 : Ce3+ */   -1.0 * 2/5/7            ,
/* 4f_2 : Pr3+ */   -1.0 * 2*2*13/3/3/5/5/11,
/* 4f_3 : Nd3+ */   -1.0 * 7/3/3/11/11      ,
/* 4f_4 : Pm3+ */    1.0 * 2*7/3/5/11/11    ,
/* 4f_5 : Sm3+ */    1.0 * 13/3/3/5/7       ,
/* 4f_6 : Eu3+ */    1.0 * 0                ,
/* 4f_7 : Gd3+ */    1.0 * 0                ,
/* 4f_8 : Tb3+ */   -1.0 * 1/3/3/11         ,
/* 4f_9 : Dy3+ */   -1.0 * 2/3/3/5/7        ,
/* 4f_10: Ho3+ */   -1.0 * 1/2/3/3/5/5      ,
/* 4f_11: Er3+ */    1.0 * 2*2/3/3/5/5/7    ,
/* 4f_12: Tm3+ */    1.0 * 1/3/3/11         ,
/* 4f_13: Yb3+ */    1.0 * 2/3/3/7          ,
/* 4f_14       */    1.0 * 0                
};
/*----------------------------------------------------------------------------
                                 beta_J
-----------------------------------------------------------------------------*/
 
DOUBLE beta_J[NR_RE]={
/* 4f_0        */    1.0 * 0                              ,
/* 4f_1 : Ce3+ */    1.0 * 2/3/3/5/7                      ,
/* 4f_2 : Pr3+ */   -1.0 * 2*2/3/3/5/11/11                ,
/* 4f_3 : Nd3+ */   -1.0 * 2*2*2*17/3/3/3/11/11/11/13     ,
/* 4f_4 : Pm3+ */    1.0 * 2*2*2*7*17/3/3/3/5/11/11/11/13 ,
/* 4f_5 : Sm3+ */    1.0 * 2*13/3/3/3/5/7/11              ,
/* 4f_6 : Eu3+ */    1.0 * 0                              ,
/* 4f_7 : Gd3+ */    1.0 * 0                              ,
/* 4f_8 : Tb3+ */    1.0 * 2/3/3/3/5/11/11                ,
/* 4f_9 : Dy3+ */   -1.0 * 2*2*2/3/3/3/5/7/11/13          ,
/* 4f_10: Ho3+ */   -1.0 * 1/2/3/5/7/11/13                ,
/* 4f_11: Er3+ */    1.0 * 2/3/3/5/7/11/13                ,
/* 4f_12: Tm3+ */    1.0 * 2*2*2/3/3/3/3/5/11/11          ,
/* 4f_13: Yb3+ */   -1.0 * 2/3/5/7/11                     ,
/* 4f_14       */    1.0 * 0                              
};
/*----------------------------------------------------------------------------
                                gamma_J
-----------------------------------------------------------------------------*/
 
DOUBLE gamma_J[NR_RE]={
/* 4f_0        */    1.0 * 0                                 ,
/* 4f_1 : Ce3+ */    1.0 * 0                                 ,
/* 4f_2 : Pr3+ */    1.0 * 2*2*2*2*17/3/3/3/3/5/7/11/11/13   ,
/* 4f_3 : Nd3+ */   -1.0 * 5*17*19/3/3/3/7/11/11/11/13/13    ,
/* 4f_4 : Pm3+ */    1.0 * 2*2*2*17*19/3/3/3/7/11/11/11/13/13,
/* 4f_5 : Sm3+ */    1.0 * 0                                 ,
/* 4f_6 : Eu3+ */    1.0 * 0                                 ,
/* 4f_7 : Gd3+ */    1.0 * 0                                 ,
/* 4f_8 : Tb3+ */   -1.0 * 1/3/3/3/3/7/11/11/13              ,
/* 4f_9 : Dy3+ */    1.0 * 2*2/3/3/3/7/11/11/13/13           ,
/* 4f_10: Ho3+ */   -1.0 * 5/3/3/3/7/11/11/13/13             ,
/* 4f_11: Er3+ */    1.0 * 2*2*2/3/3/3/7/11/11/13/13         ,
/* 4f_12: Tm3+ */   -1.0 * 5/3/3/3/3/7/11/11/13              ,
/* 4f_13: Yb3+ */    1.0 * 2*2/3/3/3/7/11/13                 ,
/* 4f_14       */    1.0 * 0                                 
};


/*----------------------------------------------------------------------------
                               info_thetakq()
-----------------------------------------------------------------------------*/
void info_thetakq()/* Liste der Stevensfaktoren alpha_J,beta_J,gamma_J zeigen*/
{
    CHAR *name = "results/thetakq.info";
 
    FILE *fopen(),*fp;
    CHAR *s,*text=" %4s  %2i     %10.4f     %10.4f      %10.4f\n";
    DOUBLE a,b,c;
    DOUBLE x=100.0;
    DOUBLE y=10000.0;
    DOUBLE z=1000000.0;
    INT i;
 
    fp = fopen_errchk(name,"w");
    clearscreen;
    printf("Information contained in the File %s .\n",name);
 
    fprintf(fp,"============================================================\n");
    fprintf(fp,"|Table of the Stevens factors alpha_J ,beta_J und gamma_J. |\n");
    fprintf(fp,"============================================================\n");
    fprintf(fp,"In the disseration of  Peter Hoffmann it is shown that     \n");
    fprintf(fp,"                                    2S+1                   \n");
    fprintf(fp,"the quantisation constant  theta   (     L   ) (k=2,4,6)   \n");
    fprintf(fp,"                                 k        J                \n");
    fprintf(fp,"of the general Stevensoperatoren with alpha_J , beta_J     \n");
    fprintf(fp,"und gamma_J can be identified with the usual               \n");
    fprintf(fp,"Stevensoperatoren.                                         \n");
    fprintf(fp,"                                                           \n");
    fprintf(fp,"The Stevens factors are rational numbers , in which        \n");
    fprintf(fp,"Abragam und Bleaney   \"Electronic Paramagnetic Resonance  \n");
    fprintf(fp,"                     of Transition Ions\",1970,Table  20   \n");
    fprintf(fp,"calculated exactly ,but here on 4 decimal places           \n");
    fprintf(fp,"are shown .                                                \n");
    fprintf(fp,"\n\n");
    fprintf(fp,"============================================================\n");
    fprintf(fp,"|                          2             4                6|\n");
    fprintf(fp,"|ion  nr.of      alpha_J*10     beta_J*10       gamma_J*10 |\n");
    fprintf(fp,"|    f-electrons                                           |\n");
    fprintf(fp,"============================================================\n");

    for (i=0;i<(int)ANZ_IONEN;++i){
                              s=IONENIMP[i].ionname;a=x*alpha_J[E4f(i)];b=y*beta_J[E4f(i)];c=z*gamma_J[E4f(i)];fprintf(fp,text,s,E4f(i),a,b,c);
                             }
    fclose(fp);
}
/*----------------------------------------------------------------------------
                                  info_rn()
-----------------------------------------------------------------------------*/
                   /*             n                   */
void info_rn()     /* Liste der <r > n=2,4,6 ausgeben */
{
    CHAR *name = "results/rn.info";
 
    FILE   *fopen(),*fp;
    DOUBLE a,b,c,a2,a4,a6;
    INT i;
 
    CHAR *s;
    CHAR *texti =" %2s       %10.4f    %10.4f    %10.4f nach U. Walter\n";
    CHAR *textfw=" %2s      %10.3f    %10.3f    %10.3f              \n";
    CHAR *sexti =" %2s      %10.4f    %10.4f    %10.4f nach U. Walter\n";
    CHAR *sextfw=" %2s      %10.4f    %10.4f    %10.4f              \n";
 
    fp = fopen_errchk(name,"w");
    clearscreen;
    printf("Information contained in the file %s.\n",name);
 
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"|             n                                           |\n");
    fprintf(fp,"|Table of <r > for n = 2 , 4 , 6 : a0=Bohr Radius         |\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"That over the wavefunction  R (r) of the 4f-Elektrons    \n");
    fprintf(fp,"                          n      4f                       \n");
    fprintf(fp,"average                 <r > after Freeman and            \n");
    fprintf(fp,"Desclaux Journal of Magnetism and Magnetic Materials       \n");
    fprintf(fp,"12 (1979) 11.                                             \n");
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"|              2    2        4    4        6    6         |\n");
    fprintf(fp,"|ion         <r >/a0       <r >/a0       <r >/a0          |\n");
    fprintf(fp,"-----------------------------------------------------------\n");

    for (i=0;i<(int)ANZ_IONEN;++i){
                              if(i==1||i==3){s=IONENIMP[i].ionname;a=r2(i);b=r4(i);c=r6(i);fprintf(fp,texti,s,a,b,c);}
                              else    {s=IONENIMP[i].ionname;a=r2(i);b=r4(i);c=r6(i);fprintf(fp,textfw,s,a,b,c);}
                             }

    fprintf(fp,"\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"|              2   o 2       4   o 4       6   o 6        |\n");
    fprintf(fp,"|RE3+        <r >/ A       <r >/ A       <r >/ A          |\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    a2 = A0_BOHR*A0_BOHR;
    a4 = a2*a2;
    a6 = a2*a4;

    for (i=0;i<(int)ANZ_IONEN;++i){
                              if(i==1||i==3){s=IONENIMP[i].ionname;a=r2(i)*a2;b=r4(i)*a4;c=r6(i)*a6;fprintf(fp,sexti,s,a,b,c);}
                              else    {s=IONENIMP[i].ionname;a=r2(i)*a2;b=r4(i)*a4;c=r6(i)*a6;fprintf(fp,sextfw,s,a,b,c);}
                             }


    fprintf(fp,"\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"|                 2                 4                  6  |\n");
    fprintf(fp,"|  a := alpha_J*10    b := beta_J*10    c := gamma_J*10   |\n");
    fprintf(fp,"|                                                         |\n");
    fprintf(fp,"|              2   o 2       4   o 4       6   o 6        |\n");
    fprintf(fp,"|RE3+      a*<r >/ A     b*<r >/ A     c*<r >/ A          |\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    a2 = A0_BOHR*A0_BOHR*100;
    a4 = a2*a2;
    a6 = a2*a4;
    for (i=0;i<(int)ANZ_IONEN;++i){
                              if(i==1||i==3){s=IONENIMP[i].ionname;a=r2(i)*a2*alpha_J[E4f(i)];b=r4(i)*a4*beta_J[E4f(i)];c=r6(i)*a6*gamma_J[E4f(i)];fprintf(fp,sexti,s,a,b,c);}
                              else          {s=IONENIMP[i].ionname;a=r2(i)*a2*alpha_J[E4f(i)];b=r4(i)*a4*beta_J[E4f(i)];c=r6(i)*a6*gamma_J[E4f(i)];fprintf(fp,sextfw,s,a,b,c);}
                             }

    fclose(fp);
}
/*------------------------------------------------------------------------------
ENDEMODUL    T H E T A     C
------------------------------------------------------------------------------*/
