 
/*-----------------------------------------------------------------------------
 
                                M   O  D  U  L
 
                                   EINGABE C
 
-------------------------------------------------------------------------------
 
Aufgabe               :  Parametersatzfiles erzeugen und lesen
 
-----------------------------------------------------------------------------*/
 
 
/*----------------------------------------------------------------------------
Includedateien holen
-----------------------------------------------------------------------------*/
#include <stdio.h>          /* damit FILE definiert wird               */
#include <stdlib.h>
#include <string.h>
#include <math.h>           /* damit sqrt in define_j.c definiert wird */
#define pi (4.0*atan(1.0))  /* atan() braucht <math.h>                 */
#include "types.c"          /* benutze Datentypen laden                */
/*----------------------------------------------------------------------------
extern definierte Dateien
-----------------------------------------------------------------------------*/
 
 
extern DOUBLE alpha_J[]; /* definiert in theta.c */
extern DOUBLE beta_J[];  /* definiert in theta.c */
extern DOUBLE gamma_J[]; /* definiert in theta.c */
 
 
extern QUADRANT *abfrage();      /* definiert in ORTHO.C    */
extern INT       null();         /* definiert in DIAHERMX.C */
extern INT      *sort();         /* definiert in DIAHERMX.C */
extern DOUBLE  accuracy();       /* definiert in DIAHERMX.c */
extern BRUCH     *is_rational(); /* definiert in STEVENS.C*/
extern LONG      ggt_l();        /* definiert in STEVENS.C*/
extern INT       is_equal();     /* definiert in DIAHERMX.C */
extern DOUBLE    exp_();             /* definiert in ORTHO.C */
 
extern VEKTOR *vr_alloc();    /* definiert in MATRIX.C */
extern MATRIX *mx_alloc();    /* definiert in MATRIX.C */
 
extern DOUBLE omegan0n();     /* definiert in MAIN.C  */
extern DOUBLE omegan1n();     /* definiert in MAIN.C  */
extern DOUBLE omegan2n();     /* definiert in MAIN.C  */
extern DOUBLE omegan3n();     /* definiert in MAIN.C  */
extern DOUBLE omegan4n();     /* definiert in MAIN.C  */
extern DOUBLE omegan5n();     /* definiert in MAIN.C  */
extern DOUBLE omegan6n();     /* definiert in MAIN.C  */
extern DOUBLE epn0n();     /* definiert in MAIN.C  */
extern DOUBLE epn1n();     /* definiert in MAIN.C  */
extern DOUBLE epn2n();     /* definiert in MAIN.C  */
extern DOUBLE epn3n();     /* definiert in MAIN.C  */
extern DOUBLE epn4n();     /* definiert in MAIN.C  */
extern DOUBLE epn5n();     /* definiert in MAIN.C  */
extern DOUBLE epn6n();     /* definiert in MAIN.C  */
 
extern INT    is_einheit_imp();  /* definiert in MAIN.C */
extern INT    isinlimits();      /* definiert in MAIN.C */
extern INT    isimplementiert(); /* definiert in MAIN.C */
extern INT    isreell();         /* definiert in MAIN.C */
extern DOUBLE a_tof();           /* definiert in MAIN.C  */
extern INT    a_toi();           /* definiert in MAIN.C */
extern CHAR  *a_tos();           /* definiert in MAIN.C */
extern INT    read_error();      /* definiert in MAIN.C */
extern INT    Bkq_error();       /* definiert in MAIN.C */
extern INT    write_title();     /* definiert in DIAHERMX.C*/
  
extern CHAR  *leftcopy();        /* definiert in MAIN.C */
 
extern IONEN     IONENIMP[];     /* definiert in MAIN.C  */
extern EINHEIT   EINHEITIMP[];   /* definiert in MAIN.C  */
extern SYMMETRIE SYMLISTE[];     /* definiert in MAIN.C  */
 
extern ITERATION *iter_alloc();  /* definiert in MAIN.C */
extern ITERATION *auswahlregel();/* definiert in MAIN.C */
extern MATRIX    *calc_Bmag();   /* definiert in MAIN.C */
extern MATRIX    *calc_Bmag_D(); /* definiert in MAIN.C */
extern MATRIX    *calcBmol();    /* definiert in MAIN.C */
extern STEVENS   *calc_Pkq();    /* definiert in STEVENS.C */

void drucke_par();
void drucke_mag();
/*INT strstr();
INT strchr();
INT strlen();
INT strncmp();
 */
/*open a file: similar fopen but with error check */
FILE * fopen_errchk(const char * filename,const char * mode)
{ FILE *file;
 
  file = fopen (filename,mode);
  if (file == NULL)
    { 
      fprintf (stderr, "Error: Couldn't open file %s - possible reason: directory results does not exist\n",filename);
      exit (EXIT_FAILURE);
     }
  return file;
} 


 
/*******************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
**********************  E I N G A B E  -  F I L E S   **************************
**********************                                **************************
**********************        E R Z E U G E N         **************************
********************************************************************************
********************************************************************************
********************************************************************************
*******************************************************************************/
 
/*------------------------------------------------------------------------------
                                create_Vkq()
------------------------------------------------------------------------------*/
void create_Vkq(einheitnr_in,einheitnr_out,ion,symmetrienr,modus,temp)
    INT    einheitnr_in,einheitnr_out;
    CHAR   *ion;
    INT    symmetrienr;
    CHAR   modus;
    DOUBLE temp; /* Temperatur der Probe */
{
 
    CHAR *name=VKQNAME;
    CHAR *einheit_in;
    CHAR *einheit_out;
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11;
    CHAR *t20,*t21,*t22,*t23,*t24,*t25;
    CHAR *t26,*t27,*t28,*t29,*t30,*t31;
    CHAR *leftcopy();
 
    FILE    *fopen(), *fp;
    TABELLE *tabelle;
    INT     dimj,ionennr;
 
    fp  = fopen_errchk(name,"w");
    write_title(fp);
 
    ionennr     = isimplementiert(ion);
    dimj        = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

    einheit_in  = EINHEITIMP[ einheitnr_in ].einheit;
    einheit_out = EINHEITIMP[ einheitnr_out].einheit;
    ion         = leftcopy(ion ,25);
 
    tabelle = TABELLE_ALLOC(1);
    t01 = "===============================================================\n";
    t02 = "|                                                             |\n";
    t03 = "|Crystal Field parameter  V in   %6s    in which           |\n";
    t04 = "|                         kq  (Vkq are NOT the same as in     |\n";
    t05 = "|                              Elliot and Stevens Notation)   |\n";
    t06 = "|        ---                ---                  *   +        |\n";
    t07 = "| H   =  >   V   P  (J)  +  >    V   P  (J)  +  V   P  (J)    |\n";
    t08 = "|  KF    ---  k0  k0        ---   kq  kq         kq  kq       |\n";
    t09 = "|         k                 k>0                               |\n";
    t10 = "|                           q>0                               |\n";
    t20 = "|                                                             |\n";
    t21 = "|                                      2S+1                   |\n";
    t22 = "| Where:   V   :=  epsilon  * theta  (      L ) * Q           |\n";
    t23 = "|           kq            kq       kq        J     kq         |\n";
    t24 = "|                                                             |\n";
    t25 = "|          Q   := Multipolmoment of the crystal Field         |\n";
    t26 = "|           kq                                                |\n";
    t27 = "|                       /                                     |\n";
    t28 = "|                      |  rho(R')       *          3          |\n";
    t29 = "|              = -|e|  | ------------- C (omega') d R'        |\n";
    t30 = "|                      |      k+1       kq                    |\n";
    t31 = "|                      /  |R'|                                |\n";
    t11 = "===============================================================\n";
    TSS = "| Energy Eigenvalues are in  : %6s                         |\n";
    T15 = "| Temperature of the probe   : %7.2f Kelvin                 |\n";
    T12 = "| Ion                        : %25s     |\n";
    T13 = "| Symmetry         : %s     Symmetry number : %1d       |\n";
    T14 = "| Magnetic field             : %s                 |\n";
    T11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T20n= "| ReV20 :                                                     |\n";
    T20r= "|*ReV20 :                                                     |\n";
    STR = "|-------------------------------------------------------------|\n";
    T21n= "| ReV21 :                      | ImV21 :                      |\n";
    T21r= "|*ReV21 :                      | ImV21 :                      |\n";
    T21c= "|*ReV21 :                      |*ImV21 :                      |\n";
 
    T22n= "| ReV22 :                      | ImV22 :                      |\n";
    T22r= "|*ReV22 :                      | ImV22 :                      |\n";
    T22c= "|*ReV22 :                      |*ImV22 :                      |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T40n= "| ReV40 :                                                     |\n";
    T40r= "|*ReV40 :                                                     |\n";
 
    T41n= "| ReV41 :                      | ImV41 :                      |\n";
    T41r= "|*ReV41 :                      | ImV41 :                      |\n";
    T41c= "|*ReV41 :                      |*ImV41 :                      |\n";
 
    T42n= "| ReV42 :                      | ImV42 :                      |\n";
    T42r= "|*ReV42 :                      | ImV42 :                      |\n";
    T42c= "|*ReV42 :                      |*ImV42 :                      |\n";
 
    T43n= "| ReV43 :                      | ImV43 :                      |\n";
    T43r= "|*ReV43 :                      | ImV43 :                      |\n";
    T43c= "|*ReV43 :                      |*ImV43 :                      |\n";
 
    T44n= "| ReV44 :                      | ImV44 :                      |\n";
    T44r= "|*ReV44 :                      | ImV44 :                      |\n";
    T44c= "|*ReV44 :                      |*ImV44 :                      |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T60n= "| ReV60 :                                                     |\n";
    T60r= "|*ReV60 :                                                     |\n";
 
    T61n= "| ReV61 :                      | ImV61 :                      |\n";
    T61r= "|*ReV61 :                      | ImV61 :                      |\n";
    T61c= "|*ReV61 :                      |*ImV61 :                      |\n";
 
    T62n= "| ReV62 :                      | ImV62 :                      |\n";
    T62r= "|*ReV62 :                      | ImV62 :                      |\n";
    T62c= "|*ReV62 :                      |*ImV62 :                      |\n";
 
    T63n= "| ReV63 :                      | ImV63 :                      |\n";
    T63r= "|*ReV63 :                      | ImV63 :                      |\n";
    T63c= "|*ReV63 :                      |*ImV63 :                      |\n";
 
    T64n= "| ReV64 :                      | ImV64 :                      |\n";
    T64r= "|*ReV64 :                      | ImV64 :                      |\n";
    T64c= "|*ReV64 :                      |*ImV64 :                      |\n";
 
    T65n= "| ReV65 :                      | ImV65 :                      |\n";
    T65r= "|*ReV65 :                      | ImV65 :                      |\n";
    T65c= "|*ReV65 :                      |*ImV65 :                      |\n";
 
    T66n= "| ReV66 :                      | ImV66 :                      |\n";
    T66r= "|*ReV66 :                      | ImV66 :                      |\n";
    T66c= "|*ReV66 :                      |*ImV66 :                      |\n";
    t11 = "===============================================================\n";
 
    fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,t03,einheit_in);fprintf(fp,"%s",t04);
    fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
    fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
    fprintf(fp,"%s",t20);
    fprintf(fp,"%s",t21);
    fprintf(fp,"%s",t22);
    fprintf(fp,"%s",t23);
    fprintf(fp,"%s",t24);
    fprintf(fp,"%s",t25);
    fprintf(fp,"%s",t26);
    fprintf(fp,"%s",t27);
    fprintf(fp,"%s",t28);
    fprintf(fp,"%s",t29);
    fprintf(fp,"%s",t30);
    fprintf(fp,"%s",t31);
    fprintf(fp,"%s",t11);
 
    drucke_par( fp,modus,dimj,tabelle,einheit_out,temp,ion,symmetrienr );
    if(modus != NOMAG )  drucke_mag( fp,modus );
 
    printf("file %s produced.\n",name);
    fclose(fp);
 
}
/*------------------------------------------------------------------------------
                                create_Dkq()
------------------------------------------------------------------------------*/
void create_Dkq(einheitnr_in,einheitnr_out,ion,symmetrienr,modus,temp)
    INT    einheitnr_in,einheitnr_out;
    CHAR   *ion;
    INT    symmetrienr;
    CHAR   modus;
    DOUBLE temp; /* Temperatur der Probe */
{
 
    CHAR *name=DKQNAME;
    CHAR *einheit_in;
    CHAR *einheit_out;
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11;
    CHAR *t20,*t21,*t22,*t23,*t24,*t25;
    CHAR *t26,*t27,*t28,*t29,*t30,*t31;
    CHAR *leftcopy();
 
    FILE    *fopen(), *fp;
    TABELLE *tabelle;
    INT     dimj,ionennr;
 
    fp  = fopen_errchk(name,"w");
    write_title(fp);
 
    ionennr     = isimplementiert(ion);
    dimj        = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

    einheit_in  = EINHEITIMP[ einheitnr_in ].einheit;
    einheit_out = EINHEITIMP[ einheitnr_out].einheit;
    ion         = leftcopy(ion ,25);
 
    tabelle = TABELLE_ALLOC(1);
    t01 = "===============================================================\n";
    t02 = "|                                                             |\n";
    t03 = "|Crystal Field parameter D  in   %6s   in which            |\n";
    t04 = "|                         kq                                  |\n";
    t05 = "|                                                             |\n";
    t06 = "|        ---             ---               *   +              |\n";
    t07 = "| H   =  >   D   C    +  >    D   C    +  D   C               |\n";
    t08 = "|  KF    ---  k0  k0     ---   kq  kq      kq  kq             |\n";
    t09 = "|         k              k>0                                  |\n";
    t10 = "|                        q>0                                  |\n";
    t20 = "|                                                             |\n";
    t21 = "|           *         q                                       |\n";
    t22 = "| mit:     D   =  (-1)   D                                    |\n";
    t23 = "|           kq            k,-q                                |\n";
    t24 = "|                                                             |\n";
    t25 = "| und                                                         |\n";
    t26 = "|           +         q                                       |\n";
    t27 = "|          C   =  (-1)   C                                    |\n";
    t28 = "|           kq            k,-q                                |\n";
    t29 = "|                                                             |\n";
    t30 = "| where the  C    T e n s o r o p e r a t o r s   are   .     |\n";
    t31 = "|             kq                                              |\n";
    t11 = "===============================================================\n";
    TSS = "| Energy Eigenvalues are in  : %6s                         |\n";
    T15 = "| Temperature of the probe   : %7.2f Kelvin                 |\n";
    T12 = "| Ion                        : %25s     |\n";
    T13 = "| Symmetry         : %s     Symmetry number : %1d       |\n";
    T14 = "| Magnetic field             : %s                 |\n";
    T11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T20n= "| ReD20 :                                                     |\n";
    T20r= "|*ReD20 :                                                     |\n";
    STR = "|-------------------------------------------------------------|\n";
    T21n= "| ReD21 :                      | ImD21 :                      |\n";
    T21r= "|*ReD21 :                      | ImD21 :                      |\n";
    T21c= "|*ReD21 :                      |*ImD21 :                      |\n";
 
    T22n= "| ReD22 :                      | ImD22 :                      |\n";
    T22r= "|*ReD22 :                      | ImD22 :                      |\n";
    T22c= "|*ReD22 :                      |*ImD22 :                      |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T40n= "| ReD40 :                                                     |\n";
    T40r= "|*ReD40 :                                                     |\n";
 
    T41n= "| ReD41 :                      | ImD41 :                      |\n";
    T41r= "|*ReD41 :                      | ImD41 :                      |\n";
    T41c= "|*ReD41 :                      |*ImD41 :                      |\n";
 
    T42n= "| ReD42 :                      | ImD42 :                      |\n";
    T42r= "|*ReD42 :                      | ImD42 :                      |\n";
    T42c= "|*ReD42 :                      |*ImD42 :                      |\n";
 
    T43n= "| ReD43 :                      | ImD43 :                      |\n";
    T43r= "|*ReD43 :                      | ImD43 :                      |\n";
    T43c= "|*ReD43 :                      |*ImD43 :                      |\n";
 
    T44n= "| ReD44 :                      | ImD44 :                      |\n";
    T44r= "|*ReD44 :                      | ImD44 :                      |\n";
    T44c= "|*ReD44 :                      |*ImD44 :                      |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T60n= "| ReD60 :                                                     |\n";
    T60r= "|*ReD60 :                                                     |\n";
 
    T61n= "| ReD61 :                      | ImD61 :                      |\n";
    T61r= "|*ReD61 :                      | ImD61 :                      |\n";
    T61c= "|*ReD61 :                      |*ImD61 :                      |\n";
 
    T62n= "| ReD62 :                      | ImD62 :                      |\n";
    T62r= "|*ReD62 :                      | ImD62 :                      |\n";
    T62c= "|*ReD62 :                      |*ImD62 :                      |\n";
 
    T63n= "| ReD63 :                      | ImD63 :                      |\n";
    T63r= "|*ReD63 :                      | ImD63 :                      |\n";
    T63c= "|*ReD63 :                      |*ImD63 :                      |\n";
 
    T64n= "| ReD64 :                      | ImD64 :                      |\n";
    T64r= "|*ReD64 :                      | ImD64 :                      |\n";
    T64c= "|*ReD64 :                      |*ImD64 :                      |\n";
 
    T65n= "| ReD65 :                      | ImD65 :                      |\n";
    T65r= "|*ReD65 :                      | ImD65 :                      |\n";
    T65c= "|*ReD65 :                      |*ImD65 :                      |\n";
 
    T66n= "| ReD66 :                      | ImD66 :                      |\n";
    T66r= "|*ReD66 :                      | ImD66 :                      |\n";
    T66c= "|*ReD66 :                      |*ImD66 :                      |\n";
    t11 = "===============================================================\n";
 
    fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,t03,einheit_in);fprintf(fp,"%s",t04);
    fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
    fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
    fprintf(fp,"%s",t20);
    fprintf(fp,"%s",t21);
    fprintf(fp,"%s",t22);
    fprintf(fp,"%s",t23);
    fprintf(fp,"%s",t24);
    fprintf(fp,"%s",t25);
    fprintf(fp,"%s",t26);
    fprintf(fp,"%s",t27);
    fprintf(fp,"%s",t28);
    fprintf(fp,"%s",t29);
    fprintf(fp,"%s",t30);
    fprintf(fp,"%s",t31);
    fprintf(fp,"%s",t11);
 
    drucke_par( fp,modus,dimj,tabelle,einheit_out,temp,ion,symmetrienr );
    if(modus != NOMAG )  drucke_mag( fp,modus );
 
    printf("file %s produced.\n",name);
    fclose(fp);
 
}
/*------------------------------------------------------------------------------
                                create_Lkq()
------------------------------------------------------------------------------*/
void create_Lkq(einheitnr_in,einheitnr_out,ion,symmetrienr,modus,temp)
    INT    einheitnr_in,einheitnr_out;
    CHAR   *ion;
    INT    symmetrienr;
    CHAR   modus;
    DOUBLE temp; /* Temperatur der Probe */
{
 
    CHAR *name=LKQNAME;
    CHAR *einheit_in;
    CHAR *einheit_out;
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09;
    CHAR *t11,*t12,*t13,*t14,*t15,*t16,*t17,*t18,*t19,*t20;
    CHAR *t21,*t22,*t23,*t24,*t25,*t26,*t27,*t28,*t29,*t30;
    CHAR *t31,*t32,*t33,*t34,*t35,*t36,*t37,*t38,*t39,*t40;
    CHAR *t41,*t42,*t43,*t44,*t45,*t46,*t47,*t48,*t49,*t50;
    CHAR *t51,*t52,*t53,*t54;
    CHAR *leftcopy();
 
    FILE    *fopen(), *fp;
    TABELLE *tabelle;
    INT     dimj,ionennr;
 
    fp  = fopen_errchk(name,"w");
    write_title(fp);
 
    ionennr     = isimplementiert(ion);
    dimj        = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

    einheit_in  = EINHEITIMP[ einheitnr_in ].einheit;
    einheit_out = EINHEITIMP[ einheitnr_out].einheit;
    ion         = leftcopy(ion ,25);
 
    tabelle = TABELLE_ALLOC(1);
    t01 = "===============================================================\n";
    t02 = "|                                                             |\n";
    t03 = "|Crystal Field paramter L    in  %6s    in which           |\n";
    t04 = "|                        kq                                   |\n";
    t05 = "|                                                             |\n";
    t06 = "|        ---            ---                                   |\n";
    t07 = "| H   =  >   L   Z   +  >    L   Z   +   L     Z              |\n";
    t08 = "|  KF    ---  k0  k0    ---   kq  kq      k,-q  k,-q          |\n";
    t09 = "|         k             k>0                                   |\n";
    t51 = "|                       q>0                                   |\n";
    t12 = "|                                                             |\n";
    t13 = "|                                                             |\n";
    t14 = "|   The L  - are defined by :                                 |\n";
    t15 = "|        kq       kq                                          |\n";
    t16 = "|                                                             |\n";
    t17 = "|                                                             |\n";
    t18 = "|          /                             q                    |\n";
    t19 = "|         |  L      :=  Re( D      ) (-1)                     |\n";
    t20 = "|         |   k,|q|          k,|q|                            |\n";
    t21 = "|         |                                                   |\n";
    t22 = "|    L  = |  L      :=      D                                 |\n";
    t23 = "|     kq  |   k,0            k,0                              |\n";
    t24 = "|         |                              q                    |\n";
    t25 = "|         |  L      :=  Im( D      ) (-1)                     |\n";
    t26 = "|         |   k,-|q|         k,|q|                            |\n";
    t27 = "|          \\                                 .                |\n";
    t28 = "|                                                             |\n";
    t29 = "|  The L-Parameter are therefore real.                        |\n";
    t30 = "|                                                             |\n";
    t31 = "|  The Tensor operator  Z   , Which are the quantised         |\n";
    t32 = "|                         kq                                  |\n";
    t33 = "|  'tesseral harmonics', are given by                         |\n";
    t34 = "|                                                             |\n";
    t35 = "|                                                             |\n";
    t36 = "|          /                             q                    |\n";
    t37 = "|         |  Z      = 2 Re( C      ) (-1)                     |\n";
    t38 = "|         |   k,|q|          k,|q|                            |\n";
    t39 = "|         |                                                   |\n";
    t40 = "|    Z  = |  Z      =       C                                 |\n";
    t41 = "|     kq  |   k,0            k,0                              |\n";
    t42 = "|         |                              q                    |\n";
    t43 = "|         |  Z      = 2 Im( C      ) (-1)                     |\n";
    t44 = "|         |   k,-|q|         k,|q|                            |\n";
    t45 = "|          \\                                                  |\n";
    t46 = "|                                                             |\n";
    t47 = "|                                        4*pi                 |\n";
    t48 = "|  with the Tensor operators  C   = sqrt(------) Y   .        |\n";
    t49 = "|                             kq         2k+1    kq           |\n";
    t50 = "|                                                             |\n";
    t52 = "|  Y   =  quantised spherical harmonics                       |\n";
    t53 = "|   kq                                                        |\n";
    t54 = "|                                                             |\n";
    t11 = "===============================================================\n";
    TSS = "| Energy eigenvalues are in  : %6s                         |\n";
    T15 = "| Temperature of the sample  : %7.2f Kelvin                 |\n";
    T12 = "| Ion                        : %25s    |\n";
    T13 = "| Symmetry  : %s            Symmetry number: %1d        |\n";
    T14 = "| Magnetic field             : %s                 |\n";
    T11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T20n= "| L 2, 0:                                                     |\n";
    T20r= "|*L 2, 0:                                                     |\n";
    STR = "|-------------------------------------------------------------|\n";
    T21n= "| L 2, 1:                      | L 2,-1:                      |\n";
    T21r= "|*L 2, 1:                      | L 2,-1:                      |\n";
    T21c= "|*L 2, 1:                      |*L 2,-1:                      |\n";
 
    T22n= "| L 2, 2:                      | L 2,-2:                      |\n";
    T22r= "|*L 2, 2:                      | L 2,-2:                      |\n";
    T22c= "|*L 2, 2:                      |*L 2,-2:                      |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T40n= "| L 4, 0:                                                     |\n";
    T40r= "|*L 4, 0:                                                     |\n";
 
    T41n= "| L 4, 1:                      | L 4,-1:                      |\n";
    T41r= "|*L 4, 1:                      | L 4,-1:                      |\n";
    T41c= "|*L 4, 1:                      |*L 4,-1:                      |\n";
 
    T42n= "| L 4, 2:                      | L 4,-2:                      |\n";
    T42r= "|*L 4, 2:                      | L 4,-2:                      |\n";
    T42c= "|*L 4, 2:                      |*L 4,-2:                      |\n";
 
    T43n= "| L 4, 3:                      | L 4,-3:                      |\n";
    T43r= "|*L 4, 3:                      | L 4,-3:                      |\n";
    T43c= "|*L 4, 3:                      |*L 4,-3:                      |\n";
 
    T44n= "| L 4, 4:                      | L 4,-4:                      |\n";
    T44r= "|*L 4, 4:                      | L 4,-4:                      |\n";
    T44c= "|*L 4, 4:                      |*L 4,-4:                      |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T60n= "| L 6,0 :                                                     |\n";
    T60r= "|*L 6,0 :                                                     |\n";
 
    T61n= "| L 6, 1:                      | L 6,-1:                      |\n";
    T61r= "|*L 6, 1:                      | L 6,-1:                      |\n";
    T61c= "|*L 6, 1:                      |*L 6,-1:                      |\n";
 
    T62n= "| L 6, 2:                      | L 6,-2:                      |\n";
    T62r= "|*L 6, 2:                      | L 6,-2:                      |\n";
    T62c= "|*L 6, 2:                      |*L 6,-2:                      |\n";
 
    T63n= "| L 6, 3:                      | L 6,-3:                      |\n";
    T63r= "|*L 6, 3:                      | L 6,-3:                      |\n";
    T63c= "|*L 6, 3:                      |*L 6,-3:                      |\n";
 
    T64n= "| L 6, 4:                      | L 6,-4:                      |\n";
    T64r= "|*L 6, 4:                      | L 6,-4:                      |\n";
    T64c= "|*L 6, 4:                      |*L 6,-4:                      |\n";
 
    T65n= "| L 6, 5:                      | L 6,-5:                      |\n";
    T65r= "|*L 6, 5:                      | L 6,-5:                      |\n";
    T65c= "|*L 6, 5:                      |*L 6,-5:                      |\n";
 
    T66n= "| L 6, 6:                      | L 6,-6:                      |\n";
    T66r= "|*L 6, 6:                      | L 6,-6:                      |\n";
    T66c= "|*L 6, 6:                      |*L 6,-6:                      |\n";
    t11 = "===============================================================\n";
 
    fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,t03,einheit_in);fprintf(fp,"%s",t04);
    fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
    fprintf(fp,"%s",t09);fprintf(fp,"%s",t51);
 
    fprintf(fp,"%s",t12);fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);
    fprintf(fp,"%s",t15);fprintf(fp,"%s",t16);fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);
    fprintf(fp,"%s",t19);fprintf(fp,"%s",t20);
 
    fprintf(fp,"%s",t21);fprintf(fp,"%s",t22);fprintf(fp,"%s",t23);fprintf(fp,"%s",t24);
    fprintf(fp,"%s",t25);fprintf(fp,"%s",t26);fprintf(fp,"%s",t27);fprintf(fp,"%s",t28);
    fprintf(fp,"%s",t29);fprintf(fp,"%s",t30);
 
    fprintf(fp,"%s",t31);fprintf(fp,"%s",t32);fprintf(fp,"%s",t33);fprintf(fp,"%s",t34);
    fprintf(fp,"%s",t35);fprintf(fp,"%s",t36);fprintf(fp,"%s",t37);fprintf(fp,"%s",t38);
    fprintf(fp,"%s",t39);fprintf(fp,"%s",t40);
 
    fprintf(fp,"%s",t41);fprintf(fp,"%s",t42);fprintf(fp,"%s",t43);fprintf(fp,"%s",t44);
    fprintf(fp,"%s",t45);fprintf(fp,"%s",t46);fprintf(fp,"%s",t47);fprintf(fp,"%s",t48);
    fprintf(fp,"%s",t49);fprintf(fp,"%s",t50);
 
    fprintf(fp,"%s",t52);fprintf(fp,"%s",t53);fprintf(fp,"%s",t54);
    fprintf(fp,"%s",t11);
 
    drucke_par( fp,modus,dimj,tabelle,einheit_out,temp,ion,symmetrienr );
    if(modus != NOMAG )  drucke_mag( fp,modus );
 
    printf("file %s produced.\n",name);
    fclose(fp);
 
}
/*------------------------------------------------------------------------------
                                create_Wkq()
------------------------------------------------------------------------------*/
void create_Wkq(einheitnr_in,einheitnr_out,ion,symmetrienr,modus,temp)
    INT    einheitnr_in,einheitnr_out;
    CHAR   *ion;
    INT    symmetrienr;
    CHAR   modus;
    DOUBLE temp; /* Temperatur der Probe */
{
 
    CHAR *name=WKQNAME;
    CHAR *einheit_in;
    CHAR *einheit_out;
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t71,*t72,*t73,*t74,*t75,*t76,*t77,*t78,*t79,*t70;
    CHAR *t81,*t82,*t83,*t84,*t85,*t86,*t87,*t88,*t80;
    CHAR *t11;
    CHAR *leftcopy();
 
    FILE    *fopen(), *fp;
    TABELLE *tabelle;
    INT     ionennr,dimj;
 
 
    fp  = fopen_errchk(name,"w");
    write_title(fp);
 
    ionennr     = isimplementiert(ion);
    dimj        = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

    einheit_in  = EINHEITIMP[ einheitnr_in ].einheit;
    einheit_out = EINHEITIMP[ einheitnr_out].einheit;
    ion         = leftcopy(ion ,25);
 
    tabelle = TABELLE_ALLOC(1);
    t01 = "===============================================================\n";
    t02 = "|                                    o**k                     |\n";
    t03 = "|Crystal Field Parameter W  in   %6s/A    in which         |\n";
    t04 = "|                         kq                                  |\n";
    t05 = "|                                                             |\n";
    t06 = "|        ---                ---                  *   +        |\n";
    t07 = "| H   =  >   W   R  (J)  +  >    W   R  (J)  +  W   R  (J)    |\n";
    t08 = "|  KF    ---  k0  k0        ---   kq  kq         kq  kq       |\n";
    t09 = "|         k                 k>0                               |\n";
    t10 = "|                           q>0                               |\n";
    t70 = "|                                                             |\n";
    t71 = "|                                                             |\n";
    t72 = "|                           with                              |\n";
    t73 = "|                                                             |\n";
    t74 = "|                                                             |\n";
    t75 = "|                          2S+1           k                   |\n";
    t76 = "|          R   :=  theta (     L   ) * < r  > * P  (J)        |\n";
    t77 = "|           kq          k        J               kq           |\n";
    t78 = "|                                                             |\n";
    t79 = "|                                                             |\n";
    t80 = "| remark : The W 's are the complex analogies of the  A 's .  |\n";
    t81 = "|              kq                                     kq      |\n";
    t82 = "|                                                             |\n";
    t83 = "|                            also                             |\n";
    t84 = "|                                                             |\n";
    t85 = "|                          2S+1          k                    |\n";
    t86 = "|            W   *  theta (    L  ) * < r  > =  V             |\n";
    t87 = "|             kq         k      J                kq           |\n";
    t88 = "|                                                             |\n";
    t11 = "===============================================================\n";
    TSS = "| Energy Eigenvalues are in  : %6s                         |\n";
    T15 = "| Temperature of the probe   : %7.2f Kelvin                 |\n";
    T12 = "| Ion                        : %25s     |\n";
    T13 = "| Symmetry         : %s     Symmetry number : %1d       |\n";
    T14 = "| Magnetic field             : %s                 |\n";
    T11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T20n= "| ReW20 :                                                     |\n";
    T20r= "|*ReW20 :                                                     |\n";
    STR = "|-------------------------------------------------------------|\n";
    T21n= "| ReW21 :                      | ImW21 :                      |\n";
    T21r= "|*ReW21 :                      | ImW21 :                      |\n";
    T21c= "|*ReW21 :                      |*ImW21 :                      |\n";
 
    T22n= "| ReW22 :                      | ImW22 :                      |\n";
    T22r= "|*ReW22 :                      | ImW22 :                      |\n";
    T22c= "|*ReW22 :                      |*ImW22 :                      |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T40n= "| ReW40 :                                                     |\n";
    T40r= "|*ReW40 :                                                     |\n";
 
    T41n= "| ReW41 :                      | ImW41 :                      |\n";
    T41r= "|*ReW41 :                      | ImW41 :                      |\n";
    T41c= "|*ReW41 :                      |*ImW41 :                      |\n";
 
    T42n= "| ReW42 :                      | ImW42 :                      |\n";
    T42r= "|*ReW42 :                      | ImW42 :                      |\n";
    T42c= "|*ReW42 :                      |*ImW42 :                      |\n";
 
    T43n= "| ReW43 :                      | ImW43 :                      |\n";
    T43r= "|*ReW43 :                      | ImW43 :                      |\n";
    T43c= "|*ReW43 :                      |*ImW43 :                      |\n";
 
    T44n= "| ReW44 :                      | ImW44 :                      |\n";
    T44r= "|*ReW44 :                      | ImW44 :                      |\n";
    T44c= "|*ReW44 :                      |*ImW44 :                      |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T60n= "| ReW60 :                                                     |\n";
    T60r= "|*ReW60 :                                                     |\n";
 
    T61n= "| ReW61 :                      | ImW61 :                      |\n";
    T61r= "|*ReW61 :                      | ImW61 :                      |\n";
    T61c= "|*ReW61 :                      |*ImW61 :                      |\n";
 
    T62n= "| ReW62 :                      | ImW62 :                      |\n";
    T62r= "|*ReW62 :                      | ImW62 :                      |\n";
    T62c= "|*ReW62 :                      |*ImW62 :                      |\n";
 
    T63n= "| ReW63 :                      | ImW63 :                      |\n";
    T63r= "|*ReW63 :                      | ImW63 :                      |\n";
    T63c= "|*ReW63 :                      |*ImW63 :                      |\n";
 
    T64n= "| ReW64 :                      | ImW64 :                      |\n";
    T64r= "|*ReW64 :                      | ImW64 :                      |\n";
    T64c= "|*ReW64 :                      |*ImW64 :                      |\n";
 
    T65n= "| ReW65 :                      | ImW65 :                      |\n";
    T65r= "|*ReW65 :                      | ImW65 :                      |\n";
    T65c= "|*ReW65 :                      |*ImW65 :                      |\n";
 
    T66n= "| ReW66 :                      | ImW66 :                      |\n";
    T66r= "|*ReW66 :                      | ImW66 :                      |\n";
    T66c= "|*ReW66 :                      |*ImW66 :                      |\n";
    t11 = "===============================================================\n";
 
    fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,t03,einheit_in);fprintf(fp,"%s",t04);
    fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
    fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
    fprintf(fp,"%s",t70);fprintf(fp,"%s",t71);fprintf(fp,"%s",t72);fprintf(fp,"%s",t73);
    fprintf(fp,"%s",t74);fprintf(fp,"%s",t75);fprintf(fp,"%s",t76);fprintf(fp,"%s",t77);
    fprintf(fp,"%s",t78);fprintf(fp,"%s",t79);fprintf(fp,"%s",t80);fprintf(fp,"%s",t81);
    fprintf(fp,"%s",t82);fprintf(fp,"%s",t83);fprintf(fp,"%s",t84);fprintf(fp,"%s",t85);
    fprintf(fp,"%s",t86);fprintf(fp,"%s",t87);fprintf(fp,"%s",t88);
    fprintf(fp,"%s",t11);
 
    drucke_par( fp,modus,dimj,tabelle,einheit_out,temp,ion,symmetrienr );
    if(modus != NOMAG )  drucke_mag( fp,modus );
 
    printf("file %s produced.\n",name);
    fclose(fp);
 
}
/*------------------------------------------------------------------------------
                                create_Akq()
------------------------------------------------------------------------------*/
void create_Akq(einheitnr_in,einheitnr_out,ion,symmetrienr,modus,temp)
    INT    einheitnr_in,einheitnr_out;
    CHAR   *ion;
    INT    symmetrienr;
    CHAR   modus;
    DOUBLE temp; /* Temperatur der Probe */
{
 
    CHAR *name=AKQNAME;
    CHAR *einheit_in;
    CHAR *einheit_out;
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t31,*t32,*t33,*t34,*t35,*t36,*t37,*t38,*t39,*t30;
    CHAR *t51,*t52,*t53,*t54,*t55,*t56,*t50;
    CHAR *t71,*t72,*t73,*t74,*t75,*t76,*t70;
    CHAR *t11;
    CHAR *leftcopy();
 
    FILE    *fopen(), *fp;
    TABELLE *tabelle;
    INT     ionennr,dimj;
 
 
    fp  = fopen_errchk(name,"w");
    write_title(fp);
 
    ionennr     = isimplementiert(ion);
    dimj        = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

    einheit_in  = EINHEITIMP[ einheitnr_in ].einheit;
    einheit_out = EINHEITIMP[ einheitnr_out].einheit;
    ion         = leftcopy(ion ,25);
 
    tabelle = TABELLE_ALLOC(1);
    t01 = "===============================================================\n";
    t02 = "|                                    **k                      |\n";
    t03 = "|Crystal Field parameter A  in   %6s/a0 and real           |\n";
    t04 = "|                        kq (compare Hutchings p255 eq 5.5)   |\n";
    t05 = "|-----------------------------------------------              |\n";
    t06 = "|                                              |              |\n";
    t07 = "|        ---        2S+1       k               |              |\n";
    t08 = "| H   =  >   theta (    L ) < r > A  STEV (J)  |              |\n";
    t09 = "|  KF    ---      k      J         kq    kq    |              |\n";
    t10 = "|        k>=0                                  |              |\n";
    t30 = "|        q>=0                                  |              |\n";
    t31 = "|-----------------------------------------------              |\n";
    t32 = "|                                                             |\n";
    t33 = "| q=0 :  STEV (J) := P  (J)                                   |\n";
    t34 = "|            k0       k0                                      |\n";
    t35 = "|                                                             |\n";
    t36 = "|                                +                            |\n";
    t37 = "| q>0 :  STEV (J) := ( P  (J) + P  (J) )/2/omega    ,V  real  |\n";
    t38 = "|            kq         kq       kq             kq    kq      |\n";
    t39 = "|            s                   +                            |\n";
    t50 = "| q>0 :  STEV (J) := ( P  (J) - P  (J) )/2/i/omega  ,V        |\n";
    t51 = "|            kq         kq       kq               kq  kq imag.|\n";
    t52 = "|-------------------------------------------------------------|\n";
    t53 = "|                 2S+1          k                             |\n";
    t54 = "| And    theta (     L  ) * < r  > * A    =   B              |\n";
    t55 = "|              k       J               kq       kq            |\n";
    t56 = "|-------------------------------------------------------------|\n";
    t70 = "| And    q=0 :  B        =  V                  | The left of  |\n";
    t71 = "|                k0          k0                | the equation |\n";
    t72 = "|        q>0 :    B       =  2 * omega  * V     | are the      |\n";
    t73 = "|               s  kq                 kq   kq   | magnitudes that|\n";
    t74 = "|        q>0 : B   / i  =  2 * omega  * V     | are real and |\n";
    t75 = "|               kq                  kq   kq   | listed       |\n";
    t76 = "|                                             | below.       |\n";
    t11 = "===============================================================\n";
    TSS = "| Energy Eigenvalues are in  : %6s                         |\n";
    T15 = "| Temperature of the probe   : %7.2f Kelvin                 |\n";
    T12 = "| Ion                        : %25s     |\n";
    T13 = "| Symmetry         : %s     Symmetry number : %1d       |\n";
    T14 = "| Magnetic field             : %s                 |\n";
    T11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T20n= "|  A20 R:                                                     |\n";
    T20r= "|* A20 R:                                                     |\n";
    STR = "|-------------------------------------------------------------|\n";
    T21n= "|  A21 R:                                                     |\n";
    T21r= "|* A21 R:                                                     |\n";
    T21c= "|* A21 R:                                                     |\n";
 
    T22n= "|  A22 R:                                                     |\n";
    T22r= "|* A22 R:                                                     |\n";
    T22c= "|* A22 R:                                                     |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T40n= "|  A40 R:                                                     |\n";
    T40r= "|* A40 R:                                                     |\n";
 
    T41n= "|  A41 R:                                                     |\n";
    T41r= "|* A41 R:                                                     |\n";
    T41c= "|* A41 R:                                                     |\n";
 
    T42n= "|  A42 R:                                                     |\n";
    T42r= "|* A42 R:                                                     |\n";
    T42c= "|* A42 R:                                                     |\n";
 
    T43n= "|  A43 R:                                                     |\n";
    T43r= "|* A43 R:                                                     |\n";
    T43c= "|* A43 R:                                                     |\n";
 
    T44n= "|  A44 R:                                                     |\n";
    T44r= "|* A44 R:                                                     |\n";
    T44c= "|* A44 R:                                                     |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T60n= "|  A60 R:                                                     |\n";
    T60r= "|* A60 R:                                                     |\n";
 
    T61n= "|  A61 R:                                                     |\n";
    T61r= "|* A61 R:                                                     |\n";
    T61c= "|* A61 R:                                                     |\n";
 
    T62n= "|  A62 R:                                                     |\n";
    T62r= "|* A62 R:                                                     |\n";
    T62c= "|* A62 R:                                                     |\n";
 
    T63n= "|  A63 R:                                                     |\n";
    T63r= "|* A63 R:                                                     |\n";
    T63c= "|* A63 R:                                                     |\n";
 
    T64n= "|  A64 R:                                                     |\n";
    T64r= "|* A64 R:                                                     |\n";
    T64c= "|* A64 R:                                                     |\n";
 
    T65n= "|  A65 R:                                                     |\n";
    T65r= "|* A65 R:                                                     |\n";
    T65c= "|* A65 R:                                                     |\n";
 
    T66n= "|  A66 R:                                                     |\n";
    T66r= "|* A66 R:                                                     |\n";
    T66c= "|* A66 R:                                                     |\n";
    t11 = "===============================================================\n";
 
    fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,t03,einheit_in);fprintf(fp,"%s",t04);
    fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
    fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
    fprintf(fp,"%s",t30);fprintf(fp,"%s",t31);fprintf(fp,"%s",t32);fprintf(fp,"%s",t33);
    fprintf(fp,"%s",t34);fprintf(fp,"%s",t35);fprintf(fp,"%s",t36);fprintf(fp,"%s",t37);
    fprintf(fp,"%s",t38);fprintf(fp,"%s",t39);
    fprintf(fp,"%s",t50);fprintf(fp,"%s",t51);fprintf(fp,"%s",t52);fprintf(fp,"%s",t53);
    fprintf(fp,"%s",t54);fprintf(fp,"%s",t55);fprintf(fp,"%s",t56);
    fprintf(fp,"%s",t70);fprintf(fp,"%s",t71);fprintf(fp,"%s",t72);fprintf(fp,"%s",t73);
    fprintf(fp,"%s",t74);fprintf(fp,"%s",t75);fprintf(fp,"%s",t76);
    fprintf(fp,"%s",t11);
 
    drucke_par( fp,modus,dimj,tabelle,einheit_out,temp,ion,symmetrienr );
    if(modus != NOMAG )  drucke_mag( fp,modus );
 
    printf("file %s produced.\n",name);
    fclose(fp);
 
 
}
/*------------------------------------------------------------------------------
                                create_Bkq()
------------------------------------------------------------------------------*/
void create_Bkq(einheitnr_in,einheitnr_out,ion,symmetrienr,modus,temp)
    INT    einheitnr_in,einheitnr_out;
    CHAR   *ion;
    INT    symmetrienr;
    CHAR   modus;
    DOUBLE temp; /* Temperatur der Probe */
{
 
    CHAR *name=BKQNAME;
    CHAR *einheit_in;
    CHAR *einheit_out;
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t31,*t32,*t33,*t34,*t35,*t36,*t37,*t38,*t39,*t30;
    CHAR *t51,*t52,*t53,*t54,*t55,*t56,*t50;
    CHAR *t11;
    CHAR *leftcopy();
 
    FILE    *fopen(), *fp;
    TABELLE *tabelle;
    INT     ionennr,dimj;
 
 
    fp  = fopen_errchk(name,"w");
    write_title(fp);
 
    ionennr     = isimplementiert(ion);
    dimj        = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

    einheit_in  = EINHEITIMP[ einheitnr_in ].einheit;
    einheit_out = EINHEITIMP[ einheitnr_out].einheit;
    ion         = leftcopy(ion ,25);
 
    tabelle = TABELLE_ALLOC(1);
    t01 = "===============================================================\n";
    t02 = "|                                                             |\n";
    t03 = "|Crystal Field parameter B  in   %6s  and real             |\n";
    t04 = "|                        kq   (compare Hutchings)             |\n";
    t05 = "|-----------------------------                                |\n";
    t06 = "|        ---    x            |                                |\n";
    t07 = "| H   =  >     B   STEV (J)  |  with  STEV (J)  hermitesch    |\n";
    t08 = "|  KF    ---    kq     kq    |           kq                   |\n";
    t09 = "|        k>=0                |                                |\n";
    t10 = "|        q>=0,x=c,s          |                                |\n";
    t30 = "|-----------------------------                                |\n";
    t31 = "|                                                             |\n";
    t32 = "| q=0 :  STEV (J) := P  (J)                                   |\n";
    t33 = "|            k0       k0                                      |\n";
    t34 = "|            c                   +                            |\n";
    t35 = "| q>0 :  STEV (J) := ( P  (J) + P  (J) )/2/omega    ,V   real |\n";
    t36 = "|            kq         kq       kq             kq    kq      |\n";
    t37 = "|            s                   +                            |\n";
    t38 = "| q>0 :  STEV (J) := ( P  (J) - P  (J) )/2/i/omega  ,V        |\n";
    t39 = "|            kq         kq       kq               kq  kq imag.|\n";
    t50 = "|                                             | The left of   |\n";
    t51 = "| und    q=0 :  B        =  V                 | the equation  |\n";
    t52 = "|                k0 c        k0               | are the       |\n";
    t53 = "|        q>0 :     B     =  2 * omega re(V  ) | magnitudesthat|\n";
    t54 = "|                s  kq               kq   kq  | are real and  |\n";
    t55 = "|        q>0 :  B        =  -2* omega*im(V  ) | are listed    |\n";
    t56 = "|                kq                  kq   kq  | below         |\n";
    t11 = "===============================================================\n";
    TSS = "| Energy Eigenvalues are in  : %6s                         |\n";
    T15 = "| Temperature of the probe   : %7.2f Kelvin                 |\n";
    T12 = "| Ion                        : %25s     |\n";
    T13 = "| Symmetry         : %s     Symmetry number : %1d       |\n";
    T14 = "| Magnetic field             : %s                 |\n";
    T11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T20n= "|  B20 c:                                                     |\n";
    T20r= "|* B20 c:                                                     |\n";
    STR = "|-------------------------------------------------------------|\n";
    T21n= "|  B21 c:                      | B21 s:                       |\n";
    T21r= "|* B21 C:                      | B21 s:                       |\n";
    T21c= "|* B21 C:                      |*B21 s:                       |\n";
 
    T22n= "|  B22 C:                      | B22 s:                       |\n";
    T22r= "|* B22 C:                      | B22 s:                       |\n";
    T22c= "|* B22 C:                      |*B22 s:                       |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T40n= "|  B40 C:                                                     |\n";
    T40r= "|* B40 C:                                                     |\n";
 
    T41n= "|  B41 C:                      | B41 s:                       |\n";
    T41r= "|* B41 C:                      | B41 s:                       |\n";
    T41c= "|* B41 C:                      |*B41 s:                       |\n";
 
    T42n= "|  B42 C:                      | B42 s:                       |\n";
    T42r= "|* B42 C:                      | B42 s:                       |\n";
    T42c= "|* B42 C:                      |*B42 s:                       |\n";
 
    T43n= "|  B43 C:                      | B43 s:                       |\n";
    T43r= "|* B43 C:                      | B43 s:                       |\n";
    T43c= "|* B43 C:                      |*B43 s:                       |\n";
 
    T44n= "|  B44 C:                      | B44 s:                       |\n";
    T44r= "|* B44 C:                      | B44 s:                       |\n";
    T44c= "|* B44 C:                      |*B44 s:                       |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    T60n= "|  B60 C:                                                     |\n";
    T60r= "|* B60 C:                                                     |\n";
 
    T61n= "|  B61 C:                      | B61 s:                       |\n";
    T61r= "|* B61 C:                      | B61 s:                       |\n";
    T61c= "|* B61 C:                      |*B61 s:                       |\n";
 
    T62n= "|  B62 C:                      | B62 s:                       |\n";
    T62r= "|* B62 C:                      | B62 s:                       |\n";
    T62c= "|* B62 C:                      |*B62 s:                       |\n";
 
    T63n= "|  B63 C:                      | B63 s:                       |\n";
    T63r= "|* B63 C:                      | B63 s:                       |\n";
    T63c= "|* B63 C:                      |*B63 s:                       |\n";
 
    T64n= "|  B64 C:                      | B64 s:                       |\n";
    T64r= "|* B64 C:                      | B64 s:                       |\n";
    T64c= "|* B64 C:                      |*B64 s:                       |\n";
 
    T65n= "|  B65 C:                      | B65 s:                       |\n";
    T65r= "|* B65 C:                      | B65 s:                       |\n";
    T65c= "|* B65 C:                      |*B65 s:                       |\n";
 
    T66n= "|  B66 C:                      | B66 s:                       |\n";
    T66r= "|* B66 C:                      | B66 s:                       |\n";
    T66c= "|* B66 C:                      |*B66 s:                       |\n";
    t11 = "===============================================================\n";
 
    fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,t03,einheit_in);fprintf(fp,"%s",t04);
    fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
    fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
    fprintf(fp,"%s",t30);fprintf(fp,"%s",t31);fprintf(fp,"%s",t32);fprintf(fp,"%s",t33);
    fprintf(fp,"%s",t34);fprintf(fp,"%s",t35);fprintf(fp,"%s",t36);fprintf(fp,"%s",t37);
    fprintf(fp,"%s",t38);fprintf(fp,"%s",t39);
    fprintf(fp,"%s",t50);fprintf(fp,"%s",t51);fprintf(fp,"%s",t52);fprintf(fp,"%s",t53);
    fprintf(fp,"%s",t54);fprintf(fp,"%s",t55);fprintf(fp,"%s",t56);
    fprintf(fp,"%s",t11);
 
    drucke_par( fp,modus,dimj,tabelle,einheit_out,temp,ion,symmetrienr );
    if(modus != NOMAG )  drucke_mag( fp,modus );
 
    printf("file %s produced.\n",name);
    fclose(fp);
 
 
}
/*------------------------------------------------------------------------------
                                create_xW()
------------------------------------------------------------------------------*/
void create_xW(einheitnr_in,einheitnr_out,ion,symmetrienr,modus,temp)
    INT    einheitnr_in,einheitnr_out;
    CHAR   *ion;
    INT    symmetrienr;
    CHAR   modus;
    DOUBLE temp; /* Temperatur der Probe */
{
 
    CHAR *name=XWNAME;
    CHAR *einheit_in;
    CHAR *einheit_out;
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11,*t12,*t13,*t14,*t15,*t16,*t17,*t18,*t19,*t20;
    CHAR *t21,*t22,*t23,*t24,*t25,*t26,*t27,*t28,*t29,*t30;
    CHAR *t31,*t32,*t33,*t34,*t35,*t36,*t37,*t38,*t39,*t40;
    CHAR *t41,*t42,*t43,*t44,*t45,*t46,*t47,*t48,*t49,*t50;
    CHAR *t51,*t52,*t53;
    CHAR *tx,*tW,*str,*tss;
    CHAR *leftcopy();
 
    FILE *fopen(), *fp;
    INT  ionennr,dimj;
 
 
    fp  = fopen_errchk(name,"w");
    write_title(fp);
 
    ionennr     = isimplementiert(ion);
    dimj        = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

    einheit_in  = EINHEITIMP[ einheitnr_in ].einheit;
    einheit_out = EINHEITIMP[ einheitnr_out].einheit;
    ion         = leftcopy(ion ,25);
 
    t01 = "===============================================================\n";
    t02 = "|              Crystal Field parameter  x,W (real)            |\n";
    t03 = "|                            in  %6s                          |\n";
    t04 = "|                    for cubic Symmetry                       |\n";
    t05 = "|                                                             |\n";
    t06 = "|                                                             |\n";
    t07 = "| Starting with                                               |\n";
    t08 = "|         ---------------------------------------             |\n";
    t09 = "|        |                                       |            |\n";
    t10 = "|        | H   =  B  *  O (J)   +   B   *  O (J) |            |\n";
    t11 = "|        |  KF     40    4           60     6    |            |\n";
    t12 = "|        |                                       |            |\n";
    t13 = "|         ---------------------------------------             |\n";
    t14 = "|         ---------------------------------------             |\n";
    t15 = "|        |                                       |            |\n";
    t16 = "|Where   |  O (J) :=  STEV (J)  +  5 * STEV (J)  |            |\n";
    t17 = "|        |   4            40               44    |            |\n";
    t18 = "|        |                                       |            |\n";
    t19 = "|        |  O (J) :=  STEV (J)  - 21 * STEV (J)  |            |\n";
    t20 = "|        |   6            60               64    |            |\n";
    t21 = "|         ---------------------------------------             |\n";
    t22 = "|                              -----------------------------  |\n";
    t23 = "|setting      ( -1<=x=<1 )    |       x * W = B  * F(4)     | |\n";
    t24 = "|                             |                40           | |\n";
    t25 = "|                             |                             | |\n";
    t26 = "|and       B                  | (1-|x|) * W = B  * F(6)     | |\n";
    t27 = "|           60                |                60           | |\n";
    t28 = "|                              -----------------------------  |\n";
    t29 = "|                                                             |\n";
    t30 = "|and     F(k) := F (J) := GGT(  { <JM'| STEV  |JM> }        ) |\n";
    t31 = "|                 k0                        k0     M.M'=-J..J |\n";
    t32 = "|                                                             |\n";
    t33 = "|In the case of J=9/2 one puts F(4)=  60  in stead of  84     |\n";
    t34 = "|                             F(6)=2520 in stead of 5040      |\n";
    t35 = "|in case of J= 8  one puts F(4)=  60 instead of  420          |\n";
    t36 = "|                                                             |\n";
    t37 = "|                                                             |\n";
    t38 = "|note :- the <JM'|STEV  |JM> are integers.                    |\n";
    t39 = "|                        k0                                   |\n";
    t40 = "|                                                             |\n";
    t41 = "|         - The B parameters contain the quantisation         |\n";
    t42 = "|                                2S+1                         |\n";
    t43 = "|           parameters   theta (     L   )  .  therfore       |\n";
    t44 = "|                             k        J                      |\n";
    t45 = "|                                                             |\n";
    t46 = "|           the x,W-Parameter and B-Parameter  follow         |\n";
    t47 = "|           th same selection rules regarding the             |\n";
    t48 = "|           quantisation parameters .                         |\n";
    t49 = "|           Especially for the ground state                   |\n";
    t50 = "|                                                             |\n";
    t51 = "|           2             3+                         O4(5/2)  |\n";
    t52 = "|            F    from  Ce   : B  = 0  =>  H  = x*W*-------   |\n";
    t53 = "|             5/2               60           kf       F(4)    |\n";
    fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,t03,einheit_in);fprintf(fp,"%s",t04);
    fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
    fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
    fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,"%s",t15);fprintf(fp,"%s",t16);
    fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);fprintf(fp,"%s",t19);fprintf(fp,"%s",t20);
    fprintf(fp,"%s",t21);fprintf(fp,"%s",t22);fprintf(fp,"%s",t23);fprintf(fp,"%s",t24);
    fprintf(fp,"%s",t25);fprintf(fp,"%s",t26);fprintf(fp,"%s",t27);fprintf(fp,"%s",t28);
    fprintf(fp,"%s",t29);fprintf(fp,"%s",t30);fprintf(fp,"%s",t31);fprintf(fp,"%s",t32);
    fprintf(fp,"%s",t33);fprintf(fp,"%s",t34);fprintf(fp,"%s",t35);fprintf(fp,"%s",t36);
    fprintf(fp,"%s",t37);fprintf(fp,"%s",t38);fprintf(fp,"%s",t39);fprintf(fp,"%s",t40);
    fprintf(fp,"%s",t41);fprintf(fp,"%s",t42);fprintf(fp,"%s",t43);fprintf(fp,"%s",t44);
    fprintf(fp,"%s",t45);fprintf(fp,"%s",t46);fprintf(fp,"%s",t47);fprintf(fp,"%s",t48);
    fprintf(fp,"%s",t49);fprintf(fp,"%s",t50);fprintf(fp,"%s",t51);fprintf(fp,"%s",t52);
    fprintf(fp,"%s",t53);
 
    t11 = "===============================================================\n";
    tss = "| Energy Eigenvalues are in  : %6s                         |\n";
    t15 = "| Temperature of the probe   : %7.2f Kelvin                 |\n";
    t12 = "| Ion                        : %25s     |\n";
    t13 = "| Symmetry         : %s     Symmetry number : %1d       |\n";
    t14 = "| Magnetic field             : %s                 |\n";
    t11 = "===============================================================\n";
    t11 = "===============================================================\n";
    tx  = "|   x  :                                                      |\n";
    str = "|-------------------------------------------------------------|\n";
    tW  = "|   W  :                                                      |\n";
    t11 = "===============================================================\n";
 
 
    fprintf(fp,"%s",t11);
    fprintf(fp,tss,einheit_out);
    fprintf(fp,t15,temp);
    fprintf(fp,t12,ion);
    fprintf(fp,t13,SYMLISTE[symmetrienr].symname,symmetrienr);
 
    switch(modus){
 
       case NOMAG :  fprintf(fp,t14,NICHTANGELEGT);
                     break;
       default    :  fprintf(fp,t14,ANGELEGT);
    }
    fprintf(fp,"%s",t11);
 
 
    fprintf(fp,"%s",t11);
    fprintf(fp,"%s",tx);
    fprintf(fp,"%s",str);
    fprintf(fp,"%s",tW);
    fprintf(fp,"%s",t11);
 
 
    if(modus != NOMAG )  drucke_mag( fp,modus );
 
    printf("file %s produced.\n",name);
    fclose(fp);
 
}
 
/*------------------------------------------------------------------------------
                                create_nn()
------------------------------------------------------------------------------*/
void create_nn(name,modus,nn,ion,temp) /* Erzeuge Eingabefile fuer Umgebungsatome */
    CHAR   *name;
    CHAR   modus;
    INT    nn;
    CHAR   *ion;
    DOUBLE temp;
{
    FILE *fopen(), *fp;
    CHAR *s_modus,*x1,*x2,*x3,*y1,*y2,*y3,*z1,*z2,*z3;
    CHAR *leftcopy();
    CHAR *rs,*rsnn1,*rsn,*rstk,*rsk,*rstl,*rstu;
    CHAR *po,*ponn1,*ponn,*potk,*pok,*potl,*potu;
    CHAR *rsm,*rsmt,*rsmtk,*rsmk,*rsmtu;
    CHAR *pom,*pomt,*pomtk,*pomk,*pomtu;
    INT  i;
 
    #define SN "%s\n"
 
 
rs   ="=======================================================================";
rsnn1="|  1    nearest neighbour                                            |";
rsn="| %3d   nearest neighbour                                             |\n";
rstk ="| Nr. |     q      | %8s | %8s | %8s | comment |\n";
rsk="| %3d |            |            |            |            |           |\n";
rstl ="------+------------+------------+------------+-------------------------";
rstu ="-----------------------------------------------------------------------";
 
rsm  = "====================================================";
rsmt = "| Magnetic field B                                     |";
rsmtk= "| %8s | %8s | %8s | comment |\n";
rsmk = "|            |            |            |           |";
rsmtu= "----------------------------------------------------";
 
po   = "==========================================================";
ponn1= "|  1    nearest neighbour                            |";
ponn = "| %3d   nearest neighbour                                |\n";
potk = "| Nr. |     q      | %8s | %8s | comment |\n";
pok  = "| %3d |            |            |            |           |\n";
potl = "------+------------+------------+------------+-----------|";
potu = "----------------------------------------------------------";
 
pom  = "=======================================";
pomt = "| Magnetic field B                        |";
pomtk= "| %8s | %8s | comment |\n";
pomk = "|            |            |           |";
pomtu= "---------------------------------------";
 
    ion = leftcopy(ion ,25);
 
    fp  = fopen_errchk(name,"w");
 
    switch(modus){
         case 'r': s_modus = "KARTESISCHEN";
                   x1      = "  x := x-koordinate   | in Angstroem";
                   x2      = "  y := y-koordinate   | in Angstroem";
                   x3      = "  z := z-koordinate   | in Angstroem";
                   y1      = "    x     ";
                   y2      = "    y     ";
                   y3      = "    z     ";
                   z1      = "    Bx    ";
                   z2      = "    By    ";
                   z3      = "    Bz    ";
                   break;
         case 's': s_modus = "SPHAERISCHEN";
                   x1      = "r     := Abstand      | in Angstroem";
                   x2      = "theta := Polarwinkel  | in Grad     ";
                   x3      = "phi   := Azimutwinkel | in Grad     ";
                   y1      = "    r     ";
                   y2      = "   phi    ";
                   y3      = "  theta   ";
                   z1      = "   |B|    ";
                   z2      = "   phi    ";
                   z3      = "  theta   ";
                   break;
         case 'p': s_modus = "POLAREN     ";
                   x1      = "r     := Abstand      | in Angstroem";
                   x2      = "phi   := Azimutwinkel | in Grad     ";
                   x3      = "                      |             ";
                   y1      = "    r     ";
                   y2      = "   phi    ";
                   y3      = "";
                   z1      = "   |B|    ";
                   z2      = "   phi    ";
                   z3      = "          ";
                   break;
         default : s_modus = "FEHLER!     ";
                   x1      = "                                    ";
                   x2      = "                                    ";
                   x3      = "                                    ";
                   y1      = "          ";
                   y2      = "          ";
                   y3      = "";
                   z1      = "          ";
                   z2      = "          ";
                   z3      = "          ";
    }
 
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"| Temperature of the sample       : %7.2f Kelvin          |\n",temp);
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"| the nearest neighbours of    : %25s |\n",ion);
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"| The co-ordinates of the surrounding ion are given in    |\n");
    fprintf(fp,"| %12s co-ordinates.                              |\n",s_modus);
    fprintf(fp,"|---------------------------------------------------------|\n");
    fprintf(fp,"|   B := Magnetic field | in Tesla                    |\n");
    fprintf(fp,"|                       |                             |\n");
    fprintf(fp,"|-----------------------+---------------------------------|\n");
    fprintf(fp,"|   q := charge of the  | in units of the elementary     |\n");
    fprintf(fp,"|        surrounding ion| charge |e|                      |\n");
    fprintf(fp,"|-----------------------+---------------------------------|\n");
    fprintf(fp,"| %36s                    |\n",x1);
    fprintf(fp,"| %36s                    |\n",x2);
    fprintf(fp,"| %36s                    |\n",x3);
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"\n");
 
    if( modus=='r' || modus=='s' ){ /* Tabelle in rechtw. oder sphaer. Koord.*/
 
         fprintf(fp,SN,rsm);
         fprintf(fp,SN,rsmt);
         fprintf(fp,SN,rsm);
         fprintf(fp,rsmtk,z1,z2,z3);
         fprintf(fp,SN,rsm);
         fprintf(fp,SN,rsmk);
         fprintf(fp,SN,rsmtu);
         fprintf(fp,"\n");
 
         fprintf(fp,SN,rs);
         if( nn==1 )
             fprintf(fp,SN,rsnn1);
         else
             fprintf(fp,rsn,nn);
         fprintf(fp,SN,rs);
         fprintf(fp,rstk,y1,y2,y3);
         fprintf(fp,SN,rs);
         for( i=1 ; i<=nn ; ++i ){
              fprintf(fp, rsk,i);
              if(i != nn )
                   fprintf(fp,SN,rstl);
              else
                   fprintf(fp,SN,rstu);
         }
    }
    else{ /* Tabelle in Polarkoordinaten */
 
         fprintf(fp,SN,pom);
         fprintf(fp,SN,pomt);
         fprintf(fp,SN,pom);
         fprintf(fp,pomtk,z1,z2);
         fprintf(fp,SN,pom);
         fprintf(fp,SN,pomk);
         fprintf(fp,SN,pomtu);
         fprintf(fp,"\n");
 
         fprintf(fp,SN,po);
         if( nn==1 )
              fprintf(fp,SN,ponn1);
         else
              fprintf(fp,ponn,nn);
         fprintf(fp,SN,po);
         fprintf(fp,potk,y1,y2,y3);
         fprintf(fp,SN,po);
         for( i=1 ; i<=nn ; ++i ){
              fprintf(fp,pok,i);
              if(i != nn )
                   fprintf(fp,SN,potl);
              else
                   fprintf(fp,SN,potu);
         }
    }
    fclose(fp);
    printf("file %s produced.\n",name);
}
 
/*******************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
**********************  E I N G A B E  -  F I L E S   **************************
**********************                                **************************
**********************            L E S E N           **************************
********************************************************************************
********************************************************************************
********************************************************************************
*******************************************************************************/
 
 
 
/*------------------------------------------------------------------------------
                          read_Vkq()
------------------------------------------------------------------------------*/
ITERATION *read_Vkq(name,vsymmetrienr_vor)  /* Vkq aus file name lesen */
    CHAR *name;
    INT  vsymmetrienr_vor;/* falls symmetrienr nicht vorgegeben  */
                          /* dann vsymmetrienr_vor <  0          */
{
    FILE      *fp,*fopen();
    INT       anz_nn,dimj/*,zwei_j*/,ionennr,symmetrienr;
    INT       buffer_size=381,einheitnr_in,einheitnr_out;
    DOUBLE    versionsnummer;
    DOUBLE    myB;
    DOUBLE    sin(),cos();
    DOUBLE    a_tof(),v40,v44,v60,v64,sqrt(),temperatur;
    CHAR  /*  *einheit_in,*einheit_out,*/modus;
    CHAR      *ion;
    CHAR      c,*string,*line,*fgets(),*a_tos();
    ITERATION *iteration,*iter_alloc();
    ITERATION *auswahlregel();
    STEVENS   *calc_Pkq();
    MATRIX    *readBmag();
 
    printf("Reading file %s ....\n",name);
    string   = STRING_ALLOC(buffer_size);
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 1.|    */
    versionsnummer =a_tof( line , 11,23);
    fclose(fp);
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    line=fgets( string , buffer_size , fp ); /* leerzeile                  */
    line=fgets( string , buffer_size , fp ); /* Kristallfeldpara....       */
    c = VALUE(line,33);
    einheitnr_in = is_einheit_imp(c);
    if( einheitnr_in == NICHTIMP )
         read_error(21,fp,name);
/*  einheit_in = EINHEITIMP[ einheitnr_in ].einheit; */
    myB        = EINHEITIMP[ einheitnr_in ].myB;
 
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    line=fgets( string , buffer_size , fp ); /* : Energieeigenwerte......*/
    c = VALUE(line,31);
    einheitnr_out= is_einheit_imp(c);
    if( einheitnr_out== NICHTIMP )
         read_error(21,fp,name);
/*  einheit_out= EINHEITIMP[ einheitnr_out].einheit; */
 
    line=fgets( string , buffer_size , fp ); /* : Temperatur der Probe :.*/
    temperatur = a_tof(line,31,61);
 
    line=fgets( string , buffer_size , fp ); /* : Betrachtetes Aufion....*/
    ion =a_tos( line , 31,56);
 
    line=fgets( string , buffer_size , fp ); /* : Symmetrie ... Symmnr...*/
    symmetrienr =a_toi( line , 53,61);
    if( vsymmetrienr_vor >= 0 ) symmetrienr = vsymmetrienr_vor;
 
    anz_nn    = 0;
    ionennr   = isimplementiert(ion);
    dimj      = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

/*  zwei_j    = dimj - 1 ; */
    iteration = iter_alloc(dimj,anz_nn);
 
 
    ANZ_NN(   iteration)     =  anz_nn;
    DIMJ(     iteration)     =  dimj;
    GJ(       iteration)     =  IONENIMP[ ionennr ].gj;
    IONNAME(  iteration)     =  IONENIMP[ ionennr ].ionname;
    EINHEITNRIN( iteration)     =  einheitnr_in;
    EINHEITNROUT(iteration)     =  einheitnr_out;
    PKQ(    iteration)     =  calc_Pkq( dimj );
    SYMMETRIENR(iteration) =  symmetrienr;
    TEMPERATUR( iteration) =  temperatur;
    EFVERSION(  iteration) =  versionsnummer;
 
    line=fgets( string , buffer_size , fp ); /* Magnetfeld: . Symmnr...  */
    c=VALUE(line,31);
    switch(c){
        case 'a':  modus = 'a';
                   break;
        case 'n':  modus = 'n';
                   break;
        default : read_error(24,fp,name); modus=0;
    }
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 3.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 4.==== */
 
    line=fgets( string , buffer_size , fp ); /*|RE V20 ............. */
    RT( V20(iteration) ) = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V21 .. |IM V21 .. */
    RT( V21(iteration) ) = a_tof(line, 9,30);
    IT( V21(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V22 .. |IM V22 .. */
    RT( V22(iteration) ) = a_tof(line, 9,30);
    IT( V22(iteration) ) = a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 5.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 6.==== */
    line=fgets( string , buffer_size , fp ); /*|RE V40 ............. */
    RT( V40(iteration) ) = v40 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V41 .. |IM V41 .. */
    RT( V41(iteration) ) = a_tof(line, 9,30);
    IT( V41(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V42 .. |IM V42 .. */
    RT( V42(iteration) ) = a_tof(line, 9,30);
    IT( V42(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V43 .. |IM V43 .. */
    RT( V43(iteration) ) = a_tof(line, 9,30);
    IT( V43(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V44 .. |IM V44 .. */
    RT( V44(iteration) ) = v44 = a_tof(line, 9,30);
    IT( V44(iteration) ) = a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 7.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 8.==== */
    line=fgets( string , buffer_size , fp ); /*|RE V60 ............. */
    RT( V60(iteration) ) = v60 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V61 .. |IM V61 .. */
    RT( V61(iteration) ) = a_tof(line, 9,30);
    IT( V61(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V62 .. |IM V62 .. */
    RT( V62(iteration) ) = a_tof(line, 9,30);
    IT( V62(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V63 .. |IM V63 .. */
    RT( V63(iteration) ) = a_tof(line, 9,30);
    IT( V63(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V64 .. |IM V64 .. */
    RT( V64(iteration) ) = v64 = a_tof(line, 9,30);
    IT( V64(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V65 .. |IM V65 .. */
    RT( V65(iteration) ) = a_tof(line, 9,30);
    IT( V65(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V66 .. |IM V66 .. */
    RT( V66(iteration) ) = a_tof(line, 9,30);
    IT( V66(iteration) ) = a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*9.==== */
 
    /* Auswahlregeln fuer Vkq beachten */
    iteration = auswahlregel(iteration,symmetrienr);
 
 
 
    if( symmetrienr==8 ){
       if(ABSD( 2.0*v44- 5.0*v40 ) > EPS1) read_error(22,fp,name);
       if(ABSD( 2.0*v64+21.0*v60 ) > EPS1) read_error(23,fp,name);
    }
 
 
    HMAG(iteration) = readBmag(fp,name,modus,myB,iteration,buffer_size,string);
 
    fclose(fp);
    return( iteration );
}
/*------------------------------------------------------------------------------
                          readBmag()
------------------------------------------------------------------------------*/
MATRIX *readBmag(fp,name,modus,myB,iteration,buffer_size,string)
  FILE      *fp;
  CHAR      *name,modus;
  DOUBLE    myB;
  ITERATION *iteration;
  INT        buffer_size;
  CHAR      *string;
 
{
  CHAR    *line,*fgets(),c;
  INT     i;
  DOUBLE  x1,x2,x3,a_tof();
  DOUBLE    h,theta,phi;
  DOUBLE    sin(),cos();
  MATRIX    *calc_Bmag();
  MATRIX    *calcBmol();
 
  QUADRANT  *q,*abfrage();
  DOUBLE    macheps,accuracy();
 
    macheps = MACHEPSFACT*accuracy();
 
    if( modus=='n' ) {      /* kein Magnetfeld angelegt */
        B1(iteration) = 0.0;
        B2(iteration) = 0.0;
        B3(iteration) = 0.0;
        HMAG(iteration) = calc_Bmag( DIMJ(iteration),GJ(iteration),myB,
                                     B1(iteration),
                                     B2(iteration),
                                     B3(iteration)  );
				     
	MODUS(iteration) = 'r';
        B1MOL(iteration) = 0.0;
        B2MOL(iteration) = 0.0;
        B3MOL(iteration) = 0.0;
 
        BMOL( iteration) = 0.0;
        PHI(  iteration) = 0.0;
        THETA(iteration) = 0.0;
 
      return( HMAG(iteration) );
    }
    /* falls ein Magnetfeld angelegt wurde : */
 
    /* externes Magnetfeld lesen  */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*1.==== */
    line=fgets( string , buffer_size , fp ); /* : die Koordinaten der ...*/
    line=fgets( string , buffer_size , fp );/* : º P... v º R... v º S...*/
    c   = VALUE(line,2);
    switch(c){
         case 'K' : MODUS(iteration) = 'r';
                    break;
         case 'S' : MODUS(iteration) = 's';
                    break;
         case 'P' : MODUS(iteration) = 'p';
                    break;
         default  : read_error(3,fp,name);
    }
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*2.==== */
 
    for( i=1 ; i<= 3 ; ++i)/* Kopf der Magnetfeldtabelle ueberlesen */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );
    line=fgets( string , buffer_size , fp );   /* : º  h1 º h2 º h3 º*/
    B1(iteration) = x1 = a_tof(line, 2,13);
    B2(iteration) = x2 = a_tof(line,15,26);
    x3=0.0;
    if(MODUS(iteration) !='p')  B3(iteration) = x3 = a_tof(line,28,39);
    isinlimits(fp,name ,0, x1,x2,x3,MODUS(iteration) );
 
 
    switch( MODUS(iteration) ){ /* transformiere B -> (Bx,By,Bz) */
 
       case 'p' :
       case 's' :  h     = B1(iteration);
                   phi   = B2(iteration)*pi/180.0;
                   theta = B3(iteration)*pi/180.0;
 
                   if( MODUS(iteration) == 'p' )  theta = pi/2;
 
                   B1(iteration) = h*sin(theta)*cos(phi);
                   B2(iteration) = h*sin(theta)*sin(phi);
                   B3(iteration) = h*cos(theta);
                   break;
 
       case 'r' :  break;
    }
 
    HMAG(iteration) = calc_Bmag( DIMJ(iteration),GJ(iteration),
                                          myB,B1(iteration),
                                              B2(iteration),
                                              B3(iteration) );
 
    /* Molekularfeld lesen  */
 if( EFVERSION(iteration) >= 2.0 ){
    for( i=1 ; i<= 3 ; ++i)/* Kopf der Magnetfeldtabelle ueberlesen */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );
    line=fgets( string , buffer_size , fp );   /* : º  h1 º h2 º h3 º*/
    B1MOL(iteration) = x1 = a_tof(line, 2,13);
    B2MOL(iteration) = x2 = a_tof(line,15,26);
    x3=0.0;
    if(MODUS(iteration) !='p') B3MOL(iteration) = x3 = a_tof(line,28,39);
    isinlimits(fp,name ,0, x1,x2,x3,MODUS(iteration) );
 
 
    switch( MODUS(iteration) ){ /* transformiere B -> (Bx,By,Bz) */
 
       case 'p' :
       case 's' :  h     = B1MOL(iteration);
                   phi   = B2MOL(iteration)*pi/180.0;
                   theta = B3MOL(iteration)*pi/180.0;
 
                   if( MODUS(iteration) == 'p' )  theta = pi/2;
 
                   B1MOL(iteration) = h*sin(theta)*cos(phi);
                   B2MOL(iteration) = h*sin(theta)*sin(phi);
                   B3MOL(iteration) = h*cos(theta);
                   BMOL( iteration) = h;
                   PHI(  iteration) = phi;
                   THETA(iteration) = theta;
                   break;
 
       case 'r' :  /* umrechnung (Bx,By,Bz) -> (B,theta,phi) */
                   q = abfrage(B2MOL(iteration),B1MOL(iteration),macheps);
                   phi = WINKEL(q);
                   h   = RADIUS(q); /* h=r*sin(theta) */
                   free_(q);
                   q = abfrage(h,B3MOL(iteration),macheps);
                   theta = WINKEL(q);
                   h     = RADIUS(q);
                   free_(q);
                   BMOL( iteration) = h;
                   PHI(  iteration) = phi;
                   THETA(iteration) = theta;
                   break;
    }
 
    HMAG(iteration) = calcBmol( DIMJ(iteration),HMAG(iteration),
                                2.0*(GJ(iteration)-1),myB,
                                B1MOL(iteration),
                                B2MOL(iteration),
                                B3MOL(iteration) );
 }
 else{
        B1MOL(iteration) = 0.0;
        B2MOL(iteration) = 0.0;
        B3MOL(iteration) = 0.0;
        BMOL( iteration) = 0.0;
        PHI(  iteration) = 0.0;
        THETA(iteration) = 0.0;
     }
   return ( HMAG(iteration) );
}
/*------------------------------------------------------------------------------
                          read_Dkq()
------------------------------------------------------------------------------*/
ITERATION *read_Dkq(name,vsymmetrienr_vor)  /* Dkq aus file name lesen */
    CHAR *name;
    INT  vsymmetrienr_vor;/* falls symmetrienr nicht vorgegeben  */
                          /* dann vsymmetrienr_vor <  0          */
{
    FILE      *fp,*fopen();
    INT       anz_nn,dimj/*,zwei_j*/,ionennr,symmetrienr,e_4f;
    INT       buffer_size=381,einheitnr_in,einheitnr_out;
    DOUBLE    versionsnummer;
    DOUBLE    myB;
   
 
    DOUBLE    sin(),cos();
    DOUBLE    a_tof(),v40,v44,v60,v64,sqrt(),temperatur;
    CHAR  /*  *einheit_in,*einheit_out,*/modus;
    CHAR      *ion;
    CHAR      c,*string,*line,*fgets(),*a_tos();
    ITERATION *iteration,*iter_alloc();
    ITERATION *auswahlregel();
    STEVENS   *calc_Pkq();
    MATRIX    *readBmag();
 
    printf("Reading file %s ....\n",name);
    string   = STRING_ALLOC(buffer_size);
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 1.==== */
    versionsnummer =a_tof( line , 11,23);
    fclose(fp);
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    line=fgets( string , buffer_size , fp ); /* leerzeile                  */
    line=fgets( string , buffer_size , fp ); /* Kristallfeldpara....       */
    c = VALUE(line,33);
    einheitnr_in = is_einheit_imp(c);
    if( einheitnr_in == NICHTIMP )
         read_error(21,fp,name);
/*  einheit_in = EINHEITIMP[ einheitnr_in ].einheit; */
    myB        = EINHEITIMP[ einheitnr_in ].myB;
 
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    line=fgets( string , buffer_size , fp ); /* : Energieeigenwerte......*/
    c = VALUE(line,31);
    einheitnr_out= is_einheit_imp(c);
    if( einheitnr_out== NICHTIMP )
         read_error(21,fp,name);
/*  einheit_out= EINHEITIMP[ einheitnr_out].einheit; */
 
    line=fgets( string , buffer_size , fp ); /* : Temperatur der Probe :.*/
    temperatur = a_tof(line,31,61);
 
    line=fgets( string , buffer_size , fp ); /* : Betrachtetes Aufion....*/
    ion =a_tos( line , 31,56);
 
    line=fgets( string , buffer_size , fp ); /* : Symmetrie ... Symmnr...*/
    symmetrienr =a_toi( line , 53,61);
    if( vsymmetrienr_vor >= 0 ) symmetrienr = vsymmetrienr_vor;
 
    anz_nn    = 0;
    ionennr   = isimplementiert(ion);
    dimj      = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

/*  zwei_j    = dimj - 1 ; */
    iteration = iter_alloc(dimj,anz_nn);
 
 
    ANZ_NN(   iteration)     =  anz_nn;
    DIMJ(     iteration)     =  dimj;
    GJ(       iteration)     =  IONENIMP[ ionennr ].gj;
    IONNAME(  iteration)     =  IONENIMP[ ionennr ].ionname;
    EINHEITNRIN( iteration)     =  einheitnr_in;
    EINHEITNROUT(iteration)     =  einheitnr_out;
    PKQ(    iteration)     =  calc_Pkq( dimj );
    SYMMETRIENR(iteration) =  symmetrienr;
    TEMPERATUR( iteration) =  temperatur;
    EFVERSION(  iteration) =  versionsnummer;
 
    line=fgets( string , buffer_size , fp ); /* Magnetfeld: . Symmnr...  */
    c=VALUE(line,31);
    switch(c){
        case 'a':  modus = 'a';
                   break;
        case 'n':  modus = 'n';
                   break;
        default : read_error(24,fp,name); modus=0;
    }
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 3.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 4.==== */
 
    line=fgets( string , buffer_size , fp ); /*|RE D20 ............. */
    RT( V20(iteration) ) = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D21 .. |IM D21 .. */
    RT( V21(iteration) ) = a_tof(line, 9,30);
    IT( V21(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D22 .. |IM D22 .. */
    RT( V22(iteration) ) = a_tof(line, 9,30);
    IT( V22(iteration) ) = a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 5.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 6.==== */
    line=fgets( string , buffer_size , fp ); /*|RE D40 ............. */
    RT( V40(iteration) ) = v40 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D41 .. |IM D41 .. */
    RT( V41(iteration) ) = a_tof(line, 9,30);
    IT( V41(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D42 .. |IM D42 .. */
    RT( V42(iteration) ) = a_tof(line, 9,30);
    IT( V42(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D43 .. |IM D43 .. */
    RT( V43(iteration) ) = a_tof(line, 9,30);
    IT( V43(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D44 .. |IM D44 .. */
    RT( V44(iteration) ) = v44 = a_tof(line, 9,30);
    IT( V44(iteration) ) = a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 7.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 8.==== */
    line=fgets( string , buffer_size , fp ); /*|RE D60 ............. */
    RT( V60(iteration) ) = v60 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D61 .. |IM D61 .. */
    RT( V61(iteration) ) = a_tof(line, 9,30);
    IT( V61(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D62 .. |IM D62 .. */
    RT( V62(iteration) ) = a_tof(line, 9,30);
    IT( V62(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D63 .. |IM D63 .. */
    RT( V63(iteration) ) = a_tof(line, 9,30);
    IT( V63(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D64 .. |IM D64 .. */
    RT( V64(iteration) ) = v64 = a_tof(line, 9,30);
    IT( V64(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D65 .. |IM D65 .. */
    RT( V65(iteration) ) = a_tof(line, 9,30);
    IT( V65(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE D66 .. |IM D66 .. */
    RT( V66(iteration) ) = a_tof(line, 9,30);
    IT( V66(iteration) ) = a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*9.==== */
 
    /* Auswahlregeln fuer Dkq beachten */
    iteration = auswahlregel(iteration,symmetrienr);
 
 
    /* Dkq auf Vkq umrechnen            */
    /*                                  */
    /* Vkq = Dkq * epsilon_kq * theta_k */
    /*                                  */
    /*                                  */
 
     e_4f = E4f( ionennr );
 
 
     RT( V20(iteration) ) *=  epn2n(0) * alpha_J[ e_4f ];
     RT( V21(iteration) ) *=  epn1n(1) * alpha_J[ e_4f ];
     RT( V22(iteration) ) *=  epn0n(2) * alpha_J[ e_4f ];
 
     RT( V40(iteration) ) *=  epn4n(0) * beta_J[ e_4f ];
     RT( V41(iteration) ) *=  epn3n(1) * beta_J[ e_4f ];
     RT( V42(iteration) ) *=  epn2n(2) * beta_J[ e_4f ];
     RT( V43(iteration) ) *=  epn1n(3) * beta_J[ e_4f ];
     RT( V44(iteration) ) *=  epn0n(4) * beta_J[ e_4f ];
 
     RT( V60(iteration) ) *=  epn6n(0) * gamma_J[ e_4f ];
     RT( V61(iteration) ) *=  epn5n(1) * gamma_J[ e_4f ];
     RT( V62(iteration) ) *=  epn4n(2) * gamma_J[ e_4f ];
     RT( V63(iteration) ) *=  epn3n(3) * gamma_J[ e_4f ];
     RT( V64(iteration) ) *=  epn2n(4) * gamma_J[ e_4f ];
     RT( V65(iteration) ) *=  epn1n(5) * gamma_J[ e_4f ];
     RT( V66(iteration) ) *=  epn0n(6) * gamma_J[ e_4f ];
 
 
 
     IT( V21(iteration) ) *=  epn1n(1) * alpha_J[ e_4f ];
     IT( V22(iteration) ) *=  epn0n(2) * alpha_J[ e_4f ];
 
     IT( V41(iteration) ) *=  epn3n(1) * beta_J[ e_4f ];
     IT( V42(iteration) ) *=  epn2n(2) * beta_J[ e_4f ];
     IT( V43(iteration) ) *=  epn1n(3) * beta_J[ e_4f ];
     IT( V44(iteration) ) *=  epn0n(4) * beta_J[ e_4f ];
 
     IT( V61(iteration) ) *=  epn5n(1) * gamma_J[ e_4f ];
     IT( V62(iteration) ) *=  epn4n(2) * gamma_J[ e_4f ];
     IT( V63(iteration) ) *=  epn3n(3) * gamma_J[ e_4f ];
     IT( V64(iteration) ) *=  epn2n(4) * gamma_J[ e_4f ];
     IT( V65(iteration) ) *=  epn1n(5) * gamma_J[ e_4f ];
     IT( V66(iteration) ) *=  epn0n(6) * gamma_J[ e_4f ];
 
 
 
    if( symmetrienr==8 ){
       if( ABSD(2.0*epn0n(4)*v44 - 5.0*epn4n(0)*v40) > EPS1)
       read_error(42,fp,name);
       if( ABSD(2.0*epn2n(4)*v64 + 21.0*epn6n(0)*v60) > EPS1)
       read_error(43,fp,name);
    }
 
 
    HMAG(iteration) = readBmag(fp,name,modus,myB,iteration,buffer_size,string);
 
    fclose(fp);
    return( iteration );
}
INT  extract(line,parnam,var,unit)
CHAR * line;
const CHAR * parnam;
DOUBLE * var;
const CHAR * unit;
{CHAR * token;
 // check if line is true comment line -> if yes return 1
 if (line[strspn(line," \t")]=='#'&&line[strspn(line," \t#")]!='!') return 0;

      if ((token = strstr (line, parnam))!=NULL)
        {          token+=strlen(parnam);
         if (strstr (token, "=")!=NULL)
           {while(strstr(token," ")==token)++token;
            if(token ==strstr (token, "="))
              {++token;while(strstr(token," ")==token)++token;
                              (*var)  = strtod (token, NULL);
               printf("%s=%g %s ",parnam,(*var),unit);
               return 1;
              }
           }
        } 
  return 0;
}

/************************************************************************************/
/*  for reading in case of opiont -L -B (read_Lkq and read_Bkq) from new file format */
/************************************************************************************/
ITERATION *read_new_format(type,iteration,name,vsymmetrienr_vor)
CHAR type;
ITERATION *iteration;
    CHAR *name;
    INT  vsymmetrienr_vor;/* falls symmetrienr nicht vorgegeben  */
                          /* dann vsymmetrienr_vor <  0          */
{   FILE      *fp,*fopen();
    INT       anz_nn,dimj/*,zwei_j*/,ionennr,symmetrienr,e_4f;
    INT       buffer_size=381,einheitnr_in,einheitnr_out,ia=0,ib=0,ic=0;
    DOUBLE    versionsnummer;
    DOUBLE    myB,d,alpha,beta,gamma;

    DOUBLE    sin(),cos();
    DOUBLE    a_tof(),v40=0,v44=0,v60=0,v64=0,sqrt(),temperatur;
    CHAR  /*  *einheit_in,*einheit_out,*/modus;
    CHAR      *ion=0,*token;
    CHAR      c,*string,*line,*fgets(),*a_tos();
    ITERATION  *iter_alloc();
    ITERATION *auswahlregel();
    STEVENS   *calc_Pkq();
    MATRIX    *readBmag();
  /* read mcphas single ion input file */
       printf("file format as single ion input module #!MODULE=so1ion or #!MODULE=cfield\n");
    string   = STRING_ALLOC(buffer_size);


    /* some fixed quantitities */
    symmetrienr =0;
     versionsnummer=VERSION;
    c = 'm'; /* unit is meV */
    einheitnr_in = is_einheit_imp(c);
    if( einheitnr_in == NICHTIMP )read_error(21,fp,name);
    myB        = EINHEITIMP[ einheitnr_in ].myB;
    c = 'm';
    einheitnr_out= is_einheit_imp(c);
    if( einheitnr_out== NICHTIMP )
         read_error(21,fp,name);
/*  einheit_out= EINHEITIMP[ einheitnr_out].einheit; */

     temperatur=10.0; /* can be modified by reading in a temperature below*/

    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
    while(feof(fp)==0&&ion==0)
    {line=fgets( string , buffer_size , fp );
     // check if line is true comment line -> if yes do not read it
     if(feof(fp)==0&&(line[strspn(line," \t")]!='#'||line[strspn(line," \t#")]=='!'))      
     {while ((token=strchr(line,'\r'))!=NULL){*token=' ';}
      /*read iontype*/
      if ((token = strstr (line, "IONTYPE"))!=NULL)
        {token+=strlen("IONTYPE");
         if(strstr(line,"perl")!=NULL){fprintf(stderr,"ERROR so1ion/cfield: perl parsing not implemented "
                                          "- try to use program singleion !\n");exit(EXIT_FAILURE);}
         if (strstr (token, "=")!=NULL)
           {token = strstr (token, "=")+1;
           while (*token==' '){++token;} /* remove starting spaces*/
           ion =a_tos( token , 0,5);
            /*strncpy(ion,token,1)*/;/*maximal 5 characters*/
            /*remove from string var all characters after delimiters*/
            strtok(ion," \n");
            printf("IONTYPE=%s ",ion);
           }
        } 

     }
    }
    fclose(fp);printf("\n");
    if(ion==0){fprintf(stderr,"ERROR so1ion/cfield:IONTYPE not found!\n");exit(EXIT_FAILURE);}

    if(strncmp(ion,"S=",2)==0)  /* J=... ion !! same as S= */
     {ion[0]='J';}    
  
    anz_nn    = 0;
    ionennr   = isimplementiert(ion);
    dimj      = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"J=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    
/*  zwei_j    = dimj - 1 ; */
 
    iteration = iter_alloc(dimj,anz_nn);
 
    ANZ_NN(   iteration)     =  anz_nn;
    DIMJ(     iteration)     =  dimj;
    GJ(       iteration)     =  IONENIMP[ ionennr ].gj;
    IONNAME(  iteration)     =  IONENIMP[ ionennr ].ionname;
    IONENNUMMER(iteration)   = ionennr;
    EINHEITNRIN( iteration)     =  einheitnr_in;
    EINHEITNROUT(iteration)     =  einheitnr_out;
    PKQ(      iteration)     =  calc_Pkq( dimj );
    SYMMETRIENR(iteration) =  symmetrienr;
    EFVERSION(  iteration) =  versionsnummer;
 
 
    modus='a';/* Magnetfeld: . Symmnr...  */
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
    while(feof(fp)==0)
    {line=fgets( string , buffer_size , fp );
     if(feof(fp)==0&&strstr (line, "#")==NULL){while ((token=strchr(line,'\r'))!=NULL){*token=' ';}
        if(extract(line,"GJ",&GJ(iteration)," "))if(fabs(IONENIMP[ ionennr ].gj-GJ(iteration))>1e-4)
           fprintf(stderr,"# Warning: taking Lande factor from sipf file which does not correspond to "
                           "internal tabulated value GJinternal=%g\n",IONENIMP[ ionennr ].gj);
        if(extract(line,"R2",&IONENIMP[ ionennr ].r2,"a0^2 a0=0.5292 A"));
        if(extract(line,"R4",&IONENIMP[ ionennr ].r4,"a0^4 a0=0.5292 A"));
        if(extract(line,"R6",&IONENIMP[ ionennr ].r6,"a0^6 a0=0.5292 A"));
        if(extract(line,"ALPHA",&alpha," "))ia=1;
        if(extract(line,"BETA",&beta," "))ib=1;
        if(extract(line,"GAMMA",&gamma," "))ic=1;
        if(extract(line,"nof_electrons",&d," "))IONENIMP[ ionennr ].elektronen_in_vier_f=(INT)d; 
   switch(type)
   {case 'L':     
     extract(line,"L20",&RT(V20(iteration)),"meV");
     extract(line,"L21",&RT(V21(iteration)),"meV");
     extract(line,"L21S",&IT(V21(iteration)),"meV");IT(V21(iteration))*=-1.0;
     extract(line,"L22",&RT(V22(iteration)),"meV");
     extract(line,"L22S",&IT(V22(iteration)),"meV");IT(V22(iteration))*=-1.0;
     extract(line,"L40",&RT(V40(iteration)),"meV");
     extract(line,"L41",&RT(V41(iteration)),"meV");
     extract(line,"L41S",&IT(V41(iteration)),"meV");IT(V41(iteration))*=-1.0;
     extract(line,"L42",&RT(V42(iteration)),"meV");
     extract(line,"L42S",&IT(V42(iteration)),"meV");IT(V42(iteration))*=-1.0;
     extract(line,"L43",&RT(V43(iteration)),"meV");
     extract(line,"L43S",&IT(V43(iteration)),"meV");IT(V43(iteration))*=-1.0;
     extract(line,"L44",&RT(V44(iteration)),"meV");
     extract(line,"L44S",&IT(V44(iteration)),"meV");IT(V44(iteration))*=-1.0;
     extract(line,"L60", &RT(V60(iteration)),"meV");
     extract(line,"L61", &RT(V61(iteration)),"meV");
     extract(line,"L61S",&IT(V61(iteration)),"meV");IT(V61(iteration))*=-1.0;
     extract(line,"L62", &RT(V62(iteration)),"meV");
     extract(line,"L62S",&IT(V62(iteration)),"meV");IT(V62(iteration))*=-1.0;
     extract(line,"L63", &RT(V63(iteration)),"meV");
     extract(line,"L63S",&IT(V63(iteration)),"meV");IT(V63(iteration))*=-1.0;
     extract(line,"L64", &RT(V64(iteration)),"meV");
     extract(line,"L64S",&IT(V64(iteration)),"meV");IT(V64(iteration))*=-1.0;
     extract(line,"L65", &RT(V65(iteration)),"meV");
     extract(line,"L65S",&IT(V65(iteration)),"meV");IT(V65(iteration))*=-1.0;
     extract(line,"L66", &RT(V66(iteration)),"meV");
     extract(line,"L66S",&IT(V66(iteration)),"meV");IT(V66(iteration))*=-1.0;
    break;
    case 'B':
     extract(line,"B20",&RT(V20(iteration)),"meV");
     extract(line,"B21",&RT(V21(iteration)),"meV");
     extract(line,"B21S",&IT(V21(iteration)),"meV");IT(V21(iteration))*=-1.0;
     extract(line,"B22",&RT(V22(iteration)),"meV");
     extract(line,"B22S",&IT(V22(iteration)),"meV");IT(V22(iteration))*=-1.0;
     extract(line,"B40",&RT(V40(iteration)),"meV");
     extract(line,"B41",&RT(V41(iteration)),"meV");
     extract(line,"B41S",&IT(V41(iteration)),"meV");IT(V41(iteration))*=-1.0;
     extract(line,"B42",&RT(V42(iteration)),"meV");
     extract(line,"B42S",&IT(V42(iteration)),"meV");IT(V42(iteration))*=-1.0;
     extract(line,"B43",&RT(V43(iteration)),"meV");
     extract(line,"B43S",&IT(V43(iteration)),"meV");IT(V43(iteration))*=-1.0;
     extract(line,"B44",&RT(V44(iteration)),"meV");
     extract(line,"B44S",&IT(V44(iteration)),"meV");IT(V44(iteration))*=-1.0;
     extract(line,"B60", &RT(V60(iteration)),"meV");
     extract(line,"B61", &RT(V61(iteration)),"meV");
     extract(line,"B61S",&IT(V61(iteration)),"meV");IT(V61(iteration))*=-1.0;
     extract(line,"B62", &RT(V62(iteration)),"meV");
     extract(line,"B62S",&IT(V62(iteration)),"meV");IT(V62(iteration))*=-1.0;
     extract(line,"B63", &RT(V63(iteration)),"meV");
     extract(line,"B63S",&IT(V63(iteration)),"meV");IT(V63(iteration))*=-1.0;
     extract(line,"B64", &RT(V64(iteration)),"meV");
     extract(line,"B64S",&IT(V64(iteration)),"meV");IT(V64(iteration))*=-1.0;
     extract(line,"B65", &RT(V65(iteration)),"meV");
     extract(line,"B65S",&IT(V65(iteration)),"meV");IT(V65(iteration))*=-1.0;
     extract(line,"B66", &RT(V66(iteration)),"meV");
     extract(line,"B66S",&IT(V66(iteration)),"meV");IT(V66(iteration))*=-1.0;
    break;
    default: printf("error - not implemented cf parameter type for this format of input file\n");exit(1);
    }
     extract(line,"Bx", &B1(iteration),"T");
     extract(line,"By", &B2(iteration),"T");
     extract(line,"Bz", &B3(iteration),"T");
      /*read temperature*/
      extract(line,"TEMP",&temperatur,"K");

     extract(line,"Dx2", &B1S(iteration),"meV");
     extract(line,"Dy2", &B2S(iteration),"meV");
     extract(line,"Dz2", &B3S(iteration),"meV"); 
     }
    }
    printf("\n");
    e_4f = E4f( ionennr );
    if(ia)alpha_J[ e_4f ]=alpha;
    if(ib)beta_J[ e_4f ]=beta;
    if(ic)gamma_J[ e_4f ]=gamma;

     /* Auswahlregeln beachten */
    iteration = auswahlregel(iteration,symmetrienr);
    

   switch(type)
   {case 'L': 
 
     /* Lkq auf Dkq umrechnen                 */
    /*                                       */
    /*             q                         */
    /* D     = (-1)  (  L      + i L      )  */
    /*  k,|q|            k,|q|      k,-|q|   */
    /*                                       */
 
     RT( V21(iteration) ) *=  -1;
 
     RT( V41(iteration) ) *=  -1;
     RT( V43(iteration) ) *=  -1;
 
     RT( V61(iteration) ) *=  -1;
     RT( V63(iteration) ) *=  -1;
     RT( V65(iteration) ) *=  -1;
 
 
     IT( V21(iteration) ) *=  -1;
 
     IT( V41(iteration) ) *=  -1;
     IT( V43(iteration) ) *=  -1;
 
     IT( V61(iteration) ) *=  -1;
     IT( V63(iteration) ) *=  -1;
     IT( V65(iteration) ) *=  -1;
 
    /* Dkq auf Vkq umrechnen            */
    /*                                  */
    /* Vkq = Dkq * epsilon_kq * theta_k */
    /*                                  */
    /*                                  */
 
     

     RT( V20(iteration) ) *=  epn2n(0) * alpha_J[ e_4f ];
     RT( V21(iteration) ) *=  epn1n(1) * alpha_J[ e_4f ];
     RT( V22(iteration) ) *=  epn0n(2) * alpha_J[ e_4f ];
 
     RT( V40(iteration) ) *=  epn4n(0) * beta_J[ e_4f ];
     RT( V41(iteration) ) *=  epn3n(1) * beta_J[ e_4f ];
     RT( V42(iteration) ) *=  epn2n(2) * beta_J[ e_4f ];
     RT( V43(iteration) ) *=  epn1n(3) * beta_J[ e_4f ];
     RT( V44(iteration) ) *=  epn0n(4) * beta_J[ e_4f ];
 
     RT( V60(iteration) ) *=  epn6n(0) * gamma_J[ e_4f ];
     RT( V61(iteration) ) *=  epn5n(1) * gamma_J[ e_4f ];
     RT( V62(iteration) ) *=  epn4n(2) * gamma_J[ e_4f ];
     RT( V63(iteration) ) *=  epn3n(3) * gamma_J[ e_4f ];
     RT( V64(iteration) ) *=  epn2n(4) * gamma_J[ e_4f ];
     RT( V65(iteration) ) *=  epn1n(5) * gamma_J[ e_4f ];
     RT( V66(iteration) ) *=  epn0n(6) * gamma_J[ e_4f ];
  
     IT( V21(iteration) ) *=  epn1n(1) * alpha_J[ e_4f ];
     IT( V22(iteration) ) *=  epn0n(2) * alpha_J[ e_4f ];
 
     IT( V41(iteration) ) *=  epn3n(1) * beta_J[ e_4f ];
     IT( V42(iteration) ) *=  epn2n(2) * beta_J[ e_4f ];
     IT( V43(iteration) ) *=  epn1n(3) * beta_J[ e_4f ];
     IT( V44(iteration) ) *=  epn0n(4) * beta_J[ e_4f ];
 
     IT( V61(iteration) ) *=  epn5n(1) * gamma_J[ e_4f ];
     IT( V62(iteration) ) *=  epn4n(2) * gamma_J[ e_4f ];
     IT( V63(iteration) ) *=  epn3n(3) * gamma_J[ e_4f ];
     IT( V64(iteration) ) *=  epn2n(4) * gamma_J[ e_4f ];
     IT( V65(iteration) ) *=  epn1n(5) * gamma_J[ e_4f ];
     IT( V66(iteration) ) *=  epn0n(6) * gamma_J[ e_4f ];
 
 
 
    if( symmetrienr==8 ){
       if( ABSD(2.0*epn0n(4)*v44 - 5.0*epn4n(0)*v40) > EPS1)
       read_error(52,fp,name);
       if( ABSD(2.0*epn2n(4)*v64 + 21.0*epn6n(0)*v60) > EPS1)
       read_error(53,fp,name);
    }
    break;
    case 'B':
  
     /* Bkq auf Vkq umrechnen       */
     /*                             */
     /*       | Bk0  ,         q =0 */
     /* Vkq = |                     */
     /*       | Bkq/2/omegakq, q >0 */
     /*                             */
 
 
     RT( V21(iteration) ) /=  2 * omegan1n(1);
     RT( V22(iteration) ) /=  2 * omegan0n(2);
 
     RT( V41(iteration) ) /=  2 * omegan3n(1);
     RT( V42(iteration) ) /=  2 * omegan2n(2);
     RT( V43(iteration) ) /=  2 * omegan1n(3);
     RT( V44(iteration) ) /=  2 * omegan0n(4);
 
     RT( V61(iteration) ) /=  2 * omegan5n(1);
     RT( V62(iteration) ) /=  2 * omegan4n(2);
     RT( V63(iteration) ) /=  2 * omegan3n(3);
     RT( V64(iteration) ) /=  2 * omegan2n(4);
     RT( V65(iteration) ) /=  2 * omegan1n(5);
     RT( V66(iteration) ) /=  2 * omegan0n(6);
 
 
     IT( V21(iteration) ) /=  2 * omegan1n(1);
     IT( V22(iteration) ) /=  2 * omegan0n(2);
 
     IT( V41(iteration) ) /=  2 * omegan3n(1);
     IT( V42(iteration) ) /=  2 * omegan2n(2);
     IT( V43(iteration) ) /=  2 * omegan1n(3);
     IT( V44(iteration) ) /=  2 * omegan0n(4);
 
     IT( V61(iteration) ) /=  2 * omegan5n(1);
     IT( V62(iteration) ) /=  2 * omegan4n(2);
     IT( V63(iteration) ) /=  2 * omegan3n(3);
     IT( V64(iteration) ) /=  2 * omegan2n(4);
     IT( V65(iteration) ) /=  2 * omegan1n(5);
     IT( V66(iteration) ) /=  2 * omegan0n(6);
    break;
    default: read_error(53,fp,name);    
   }

     MODUS(iteration) = 'r';
 
    TEMPERATUR( iteration) =  temperatur;
    /* falls ein Magnetfeld angelegt wurde : */
    HMAG(iteration) = calc_Bmag( DIMJ(iteration),GJ(iteration),
                                          myB,B1(iteration),
                                              B2(iteration),
                                              B3(iteration) );
 
        B1MOL(iteration) = 0.0;
        B2MOL(iteration) = 0.0;
        B3MOL(iteration) = 0.0;
        BMOL( iteration) = 0.0;
        PHI(  iteration) = 0.0;
        THETA(iteration) = 0.0;
 
       fclose(fp);
 return (iteration);
}
/*------------------------------------------------------------------------------
                          read_Lkq()
------------------------------------------------------------------------------*/
ITERATION *read_Lkq(name,vsymmetrienr_vor)  /* Lkq aus file name lesen */
    CHAR *name;
    INT  vsymmetrienr_vor;/* falls symmetrienr nicht vorgegeben  */
                          /* dann vsymmetrienr_vor <  0          */
{
    FILE      *fp,*fopen();
    INT       anz_nn,dimj/*,zwei_j*/,ionennr,symmetrienr,e_4f;
    INT       buffer_size=381,einheitnr_in,einheitnr_out;
    DOUBLE    versionsnummer;
    DOUBLE    myB;

    DOUBLE    sin(),cos();
    DOUBLE    a_tof(),v40=0,v44=0,v60=0,v64=0,sqrt(),temperatur;
    CHAR  /*  *einheit_in,*einheit_out,*/modus;
    CHAR      *ion=0;
    CHAR      c,*string,*line,*fgets(),*a_tos();
    ITERATION *iteration,*iter_alloc();
    ITERATION *auswahlregel();
    STEVENS   *calc_Pkq();
    MATRIX    *readBmag();
 
    printf("Reading file %s ....\n",name);
    string   = STRING_ALLOC(buffer_size);
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
    line=fgets( string , buffer_size , fp );
    fclose(fp);
    if(strncmp(line,"#!cfield",8)==0||strncmp(line,"#!so1ion",8)==0||strncmp(line,"#!MODULE",8)==0)
   {iteration=read_new_format('L',iteration,name,vsymmetrienr_vor);
   }
    else
   { /*read classical input file bkq.parmeter*/
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 1.==== */
    versionsnummer =a_tof( line , 11,23);
    fclose(fp);
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    line=fgets( string , buffer_size , fp ); /* leerzeile                  */
    line=fgets( string , buffer_size , fp ); /* Kristallfeldpara....       */
    c = VALUE(line,33);
    einheitnr_in = is_einheit_imp(c);
    if( einheitnr_in == NICHTIMP )
         read_error(21,fp,name);
/*  einheit_in = EINHEITIMP[ einheitnr_in ].einheit; */
    myB        = EINHEITIMP[ einheitnr_in ].myB;
 
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    line=fgets( string , buffer_size , fp ); /* : Energieeigenwerte......*/
    c = VALUE(line,31);
    einheitnr_out= is_einheit_imp(c);
    if( einheitnr_out== NICHTIMP )
         read_error(21,fp,name);
/*  einheit_out= EINHEITIMP[ einheitnr_out].einheit; */
 
    line=fgets( string , buffer_size , fp ); /* : Temperatur der Probe :.*/
    temperatur = a_tof(line,31,61);
 
    line=fgets( string , buffer_size , fp ); /* : Betrachtetes Aufion....*/
    ion =a_tos( line , 31,56);
 
    line=fgets( string , buffer_size , fp ); /* : Symmetrie ... Symmnr...*/
    symmetrienr =a_toi( line , 53,61);
    if( vsymmetrienr_vor >= 0 ) symmetrienr = vsymmetrienr_vor;
 
    anz_nn    = 0;
    ionennr   = isimplementiert(ion);
    dimj      = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

/*  zwei_j    = dimj - 1 ; */
    iteration = iter_alloc(dimj,anz_nn);
 
 
    ANZ_NN(   iteration)     =  anz_nn;
    DIMJ(     iteration)     =  dimj;
    GJ(       iteration)     =  IONENIMP[ ionennr ].gj;
    IONNAME(  iteration)     =  IONENIMP[ ionennr ].ionname;
    EINHEITNRIN( iteration)     =  einheitnr_in;
    EINHEITNROUT(iteration)     =  einheitnr_out;
    PKQ(    iteration)     =  calc_Pkq( dimj );
    SYMMETRIENR(iteration) =  symmetrienr;
    TEMPERATUR( iteration) =  temperatur;
    EFVERSION(  iteration) =  versionsnummer;
 
    line=fgets( string , buffer_size , fp ); /* Magnetfeld: . Symmnr...  */
    c=VALUE(line,31);
    switch(c){
        case 'a':  modus = 'a';
                   break;
        case 'n':  modus = 'n';
                   break;
        default : read_error(24,fp,name); modus=0;
    }
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 3.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 4.==== */
 
    line=fgets( string , buffer_size , fp ); /*| L 2,0 ............. */
    RT( V20(iteration) ) = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*| L 2,1 .. | L 2,-1.. */
    RT( V21(iteration) ) = a_tof(line, 9,30);
    IT( V21(iteration) ) = -a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*| L 2,2 .. | L 2,-2.. */
    RT( V22(iteration) ) = a_tof(line, 9,30);
    IT( V22(iteration) ) = -a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 5.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 6.==== */
    line=fgets( string , buffer_size , fp ); /*|L 4,0  ............. */
    RT( V40(iteration) ) = v40 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|L 4,1  .. |L 4,-1 .. */
    RT( V41(iteration) ) = a_tof(line, 9,30);
    IT( V41(iteration) ) = -a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|L 4,2  .. |L 4,-2 .. */
    RT( V42(iteration) ) = a_tof(line, 9,30);
    IT( V42(iteration) ) = -a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|L 4,3  .. |L 4,-3 .. */
    RT( V43(iteration) ) = a_tof(line, 9,30);
    IT( V43(iteration) ) = -a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|L 4,4  .. |L 4,-4 .. */
    RT( V44(iteration) ) = v44 = a_tof(line, 9,30);
    IT( V44(iteration) ) = -a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 7.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 8.==== */
    line=fgets( string , buffer_size , fp ); /*|L 6,0  ............. */
    RT( V60(iteration) ) = v60 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|L 6,1  .. |L 6,-1 .. */
    RT( V61(iteration) ) = a_tof(line, 9,30);
    IT( V61(iteration) ) = -a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|L 6,2  .. |L 6,-2 .. */
    RT( V62(iteration) ) = a_tof(line, 9,30);
    IT( V62(iteration) ) = -a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|L 6,3  .. |L 6,-3 .. */
    RT( V63(iteration) ) = a_tof(line, 9,30);
    IT( V63(iteration) ) = -a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|L 6,4  .. |L 6,-4 .. */
    RT( V64(iteration) ) = v64 = a_tof(line, 9,30);
    IT( V64(iteration) ) = -a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|L 6,5  .. |L 6,-5 .. */
    RT( V65(iteration) ) = a_tof(line, 9,30);
    IT( V65(iteration) ) = -a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|L 6,6  .. |L 6,-6 .. */
    RT( V66(iteration) ) = a_tof(line, 9,30);
    IT( V66(iteration) ) = -a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*9.==== */
 
    /* Auswahlregeln fuer Dkq beachten */
    iteration = auswahlregel(iteration,symmetrienr);
 
 
    /* Lkq auf Dkq umrechnen                 */
    /*                                       */
    /*             q                         */
    /* D     = (-1)  (  L      + i L      )  */
    /*  k,|q|            k,|q|      k,-|q|   */
    /*                                       */
 
     RT( V21(iteration) ) *=  -1;
 
     RT( V41(iteration) ) *=  -1;
     RT( V43(iteration) ) *=  -1;
 
     RT( V61(iteration) ) *=  -1;
     RT( V63(iteration) ) *=  -1;
     RT( V65(iteration) ) *=  -1;
 
 
     IT( V21(iteration) ) *=  -1;
 
     IT( V41(iteration) ) *=  -1;
     IT( V43(iteration) ) *=  -1;
 
     IT( V61(iteration) ) *=  -1;
     IT( V63(iteration) ) *=  -1;
     IT( V65(iteration) ) *=  -1;
 
 
 
 
    /* Dkq auf Vkq umrechnen            */
    /*                                  */
    /* Vkq = Dkq * epsilon_kq * theta_k */
    /*                                  */
    /*                                  */
 
     e_4f = E4f( ionennr );
 
 
     RT( V20(iteration) ) *=  epn2n(0) * alpha_J[ e_4f ];
     RT( V21(iteration) ) *=  epn1n(1) * alpha_J[ e_4f ];
     RT( V22(iteration) ) *=  epn0n(2) * alpha_J[ e_4f ];
 
     RT( V40(iteration) ) *=  epn4n(0) * beta_J[ e_4f ];
     RT( V41(iteration) ) *=  epn3n(1) * beta_J[ e_4f ];
     RT( V42(iteration) ) *=  epn2n(2) * beta_J[ e_4f ];
     RT( V43(iteration) ) *=  epn1n(3) * beta_J[ e_4f ];
     RT( V44(iteration) ) *=  epn0n(4) * beta_J[ e_4f ];
 
     RT( V60(iteration) ) *=  epn6n(0) * gamma_J[ e_4f ];
     RT( V61(iteration) ) *=  epn5n(1) * gamma_J[ e_4f ];
     RT( V62(iteration) ) *=  epn4n(2) * gamma_J[ e_4f ];
     RT( V63(iteration) ) *=  epn3n(3) * gamma_J[ e_4f ];
     RT( V64(iteration) ) *=  epn2n(4) * gamma_J[ e_4f ];
     RT( V65(iteration) ) *=  epn1n(5) * gamma_J[ e_4f ];
     RT( V66(iteration) ) *=  epn0n(6) * gamma_J[ e_4f ];
 
 
 
     IT( V21(iteration) ) *=  epn1n(1) * alpha_J[ e_4f ];
     IT( V22(iteration) ) *=  epn0n(2) * alpha_J[ e_4f ];
 
     IT( V41(iteration) ) *=  epn3n(1) * beta_J[ e_4f ];
     IT( V42(iteration) ) *=  epn2n(2) * beta_J[ e_4f ];
     IT( V43(iteration) ) *=  epn1n(3) * beta_J[ e_4f ];
     IT( V44(iteration) ) *=  epn0n(4) * beta_J[ e_4f ];
 
     IT( V61(iteration) ) *=  epn5n(1) * gamma_J[ e_4f ];
     IT( V62(iteration) ) *=  epn4n(2) * gamma_J[ e_4f ];
     IT( V63(iteration) ) *=  epn3n(3) * gamma_J[ e_4f ];
     IT( V64(iteration) ) *=  epn2n(4) * gamma_J[ e_4f ];
     IT( V65(iteration) ) *=  epn1n(5) * gamma_J[ e_4f ];
     IT( V66(iteration) ) *=  epn0n(6) * gamma_J[ e_4f ];
 
 
 
    if( symmetrienr==8 ){
       if( ABSD(2.0*epn0n(4)*v44 - 5.0*epn4n(0)*v40) > EPS1)
       read_error(52,fp,name);
       if( ABSD(2.0*epn2n(4)*v64 + 21.0*epn6n(0)*v60) > EPS1)
       read_error(53,fp,name);
    }
 
    HMAG(iteration) = readBmag(fp,name,modus,myB,iteration,buffer_size,string);
    fclose(fp);
  } 
 
    return( iteration );
}
/*------------------------------------------------------------------------------
                          read_Wkq()
------------------------------------------------------------------------------*/
ITERATION *read_Wkq(name,vsymmetrienr_vor)  /* Wkq aus file name lesen */
    CHAR *name;
    INT  vsymmetrienr_vor;/* falls symmetrienr nicht vorgegeben  */
                          /* dann vsymmetrienr_vor <  0          */
{
    FILE      *fp,*fopen();
    INT       anz_nn,dimj/*,zwei_j*/,ionennr,symmetrienr,e_4f;
    INT       buffer_size=381,einheitnr_in,einheitnr_out;
    DOUBLE    versionsnummer;
    DOUBLE    myB,f2,f4,f6;
    
    DOUBLE    sin(),cos();
    DOUBLE    a_tof(),w40,w44,w60,w64,sqrt(),temperatur;
    CHAR  /*  *einheit_in,*einheit_out,*/modus;
    CHAR      *ion;
    CHAR      c,*string,*line,*fgets(),*a_tos();
    ITERATION *iteration,*iter_alloc();
    ITERATION *auswahlregel();
    STEVENS   *calc_Pkq();
    MATRIX    *readBmag();
 
    printf("Reading file %s ....\n",name);
    string   = STRING_ALLOC(buffer_size);
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 1.==== */
    versionsnummer =a_tof( line , 11,23);
    fclose(fp);
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    line=fgets( string , buffer_size , fp ); /* leerzeile                  */
    line=fgets( string , buffer_size , fp ); /* Kristallfeldpara....       */
    c = VALUE(line,33);
    einheitnr_in = is_einheit_imp(c);
    if( einheitnr_in == NICHTIMP )
         read_error(21,fp,name);
/*  einheit_in = EINHEITIMP[ einheitnr_in ].einheit; */
    myB        = EINHEITIMP[ einheitnr_in ].myB;
 
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    line=fgets( string , buffer_size , fp ); /* : Energieeigenwerte......*/
    c = VALUE(line,31);
    einheitnr_out= is_einheit_imp(c);
    if( einheitnr_out== NICHTIMP )
         read_error(21,fp,name);
/*  einheit_out= EINHEITIMP[ einheitnr_out].einheit; */
 
    line=fgets( string , buffer_size , fp ); /* : Temperatur der Probe :.*/
    temperatur = a_tof(line,31,61);
 
    line=fgets( string , buffer_size , fp ); /* : Betrachtetes Aufion....*/
    ion =a_tos( line , 31,56);
 
    line=fgets( string , buffer_size , fp ); /* : Symmetrie ... Symmnr...*/
    symmetrienr =a_toi( line , 53,61);
    if( vsymmetrienr_vor >= 0 ) symmetrienr = vsymmetrienr_vor;
 
    anz_nn    = 0;
    ionennr   = isimplementiert(ion);
    dimj      = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

/*  zwei_j    = dimj - 1 ; */
    iteration = iter_alloc(dimj,anz_nn);
 
 
    ANZ_NN(   iteration)     =  anz_nn;
    DIMJ(     iteration)     =  dimj;
    GJ(       iteration)     =  IONENIMP[ ionennr ].gj;
    IONNAME(  iteration)     =  IONENIMP[ ionennr ].ionname;
    EINHEITNRIN( iteration)     =  einheitnr_in;
    EINHEITNROUT(iteration)     =  einheitnr_out;
    PKQ(    iteration)     =  calc_Pkq( dimj );
    SYMMETRIENR(iteration) =  symmetrienr;
    TEMPERATUR( iteration) =  temperatur;
    EFVERSION(  iteration) =  versionsnummer;
 
    line=fgets( string , buffer_size , fp ); /* Magnetfeld: . Symmnr...  */
    c=VALUE(line,31);
    switch(c){
        case 'a':  modus = 'a';
                   break;
        case 'n':  modus = 'n';
                   break;
        default : read_error(24,fp,name); modus=0;
    }
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 3.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 4.==== */
 
    line=fgets( string , buffer_size , fp ); /*|RE W20 ............. */
    RT( V20(iteration) ) = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W21 .. |IM W21 .. */
    RT( V21(iteration) ) = a_tof(line, 9,30);
    IT( V21(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W22 .. |IM W22 .. */
    RT( V22(iteration) ) = a_tof(line, 9,30);
    IT( V22(iteration) ) = a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 5.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 6.==== */
    line=fgets( string , buffer_size , fp ); /*|RE W40 ............. */
    RT( V40(iteration) ) = w40 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W41 .. |IM W41 .. */
    RT( V41(iteration) ) = a_tof(line, 9,30);
    IT( V41(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE V42 .. |IM W42 .. */
    RT( V42(iteration) ) = a_tof(line, 9,30);
    IT( V42(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W43 .. |IM W43 .. */
    RT( V43(iteration) ) = a_tof(line, 9,30);
    IT( V43(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W44 .. |IM W44 .. */
    RT( V44(iteration) ) = w44 = a_tof(line, 9,30);
    IT( V44(iteration) ) = a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 7.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 8.==== */
    line=fgets( string , buffer_size , fp ); /*|RE W60 ............. */
    RT( V60(iteration) ) = w60 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W61 .. |IM W61 .. */
    RT( V61(iteration) ) = a_tof(line, 9,30);
    IT( V61(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W62 .. |IM W62 .. */
    RT( V62(iteration) ) = a_tof(line, 9,30);
    IT( V62(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W63 .. |IM W63 .. */
    RT( V63(iteration) ) = a_tof(line, 9,30);
    IT( V63(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W64 .. |IM W64 .. */
    RT( V64(iteration) ) = w64 = a_tof(line, 9,30);
    IT( V64(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W65 .. |IM W65 .. */
    RT( V65(iteration) ) = a_tof(line, 9,30);
    IT( V65(iteration) ) = a_tof(line,40,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE W66 .. |IM W66 .. */
    RT( V66(iteration) ) = a_tof(line, 9,30);
    IT( V66(iteration) ) = a_tof(line,40,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*9.==== */
 
    /* Auswahlregeln fuer Vkq beachten */
    iteration = auswahlregel(iteration,symmetrienr);
 
 
 
     /* Wkq auf Vkq umrechnen                */
     /*                                      */
     /*              2S+1         k          */
     /* Vkq = alpha (    L  ) * <r > * W     */
     /*            k      J             kq   */
     /*                                      */
 
     e_4f = E4f( ionennr );
 
     f2= A0_BOHR*A0_BOHR;
     RT( V20(iteration) ) *=  alpha_J[ e_4f ] * r2(ionennr)*f2;
     RT( V21(iteration) ) *=  alpha_J[ e_4f ] * r2(ionennr)*f2;
     RT( V22(iteration) ) *=  alpha_J[ e_4f ] * r2(ionennr)*f2;
 
     f4= f2 * f2;
     RT( V40(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr)*f4;
     RT( V41(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr)*f4;
     RT( V42(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr)*f4;
     RT( V43(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr)*f4;
     RT( V44(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr)*f4;
 
     f6= f4 * f2;
     RT( V60(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     RT( V61(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     RT( V62(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     RT( V63(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     RT( V64(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     RT( V65(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     RT( V66(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
 
 
 
     f2= A0_BOHR*A0_BOHR;
     IT( V21(iteration) ) *=  alpha_J[ e_4f ] * r2(ionennr)*f2;
     IT( V22(iteration) ) *=  alpha_J[ e_4f ] * r2(ionennr)*f2;
 
     f4= f2 * f2;
     IT( V41(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr)*f4;
     IT( V42(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr)*f4;
     IT( V43(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr)*f4;
     IT( V44(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr)*f4;
 
     f6= f4 * f2;
     IT( V61(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     IT( V62(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     IT( V63(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     IT( V64(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     IT( V65(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
     IT( V66(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr)*f6;
 
 
 
    if( symmetrienr==8 ){
       if(ABSD( 2.0*w44- 5.0*w40 ) > EPS1) read_error(28,fp,name);
       if(ABSD( 2.0*w64+21.0*w60 ) > EPS1) read_error(29,fp,name);
    }
 
 
    HMAG(iteration) = readBmag(fp,name,modus,myB,iteration,buffer_size,string);
 
    fclose(fp);
    return( iteration );
}
/*------------------------------------------------------------------------------
                          read_Akq()
------------------------------------------------------------------------------*/
ITERATION *read_Akq(name,vsymmetrienr_vor)  /* Akq aus file name lesen */
    CHAR *name;
    INT  vsymmetrienr_vor;/* falls symmetrienr nicht vorgegeben  */
                          /* dann vsymmetrienr_vor <  0          */
{
    FILE      *fp,*fopen();
    INT       anz_nn,dimj/*,zwei_j*/,ionennr,symmetrienr,e_4f;
    INT       buffer_size=381,einheitnr_in,einheitnr_out;
    DOUBLE    versionsnummer;
    DOUBLE    myB/*,f2,f4,f6*/; 
    
    DOUBLE    sin(),cos();
    DOUBLE    omegan0n(),omegan1n(),omegan2n();
    DOUBLE    omegan3n(),omegan4n(),omegan5n();
    DOUBLE    omegan6n();
    DOUBLE    a_tof(),a40,a44,a60,a64,sqrt(),temperatur;
    CHAR  /*  *einheit_in,*einheit_out,*/modus;
    CHAR      *ion;
    CHAR      c,*string,*line,*fgets(),*a_tos();
    ITERATION *iteration,*iter_alloc();
    ITERATION *auswahlregel();
    STEVENS   *calc_Pkq();
    MATRIX    *readBmag();
 
    printf("Reading file %s ....\n",name);
    string   = STRING_ALLOC(buffer_size);
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 1.==== */
    versionsnummer =a_tof( line , 11,23);
    fclose(fp);
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    line=fgets( string , buffer_size , fp ); /* leerzeile                  */
    line=fgets( string , buffer_size , fp ); /* Kristallfeldpara....       */
    c = VALUE(line,33 );
    einheitnr_in = is_einheit_imp(c);
    if( einheitnr_in == NICHTIMP )
         read_error(21,fp,name);
/*  einheit_in = EINHEITIMP[ einheitnr_in ].einheit; */
    myB        = EINHEITIMP[ einheitnr_in ].myB;
 
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    line=fgets( string , buffer_size , fp ); /* : Energieeigenwerte......*/
    c = VALUE(line,31);
    einheitnr_out= is_einheit_imp(c);
    if( einheitnr_out== NICHTIMP )
         read_error(21,fp,name);
/*  einheit_out= EINHEITIMP[ einheitnr_out].einheit; */
 
    line=fgets( string , buffer_size , fp ); /* : Temperatur der Probe :.*/
    temperatur = a_tof(line,31,61);
 
    line=fgets( string , buffer_size , fp ); /* : Betrachtetes Aufion....*/
    ion =a_tos( line , 31,56);
 
    line=fgets( string , buffer_size , fp ); /* : Symmetrie ... Symmnr...*/
    symmetrienr =a_toi( line , 53,61);
    if( vsymmetrienr_vor >= 0 )            symmetrienr = vsymmetrienr_vor;
    if( ! isreell( symmetrienr ,ion ) )  { fclose(fp);
                                           Bkq_error( AKQNAME,ion,symmetrienr);}
 
 
    anz_nn    = 0;
    ionennr   = isimplementiert(ion);
    dimj      = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

/*  zwei_j    = dimj - 1 ; */
    iteration = iter_alloc(dimj,anz_nn);
 
    ANZ_NN(   iteration)     =  anz_nn;
    DIMJ(     iteration)     =  dimj;
    GJ(       iteration)     =  IONENIMP[ ionennr ].gj;
    IONNAME(  iteration)     =  IONENIMP[ ionennr ].ionname;
    EINHEITNRIN( iteration)     =  einheitnr_in;
    EINHEITNROUT(iteration)     =  einheitnr_out;
    PKQ(      iteration)     =  calc_Pkq( dimj );
    SYMMETRIENR(iteration) =  symmetrienr;
    TEMPERATUR( iteration) =  temperatur;
    EFVERSION(  iteration) =  versionsnummer;
 
 
    line=fgets( string , buffer_size , fp ); /* Magnetfeld: . Symmnr...  */
    c=VALUE(line,31);
    switch(c){
        case 'a':  modus = 'a';
                   break;
        case 'n':  modus = 'n';
                   break;
        default : read_error(24,fp,name); modus=0;
    }
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 3.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 4.==== */
 
    line=fgets( string , buffer_size , fp ); /*|RE A20 ............. */
    RT( V20(iteration) ) = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A21 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR21(iteration) = -1.0;
    RT( V21(iteration) ) = a_tof(line, 9,61);
    IT( V21(iteration) ) = 0.0;
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A22 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR22(iteration) = -1.0;
    RT( V22(iteration) ) = a_tof(line, 9,61);
    IT( V22(iteration) ) = 0.0;
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 5.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 6.==== */
    line=fgets( string , buffer_size , fp ); /*|RE A40 ............. */
    RT( V40(iteration) ) = a40 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A41 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR41(iteration) = -1.0;
    RT( V41(iteration) ) = a_tof(line, 9,61);
    IT( V41(iteration) ) = 0.0;
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A42 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR42(iteration) = -1.0;
    RT( V42(iteration) ) = a_tof(line, 9,61);
    IT( V42(iteration) ) = 0.0;
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A43 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR43(iteration) = -1.0;
    RT( V43(iteration) ) = a_tof(line, 9,61);
    IT( V43(iteration) ) = 0.0;
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A44 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR44(iteration) = -1.0;
    RT( V44(iteration) ) = a44 = a_tof(line, 9,61);
    IT( V44(iteration) ) = 0.0;
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 7.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 8.==== */
    line=fgets( string , buffer_size , fp ); /*|RE A60 ............. */
    RT( V60(iteration) ) = a60 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A61 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR61(iteration) = -1.0;
    RT( V61(iteration) ) = a_tof(line, 9,61);
    IT( V61(iteration) ) = 0.0;
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A62 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR62(iteration) = -1.0;
    RT( V62(iteration) ) = a_tof(line, 9,61);
    IT( V62(iteration) ) = 0.0;
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A63 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR63(iteration) = -1.0;
    RT( V63(iteration) ) = a_tof(line, 9,61);
    IT( V63(iteration) ) = 0.0;
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A64 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR64(iteration) = -1.0;
    RT( V64(iteration) ) = a64 = a_tof(line, 9,61);
    IT( V64(iteration) ) = 0.0;
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A65 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR65(iteration) = -1.0;
    RT( V65(iteration) ) = a_tof(line, 9,61);
    IT( V65(iteration) ) = 0.0;
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE A66 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR66(iteration) = -1.0;
    RT( V66(iteration) ) = a_tof(line, 9,61);
    IT( V66(iteration) ) = 0.0;
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*9.==== */
 
    /* Auswahlregeln fuer Akq beachten */
    iteration = auswahlregel(iteration,symmetrienr);
 
 
 
     /* Akq auf Bkq umrechnen                */
     /*                                      */
     /*              2S+1         k          */
     /* Bkq = theta (    L  ) * <r > * A     */
     /*            k      J             kq   */
     /*                                      */
 
     e_4f = E4f( ionennr );
 
/*   f2= A0_BOHR*A0_BOHR; */
     RT( V20(iteration) ) *=  alpha_J[ e_4f ] * r2(ionennr);
     RT( V21(iteration) ) *=  alpha_J[ e_4f ] * r2(ionennr);
     RT( V22(iteration) ) *=  alpha_J[ e_4f ] * r2(ionennr);
 
/*   f4= f2 * f2; */
     RT( V40(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr);
     RT( V41(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr);
     RT( V42(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr);
     RT( V43(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr);
     RT( V44(iteration) ) *=   beta_J[ e_4f ] * r4(ionennr);
 
/*   f6= f4 * f2; */
     RT( V60(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr);
     RT( V61(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr);
     RT( V62(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr);
     RT( V63(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr);
     RT( V64(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr);
     RT( V65(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr);
     RT( V66(iteration) ) *=  gamma_J[ e_4f ] * r6(ionennr);
 
 
    if( symmetrienr==8 ){
       if(ABSD( a44- 5.0*a40 ) > EPS1) read_error(25,fp,name);
       if(ABSD( a64+21.0*a60 ) > EPS1) read_error(26,fp,name);
    }
 
 
     /* Bkq auf Vkq umrechnen       */
     /*                             */
     /*       | Bk0  ,         q =0 */
     /* Vkq = |                     */
     /*       | Bkq/2/omega_kq ,q >0*/
     /*                             */
 
 
     RT( V21(iteration) ) /=  2 * omegan1n(1);
     RT( V22(iteration) ) /=  2 * omegan0n(2);
 
     RT( V41(iteration) ) /=  2 * omegan3n(1);
     RT( V42(iteration) ) /=  2 * omegan2n(2);
     RT( V43(iteration) ) /=  2 * omegan1n(3);
     RT( V44(iteration) ) /=  2 * omegan0n(4);
 
     RT( V61(iteration) ) /=  2 * omegan5n(1);
     RT( V62(iteration) ) /=  2 * omegan4n(2);
     RT( V63(iteration) ) /=  2 * omegan3n(3);
     RT( V64(iteration) ) /=  2 * omegan2n(4);
     RT( V65(iteration) ) /=  2 * omegan1n(5);
     RT( V66(iteration) ) /=  2 * omegan0n(6);
 
 
 
 
 
 
    HMAG(iteration) = readBmag(fp,name,modus,myB,iteration,buffer_size,string);
 
    fclose(fp);
    return( iteration );
}
/*------------------------------------------------------------------------------
                          read_Bkqnew()
 ... routine without file reading
------------------------------------------------------------------------------*/
ITERATION *read_Bkqnew(ion)  /* Vkq aus file name lesen */
    CHAR *ion;
{
    FILE      *fopen();
    INT       anz_nn,dimj/*,zwei_j*/,ionennr,symmetrienr;
    INT    /* buffer_size=381,*/einheitnr_in,einheitnr_out;
    DOUBLE    versionsnummer;
    DOUBLE    myB;
   
    DOUBLE    sin(),cos();
    DOUBLE    omegan0n(),omegan1n(),omegan2n();
    DOUBLE    omegan3n(),omegan4n(),omegan5n();
    DOUBLE    omegan6n();
    DOUBLE    a_tof(),b40,b44,b60,b64,sqrt();
/*  CHAR      *einheit_in,*einheit_out; */
    CHAR      c/*,*string*/,*fgets(),*a_tos();
    ITERATION *iteration,*iter_alloc();
    ITERATION *auswahlregel();
    STEVENS   *calc_Pkq();
    MATRIX    *readBmag();
 
/*  string   = STRING_ALLOC(buffer_size); */
    versionsnummer =VERSION;
    c = 'm';
 /* units are in meV */
    einheitnr_in = is_einheit_imp(c);
/*  einheit_in = EINHEITIMP[ einheitnr_in ].einheit; */
    myB        = EINHEITIMP[ einheitnr_in ].myB;
 
    einheitnr_out= einheitnr_in;
/*  einheit_out= EINHEITIMP[ einheitnr_out].einheit; */
 
    symmetrienr =0;
 /*triclinic */
 
  if(strncmp(ion,"S=",2)==0) /* S=... ion !! extract dimj from string ion  equivalent to J=*/
     {ion[0]='J';     }    
    anz_nn    = 0;
    ionennr   = isimplementiert(ion);
    dimj      = IONENIMP[ ionennr ].dimj;
  if(strncmp(ion,"J=",2)==0) /* J=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

/*  zwei_j    = dimj - 1 ; */
 
    iteration = iter_alloc(dimj,anz_nn);
 
    ANZ_NN(   iteration)     =  anz_nn;
    DIMJ(     iteration)     =  dimj;
    GJ(       iteration)     =  IONENIMP[ ionennr ].gj;
    IONNAME(  iteration)     =  IONENIMP[ ionennr ].ionname;
    IONENNUMMER(iteration)   = ionennr;
    EINHEITNRIN( iteration)     =  einheitnr_in;
    EINHEITNROUT(iteration)     =  einheitnr_out;
    PKQ(      iteration)     =  calc_Pkq( dimj );
    SYMMETRIENR(iteration) =  symmetrienr;
    TEMPERATUR( iteration) =  10;
    EFVERSION(  iteration) =  versionsnummer;
 
 
    RT( V20(iteration) ) = 0;
    RT( V21(iteration) ) = 0;
    IT( V21(iteration) ) = 0.0;
    RT( V22(iteration) ) = 0;
    IT( V22(iteration) ) = 0.0;
    RT( V40(iteration) ) = b40 = 0;
    RT( V41(iteration) ) = 0;
    IT( V41(iteration) ) = 0.0;
    RT( V42(iteration) ) =0;
    IT( V42(iteration) ) = 0.0;
    RT( V43(iteration) ) = 0.0;
    IT( V43(iteration) ) = 0.0;
    RT( V44(iteration) ) = b44 = 0;
    IT( V44(iteration) ) = 0.0;
    RT( V60(iteration) ) = b60 = 0;
    RT( V61(iteration) ) = 0.0;
    IT( V61(iteration) ) = 0.0;
 
    RT( V62(iteration) ) = 0;
    IT( V62(iteration) ) = 0.0;
    RT( V63(iteration) ) = 0;
    IT( V63(iteration) ) = 0.0;
    RT( V64(iteration) ) = b64 = 0;
    IT( V64(iteration) ) = 0.0;
 
    RT( V65(iteration) ) = 0;
    IT( V65(iteration) ) = 0.0;
    RT( V66(iteration) ) = 0;
    IT( V66(iteration) ) = 0.0;
    /* Auswahlregeln fuer Bkq beachten */
    iteration = auswahlregel(iteration,symmetrienr);
 
       B1(iteration) = 0.0;
        B2(iteration) = 0.0;
        B3(iteration) = 0.0;
       B1S(iteration) = 0.0;
        B2S(iteration) = 0.0;
        B3S(iteration) = 0.0;
        HMAG(iteration) = calc_Bmag( DIMJ(iteration),GJ(iteration),myB,
                                     B1(iteration),
                                     B2(iteration),
                                     B3(iteration)  );
 
        B1MOL(iteration) = 0.0;
        B2MOL(iteration) = 0.0;
        B3MOL(iteration) = 0.0;
 
        BMOL( iteration) = 0.0;
        PHI(  iteration) = 0.0;
        THETA(iteration) = 0.0;
 
    return( iteration );
}
/*------------------------------------------------------------------------------
                          read_Bkq()
------------------------------------------------------------------------------*/

ITERATION *read_Bkq(name,vsymmetrienr_vor)  /* Vkq aus file name lesen */
    CHAR *name;
    INT  vsymmetrienr_vor;/* falls symmetrienr nicht vorgegeben  */
                          /* dann vsymmetrienr_vor <  0          */
{
    FILE      *fp,*fopen();
    INT       anz_nn,dimj/*,zwei_j*/,ionennr,symmetrienr;
    INT       buffer_size=381,einheitnr_in,einheitnr_out;
    DOUBLE    versionsnummer;
    DOUBLE    myB;
    DOUBLE    sin(),cos();
    DOUBLE    omegan0n(),omegan1n(),omegan2n();
    DOUBLE    omegan3n(),omegan4n(),omegan5n();
    DOUBLE    omegan6n();
    DOUBLE    a_tof(),b40,b44,b60,b64,sqrt(),temperatur;
    CHAR  /*  *einheit_in,*einheit_out,*/modus;
    CHAR      *ion=0;
    CHAR      c,*string,*line,*fgets(),*a_tos();
    ITERATION *iteration,*iter_alloc();
    ITERATION *auswahlregel();
    STEVENS   *calc_Pkq();
    MATRIX    *readBmag();
 
    printf("Reading file %s ....\n",name);
    string   = STRING_ALLOC(buffer_size);
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
    line=fgets( string , buffer_size , fp );
    fclose(fp);
    if(strncmp(line,"#!cfield",8)==0||strncmp(line,"#!so1ion",8)==0||strncmp(line,"#!MODULE",8)==0)
   {iteration=read_new_format('B',iteration,name,vsymmetrienr_vor);
    }
    else
   { /*read classical input file bkq.parmeter*/
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 1.==== */
    versionsnummer =a_tof( line , 11,23);
    fclose(fp);
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    line=fgets( string , buffer_size , fp ); /* leerzeile                  */
    line=fgets( string , buffer_size , fp ); /* Kristallfeldpara....       */
    c = VALUE(line,33 );
    einheitnr_in = is_einheit_imp(c);
    if( einheitnr_in == NICHTIMP )
         read_error(21,fp,name);
/*  einheit_in = EINHEITIMP[ einheitnr_in ].einheit; */
    myB        = EINHEITIMP[ einheitnr_in ].myB;
 
 
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    line=fgets( string , buffer_size , fp ); /* : Energieeigenwerte......*/
    c = VALUE(line,31);
    einheitnr_out= is_einheit_imp(c);
    if( einheitnr_out== NICHTIMP )
         read_error(21,fp,name);
/*  einheit_out= EINHEITIMP[ einheitnr_out].einheit; */
 
    line=fgets( string , buffer_size , fp ); /* : Temperatur der Probe :.*/
    temperatur = a_tof(line,31,61);
 
    line=fgets( string , buffer_size , fp ); /* : Betrachtetes Aufion....*/
    ion =a_tos( line , 31,56);
 
    line=fgets( string , buffer_size , fp ); /* : Symmetrie ... Symmnr...*/
    symmetrienr =a_toi( line , 53,61);
    if( vsymmetrienr_vor >= 0 )            symmetrienr = vsymmetrienr_vor;
/*  if( ! isreell( symmetrienr ,ion ) )  { fclose(fp);
                                           Bkq_error( BKQNAME,ion,symmetrienr);}
*/ 
 
    anz_nn    = 0;
    ionennr   = isimplementiert(ion);
    dimj      = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    
/*  zwei_j    = dimj - 1 ; */
 
    iteration = iter_alloc(dimj,anz_nn);
 
    ANZ_NN(   iteration)     =  anz_nn;
    DIMJ(     iteration)     =  dimj;
    GJ(       iteration)     =  IONENIMP[ ionennr ].gj;
    IONNAME(  iteration)     =  IONENIMP[ ionennr ].ionname;
    IONENNUMMER(iteration)   = ionennr;
    EINHEITNRIN( iteration)     =  einheitnr_in;
    EINHEITNROUT(iteration)     =  einheitnr_out;
    PKQ(      iteration)     =  calc_Pkq( dimj );
    SYMMETRIENR(iteration) =  symmetrienr;
    TEMPERATUR( iteration) =  temperatur;
    EFVERSION(  iteration) =  versionsnummer;
 
 
    line=fgets( string , buffer_size , fp ); /* Magnetfeld: . Symmnr...  */
    c=VALUE(line,31);
    switch(c){
        case 'a':  modus = 'a';
                   break;
        case 'n':  modus = 'n';
                   break;
        default : read_error(24,fp,name); modus=0;
    }
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 3.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 4.==== */
 
    line=fgets( string , buffer_size , fp ); /*|RE B20 ............. */
    RT( V20(iteration) ) = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B21 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR21(iteration) = -1.0;
    RT( V21(iteration) ) = a_tof(line, 9,30);
    IT( V21(iteration) ) = -a_tof(line, 39,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B22 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR22(iteration) = -1.0;
    RT( V22(iteration) ) = a_tof(line, 9,30);
    IT( V22(iteration) ) = -a_tof(line, 39,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 5.==== */

    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 6.==== */
    line=fgets( string , buffer_size , fp ); /*|RE B40 ............. */
    RT( V40(iteration) ) = b40 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B41 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR41(iteration) = -1.0;
    RT( V41(iteration) ) = a_tof(line, 9,30);
    IT( V41(iteration) ) = -a_tof(line, 39,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B42 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR42(iteration) = -1.0;
    RT( V42(iteration) ) = a_tof(line, 9,30);
    IT( V42(iteration) ) = -a_tof(line, 39,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B43 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR43(iteration) = -1.0;
    RT( V43(iteration) ) = a_tof(line, 9,30);
    IT( V43(iteration) ) = -a_tof(line, 39,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B44 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR44(iteration) = -1.0;
    RT( V44(iteration) ) = b44 = a_tof(line, 9,30);
    IT( V44(iteration) ) = -a_tof(line, 39,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 7.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 8.==== */
    line=fgets( string , buffer_size , fp ); /*|RE B60 ............. */
    RT( V60(iteration) ) = b60 = a_tof(line,9,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B61 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR61(iteration) = -1.0;
    RT( V61(iteration) ) = a_tof(line, 9,30);
    IT( V61(iteration) ) = -a_tof(line, 39,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B62 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR62(iteration) = -1.0;
    RT( V62(iteration) ) = a_tof(line, 9,30);
    IT( V62(iteration) ) = -a_tof(line, 39,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B63 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR63(iteration) = -1.0;
    RT( V63(iteration) ) = a_tof(line, 9,30);
    IT( V63(iteration) ) = -a_tof(line, 39,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B64 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR64(iteration) = -1.0;
    RT( V64(iteration) ) = b64 = a_tof(line, 9,30);
    IT( V64(iteration) ) = -a_tof(line, 39,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B65 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR65(iteration) = -1.0;
    RT( V65(iteration) ) = a_tof(line, 9,30);
    IT( V65(iteration) ) = -a_tof(line, 39,61);
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|RE B66 ............. */
    c = VALUE(line,7); if( c=='i'||c=='I' ) VOR66(iteration) = -1.0;
    RT( V66(iteration) ) = a_tof(line, 9,30);
    IT( V66(iteration) ) = -a_tof(line, 39,61);
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*9.==== */

    /* Auswahlregeln fuer Bkq beachten */
    iteration = auswahlregel(iteration,symmetrienr);

 
 
    if( symmetrienr==8 ){
       if(ABSD( b44- 5.0*b40 ) > EPS1) read_error(25,fp,name);
       if(ABSD( b64+21.0*b60 ) > EPS1) read_error(26,fp,name);
    }
 
 
     /* Bkq auf Vkq umrechnen       */
     /*                             */
     /*       | Bk0  ,         q =0 */
     /* Vkq = |                     */
     /*       | Bkq/2/omegakq, q >0 */
     /*                             */
 
 
     RT( V21(iteration) ) /=  2 * omegan1n(1);
     RT( V22(iteration) ) /=  2 * omegan0n(2);
 
     RT( V41(iteration) ) /=  2 * omegan3n(1);
     RT( V42(iteration) ) /=  2 * omegan2n(2);
     RT( V43(iteration) ) /=  2 * omegan1n(3);
     RT( V44(iteration) ) /=  2 * omegan0n(4);
 
     RT( V61(iteration) ) /=  2 * omegan5n(1);
     RT( V62(iteration) ) /=  2 * omegan4n(2);
     RT( V63(iteration) ) /=  2 * omegan3n(3);
     RT( V64(iteration) ) /=  2 * omegan2n(4);
     RT( V65(iteration) ) /=  2 * omegan1n(5);
     RT( V66(iteration) ) /=  2 * omegan0n(6);
 
 
     IT( V21(iteration) ) /=  2 * omegan1n(1);
     IT( V22(iteration) ) /=  2 * omegan0n(2);
 
     IT( V41(iteration) ) /=  2 * omegan3n(1);
     IT( V42(iteration) ) /=  2 * omegan2n(2);
     IT( V43(iteration) ) /=  2 * omegan1n(3);
     IT( V44(iteration) ) /=  2 * omegan0n(4);
 
     IT( V61(iteration) ) /=  2 * omegan5n(1);
     IT( V62(iteration) ) /=  2 * omegan4n(2);
     IT( V63(iteration) ) /=  2 * omegan3n(3);
     IT( V64(iteration) ) /=  2 * omegan2n(4);
     IT( V65(iteration) ) /=  2 * omegan1n(5);
     IT( V66(iteration) ) /=  2 * omegan0n(6);
 
 

 
 
    HMAG(iteration) = readBmag(fp,name,modus,myB,iteration,buffer_size,string);
    fclose(fp);
   }
    return( iteration );
}
/*------------------------------------------------------------------------------
                          read_xW()
------------------------------------------------------------------------------*/
ITERATION *read_xW(name,vsymmetrienr_vor)  /* x,W aus file name lesen */
    CHAR *name;
    INT  vsymmetrienr_vor;/* falls symmetrienr nicht vorgegeben  */
                          /* dann vsymmetrienr_vor <  0          */
{
    UNUSED_PARAMETER(vsymmetrienr_vor);

    FILE      *fp,*fopen();
    INT       anz_nn,dimj/*,zwei_j*/,ionennr,symmetrienr;
    INT       buffer_size=381,einheitnr_in,einheitnr_out;
    DOUBLE    versionsnummer;
    DOUBLE    myB;
    DOUBLE    x,W,f4,f6;
   
    DOUBLE    sin(),cos();
    DOUBLE    omegan0n(),omegan1n(),omegan2n();
    DOUBLE    omegan3n(),omegan4n(),omegan5n();
    DOUBLE    omegan6n();
    DOUBLE    a_tof(),sqrt(),temperatur;
    CHAR  /*  *einheit_in,*einheit_out,*/modus;
    CHAR      *ion;
    CHAR      c,*string,*line,*fgets(),*a_tos();
    ITERATION *iteration,*iter_alloc();
    ITERATION *auswahlregel();
    STEVENS   *calc_Pkq();
    MATRIX    *readBmag();
 
    printf("lese file %s ....\n",name);
    string   = STRING_ALLOC(buffer_size);
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 1.==== */
    versionsnummer =a_tof( line , 11,23);
    fclose(fp);
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    line=fgets( string , buffer_size , fp ); /* leerzeile                  */
    line=fgets( string , buffer_size , fp ); /* Kristallfeldpara....       */
    c = VALUE(line,33 );
    einheitnr_in = is_einheit_imp(c);
    if( einheitnr_in == NICHTIMP )
         read_error(21,fp,name);
/*  einheit_in = EINHEITIMP[ einheitnr_in ].einheit; */
    myB        = EINHEITIMP[ einheitnr_in ].myB;
 
 
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    line=fgets( string , buffer_size , fp ); /* : Energieeigenwerte......*/
    c = VALUE(line,31);
    einheitnr_out= is_einheit_imp(c);
    if( einheitnr_out== NICHTIMP )
         read_error(21,fp,name);
/*  einheit_out= EINHEITIMP[ einheitnr_out].einheit; */
 
    line=fgets( string , buffer_size , fp ); /* : Temperatur der Probe :.*/
    temperatur = a_tof(line,31,61);
 
    line=fgets( string , buffer_size , fp ); /* : Betrachtetes Aufion....*/
    ion =a_tos( line , 31,56);
 
    line=fgets( string , buffer_size , fp ); /* : Symmetrie ... Symmnr...*/
    symmetrienr =a_toi( line , 53,61);
    if( symmetrienr != 8 ) symmetrienr = 8;
 
    anz_nn    = 0;
    ionennr   = isimplementiert(ion);
    dimj      = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

/*  zwei_j    = dimj - 1 ; */
    iteration = iter_alloc(dimj,anz_nn);
 
 
    ANZ_NN(   iteration)     =  anz_nn;
    DIMJ(     iteration)     =  dimj;
    GJ(       iteration)     =  IONENIMP[ ionennr ].gj;
    IONNAME(  iteration)     =  IONENIMP[ ionennr ].ionname;
    EINHEITNRIN( iteration)     =  einheitnr_in;
    EINHEITNROUT(iteration)     =  einheitnr_out;
    PKQ(      iteration)     =  calc_Pkq( dimj );
    SYMMETRIENR(iteration) =  symmetrienr;
    TEMPERATUR( iteration) =  temperatur;
    EFVERSION(  iteration) =  versionsnummer;
 
    line=fgets( string , buffer_size , fp ); /* Magnetfeld: . Symmnr...  */
    c=VALUE(line,31);
    switch(c){
        case 'a':  modus = 'a';
                   break;
        case 'n':  modus = 'n';
                   break;
        default : read_error(24,fp,name); modus = 0;
    }
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 3.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 4.==== */
 
    line=fgets( string , buffer_size , fp ); /*|   x   :............ */
    x = a_tof(line,9,61);
    if( !( ABSD(x)<= 1.0) ){
        printf("\n Fehler :  Der Betrag von x ist nicht kleiner gleich 1 !\n");
        exit(1);
    }
 
    line=fgets( string , buffer_size , fp ); /*--------------------- */
 
    line=fgets( string , buffer_size , fp ); /*|   W   :............ */
    W = a_tof(line, 9,61);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 5.==== */
 
    f4   = F4(  ionennr );
    RT( V40(iteration) ) = x*W/f4;
    RT( V44(iteration) ) =   5 * RT( V40(iteration) );
 
    f6   = F6(  ionennr );
    RT( V60(iteration) ) =  ( 1.0 - ABSD(x) )*W/f6;
    RT( V64(iteration) ) = -21 * RT( V60(iteration) );
 
     /* Bkq auf Vkq umrechnen       */
     /*                             */
     /*       | Bk0  ,         q =0 */
     /* Vkq = |                     */
     /*       | Bkq/2/omegakq ,q >0 */
     /*                             */
 
     RT( V44(iteration) ) /=  2 * omegan0n(4);
     RT( V64(iteration) ) /=  2 * omegan2n(4);
 
 
    /* Auswahlregeln fuer Vkq beachten */
    iteration = auswahlregel(iteration,symmetrienr);
 
 
 
 
 
    HMAG(iteration) = readBmag(fp,name,modus,myB,iteration,buffer_size,string);
 
    fclose(fp);
    return( iteration );
}
 
/*------------------------------------------------------------------------------
                                  read_nn()
------------------------------------------------------------------------------*/
UMGEBUNG *read_nn(name) /* Lese Eingabefile name der Umgebungsatome */
    CHAR *name;
{
    UMGEBUNG *umgebung;
    DOUBLE   a_tof(),x1,x2,x3=0.,temperatur;
    FILE     *fp,*fopen();
    INT      buffer_size = 81;
    INT      i,anz_nn,nummer;
    CHAR     c,*string,*line,*fgets(),*ion,*a_tos();
 
    printf("lese file %s ....\n",name);
    string   = STRING_ALLOC(buffer_size);
    umgebung = UMGEBUNG_ALLOC(1); /* Speicher fuer Struktur UMGEBUNG holen  */
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    line=fgets( string , buffer_size , fp ); /* : Temperatur der Probe :.*/
    temperatur = a_tof(line,31,61);
    TEMPERATUR(umgebung) = temperatur;
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    line=fgets( string , buffer_size , fp ); /* : die naechsten Nachb ...*/
    ion =a_tos( line , 31,56);
 
    IONENNR(umgebung) = isimplementiert(ion); /* ion implementiert ? */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 3.==== */
 
    line=fgets( string , buffer_size , fp ); /* : die Koordinaten der ...*/
    line=fgets( string , buffer_size , fp );/* : º P... v º R... v º S...*/
    c   = VALUE(line,2);
    switch(c){
         case 'K' : MODUS(umgebung) = 'r';
                    break;
         case 'S' : MODUS(umgebung) = 's';
                    break;
         case 'P' : MODUS(umgebung) = 'p';
                    break;
         default  : read_error(3,fp,name);
    }
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*4.==== */
 
    for( i=1 ; i<= 3 ; ++i)/* Kopf der Magnetfeldtabelle ueberlesen */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );
    line=fgets( string , buffer_size , fp );   /* : º  h1 º h2 º h3 º*/
    B1(umgebung) = x1 = a_tof(line, 2,13);
    B2(umgebung) = x2 = a_tof(line,15,26);
    if(MODUS(umgebung) !='p')  B3(umgebung) = x3 = a_tof(line,28,39);
    isinlimits(fp,name ,0, x1,x2,x3,MODUS(umgebung) );
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*8.==== */
    line=fgets( string , buffer_size , fp );   /* : º_nnn__  */
 
    ANZ_NN(umgebung) = anz_nn = a_toi( line , 1 , 6 );
    Q_P(umgebung)    = DOUBLE_ALLOC(anz_nn);
    X1_P(umgebung)   = DOUBLE_ALLOC(anz_nn);
    X2_P(umgebung)   = DOUBLE_ALLOC(anz_nn);
    X3_P(umgebung)   = DOUBLE_ALLOC(anz_nn);
 
 
    for( i=1 ; i<=2 ; ++i )  /* 2.+ 3.== der kleinen Tabelle ueberlesen*/
         while(  *(line=fgets( string , buffer_size , fp )) != '='  );
 
    for( i=1 ; i<= anz_nn ; ++i ){/* Ort und Ladung lesen */
 
         line=fgets( string , buffer_size , fp ); /* º nnn º ... */
 
         if(  (nummer = a_toi(line,1,5)) != i) read_error(4,fp,name)     ;
 
         Q(umgebung,i)  = a_tof(line, 7,18);
         X1(umgebung,i) = x1 = a_tof(line,20,31);
         X2(umgebung,i) = x2 = a_tof(line,33,44);
         if(MODUS(umgebung) !='p')  X3(umgebung,i) = x3 = a_tof(line,46,57);
         isinlimits(fp,name ,i, x1,x2,x3,MODUS(umgebung) );
         if( i==anz_nn ) break;
         while(  *(line=fgets( string , buffer_size , fp )) != '-'  );
    }
 
    fclose(fp);
    return( umgebung );
}
/*------------------------------------------------------------------------------
                                stern_setzen()
------------------------------------------------------------------------------*/
TABELLE *stern_setzen(tabelle,symmetrienr,dimj)
    TABELLE *tabelle;
    INT     symmetrienr;
    INT     dimj;
{
    INT zwei_j;
 
    zwei_j = dimj - 1;
 
    switch(symmetrienr){
 
       case  0  :  switch( zwei_j - 2 >= 0 ){
                        case JA : T20 = T20r;
                                  T21 = T21c;
                                  T22 = T22c;
                                  break;
                        default : T20 = T20n;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                   }
                   switch( zwei_j - 4 >= 0 ){
                        case JA : T40 = T40r;
                                  T41 = T41c;
                                  T42 = T42c;
                                  T43 = T43c;
                                  T44 = T44c;
                                  break;
                        default : T40 = T40n;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44n;
                                  break;
                   }
                   switch( zwei_j - 6 >= 0 ){
                        case JA : T60 = T60r;
                                  T61 = T61c;
                                  T62 = T62c;
                                  T63 = T63c;
                                  T64 = T64c;
                                  T65 = T65c;
                                  T66 = T66c;
                                  break;
                        default : T60 = T60n;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                   }
                   break;
 
       case  1  :  switch( zwei_j - 2 >= 0 ){
                        case JA : T20 = T20r;
                                  T21 = T21n;
                                  T22 = T22c;
                                  break;
                        default : T20 = T20n;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                   }
                   switch( zwei_j - 4 >= 0 ){
                        case JA : T40 = T40r;
                                  T41 = T41n;
                                  T42 = T42c;
                                  T43 = T43n;
                                  T44 = T44c;
                                  break;
                        default : T40 = T40n;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44n;
                                  break;
                   }
                   switch( zwei_j - 6 >= 0 ){
                        case JA : T60 = T60r;
                                  T61 = T61n;
                                  T62 = T62c;
                                  T63 = T63n;
                                  T64 = T64c;
                                  T65 = T65n;
                                  T66 = T66c;
                                  break;
                        default : T60 = T60n;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                   }
                   break;
 
       case  2  :  switch( zwei_j - 2 >= 0 ){
                        case JA : T20 = T20r;
                                  T21 = T21n;
                                  T22 = T22c;
                                  break;
                        default : T20 = T20n;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                   }
                   switch( zwei_j - 4 >= 0 ){
                        case JA : T40 = T40r;
                                  T41 = T41n;
                                  T42 = T42c;
                                  T43 = T43n;
                                  T44 = T44c;
                                  break;
                        default : T40 = T40n;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44n;
                                  break;
                   }
                   switch( zwei_j - 6 >= 0 ){
                        case JA : T60 = T60r;
                                  T61 = T61n;
                                  T62 = T62c;
                                  T63 = T63n;
                                  T64 = T64c;
                                  T65 = T65n;
                                  T66 = T66c;
                                  break;
                        default : T60 = T60n;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                   }
                   break;
 
       case  3  :  switch( zwei_j - 2 >= 0 ){
                        case JA : T20 = T20r;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                        default : T20 = T20n;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                   }
                   switch( zwei_j - 4 >= 0 ){
                        case JA : T40 = T40r;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44c;
                                  break;
                        default : T40 = T40n;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44n;
                                  break;
                   }
                   switch( zwei_j - 6 >= 0 ){
                        case JA : T60 = T60r;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64c;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                        default : T60 = T60n;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                   }
                   break;
 
       case  4  :  switch( zwei_j - 2 >= 0 ){
                        case JA : T20 = T20r;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                        default : T20 = T20n;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                   }
                   switch( zwei_j - 4 >= 0 ){
                        case JA : T40 = T40r;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44c;
                                  break;
                        default : T40 = T40n;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44n;
                                  break;
                   }
                   switch( zwei_j - 6 >= 0 ){
                        case JA : T60 = T60r;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64c;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                        default : T60 = T60n;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                   }
                   break;
 
       case  5  :  switch( zwei_j - 2 >= 0 ){
                        case JA : T20 = T20r;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                        default : T20 = T20n;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                   }
                   switch( zwei_j - 4 >= 0 ){
                        case JA : T40 = T40r;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43c;
                                  T44 = T44n;
                                  break;
                        default : T40 = T40n;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44n;
                                  break;
                   }
                   switch( zwei_j - 6 >= 0 ){
                        case JA : T60 = T60r;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63c;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66c;
                                  break;
                        default : T60 = T60n;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                   }
                   break;
 
       case  6  :  switch( zwei_j - 2 >= 0 ){
                        case JA : T20 = T20r;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                        default : T20 = T20n;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                   }
                   switch( zwei_j - 4 >= 0 ){
                        case JA : T40 = T40r;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43c;
                                  T44 = T44n;
                                  break;
                        default : T40 = T40n;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44n;
                                  break;
                   }
                   switch( zwei_j - 6 >= 0 ){
                        case JA : T60 = T60r;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63c;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66c;
                                  break;
                        default : T60 = T60n;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                   }
                   break;
 
       case  7  :  switch( zwei_j - 2 >= 0 ){
                        case JA : T20 = T20r;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                        default : T20 = T20n;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                   }
                   switch( zwei_j - 4 >= 0 ){
                        case JA : T40 = T40r;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44n;
                                  break;
                        default : T40 = T40n;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44n;
                                  break;
                   }
                   switch( zwei_j - 6 >= 0 ){
                        case JA : T60 = T60r;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66c;
                                  break;
                        default : T60 = T60n;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                   }
                   break;
 
       case  8  :  switch( zwei_j - 2 >= 0 ){
                        case JA : T20 = T20n;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                        default : T20 = T20n;
                                  T21 = T21n;
                                  T22 = T22n;
                                  break;
                   }
                   switch( zwei_j - 4 >= 0 ){
                        case JA : T40 = T40r;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44r;
                                  break;
                        default : T40 = T40n;
                                  T41 = T41n;
                                  T42 = T42n;
                                  T43 = T43n;
                                  T44 = T44n;
                                  break;
                   }
                   switch( zwei_j - 6 >= 0 ){
                        case JA : T60 = T60r;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64r;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                        default : T60 = T60n;
                                  T61 = T61n;
                                  T62 = T62n;
                                  T63 = T63n;
                                  T64 = T64n;
                                  T65 = T65n;
                                  T66 = T66n;
                                  break;
                   }
                   break;
 
    }
    return( tabelle );
}
/*------------------------------------------------------------------------------
                                drucke_par();
------------------------------------------------------------------------------*/
void drucke_par( fp,modus,dimj,tabelle,einheit_out,temp,ion,symmetrienr )
  FILE    *fp;
  CHAR    modus;
  INT     dimj;
  TABELLE *tabelle;
  CHAR    *einheit_out;
  DOUBLE  temp;
  CHAR    *ion;
  INT     symmetrienr;
{
    TABELLE *stern_setzen();
 
    fprintf(fp,TSS,einheit_out);
    fprintf(fp,T15,temp);
    fprintf(fp,T12,ion);
    fprintf(fp,T13,SYMLISTE[symmetrienr].symname,symmetrienr);
 
    switch(modus){
 
       case NOMAG :  fprintf(fp,T14,NICHTANGELEGT);
                     break;
       default    :  fprintf(fp,T14,ANGELEGT);
    }
    fprintf(fp,"%s",T11);
 
    tabelle = stern_setzen( tabelle,symmetrienr,dimj  );
 
    fprintf(fp,"%s",T11);
    fprintf(fp,"%s",T20);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T21);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T22);fprintf(fp,"%s",T11);
 
    fprintf(fp,"%s",T11);
    fprintf(fp,"%s",T40);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T41);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T42);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T43);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T44);fprintf(fp,"%s",T11);
 
    fprintf(fp,"%s",T11);
    fprintf(fp,"%s",T60);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T61);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T62);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T63);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T64);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T65);fprintf(fp,"%s",STR);
    fprintf(fp,"%s",T66);fprintf(fp,"%s",T11);
 
}
/*------------------------------------------------------------------------------
                                drucke_mag()
------------------------------------------------------------------------------*/
void drucke_mag( fp,modus ) /* Tabelle fuer Magnetfeld ausgeben */
  FILE *fp;
  CHAR modus;
{
    CHAR *s_modus,*z1,*z2,*z3;
    CHAR *rsm,*rsmt,*rsmtk,*rsmk,*rsmtu;
    CHAR *pom,*pomt,*pomtk,*pomk,*pomtu;
 
    #define SN "%s\n"
    rsm  = "====================================================";
    rsmt = "º External magnetic field  B (scaled with g_J)     º";
    rsmtk= "º %8s º %8s º %8s º comment   º\n";
    rsmk = "º            º            º            º           º";
    rsmtu= "----------------------------------------------------";
 
    pom  = "=======================================";
    pomt = "º External Magnetic field B           º";
    pomtk= "º %8s º %8s º comment º\n";
    pomk = "º            º            º           º";
    pomtu= "---------------------------------------";
 
    switch(modus){
         case 'r': s_modus = "KARTESISCHEN";
                   z1      = "    Bx    ";
                   z2      = "    By    ";
                   z3      = "    Bz    ";
                   break;
         case 's': s_modus = "SPHAERISCHEN";
                   z1      = "   |B|    ";
                   z2      = "   phi    ";
                   z3      = "  theta   ";
                   break;
         case 'p': s_modus = "POLAREN     ";
                   z1      = "   |B|    ";
                   z2      = "   phi    ";
                   z3      = "          ";
                   break;
         default : s_modus = "FEHLER!     ";
                   z1      = "          ";
                   z2      = "          ";
                   z3      = "          ";
    }
 
    fprintf(fp,"\n");
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"º The co-ordinates of the magnetic field are given in     º\n");
    fprintf(fp,"º %12s co-ordinates.                              º\n",s_modus);
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"º   B := Magnetic field     º in Tesla                    º\n");
    fprintf(fp,"º                           º                             º\n");
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"\n");
 
    if( modus=='r' || modus=='s' ){ /* Tabelle in rechtw. oder sphaer. Koord.*/
 
         fprintf(fp,SN,rsm);
         fprintf(fp,SN,rsmt);
         fprintf(fp,SN,rsm);
         fprintf(fp,rsmtk,z1,z2,z3);
         fprintf(fp,SN,rsm);
         fprintf(fp,SN,rsmk);
         fprintf(fp,SN,rsmtu);
         fprintf(fp,"\n");
 
    }
    else{ /* Tabelle in Polarkoordinaten */
 
         fprintf(fp,SN,pom);
         fprintf(fp,SN,pomt);
         fprintf(fp,SN,pom);
         fprintf(fp,pomtk,z1,z2);
         fprintf(fp,SN,pom);
         fprintf(fp,SN,pomk);
         fprintf(fp,SN,pomtu);
         fprintf(fp,"\n");
 
        }
 
    rsm  = "====================================================";
    rsmt = "º Molecular field B_mol (scaled with  2[g_J -1])   º";
    rsmtk= "º %8s º %8s º %8s º comment   º\n";
    rsmk = "º            º            º            º           º";
    rsmtu= "----------------------------------------------------";
 
    pom  = "=======================================";
    pomt = "º molecular field B_mol                 º";
    pomtk= "º %8s º %8s º comment º\n";
    pomk = "º            º            º           º";
    pomtu= "---------------------------------------";
    if( modus=='r' || modus=='s' ){ /* Tabelle in rechtw. oder sphaer. Koord.*/
 
         fprintf(fp,SN,rsm);
         fprintf(fp,SN,rsmt);
         fprintf(fp,SN,rsm);
         fprintf(fp,rsmtk,z1,z2,z3);
         fprintf(fp,SN,rsm);
         fprintf(fp,SN,rsmk);
         fprintf(fp,SN,rsmtu);
         fprintf(fp,"\n");
 
    }
    else{ /* Tabelle in Polarkoordinaten */
 
         fprintf(fp,SN,pom);
         fprintf(fp,SN,pomt);
         fprintf(fp,SN,pom);
         fprintf(fp,pomtk,z1,z2);
         fprintf(fp,SN,pom);
         fprintf(fp,SN,pomk);
         fprintf(fp,SN,pomtu);
         fprintf(fp,"\n");
 
        }
}
 
 
/*------------------------------------------------------------------------------
                                neben_create()
------------------------------------------------------------------------------*/
void neben_create( modus,name_chi2 ) /* Eingabefile Chi2 erzeugen  */
  CHAR modus,*name_chi2;
{
   READ *read_einheit(),*read;
   CHAR *einheit_out,*einheit_in,*name_par;
   FILE    *fopen(), *fp;
   CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10,*t11,*t12;
   CHAR *t13,*t14,*t15,*t16,*t17,*t18,*t19,*t20,*t21,*t22;
   CHAR *t23,*t24,*t25,*t26,*t27,*t28,*t29,*t30,*t31,*t32;
   CHAR *t33,*t34,*t35,*t36,*t37,*t38,*t39,*t40,*t41,*t42;
   CHAR *t43,*t44,*t45;
   CHAR *tl,*t[14],*td,*ionname,is_feld;
   INT  i,ionennr;
   DOUBLE zwei_j,temperatur,b1,b2,b3;
   DOUBLE b_norm,sqrt();
   BRUCH  *is_rational(),*z1,*z2,*z3;
   LONG hauptnenner,ggt_l(),r1,r2,r3;
 
   switch(modus){
     case AKQ  : name_par = AKQNAME; break;
     case BKQ  : name_par = BKQNAME; break;
     case DKQ  : name_par = DKQNAME; break;
     case LKQ  : name_par = LKQNAME; break;
     case VKQ  : name_par = VKQNAME; break;
     case WKQ  : name_par = WKQNAME; break;
     case XW   : name_par = XWNAME ; break;
     default   : name_par = "fehler      ";
   }
 
   read        = read_einheit(name_par,modus);
 
   einheit_out = read ->einheit_out;
   einheit_in  = read ->einheit_in ;
   zwei_j      = read ->zwei_j;
   ionennr     = read ->ionennr;
   temperatur  = read ->temperatur;
   is_feld     = read ->is_feld;
   b1          = B1(read);
   b2          = B2(read);
   b3          = B3(read);
 
   free_(read);
   ionname     = IONENIMP[ ionennr ].ionname;
 
   fp  = fopen_errchk(name_chi2,"w");
   write_title(fp);
 
   fprintf(fp,"\n\n");
 
t01 =" -------------------------------------------------------- \n";
t02 ="|          constraints for a (ortho)rhombic              |\n";
t03 ="|                 crsytal field fit                      |\n";
t04 ="|                     for   %4s                          |\n";
t05 =" -------------------------------------------------------- \n";
  fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);
  fprintf(fp,t04,ionname);
  fprintf(fp,"%s",t05);
  fprintf(fp,"\n\n");
 
t01 ="==========================================================\n";
t02 ="| Selection of the transition units                      |\n";
t03 ="==========================================================\n";
t04 ="| If no mark is given in the lower field, or             |\n";
t05 ="| is put to 0, then the transition intensities are       |\n";
t06 ="| given in barn to the temperature listed below          |\n";
t07 ="| otherwise,                                             |\n";
t08 ="| the transition elements themselves are the intensities |\n";
t09 ="| (note of the sum rules!).                              |\n";
t10 ="|                                                        |\n";
t11 ="|-------------------------===============================|\n";
t12 ="| Please input a number: 0                               |\n";
t13 ="|-------------------------===============================|\n";
t14 ="| Default                 : Intensities in barns         |\n";
t15 ="|                         : Temperature : %8.2f Kelvin  |\n";
t16 ="----------------------------------------------------------\n";
  fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
  fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
  fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
  fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,t15,temperatur);
  fprintf(fp,"%s",t16);
  fprintf(fp,"\n\n");
 
 
t01 ="==============================    ===================================\n";
t02 ="| Energy values   Ei and     |    | transition     U(i->k) and       |\n";
t03 ="| its error  DEi             |    | its error     DU(i->k)           |\n";
t04 ="==============================    ===================================\n";
t05 ="|  i | Ei/%6s |DEi/%6s |    | i -> k | U(i->k)   | DU(i->k)   |\n";
t06 ="|----|-----------|-----------|    |--------|-----------|------------|\n";
t07 ="| %2d |           |           |    |   ->   |           |            |\n";
t09 ="| %2d |     0     |     0     |    |   ->   |           |            |\n";
t08 ="------------------------------    -----------------------------------\n";
 
  fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
  fprintf(fp,t05,einheit_out,einheit_out);
  fprintf(fp,"%s",t06);
  for(i=1;i<=ANZ_NIVEAUS;++i){
     if( i==1)  fprintf(fp,t09,i );
     else       fprintf(fp,t07,i );
  }
  fprintf(fp,"%s",t08);
  fprintf(fp,"\n\n");
 
t01 ="==========================================================\n";
t02 ="| Approximate positions and intensities?                 |\n";
t03 ="==========================================================\n";
t04 ="| In this case the temperature dependant characteristics |\n";
t05 ="| of CF transitions were determined out of the neutron   |\n";
t06 ="| spectra. with the line position and the transition     |\n";
t07 ="| intensities.                                           |\n";
t08 ="| If another number other than 0 is given below then     |\n";
t09 ="| the intensities and positions are read from a file     |\n";
t10 ="| <filename>, other wise not.                            |\n";
t11 ="|                                                        |\n";
t12 ="|----------------------------============================|\n";
t13 ="| Please give a number   : 0                             |\n";
t14 ="|----------------------------============================|\n";
t15 ="| Default                    : no positions or          |\n";
t16 ="|                            : Intensities given.        |\n";
t17 ="|----------------------------============================|\n";
t18 ="| Filename  (variable)       : Posint.Chi2               |\n";
t19 ="|----------------------------============================|\n";
t20 ="|            |1 Header      : I5I5F12F12F12Text         |\n";
t21 ="|I5:         |n Data lines  ...                         |\n";
t22 ="|Record      |-------------------------------------------|\n";
t23 ="|            |                                           |\n";
t24 ="|            |If the data record numbers are positive    |\n";
t25 ="|I5:         |then first all the x-values are read and   |\n";
t26 ="|number of   |then all the y-values.                     |\n";
t27 ="|Records     |                                           |\n";
t28 ="|            |                 =====      =====          |\n";
t29 ="|F12:        |per data line :  |  7| times| 11|          |\n";
t30 ="|Temperature |                 =====      =====          |\n";
t31 ="|            |-------------------------------------------|\n";
t32 ="|F12:        |                                           |\n";
t33 ="|scaling     |if the data record number is negative      |\n";
t34 ="|factor (SF) |then the x and y-values                    |\n";
t35 ="|inelastic   |are next to eachother.                     |\n";
t36 ="|scattering  |                                           |\n";
t37 ="|intensity   |                 =====      =====          |\n";
t38 ="|            |per dataline   : |  2|times | 12|          |\n";
t39 ="|F12:        |                 =====      =====          |\n";
t40 ="|SF for  QE  |                                           |\n";
t41 ="|------------------====---====---------------------------|\n";
t42 ="| Record number : 00 until 00 which the fit uses.        |\n";
t43 ="----------------------------------------------------------\n";
  fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
  fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
  fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
  fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,"%s",t15);fprintf(fp,"%s",t16);
  fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);fprintf(fp,"%s",t19);fprintf(fp,"%s",t20);
  fprintf(fp,"%s",t21);fprintf(fp,"%s",t22);fprintf(fp,"%s",t23);fprintf(fp,"%s",t24);
  fprintf(fp,"%s",t25);fprintf(fp,"%s",t26);fprintf(fp,"%s",t27);fprintf(fp,"%s",t28);
  fprintf(fp,"%s",t29);fprintf(fp,"%s",t30);fprintf(fp,"%s",t31);fprintf(fp,"%s",t32);
  fprintf(fp,"%s",t33);fprintf(fp,"%s",t34);fprintf(fp,"%s",t35);fprintf(fp,"%s",t36);
  fprintf(fp,"%s",t37);fprintf(fp,"%s",t38);fprintf(fp,"%s",t39);fprintf(fp,"%s",t40);
  fprintf(fp,"%s",t41);fprintf(fp,"%s",t42);fprintf(fp,"%s",t43);
  fprintf(fp,"\n\n");
 
 
 
t01 ="==========================================================\n";
t02 ="| Inverse Susceptibility approximation  ?                |\n";
t03 ="==========================================================\n";
t04 ="| If there is a number other than zero given below       |\n";
t05 ="| then the inverse susceptability curve is read from a   |\n";
t06 ="| file <filename> other wise not                         |\n";
t07 ="|                                                        |\n";
t08 ="|----------------------------============================|\n";
t09 ="| Please give a number       : 0                         |\n";
t10 ="|----------------------------============================|\n";
t11 ="| Default                    : no inverse Suscepti-      |\n";
t12 ="|                            : bility given              |\n";
t13 ="|----------------------------============================|\n";
t14 ="| Filename  (variable)       : Suszept.Chi2              |\n";
t15 ="|----------------------------============================|\n";
t16 ="|            |1 Header       : I5I5Text                  |\n";
t17 ="|I5:         |n Records      ...                         |\n";
t18 ="|Record-     |-------------------------------------------|\n";
t19 ="|number      |                                           |\n";
t20 ="|            |If the data record numbers are positive    |\n";
t21 ="|I5:         |then first all the x-values are read and   |\n";
t22 ="|number of   |then all the y-values.                     |\n";
t23 ="|Records     |                                           |\n";
t24 ="|            |                 =====      =====          |\n";
t25 ="|Text:       |per dataline    :|  7| times| 11|          |\n";
t26 ="|comment     |                 =====      =====          |\n";
t27 ="|            |-------------------------------------------|\n";
t28 ="|            |                                           |\n";
t29 ="|            |if the data record number is negative      |\n";
t30 ="|            |then the x and y-values                    |\n";
t31 ="|            |are next to eachother.                     |\n";
t32 ="|            |                                           |\n";
t33 ="|            |                 =====      =====          |\n";
t34 ="|            |per dataline  :  |  2| times| 12|          |\n";
t35 ="|            |                 =====      =====          |\n";
t36 ="|            |                                           |\n";
t37 ="|---------------------------------------------------=====|\n";
t38 ="| Recordnr.    the inverse Suscept. in [100]    : 00     |\n";
t39 ="| Recordnr.    the inverse Suscept. in [010]    : 00     |\n";
t40 ="| Recordnr.    the inverse Suscept. in [001]    : 00     |\n";
t41 ="| Recordnr   .    der inverse Suscept. Polycrystal: 00   |\n";
t42 ="| Molecular field constant Lambda (in  mol/emu ): 0000.00|\n";
t43 ="| Filename (var.) Theta(T)   : Theta.chi2                |\n";
t44 ="------------------------------===========================-\n";
  fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
  fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
  fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
  fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,"%s",t15);fprintf(fp,"%s",t16);
  fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);fprintf(fp,"%s",t19);fprintf(fp,"%s",t20);
  fprintf(fp,"%s",t21);fprintf(fp,"%s",t22);fprintf(fp,"%s",t23);fprintf(fp,"%s",t24);
  fprintf(fp,"%s",t25);fprintf(fp,"%s",t26);fprintf(fp,"%s",t27);fprintf(fp,"%s",t28);
  fprintf(fp,"%s",t29);fprintf(fp,"%s",t30);fprintf(fp,"%s",t31);fprintf(fp,"%s",t32);
  fprintf(fp,"%s",t33);fprintf(fp,"%s",t34);fprintf(fp,"%s",t35);fprintf(fp,"%s",t36);
  fprintf(fp,"%s",t37);fprintf(fp,"%s",t38);fprintf(fp,"%s",t39);fprintf(fp,"%s",t40);
  fprintf(fp,"%s",t41);fprintf(fp,"%s",t42);fprintf(fp,"%s",t43);fprintf(fp,"%s",t44);
  fprintf(fp,"\n\n");
 
/* Magnetisierungskurven fitten?*/
  b_norm = sqrt( b1*b1 + b2*b2 + b3*b3 );
  if( !null(b_norm,MACHEPSFACT*accuracy())  ){
        b1  /=  b_norm;
        b2  /=  b_norm;
        b3  /=  b_norm;
  }
  else{
        b1  =  0.0;
        b2  =  0.0;
        b3  =  1.0;
      }
 
 /* Richtung von B ausrechnen  */
    z1 = is_rational(b1*b1);
    hauptnenner = ggt_l(ZAEHLER(z1),NENNER(z1));
    ZAEHLER(z1) /= hauptnenner;
    NENNER( z1) /= hauptnenner;
 
    z2 = is_rational(b2*b2);
    hauptnenner = ggt_l(ZAEHLER(z2),NENNER(z2));
    ZAEHLER(z2) /= hauptnenner;
    NENNER( z2) /= hauptnenner;
 
    z3 = is_rational(b3*b3);
    hauptnenner = ggt_l(ZAEHLER(z3),NENNER(z3));
    ZAEHLER(z3) /= hauptnenner;
    NENNER( z3) /= hauptnenner;
 
    hauptnenner = 1;
    if( NENNER(z1) > hauptnenner )  hauptnenner = NENNER(z1);
    if( NENNER(z2) > hauptnenner )  hauptnenner = NENNER(z2);
    if( NENNER(z3) > hauptnenner )  hauptnenner = NENNER(z3);
 
    ZAEHLER(z1) *= (hauptnenner/NENNER(z1));
    ZAEHLER(z2) *= (hauptnenner/NENNER(z2));
    ZAEHLER(z2) *= (hauptnenner/NENNER(z2));
 
    r1 = DSIGN(b1)*(LONG)sqrt( (DOUBLE)ZAEHLER(z1) );
    r2 = DSIGN(b2)*(LONG)sqrt( (DOUBLE)ZAEHLER(z2) );
    r3 = DSIGN(b3)*(LONG)sqrt( (DOUBLE)ZAEHLER(z3) );
    free_(z1);free_(z2);free_(z3);
 
t01 ="==========================================================\n";
t02 ="| Magnetisation curve approximation ?                    |\n";
t03 ="==========================================================\n";
t04 ="| If there is a number other than zero given below       |\n";
t05 ="| then the magnetisation curve is read from a file       |\n";
t06 ="| <filename> other wise not                              |\n";
t07 ="|                                                        |\n";
t08 ="|----------------------------============================|\n";
t09 ="| Please enter a number      : 0                         |\n";
t10 ="|----------------------------============================|\n";
t11 ="| Default                    : no magnetisation curve    |\n";
t12 ="|                            : given                     |\n";
t13 ="|----------------------------============================|\n";
t14 ="| Filename  (variable)       : Moment.Chi2               |\n";
t15 ="|----------------------------============================|\n";
t16 ="|            |1 header       : I5I5Text                  |\n";
t17 ="|I5:         |n records      ...                         |\n";
t18 ="|record-     |-------------------------------------------|\n";
t19 ="|number      |                                           |\n";
t20 ="|            |If the data record numbers are positive    |\n";
t21 ="|I5:         |then first all the x-values are read and   |\n";
t22 ="|number of   |then all the y-values.                     |\n";
t23 ="|data lines  |                                           |\n";
t24 ="|            |                 =====      =====          |\n";
t25 ="|Text:       |per dataline  :  |  7|times | 11|          |\n";
t26 ="|comment     |                 =====      =====          |\n";
t27 ="|            |-------------------------------------------|\n";
t28 ="|            |                                           |\n";
t29 ="|            |if the data record number is negative      |\n";
t30 ="|            |then the x and y-values		       |\n";
t31 ="|            |are next to eachother.                     |\n";
t32 ="|            |                                           |\n";
t33 ="|            |                 =====      =====          |\n";
t34 ="|            |per dataline :   |  2| times | 12|          |\n";
t35 ="|            |                 =====      =====          |\n";
t36 ="|            |                                           |\n";
t37 ="|--------------------------------------------------------|\n";
t38 ="|  B_ex := external Field                                |\n";
t39 ="|                                                        |\n";
t40 ="|                                   ======    ========== |\n";
t41 ="| Datnr of m for B_ex in     [ 1 0 0] : 00 |Temp: %7.2f K|\n";
t42 ="| Datnr of m for B_ex in     [ 0 1 0] : 00 |Temp: %7.2f K|\n";
t43 ="| Datnr of m for B_ex in     [ 0 0 1] : 00 |Temp: %7.2f K|\n";
t44 ="| Datnr of m for B_ex in  [%2d%2d%2d] : 00 |Temp: %7.2f K|\n";
t45 ="------------------------------------======----==========--\n";
  fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
  fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
  fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
  fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,"%s",t15);fprintf(fp,"%s",t16);
  fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);fprintf(fp,"%s",t19);fprintf(fp,"%s",t20);
  fprintf(fp,"%s",t21);fprintf(fp,"%s",t22);fprintf(fp,"%s",t23);fprintf(fp,"%s",t24);
  fprintf(fp,"%s",t25);fprintf(fp,"%s",t26);fprintf(fp,"%s",t27);fprintf(fp,"%s",t28);
  fprintf(fp,"%s",t29);fprintf(fp,"%s",t30);fprintf(fp,"%s",t31);fprintf(fp,"%s",t32);
  fprintf(fp,"%s",t33);fprintf(fp,"%s",t34);fprintf(fp,"%s",t35);fprintf(fp,"%s",t36);
  fprintf(fp,"%s",t37);fprintf(fp,"%s",t38);fprintf(fp,"%s",t39);fprintf(fp,"%s",t40);
  fprintf(fp,t41,temperatur);
  fprintf(fp,t42,temperatur);
  fprintf(fp,t43,temperatur);
  fprintf(fp,t44,r1,r2,r3,temperatur);
  fprintf(fp,"%s",t45);
  fprintf(fp,"\n\n");
 
 
 
t01 ="==========================================================\n";
t02 ="|Weigting of the different conditions	               |\n";
t03 ="==========================================================\n";
t04 ="| Energy Eigenvalue                 : 1                  |\n";
t05 ="| Position + Intensity              : 1                  |\n";
t06 ="| Intensity wrt. Matrix elements    : 1                  |\n";
t07 ="| Susceptibility curve              : 1                  |\n";
t08 ="| Magnetisation curve               : 1                  |\n";
t09 ="----------------------------------------------------------\n";
  fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
  fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
  fprintf(fp,"%s",t09);
  fprintf(fp,"\n\n");
 
 
t01 ="==========================================================\n";
t02 ="| number of Fitroutine                                   |\n";
t03 ="==========================================================\n";
t04 ="| 1 | DOWNHILL Simplex Methode in N - Dimensions         |\n";
t05 ="|   |  (from Numerical Recipes, Seite 292-293)           |\n";
t06 ="| 2 | POWELLs  Method         in N - Dimensions          |\n";
t07 ="|   |  (from Numerical Recipes, Seite 299-301)           |\n";
t08 ="| 3 | GRADIENTen Method       in N - Dimensions          |\n";
t09 ="|   |  (from Numerical Recipes, Seite 301-306)           |\n";
t10 ="| 4 | VA05A                                              |\n";
t11 ="|   |  (Harwell Subroutine Library,in C) 		       |\n";
t12 ="| 5 | DUMMY                                              |\n";
t13 ="|   |  (not implemented    )                             |\n";
t14 ="|   |                                                    |\n";
t15 ="|   | Numerical Recipes, Cambridge University Press      |\n";
t16 ="|   |  (ISBN 0-521-30811-9)                              |\n";
t17 ="|-------------------------------=========================|\n";
t18 ="| please give  1,2,3,4,5        : 1                      |\n";
t19 ="|-------------------------------=========================|\n";
t20 ="| Default                       : Fitroutine DOWNHILL    |\n";
t21 ="|-------------------------------=========================|\n";
t22 ="| Stop in the interation step   : %12d           |\n";
t23 ="|-------------------------------=========================|\n";
t24 ="| Default                       : DOWNHILL  %12d |\n";
t25 ="|                               : POWELL    %12d |\n";
t26 ="|                               : GRADIENT  %12d |\n";
t27 ="|                               : VA05A     %12d |\n";
t28 ="|                               : DUMMY     %12d |\n";
t29 ="----------------------------------------------------------\n";
  fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
  fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
  fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
  fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,"%s",t15);fprintf(fp,"%s",t16);
  fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);fprintf(fp,"%s",t19);fprintf(fp,"%s",t20);
  fprintf(fp,"%s",t21);fprintf(fp,t22,MAXDOWNHILL);
  fprintf(fp,"%s",t23);
  fprintf(fp,t24,MAXDOWNHILL);
  fprintf(fp,t25,MAXPOWELL  );
  fprintf(fp,t26,MAXGRADIENT);
  fprintf(fp,t27,MAXVA05A   );
  fprintf(fp,t28,MAXDUMMY   );
  fprintf(fp,"%s",t29);
  fprintf(fp,"\n\n");
 
 
t01  ="===========================================================\n";
t02  ="| Operator of the   | should the parameter to the Operator|\n";
t03  ="|Hamiltonianoperator| be fixed?   			 |\n";
t04  ="|                   | (Give a number)		         |\n";
t05  ="|                   | ( 0 = no otherwise yes)             |\n";
t06  ="===========================================================\n";
t[1] ="| O40+ 5*O44        |                                     |\n";
t[2] ="| O60-21*O64        |                                     |\n";
t[3] ="| O40- 5*O44        |                                     |\n";
t[4] ="| O60+21*O64        |                                     |\n";
t[5] ="| O20               |                                     |\n";
t[6] ="| O22               |                                     |\n";
t[7] ="| O42               |                                     |\n";
t[8] ="| O62               |                                     |\n";
t[9] ="| O66               |                                     |\n";
t[10]="| J : B_mol         | 1                                   |\n";
t[11]="| J : phi           | 1                                   |\n";
t[12]="| J : theta         | 1                                   |\n";
t07 ="-----------------------------------------------------------\n";
tl  ="|                   |                                     |\n";
 
  fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
  fprintf(fp,t05,einheit_in);
  fprintf(fp,"%s",t06);
  for( i=1; i<= 9; ++i){
       switch(i){
          case 1:
          case 3:
          case 7:  if( zwei_j < 4.0) td = tl;
                   else              td = t[i];
                   break;
          case 2:
          case 4:
          case 8:
          case 9:  if( zwei_j < 6.0) td = tl;
                   else              td = t[i];
                   break;
          case 5:
          case 6:
          default: if( zwei_j < 2.0) td = tl;
                   else              td = t[i];
                   break;
       }
       fprintf(fp,"%s",td);
  }
  if(is_feld=='a')
       for(i=10; i<=12; ++i)
           fprintf(fp,"%s",t[i]);
  else for(i=10; i<=12; ++i)
           fprintf(fp,"%s",tl);
  fprintf(fp,"%s",t07);
  fprintf(fp,"\n\n");
 
 
t01 ="==========================================================\n";
t02 ="| Synopitcal table during Fits given ?    	       |\n";
t03 ="==========================================================\n";
t04 ="| If there is a zero or nothing input below              |\n";
t05 ="| then there is no table produced                        |\n";
t06 ="| 					               |\n";
t07 ="|-------------------------===============================|\n";
t08 ="| Please enter a number : 1                              |\n";
t09 ="|-------------------------===============================|\n";
t10 ="| Default                 : table                        |\n";
t11 ="|                         : given                        |\n";
t12 ="|-------------------------===============================|\n";
t13 ="| Modulo                  : 5                            |\n";
t14 ="|-------------------------===============================|\n";
t15 ="| <Modulo> determines in which Iteration step             |\n";
t16 ="| the synoptical table is produced. This follows from    |\n";
t17 ="| :				                       |\n";
t18 ="| Iteration step modulo  <Modulo> is zero .              |\n";
t19 ="| example: <modulo>=5 ==> displays the table             |\n";
t20 ="| after every 5th Iteration step.                        |\n";
t21 ="----------------------------------------------------------\n";
  fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
  fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
  fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
  fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,"%s",t15);fprintf(fp,"%s",t16);
  fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);fprintf(fp,"%s",t19);
  fprintf(fp,"%s",t20);fprintf(fp,"%s",t21);
  fprintf(fp,"\n\n");
 
 
  fclose(fp);
  printf("file %s produced.\n",name_chi2);
 
}
/*------------------------------------------------------------------------------
                          read_einheit()
------------------------------------------------------------------------------*/
READ *read_einheit(name,art)
    CHAR *name;
    CHAR art;
{
    READ      *read;
    FILE      *fp,*fopen();
    INT   /*  anz_nn,*/dimj,zwei_j,ionennr/*,symmetrienr*/;
    INT       buffer_size=381,i,einheitnr_in,einheitnr_out;
/*  DOUBLE    versionsnummer; */
    DOUBLE    x1,x2,x3/*,myB*/;
    DOUBLE    h,theta,phi;
    DOUBLE    sin(),cos();
    DOUBLE    a_tof(),sqrt(),temperatur;
    CHAR      *einheit_in,*einheit_out,modus;
    CHAR      *ion;
    CHAR      c,*string,*line,*fgets(),*a_tos();
 
    printf("Reading file %s ....\n",name);
    string   = STRING_ALLOC(buffer_size);
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 1.==== */
/*  versionsnummer =a_tof( line , 11,23); */
    fclose(fp);
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    line=fgets( string , buffer_size , fp ); /* leerzeile                  */
    line=fgets( string , buffer_size , fp ); /* Kristallfeldpara....       */
    c = VALUE(line,33 );
    einheitnr_in = is_einheit_imp(c);
    if( einheitnr_in == NICHTIMP )
         read_error(21,fp,name);
    einheit_in = EINHEITIMP[ einheitnr_in ].einheit;
/*  myB        = EINHEITIMP[ einheitnr_in ].myB; */
 
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    line=fgets( string , buffer_size , fp ); /* : Energieeigenwerte......*/
    c = VALUE(line,31);
    einheitnr_out= is_einheit_imp(c);
    if( einheitnr_out== NICHTIMP )
         read_error(21,fp,name);
    einheit_out= EINHEITIMP[ einheitnr_out].einheit;
 
    line=fgets( string , buffer_size , fp ); /* : Temperatur der Probe :.*/
    temperatur = a_tof(line,31,61);
 
    line=fgets( string , buffer_size , fp ); /* : Betrachtetes Aufion....*/
    ion =a_tos( line , 31,56);
 
    line=fgets( string , buffer_size , fp ); /* : Symmetrie ... Symmnr...*/
/*  symmetrienr =a_toi( line , 53,61); */
 
 
/*  anz_nn    = 0; */
    ionennr   = isimplementiert(ion);
    dimj      = IONENIMP[ ionennr ].dimj;
    if(strncmp(ion,"S=",2)==0)  /* S=... ion !! extract dimj from string ion */
     {dimj =(int)( 2 * strtod (ion+2, NULL)+1);
      IONENIMP[ ionennr ].dimj=dimj;
     }    

    zwei_j    = dimj - 1 ;
 
 
 
    line=fgets( string , buffer_size , fp ); /* Magnetfeld: . Symmnr...  */
    c=VALUE(line,31);
    switch(c){
        case 'a':  modus = 'a';
                   break;
        case 'n':  modus = 'n';
                   break;
        default : read_error(24,fp,name); modus=0;
    }
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 3.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 4.==== */
 
    if( art != XW ){
        while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 5.==== */
        while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 6.==== */
 
        while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 7.==== */
        while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 8.==== */
    }
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*9.==== */
 
 
    if( modus=='n' ) {      /* kein Magnetfeld angelegt */
        read = READ_ALLOC(1);
        B1(   read ) = 0.0;
        B2(   read ) = 0.0;
        B3(   read ) = 1.0;
        read -> einheit_in = einheit_in;
        read -> einheit_out= einheit_out;
        read -> zwei_j     = zwei_j;
        read -> ionennr    = ionennr;
        read -> temperatur = temperatur;
        read -> is_feld    = modus;
        fclose(fp);
        return( read );
 
    }
    /* falls ein Magnetfeld angelegt wurde : */
    read = READ_ALLOC(1);
 
    /* externes Magnetfeld lesen  */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*1.==== */
    line=fgets( string , buffer_size , fp ); /* : die Koordinaten der ...*/
    line=fgets( string , buffer_size , fp );/* : º P... v º R... v º S...*/
    c   = VALUE(line,2);
    switch(c){
         case 'K' : MODUS(read) = 'r';
                    break;
         case 'S' : MODUS(read) = 's';
                    break;
         case 'P' : MODUS(read) = 'p';
                    break;
         default  : read_error(3,fp,name);
    }
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/*2.==== */
 
    for( i=1 ; i<= 3 ; ++i)/* Kopf der Magnetfeldtabelle ueberlesen */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );
    line=fgets( string , buffer_size , fp );   /* : º  h1 º h2 º h3 º*/
    B1(read) = x1 = a_tof(line, 2,13);
    B2(read) = x2 = a_tof(line,15,26);
    x3=0.0;
    if(MODUS(read) !='p')  B3(read) = x3 = a_tof(line,28,39);
 
    isinlimits(fp,name ,0, x1,x2,x3,MODUS(read) );
 
 
    switch( MODUS(read) ){ /* transformiere B -> (Bx,By,Bz) */
 
       case 'p' :
       case 's' :  h     = B1(read);
                   phi   = B2(read)*pi/180.0;
                   theta = B3(read)*pi/180.0;
 
                   if( MODUS(read) == 'p' )  theta = pi/2;
 
                   B1(read) = h*sin(theta)*cos(phi);
                   B2(read) = h*sin(theta)*sin(phi);
                   B3(read) = h*cos(theta);
                   break;
 
       case 'r' :  break;
    }
 
 
    read -> einheit_in = einheit_in;
    read -> einheit_out= einheit_out;
    read -> zwei_j     = zwei_j;
    read -> ionennr    = ionennr;
    read -> temperatur = temperatur;
    read -> is_feld    = modus;
    fclose(fp);
    return( read );
}
 
/*------------------------------------------------------------------------------
                          is_sbekannt()
------------------------------------------------------------------------------*/
INT is_sbekannt(datnr,iteration)
    INT       datnr;
    ITERATION *iteration;
{
    INT d,a,b,c,p;
 
    d = ABS(datnr);
    a = ABS(NRSUSA(iteration));
    b = ABS(NRSUSB(iteration));
    c = ABS(NRSUSC(iteration));
    p = ABS(NRSUSP(iteration));
 
    if( (d==a||d==b||d==c||d==p) && d!=0 )  return(JA);
    else                                    return(NEIN);
}
/*------------------------------------------------------------------------------
                          is_mbekannt()
------------------------------------------------------------------------------*/
INT is_mbekannt(datnr,iteration)
    INT       datnr;
    ITERATION *iteration;
{
    INT d,a,b,c,p;
 
    d = ABS(datnr);
    a = ABS(NRMAGA(iteration));
    b = ABS(NRMAGB(iteration));
    c = ABS(NRMAGC(iteration));
    p = ABS(NRMAGP(iteration));
 
    if( (d==a||d==b||d==c||d==p) && d!=0 )  return(JA);
    else                                    return(NEIN);
 
}
/*------------------------------------------------------------------------------
                          is_pbekannt()
------------------------------------------------------------------------------*/
INT is_pbekannt(datnr,iteration)
    INT       datnr;
    ITERATION *iteration;
{
    INT d,a,e;
 
    d = ABS(datnr);
    a = ABS(NRPOSA(iteration));
    e = ABS(NRPOSE(iteration));
 
    if( (d>=a && d<=e) && d!=0 )  return(JA);
    else                          return(NEIN);
 
}
/*------------------------------------------------------------------------------
                          lesedaten()
------------------------------------------------------------------------------*/
void lesedaten(fp,x,f,n,name,nummerierung,ip,imax,bekannt,pdatnr,
              iteration,anzdatnr)
    FILE      *fp;
    DOUBLE    *x,*f;
    INT       *n,(*bekannt)(),*pdatnr;
    INT       nummerierung;
    INT       *ip,imax;
    CHAR      *name;
    ITERATION *iteration;
    INT       *anzdatnr;
{
    INT      buffer_size=381,i,k,j,z,diff;
    INT      anz_dat,a_toi(),datnr;
    INT      fpa,fpzl,fna,fnzl,loopmax;
    CHAR     *string,*line,*fgets();
    FILE     *fopen();
 
    string   = STRING_ALLOC(buffer_size);
 
    fpa  = FORMATPA(iteration);   /* 7 mal 11 */
    fpzl = FORMATPB(iteration);
    fna  = FORMATNA(iteration);   /* 2 mal 12 */
    fnzl = FORMATNB(iteration);
 
    i     = *ip;
    datnr = *pdatnr;
    diff  = imax-i;
 
     printf("lese file %s ....\n",name);
     if( (fp=fopen(name,"rb"))==(FILE*)0 )
          read_error(2,fp,name);
 
     *anzdatnr=0;
 
 
     loopmax = LOOPMAX+1;
     loop:
        --loopmax;
       if(loopmax==0){
             printf("Fehler : kann file %s nicht lesen!\n",name);
             printf("Ursache: Datenformat falsch oder Datensatznummer nicht vorhanden.\n");
             fclose(fp);
             exit(1);
        }
        line=fgets( string , buffer_size , fp ); /*Kopfzeile*/
        datnr    = a_toi(line,0,4);
        anz_dat  = a_toi(line,5,9);
        if(nummerierung!=JA){
          *anzdatnr += datnr;
          datnr      = *anzdatnr;
        }
        if(!(*bekannt)(datnr,iteration) ){/* ueberlesen*/
           if( datnr>0 ) anz_dat = 2*(anz_dat/fpa+1);
           for(k=1; k<= anz_dat; ++k)
              line=fgets( string , buffer_size , fp ); /*datenzeile*/
           goto loop;
        }
        if(diff!=0){
           --diff;
           if( datnr>0 ) anz_dat = 2*(anz_dat/fpa+1);
           for(k=1; k<= anz_dat; ++k)
              line=fgets( string , buffer_size , fp ); /*datenzeile*/
           goto loop;
        }
 
 
 
        if( anz_dat <= 0 ) goto loop;
        --i;
 
        if( *n < anz_dat ) anz_dat = *n;
 
        if(datnr>0){
          j=1;
          for(k=0; k<= anz_dat/fpa; ++k){
              line=fgets( string , buffer_size , fp ); /*datenzeile*/
              for(z=0; z<=fpa-1 && j<=anz_dat; ++z)
                  VALUE(x,j++)=a_tof(line,z*fpzl,fpzl-1+z*fpzl);
          }
          j=1;
          for(k=0; k<= anz_dat/fpa; ++k){
              line=fgets( string , buffer_size , fp ); /*datenzeile*/
              for(z=0; z<=fpa-1 && j<=anz_dat; ++z)
                  VALUE(f,j++)=a_tof(line,z*fpzl,fpzl-1+z*fpzl);
          }
       }
       else if(datnr<0)
          for(k=1; k<= anz_dat; ++k){
              line=fgets( string , buffer_size , fp ); /*Datenzeile*/
              for(z=0; z<=fna-1; z+=2){
               VALUE(x,k )=a_tof(line,z*fnzl,fnzl-1+z*fnzl);
               VALUE(f,k )=a_tof(line,(z+1)*fnzl,fnzl-1+(z+1)*fnzl);
              }
          }
 
 
       *n      = anz_dat;
       *ip     = i;
       *pdatnr = datnr;
 
}
/*------------------------------------------------------------------------------
                          neben_read()
------------------------------------------------------------------------------*/
 
NEBENBEDINGUNG *neben_read(setup,name,kristallfeld,einheitnr_out,macheps)
    SETUP     *setup;
    CHAR      *name;
    KRISTALLFELD *kristallfeld;
    INT       einheitnr_out;
    DOUBLE    macheps;
{
    ITERATION *iteration;
    NEBENBEDINGUNG *neben;
    MATRIX    *intensit,*d_intensit,*mx_alloc();
    VEKTOR    *ew      ,*d_ew      ,*vr_alloc(), *fix;
    FILE      *fp,*fopen(),*fp_sus=0,*fp_theta;
    INT       buffer_size=381,i,k,einheitnr_in,a_toi(),zeile,anzahl=0,j;
    INT       anz_lines,anz_par=0,anz_var=0,spalte,datnr,anz_dat;
    INT       fpzl, fnzl, fpa, fna, z, flag,imax,ipos;
    INT       anzdatnr=0,loopmax,anz;
    DOUBLE    versionsnummer,*datent,*datensus;
    DOUBLE    a_tof(),sqrt(),temperatur;
    DOUBLE    wew,wint,wmat,wsus,wmag,wsum,wpos;
    DOUBLE    pos_t,pos_icin,pos_icqe,*datene,*dateni,*xx,*ff;
    CHAR      c,*string,*line,*fgets(),*susname,*magname,*posname;
    CHAR      *thetaname;
    CHAR      *file_name();
    INT       *sort(),*nummer;
    DOUBLE    *werte,*werti;
    INT       is_sbekannt(), is_pbekannt(), is_mbekannt();
    void       lesedaten();
 
    iteration  = ITERATION(kristallfeld);
    temperatur = TEMPERATUR(iteration);
 
    printf("lese file %s ....\n",name);
    string   = STRING_ALLOC(buffer_size);
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
 
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 1.|  = */
    c = VALUE(line,2);
    versionsnummer =a_tof(line,11,23);
    if( c != 'V' || versionsnummer < VERSION  ){
         printf("file %s has a OLD Version number (%6.2f).\n",name,versionsnummer);
         printf("file %s wird daher NEU erzeugt.\n",name);
         neben_create(EINGABEPARAMETERART(kristallfeld),name);
         exit(1);
    }
    fclose(fp);
 
    if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
/**********************/
/* Art des Uebergangs */
/**********************/
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    line=fgets( string , buffer_size , fp ); /* | i | Ei/Einheit ...       */
    c = VALUE(line,10 );
    einheitnr_in = is_einheit_imp(c);
    if( einheitnr_in == NICHTIMP )
         read_error(21,fp,name);
 
    line=fgets( string , buffer_size , fp ); /* |----|--------------       */
    anz_lines = 0;
    while(JA){
       line=fgets( string , buffer_size , fp ); /* |  1 |----------        */
       i = a_toi(line,1,4);
       if( i!=0 ) ++anz_lines;
       else       break;
    }
    fclose(fp);
    neben      =  NEBEN_ALLOC(1);
      ew       =  vr_alloc(anz_lines);
    d_ew       =  vr_alloc(anz_lines);
      intensit =  mx_alloc(anz_lines,anz_lines);
    d_intensit =  mx_alloc(anz_lines,anz_lines);
 
 
 
 
 if( (fp=fopen(name,"rb"))==(FILE*)0 )  read_error(2,fp,name);
/**********************/
/* Art des Uebergangs */
/**********************/
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    for(i=1; i<=9; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    i = a_toi(line,27,56);
    IS_INTENSITAET(  iteration)=NEIN;
    IS_MATRIXELEMENT(iteration)=NEIN;
    IS_EIGENWERT(    iteration)=NEIN;
    if( i==0      ) IS_INTENSITAET(  iteration) = JA;
    else            IS_MATRIXELEMENT(iteration) = JA;
 
    for(i=1; i<=3; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    INTTEMP(iteration) = a_tof(line,40,50);
    if(INTTEMP(iteration)<=0.0){
        printf("Error : in file %s !\n",name);
        printf("Cause: The temperature to the given intensities");
        printf(" is zero!\n");
        fclose(fp);
        exit(1);
    }
 
 while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.===   = */
 while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.===   = */
 line=fgets( string , buffer_size , fp ); /* | i | Ei/Einheit ...          */
 line=fgets( string , buffer_size , fp ); /* |-------------------          */
 
 for(zeile=1; zeile<=anz_lines; ++zeile){
    line=fgets( string , buffer_size , fp ); /* | 1  |--------       ....*/
 
    i = a_toi(line,1,4);
    RV(  ew,i) = a_tof( line , 6,16)*EINHEITIMP[einheitnr_in].fek
                                    *EINHEITIMP[einheitnr_out].fke;
    if( !is_equal( RV(ew,i),0.0,macheps) ){
        IS_EIGENWERT(iteration)=JA;
        ++anzahl;
    }
    RV(d_ew,i) = a_tof( line ,18,28)*EINHEITIMP[einheitnr_in].fek
                                    *EINHEITIMP[einheitnr_out].fke;
    if( is_equal( RV(d_ew,i),0.0,macheps) )
         RV(d_ew,i) = 0.1 * RV(ew,i);
 
    if( !is_equal( RV(d_ew,1),0.0,macheps) ) RV(d_ew,1)=0.0;
    if( !is_equal( RV(  ew,1),0.0,macheps) ) RV(  ew,1)=0.0;
 
 
    i = a_toi(line,35,37);
    k = a_toi(line,40,42);
    if( i!=0 && k!=0 ) ++anzahl;
 
    R(  intensit,i,k) = a_tof( line ,44,54);
    if( i!=0 && k!=0 && is_equal(R(intensit,i,k),0.0,macheps) )
            R(intensit,i,k) = 0.01;
 
    R(d_intensit,i,k) = a_tof( line ,56,67);
    if( is_equal( R(d_intensit,i,k),0.0,macheps) )
         R(d_intensit,i,k) = 0.1*R(intensit,i,k);
 
 }
 
 
   flag=0;
   for(i=1; i<=MXDIM(intensit); ++i)
       for(k=1; k<=MXDIM(intensit); ++k)
          if( ! is_equal( R(d_intensit,i,k),0.0,macheps) ) ++flag;
   if( flag==0 ){
       IS_INTENSITAET(  iteration)=NEIN;
       IS_MATRIXELEMENT(iteration)=NEIN;
   }
 
 /***********************************************/
 /* detailiertes gleichgewicht beruecksichtigen */
 /***********************************************/
if( IS_INTENSITAET(iteration) ){
 for(zeile=1; zeile<=anz_lines; ++zeile) /* zeile -> spalte */
   for(spalte=zeile;spalte<=anz_lines; ++spalte)
     if( spalte != zeile && !is_equal(R(intensit,zeile,spalte),0.0,macheps) )
       if( RV(ew,spalte)!=0 || RV(ew,zeile)!=0 || zeile==1){
             R(intensit,spalte,zeile) = R(intensit,zeile,spalte)
           *exp_( -RV(ew,spalte)/temperatur*EINHEITIMP[einheitnr_out].fek)
           *exp_(  RV(ew,zeile )/temperatur*EINHEITIMP[einheitnr_out].fek);
            if( ! is_equal(R(intensit,spalte,zeile),0.0,macheps) ){
               R(d_intensit,spalte,zeile)=0.1* R(intensit,spalte,zeile);
               ++anzahl;
            }
            else{
                  R(  intensit,spalte,zeile)=0.0;
                  R(d_intensit,spalte,zeile)=0.0;
                }
         }
}
else if( IS_MATRIXELEMENT(iteration) ){
 for(zeile=1; zeile<=anz_lines; ++zeile) /* zeile -> spalte */
   for(spalte=zeile;spalte<=anz_lines; ++spalte)
     if( spalte != zeile && !is_equal(R(intensit,zeile,spalte),0.0,macheps) ){
            R(intensit,spalte,zeile) = R(intensit,zeile,spalte);
            ++anzahl;
            R(d_intensit,spalte,zeile) = R(d_intensit,zeile,spalte);
   }
}
    /************************************************/
    /* sind Positionen und Intensitaeten zu fitten? */
    /************************************************/
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    for(i=1; i<=10; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    i = a_toi(line,30,56);
    IS_POSFIT(iteration)=NEIN;
    if( i!=0      ) IS_POSFIT(iteration) = JA;
 
    for(i=1; i<= 5; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    posname = file_name(line,30,56);
 
    for(i=1; i<=11; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    fpa  = a_toi(line,32,34);  fpzl = a_toi(line,43,45);
 
    for(i=1; i<= 9; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    fna  = a_toi(line,32,34);  fnzl = a_toi(line,43,45);
 
    for(i=1; i<= 3; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    /* Anzahl der datensaetze fuer fit: nrpose-nrposa+1 */
    NRPOSA(iteration) = a_toi(line,20,23); /* Datensatznr. anfang */
    NRPOSE(iteration) = a_toi(line,27,30); /* datensatznr. ende */
    if( ABS(NRPOSA(iteration)) > ABS(NRPOSE(iteration)) ){
        i                 = NRPOSA(iteration);
        NRPOSA(iteration) = NRPOSE(iteration);
        NRPOSE(iteration) = i;
    }
 
    POSDATANZ(iteration)= ABS(ABS(NRPOSE(iteration))-ABS(NRPOSA(iteration))+1);
    POSNAME(  iteration) = posname;
    FORMATPA( iteration) = fpa;   /* 7 mal 11 */
    FORMATPB( iteration) = fpzl;
    FORMATNA( iteration) = fna;   /* 2 mal 12 */
    FORMATNB( iteration) = fnzl;
 
if( IS_POSFIT(iteration)==JA  &&  POSDATANZ(iteration) != 0 ){
 
    i                  = POSDATANZ(iteration);
    POST(   iteration) = DOUBLE_ALLOC( i );
    POSANZ( iteration) = INT_ALLOC(    i );
    POSICIN(iteration) = DOUBLE_ALLOC( i );
    POSICQE(iteration) = DOUBLE_ALLOC( i );
    POSE(   iteration) = DOUBLEP_ALLOC(i );
    POSI(   iteration) = DOUBLEP_ALLOC(i );
    imax               = i;
    while( i>0 ){
     if( i==imax ){
        printf("reading file %s ....\n",POSNAME(iteration));
 
        if( (fp_sus=fopen(POSNAME(iteration),"rb"))==(FILE*)0 )
             read_error(2,fp_sus,POSNAME(iteration));
 
        anzdatnr=0;
 
     }
     loopmax = LOOPMAX+1;
     looppos:
        --loopmax;
       if(loopmax==0){
             printf("Error : cannot read %s !\n",POSNAME(iteration));
             printf("cause: wrong dataformat or record numbers not available.\n");
             fclose(fp_sus);
             exit(1);
        }
 
        line=fgets( string , buffer_size , fp_sus ); /*Kopfzeile*/
        datnr    = a_toi(line, 0, 4);
        anz_dat  = a_toi(line, 5, 9);
        pos_t    = a_tof(line,10,21);
        pos_icin = a_tof(line,22,33);
        pos_icqe = a_tof(line,34,45);
        if( pos_icin == 0.0 ) pos_icin = 1.0;
        if( pos_icqe == 0.0 ) pos_icqe = 1.0;
 
        if(NUMMERIERUNG(setup)!=JA){
          anzdatnr += datnr;
          datnr     = anzdatnr;
        }
        if(!is_pbekannt(datnr,iteration)){/* ueberlesen*/
           if( datnr>0 ) anz_dat = 2*(anz_dat/fpa+1);
           for(k=1; k<= anz_dat; ++k)
              line=fgets( string , buffer_size , fp_sus ); /*datenzeile*/
           goto looppos;
        }
        if( anz_dat <= 0 ) goto looppos;
        datene   = DOUBLE_ALLOC(anz_dat);
        dateni   = DOUBLE_ALLOC(anz_dat);
        if(datnr>0){
          j=1;
          for(k=0; k<= anz_dat/fpa; ++k){
              line=fgets( string , buffer_size , fp_sus ); /*datenzeile*/
              for(z=0; z<=fpa-1 && j<=anz_dat; ++z)
                  VALUE(datene,j++)=a_tof(line,z*fpzl,fpzl-1+z*fpzl);
          }
          j=1;
          for(k=0; k<= anz_dat/fpa; ++k){
              line=fgets( string , buffer_size , fp_sus ); /*datenzeile*/
              for(z=0; z<=fpa-1 && j<=anz_dat; ++z)
                  VALUE(dateni  ,j++)=a_tof(line,z*fpzl,fpzl-1+z*fpzl);
          }
       }
       else if(datnr<0)
          for(k=1; k<= anz_dat; ++k){
              line=fgets( string , buffer_size , fp_sus ); /*Datenzeile*/
              for(z=0; z<=fna-1; z+=2){
               VALUE(datene  ,k )=a_tof(line,z*fnzl,fnzl-1+z*fnzl);
               VALUE(dateni  ,k )=a_tof(line,(z+1)*fnzl,fnzl-1+(z+1)*fnzl);
              }
          }
 
           ipos = imax - i + 1;
           VALUE(POSANZ(iteration) ,ipos) = anz_dat;
           anzahl                        += anz_dat;
           VALUE(POST(  iteration) ,ipos) = pos_t;
           VALUE(POSICIN(iteration),ipos) = pos_icin;
           VALUE(POSICQE(iteration),ipos) = pos_icqe;
   /* daten nach aufsteigender energie ordnen */
           werte  = DOUBLE_ALLOC(anz_dat);
           werti  = DOUBLE_ALLOC(anz_dat);
           nummer = INT_ALLOC(   anz_dat);
           for(k=1; k<=anz_dat; ++k ){
               VALUE(werte,k) = VALUE(datene,k);
               VALUE(werti,k) = VALUE(dateni,k);
           }
           nummer  = sort( werte,nummer ,anz_dat);
           for(k=1; k<=anz_dat; ++k ){
               VALUE(datene,k) = VALUE(werte,VALUE(nummer,k));
               VALUE(dateni,k) = VALUE(werti,VALUE(nummer,k));
           }
           free_(werte);
           free_(werti);
           free_(nummer);
 
           VALUE(POSE(   iteration),ipos) = datene;
           VALUE(POSI(   iteration),ipos) = dateni;
 
           if( datene == (DOUBLE*)0 || dateni ==(DOUBLE*)0 ) {
               if( datene  != (DOUBLE*)0) free_(datene);
               if( dateni  != (DOUBLE*)0) free_(dateni);
           }
 
        --i;  /* beachte i:=anz_datensaetze wird hier um 1 erniedrigt ! */
 
    }/* end while(i>0) */
 }/* end if */
 else IS_POSFIT(iteration) = NEIN;
 
    /**********************************************/
    /* sind inverse Suszeptibilitaeten zu fitten ?*/
    /**********************************************/
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    for(i=1; i<=6; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    i = a_toi(line,30,56);
    IS_SUSFIT(iteration)=NEIN;
    if( i!=0      ) IS_SUSFIT(  iteration) = JA;
 
    for(i=1; i<= 5; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    susname = file_name(line,30,56);
 
    for(i=1; i<=11; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    fpa  = a_toi(line,32,34);  fpzl = a_toi(line,43,45);
 
    for(i=1; i<= 9; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    fna  = a_toi(line,32,34);  fnzl = a_toi(line,43,45);
 
    for(i=1; i<= 3; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    NRSUSA(iteration) = a_toi(line,53,56);
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    NRSUSB(iteration) = a_toi(line,53,56);
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    NRSUSC(iteration) = a_toi(line,53,56);
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    NRSUSP(iteration) = a_toi(line,53,56);
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    LAMBDA(iteration) = a_tof(line,48,56);
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    thetaname = file_name(line,30,56);
    fp_theta = fopen_errchk(thetaname,"rb");
    NAMETHETAFILE(iteration) = thetaname;
    LESETHETAFILE(iteration) = NEIN;
    THETAA(       iteration) = 0.0;
    if( fp_theta != (FILE*)0 )     LESETHETAFILE(iteration) = JA;
    fclose(fp_theta);
    if(LESETHETAFILE(iteration)==JA){
       NRPOSA(iteration) = i= 1;
       NRPOSE(iteration) = 1;
       FORMATPA(iteration)=7;   /* 7 mal 11 */
       FORMATPB(iteration)=11;
       FORMATNA(iteration)=2;   /* 2 mal 12 */
       FORMATNB(iteration)=12;
       anz = 200;
       xx  = DOUBLE_ALLOC(anz);
       ff  = DOUBLE_ALLOC(anz);
       lesedaten(fp_theta,xx,ff,&anz,thetaname,NUMMERIERUNG(setup),
                 &i,1,is_pbekannt,&datnr,iteration,&anzdatnr);
       fclose(fp_theta);
       THETAX(  iteration) = xx;
       THETAF(  iteration) = ff;
       THETAANZ(iteration) = anz;
    }
 
    SUSNAME(iteration)  = susname;
    FORMATPA(iteration) = fpa;   /* 7 mal 11 */
    FORMATPB(iteration) = fpzl;
    FORMATNA(iteration) = fna;   /* 2 mal 12 */
    FORMATNB(iteration) = fnzl;
 
    SUSAT(iteration) = (DOUBLE*)0;
    SUSBT(iteration) = (DOUBLE*)0;
    SUSCT(iteration) = (DOUBLE*)0;
    SUSPT(iteration) = (DOUBLE*)0;
    SUSASUS(iteration) = (DOUBLE*)0;
    SUSBSUS(iteration) = (DOUBLE*)0;
    SUSCSUS(iteration) = (DOUBLE*)0;
    SUSPSUS(iteration) = (DOUBLE*)0;
    IS_SUSA(iteration) =NEIN;
    IS_SUSB(iteration) =NEIN;
    IS_SUSC(iteration) =NEIN;
    IS_SUSP(iteration) =NEIN;
 
if( IS_SUSFIT(iteration)==JA &&
        (NRSUSA(iteration)!=0||NRSUSB(iteration)!=0||
         NRSUSC(iteration)!=0||NRSUSP(iteration)!=0)   ){
        i=0;
        if( NRSUSA(iteration) != 0 ) ++i;
 
        if( ABS(NRSUSB(iteration)) != 0 &&
            ABS(NRSUSB(iteration)) != ABS(NRSUSA(iteration)) ) ++i;
 
        if( ABS(NRSUSC(iteration)) != 0 &&
            ABS(NRSUSC(iteration)) != ABS(NRSUSA(iteration)) &&
            ABS(NRSUSC(iteration)) != ABS(NRSUSB(iteration)) ) ++i;
 
        if( ABS(NRSUSP(iteration)) != 0 &&
            ABS(NRSUSP(iteration)) != ABS(NRSUSA(iteration)) &&
            ABS(NRSUSP(iteration)) != ABS(NRSUSB(iteration)) &&
            ABS(NRSUSP(iteration)) != ABS(NRSUSC(iteration)) ) ++i;
 
         imax = i;
 
    while( i>0 ){
        anz_dat  = 200;
        datent   = DOUBLE_ALLOC(anz_dat);
        datensus = DOUBLE_ALLOC(anz_dat);
        lesedaten(fp_sus,datent,datensus,&anz_dat,
                  SUSNAME(iteration),NUMMERIERUNG(setup),
                  &i,imax,is_sbekannt,&datnr,iteration,&anzdatnr);
        if(ABS(datnr)==ABS(NRSUSA(iteration)) ){
           SUSAT(  iteration) = datent;
           SUSASUS(iteration) = datensus;
           ANZSA(  iteration) = anz_dat;
           anzahl            += anz_dat;
           if( datent!= (DOUBLE*)0 && datensus!=(DOUBLE*)0 )
              IS_SUSA(iteration) = JA;
           else{ if( datent  != (DOUBLE*)0) free_(datent);
                 if( datensus!= (DOUBLE*)0) free_(datensus); }
        }
        if(ABS(datnr)==ABS(NRSUSB(iteration)) ){
           SUSBT(  iteration) = datent;
           SUSBSUS(iteration) = datensus;
           ANZSB(  iteration) = anz_dat;
           anzahl            += anz_dat;
           if( datent!= (DOUBLE*)0 && datensus!=(DOUBLE*)0 )
              IS_SUSB(iteration) = JA;
           else{ if( datent  != (DOUBLE*)0) free_(datent);
                 if( datensus!= (DOUBLE*)0) free_(datensus); }
        }
        if(ABS(datnr)==ABS(NRSUSC(iteration)) ){
           SUSCT(  iteration) = datent;
           SUSCSUS(iteration) = datensus;
           ANZSC(  iteration) = anz_dat;
           anzahl            += anz_dat;
           if( datent!= (DOUBLE*)0 && datensus!=(DOUBLE*)0 )
              IS_SUSC(iteration) = JA;
           else{ if( datent  != (DOUBLE*)0) free_(datent);
                 if( datensus!= (DOUBLE*)0) free_(datensus); }
        }
        if(ABS(datnr)==ABS(NRSUSP(iteration)) ){
           SUSPT(  iteration) = datent;
           SUSPSUS(iteration) = datensus;
           ANZSP(  iteration) = anz_dat;
           anzahl            += anz_dat;
           if( datent!= (DOUBLE*)0 && datensus!=(DOUBLE*)0 )
              IS_SUSP(iteration) = JA;
           else{ if( datent  != (DOUBLE*)0) free_(datent);
                 if( datensus!= (DOUBLE*)0) free_(datensus); }
        }
    }/* end while(i>0) */
 }/* end if */
 else IS_SUSFIT(iteration) = NEIN;
 
    /*****************************************/
    /* sind Magnetisierungskurven zu fitten ?*/
    /*****************************************/
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    for(i=1; i<=6; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    i = a_toi(line,30,56);
    IS_MAGFIT(iteration)=NEIN;
    if( i!=0      ) IS_MAGFIT(  iteration) = JA;
 
    for(i=1; i<= 5; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    magname = file_name(line,30,56);
 
    for(i=1; i<=11; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    fpa  = a_toi(line,32,34);  fpzl = a_toi(line,43,45);
 
    for(i=1; i<= 9; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    fna  = a_toi(line,32,34);  fnzl = a_toi(line,43,45);
 
    for(i=1; i<= 6; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
 
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    NRMAGA(    iteration) = a_toi(line,37,40);
      MAGATEMP(iteration) = a_tof(line,47,55);
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    NRMAGB(    iteration) = a_toi(line,37,40);
      MAGBTEMP(iteration) = a_tof(line,47,55);
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    NRMAGC(    iteration) = a_toi(line,37,40);
      MAGCTEMP(iteration) = a_tof(line,47,55);
    while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    NRMAGP(iteration) = a_toi(line,37,40);
      MAGPTEMP(iteration) = a_tof(line,47,55);
 
if( IS_MAGFIT(iteration) == JA  ){
    flag = 0;
    if(NRMAGA(iteration)!=0 && MAGATEMP(iteration)<=0.0) flag=1;
    if(NRMAGB(iteration)!=0 && MAGBTEMP(iteration)<=0.0) flag=2;
    if(NRMAGC(iteration)!=0 && MAGCTEMP(iteration)<=0.0) flag=3;
    if(NRMAGP(iteration)!=0 && MAGPTEMP(iteration)<=0.0) flag=4;
    if(flag!=0){
        printf("Error : in file %s !\n",name);
        printf("Cause: the temperature to the %d.ten ",flag);
        printf("magnetisation curve is zero!\n");
        fclose(fp);
        exit(1);
    }
}
 
    MAGNAME(iteration)  = magname;
    FORMATPA(iteration) = fpa;   /* 7 mal 11 */
    FORMATPB(iteration) = fpzl;
    FORMATNA(iteration) = fna;   /* 2 mal 12 */
    FORMATNB(iteration) = fnzl;
 
    MAGAB(iteration) = (DOUBLE*)0;
    MAGBB(iteration) = (DOUBLE*)0;
    MAGCB(iteration) = (DOUBLE*)0;
    MAGPB(iteration) = (DOUBLE*)0;
    MAGAMAG(iteration) = (DOUBLE*)0;
    MAGBMAG(iteration) = (DOUBLE*)0;
    MAGCMAG(iteration) = (DOUBLE*)0;
    MAGPMAG(iteration) = (DOUBLE*)0;
    IS_MAGA(iteration) =NEIN;
    IS_MAGB(iteration) =NEIN;
    IS_MAGC(iteration) =NEIN;
    IS_MAGP(iteration) =NEIN;
 
if( IS_MAGFIT(iteration)==JA &&
        (NRMAGA(iteration)!=0||NRMAGB(iteration)!=0||
         NRMAGC(iteration)!=0||NRMAGP(iteration)!=0)   ){
        i=0;
        if( NRMAGA(iteration) != 0 ) ++i;
 
        if( ABS(NRMAGB(iteration)) != 0 &&
            ABS(NRMAGB(iteration)) != ABS(NRMAGA(iteration)) ) ++i;
 
        if( ABS(NRMAGC(iteration)) != 0 &&
            ABS(NRMAGC(iteration)) != ABS(NRMAGA(iteration)) &&
            ABS(NRMAGC(iteration)) != ABS(NRMAGB(iteration)) ) ++i;
 
        if( ABS(NRMAGP(iteration)) != 0 &&
            ABS(NRMAGP(iteration)) != ABS(NRMAGA(iteration)) &&
            ABS(NRMAGP(iteration)) != ABS(NRMAGB(iteration)) &&
            ABS(NRMAGP(iteration)) != ABS(NRMAGC(iteration)) ) ++i;
 
         imax = i;
 
    while( i>0 ){
        anz_dat  = 200;
        datent   = DOUBLE_ALLOC(anz_dat);
        datensus = DOUBLE_ALLOC(anz_dat);
        lesedaten(fp_sus,datent,datensus,&anz_dat,
                  MAGNAME(iteration),NUMMERIERUNG(setup),
                  &i,imax,is_mbekannt,&datnr,iteration,&anzdatnr);
        if(ABS(datnr)==ABS(NRMAGA(iteration)) ){
           MAGAB(  iteration) = datent;
           MAGAMAG(iteration) = datensus;
           ANZMA(  iteration) = anz_dat;
           anzahl            += anz_dat;
           if( datent!= (DOUBLE*)0 && datensus!=(DOUBLE*)0 )
              IS_MAGA(iteration) = JA;
           else{ if( datent  != (DOUBLE*)0) free_(datent);
                 if( datensus!= (DOUBLE*)0) free_(datensus); }
        }
        if(ABS(datnr)==ABS(NRMAGB(iteration)) ){
           MAGBB(  iteration) = datent;
           MAGBMAG(iteration) = datensus;
           ANZMB(  iteration) = anz_dat;
           anzahl            += anz_dat;
           if( datent!= (DOUBLE*)0 && datensus!=(DOUBLE*)0 )
              IS_MAGB(iteration) = JA;
           else{ if( datent  != (DOUBLE*)0) free_(datent);
                 if( datensus!= (DOUBLE*)0) free_(datensus); }
        }
        if(ABS(datnr)==ABS(NRMAGC(iteration)) ){
           MAGCB(  iteration) = datent;
           MAGCMAG(iteration) = datensus;
           ANZMC(  iteration) = anz_dat;
           anzahl            += anz_dat;
           if( datent!= (DOUBLE*)0 && datensus!=(DOUBLE*)0 )
              IS_MAGC(iteration) = JA;
           else{ if( datent  != (DOUBLE*)0) free_(datent);
                 if( datensus!= (DOUBLE*)0) free_(datensus); }
        }
        if(ABS(datnr)==ABS(NRMAGP(iteration)) ){
           MAGPB(  iteration) = datent;
           MAGPMAG(iteration) = datensus;
           ANZMP(  iteration) = anz_dat;
           anzahl            += anz_dat;
           if( datent!= (DOUBLE*)0 && datensus!=(DOUBLE*)0 )
              IS_MAGP(iteration) = JA;
           else{ if( datent  != (DOUBLE*)0) free_(datent);
                 if( datensus!= (DOUBLE*)0) free_(datensus); }
        }
    }/* end while(i>0) */
 }/* end if */
 else IS_MAGFIT(iteration) = NEIN;
 
    /*********************************/
    /* Wichtung der Nebenbedingungen */
    /*********************************/
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
 
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 2.==== */
    wew  = a_tof(line,37,56);
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 2.==== */
    wpos = a_tof(line,37,56);
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 2.==== */
    wint = a_tof(line,37,56);
    wmat = wint;
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 2.==== */
    wsus = a_tof(line,37,56);
    while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 2.==== */
    wmag = a_tof(line,37,56);
 
    wew  *= IS_EIGENWERT(    iteration);
    wpos *= IS_POSFIT(       iteration);
    wint *= IS_INTENSITAET(  iteration);
    wmat *= IS_MATRIXELEMENT(iteration);
    wsus *= IS_SUSFIT(       iteration);
    wmag *= IS_MAGFIT(       iteration);
 
    wew =ABSD(wew ); wpos=ABSD(wpos); wint=ABSD(wint);
    wmat=ABSD(wmat); wsus=ABSD(wsus); wmag=ABSD(wmag);
    wsum=wew+wpos+wint+wmat+wsus+wmag;
 
    if( wsum==0.0 ){
       wew=wpos=wint=wmat=wsus=wmag=1.0;
       wsum=wew+wpos+wint+wmat+wsus+wmag;
    }
 
    WEW( iteration) = wew /wsum;
    WPOS(iteration) = wpos/wsum;
    WMAT(iteration) = wmat/wsum;
    WINT(iteration) = wint/wsum;
    WSUS(iteration) = wsus/wsum;
    WMAG(iteration) = wmag/wsum;
 
 
    /*************************/
    /* Fitroutine auswaehlen */
    /*************************/
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    for(i=1; i<=15; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    i = a_toi(line,33,56);
    if(i>5 || i<1 ) i = FITNR;
    FITROUTINENNR(iteration) = i;
    for(i=1; i<=4; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    i = a_toi(line,33,56);
    if( i==0 ){
        if( FITROUTINENNR(iteration)==1 ) FITMAX(iteration)=MAXDOWNHILL;
        if( FITROUTINENNR(iteration)==2 ) FITMAX(iteration)=MAXPOWELL;
        if( FITROUTINENNR(iteration)==3 ) FITMAX(iteration)=MAXGRADIENT;
        if( FITROUTINENNR(iteration)==4 ) FITMAX(iteration)=MAXVA05A;
        if( FITROUTINENNR(iteration)==5 ) FITMAX(iteration)=MAXDUMMY;
    }
    else FITMAX( iteration ) =i;
 
 
    /**********************************/
    /*Welche CF-Operatoren sind frei? */
    /**********************************/
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
 
    fix = vr_alloc(12);
 
    for(k=1; k<=9; ++k){
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
       c    = VALUE(line,2);
       RV(fix,k) = 1.0;    /* != bedeutet Parameter wird festgehalten */
       if( c=='O' ){
           ++anz_par;
           i = a_toi(line,21,57);
           if( i==0 ) {++anz_var; RV(fix,k) = 0.0;}
       }
    }
    for(k=10; k<=12; ++k){
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
       c    = VALUE(line,2);
       RV(fix,k) = 1.0;    /* != bedeutet Parameter wird festgehalten */
       if( c=='J' ){
           ++anz_par;
           i = a_toi(line,21,57);
           if( i==0 ) {++anz_var; RV(fix,k) = 0.0;}
       }
    }
 
 
/****************************************************/
/* Soll Menue waehrend des Fits ausgegeben werden ? */
/****************************************************/
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
    while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
    for(i=1; i<=5; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    i = a_toi(line,27,56);
    if( i==0 ) IS_MENUE(iteration) = NEIN;
    else       IS_MENUE(iteration) = JA;
    for(i=1; i<=5; ++i)
      while(  *(line=fgets( string , buffer_size , fp )) != '|');/* 1.|*/
    i = a_toi(line,27,56);
    if( i==0 ){
        if( FITROUTINENNR(iteration)==1 ) MENUE(iteration)=MDOWNHILL;
        if( FITROUTINENNR(iteration)==2 ) MENUE(iteration)=MPOWELL;
        if( FITROUTINENNR(iteration)==3 ) MENUE(iteration)=MGRADIENT;
        if( FITROUTINENNR(iteration)==4 ) MENUE(iteration)=MVA05A;
        if( FITROUTINENNR(iteration)==5 ) MENUE(iteration)=MDUMMY;
    }
    if( i!=0 && FITMAX(iteration)/i==0  )
       IS_MENUE(iteration) = NEIN;
    if( i!=0 && FITMAX(iteration)/i!=0  )
        MENUE(iteration) = i;
 
 
 
    fclose(fp);
 
      ANZAHL(neben)        =  anzahl;
      ANZ_PAR(neben)       =  anz_par; /* wieviele operatoren */
      ANZ_VAR(neben)       =  anz_var; /* wieviele davon variabel */
      FIX(    neben)       =  fix;
      EIGENWERTE(neben)    =   ew;
    D_EIGENWERTE(neben)    = d_ew;
      INTENSITAETEN(neben) =   intensit;
    D_INTENSITAETEN(neben) = d_intensit;
 
 
    return( neben );
}
/*------------------------------------------------------------------------------
                                  file_name()
------------------------------------------------------------------------------*/
CHAR  *file_name(string,anfang,ende)
    CHAR *string;
    INT  anfang,ende;
{
    CHAR   *buffer, *neu_buffer;
    INT    i,len,neu_anfang,neu_ende,dummy;
 
    /* buffer hat die laenge ende-anfang+2         */
    /* wobei *(buffer+ende-anfang + 2 ) = '\0' ist */
    buffer     = a_tos(string,anfang,ende);
    len        = ende-anfang+1;
 
    /* anzahl der fuehrenden leerstellen bestimmen */
    neu_anfang = 0;
    for(i=0; i<=len-1; ++i)
        if( *(buffer+i) != ' ')  break;
    neu_anfang = i;
 
    /* anzahl der folgenden  leerstellen bestimmen */
    neu_ende   = len-1;
    for(i=len-1; i>=0; --i)
        if( *(buffer+i) != ' ')  break;
    neu_ende   = i;
 
    if( neu_anfang > neu_ende ){
        dummy      = neu_anfang;
        neu_anfang = neu_ende;
        neu_ende   = dummy;
    }
 
    len        = neu_ende-neu_anfang+1;
    neu_buffer = STRING_ALLOC(len+1);
    for(i=neu_anfang; i<=neu_ende; ++i)
        VALUE(neu_buffer,i-neu_anfang) = VALUE(buffer,i);
    VALUE(neu_buffer,len+1) = '\0';
    free_(buffer);
 
    return( neu_buffer );
}
/*------------------------------------------------------------------------------
ENDEMODUL    E I N G A B E    C
------------------------------------------------------------------------------*/
