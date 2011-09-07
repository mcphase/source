 
/*-----------------------------------------------------------------------------
 
                                M   O  D  U  L
 
                                 INTENSITY C
 
-------------------------------------------------------------------------------
Benutzte Literatur    :  J.Phys.Chem. Solids 1972 Vol.33 pp 59-68
                         "Birgeneau"
-----------------------------------------------------------------------------*/
 
 
/*----------------------------------------------------------------------------
Includedateien holen
-----------------------------------------------------------------------------*/
#include <stdio.h>          /* damit FILE definiert wird               */
#include <stdlib.h>
#include <time.h>
#include <math.h>           /* damit sqrt in define_j.c definiert wird */
#define pi (4.0*atan(1.0))  /* atan() braucht <math.h>                 */
#include "types.c"          /* benutze Datentypen laden                */
 
 
#define r0    (-1.91*_e*_e/_m )/*  * 10**(-12) cm  = -0.54*10(-12)cm */
#define const (r0*r0*pi)       /* in barn            2   */
                               /* spaeter noch mit gj  multiplizieren*/
/*----------------------------------------------------------------------------
extern definierte Funktionen
-----------------------------------------------------------------------------*/
extern DOUBLE  exp_();           /* definiert in ORTHO.C */
extern void    fit_ortho();      /* definiert in ORTHO.C */
extern INT     is_equal();       /* definiert in DIAHERMX.C */
extern INT     isimplementiert();/* definiert in CFIELD.C.*/
extern DOUBLE  is_null();        /* definiert in DIAHERMX.C */
extern INT       null();         /* definiert in DIAHERMX.C */
extern DOUBLE  einh_Kelvin();    /* definiert in CFIELD.C */
extern DOUBLE  accuracy();       /* definiert in DIAHERMX.c */
extern INT     is_einheit_imp(); /* definiert in CFIELD.C */
extern MATRIX  *mx_alloc();      /* definiert in MATRIX.C */
extern VEKTOR  *vr_alloc();      /* definiert in MATRIX.C */
extern void     free_vr();       /* definiert in MATRIX.C */
extern void     free_mx();       /* definiert in MATRIX.C */
extern KOMPLEX *ckon();          /* definiert in KOMPLEX.C*/
extern KOMPLEX *cadd();          /* definiert in KOMPLEX.C*/
extern KOMPLEX *csub();          /* definiert in KOMPLEX.C*/
extern KOMPLEX *cmult();         /* definiert in KOMPLEX.C*/
extern KOMPLEX *vskalar();       /* definiert in KOMPLEX.C*/
extern INT     write_titlecom(); /* definiert in DIAHERMX.C */
extern INT     test_nullmatrix();/* definiert in DIAHERMX.C */
extern CHAR    up();             /* definiert in CFIELD.C */
extern INT     isreell();        /* definiert in CFIELD.C */
extern MATRIX    *calc_iBmag();  /* definiert in CFIELD.C */
extern EWPROBLEM *solve();       /* definiert in CFIELD.C */
extern EWPROBLEM *set_ewproblem();/* definiert in DIAHERMX.c */
extern BRUCH     *is_rational(); /* definiert in STEVENS.C*/
extern LONG      ggt_l();        /* definiert in STEVENS.C*/
extern FILE     *fopen();        /* definiert in stdio.h  */
 
 
extern IONEN     IONENIMP[];     /* definiert in CFIELD.C*/
extern SYMMETRIE SYMLISTE[];     /* definiert in CFIELD.C*/
extern EINHEIT   EINHEITIMP[];   /* definiert in CFIELD.C*/
extern DOUBLE    alpha_J[];      /* definiert in THETA.C */
extern DOUBLE     beta_J[];      /* definiert in THETA.C */
extern DOUBLE    gamma_J[];      /* definiert in THETA.C */
extern DOUBLE    omegan0n();      /* definiert in CFIELD.C*/
extern DOUBLE    omegan1n();      /* definiert in CFIELD.C*/
extern DOUBLE    omegan2n();      /* definiert in CFIELD.C*/
extern DOUBLE    omegan3n();      /* definiert in CFIELD.C*/
extern DOUBLE    omegan4n();      /* definiert in CFIELD.C*/
extern DOUBLE    omegan5n();      /* definiert in CFIELD.C*/
extern DOUBLE    omegan6n();      /* definiert in CFIELD.C*/
extern DOUBLE    epn0n();       /* definiert in CFIELD.C*/
extern DOUBLE    epn1n();       /* definiert in CFIELD.C*/
extern DOUBLE    epn2n();       /* definiert in CFIELD.C*/
extern DOUBLE    epn3n();       /* definiert in CFIELD.C*/
extern DOUBLE    epn4n();       /* definiert in CFIELD.C*/
extern DOUBLE    epn5n();       /* definiert in CFIELD.C*/
extern DOUBLE    epn6n();       /* definiert in CFIELD.C*/
 
extern DOUBLE    spline();      /* definiert in SPLINE.C*/
extern INT       lesedaten();   /* definiert in EINGABE.C*/
extern INT       is_pbekannt(); /* definiert in EINGABE.C*/
extern FILE *fopen_errchk();         /* definiert in EINGABE.C*/
 
/*----------------------------------------------------------------------------
   Internal function declarations
-----------------------------------------------------------------------------*/
INT is_parametersatz_null(ITERATION *iter, INT symmetrienr, DOUBLE macheps);
void parametersatz(FILE*, CHAR, KRISTALLFELD*, INT, CHAR*, CHAR);
void richtung(ITERATION *iteration);
void raus_kpoly(SETUP*, EWPROBLEM*, KRISTALLFELD*, DOUBLE, DOUBLE, DOUBLE, CHAR*);
void raus_suszept(SETUP *setup, ITERATION *iteration, EWPROBLEM *ewproblem, KRISTALLFELD *kristallfeld);
void raus_magnetm(SETUP*, EWPROBLEM*, KRISTALLFELD*, DOUBLE, DOUBLE, DOUBLE, CHAR*);
void tabelle(FILE *fp, INT anz_niveaus, MATRIX *inten);
void mittlung(FILE*, INT, MATRIX *inten, VEKTOR*, MATRIX*, INT, INT, CHAR*);
void raus_mkommentar(FILE*, ITERATION*, DOUBLE, DOUBLE, DOUBLE, INT, INT, INT, CHAR*);
void raus_kkommentar(FILE*, ITERATION*, DOUBLE, DOUBLE, DOUBLE, CHAR*);
void raus_kommentar(FILE*, DOUBLE, DOUBLE, DOUBLE, DOUBLE, INT, CHAR*);
void kopf(FILE *fp, INT anz_niveaus, INT flag);

/*----------------------------------------------------------------------------
                              output()
-----------------------------------------------------------------------------*/
void output( setup,ewproblem,kristallfeld,modus )
   SETUP        *setup;
   EWPROBLEM    *ewproblem;
   KRISTALLFELD *kristallfeld;
   CHAR         modus;         /* 'A','a'   'B','b'  'V','v'  */
{
    CHAR art,*t13s,*dummy;
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11,*t12,*t13,*t14,*t15,*t16,*t17,*t18,*t19,*t20;
    CHAR *t21,*t22,*t23,*t24,*t25,*t26,*t27,*t28,*t29,*t30;
    CHAR *t31,*t32,*t33,*t34,*t35,*t36,*t37,*t38,*t39,*t40;
    CHAR *t41,*t42,*t43,*t44,*t45,*t46,*t47,*t48,*t49,*t50;
    CHAR *t51,*t52,*t53,*t54,*t55,*t56,*t57/*,*t58,*t59,*t60*/;
 
    CHAR *t31p,*t31m,*t32p,*t32m,*t32mm,*t38p,*t38m,*t39p,*t39m;
    MATRIX    *entartung,*ev;
    MATRIX    *aJtb2,*aJtb_2();
    VEKTOR    *ew,*v,*ev_ir,*ev_ks;
    ITERATION *iteration;
    FILE      *fopen(),*fp,*ewev;
    DOUBLE    temperatur,re,im,energie,gesamte_intensitaet;
    DOUBLE    gj,shift,j,mj,zu_summe,zustandssumme(),faktor;
    DOUBLE    Bx,By,Bz/*,anf_temp,end_temp,lambda*/;
    DOUBLE    macheps/*,Bx_ex,By_ex,Bz_ex,Bx_mol,By_mol,Bz_mol*/;
    DOUBLE    Bx2,By2,Bz2;
    DOUBLE /* aJxb_2,aJyb_2,aJzb_2,*/sumx,sumy,sumz;
    DOUBLE    mat_Jx2(),mat_Jy2(),mat_Jz2();
    DOUBLE    anf_feld,end_feld,temp;
    CHAR      *ionname/*,*symname*/,*einheit_in,*einheit_out,*f_3s();
    INT       symmetrienr,*gi/*,gi_ze,f*/;
    INT       einheitnr_in,einheitnr_out;
    INT       ionennr,elektronen4f/*,e_4f*/,anz_niveaus;
    INT       dimj,zeile,spalte,i,k,r,s/*,q*/,ze,sp,first_line,ps_null;
/*  DOUBLE    theta; */
    KOMPLEX *mat_Jx(),*mat_Jy(),*mat_Jz();
    DOUBLE magnetm();
 
    if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ){
        fp=fopen_errchk(FILENAME(kristallfeld),"w");
	
        write_titlecom(fp);
t01="#-------------------------------------------------------------- \n";
t02="#    O U T P U T         ";

t03="#-------------------------------------------------------------- \n";
t04="#\n#\n";
/*        fprintf(fp,"%s",t01);*/fprintf(fp,"%s",t02);/*fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);*/
  time_t curtime;
  struct tm *loctime;
   curtime=time(NULL); loctime=localtime(&curtime); fprintf(fp,"%s",asctime(loctime));
 
    }
    else      fp=fopen_errchk(FILENAME(kristallfeld),"a");
 
 
    iteration   = ITERATION(   kristallfeld );
 
 
t01="#-------------------------------------------------------------- \n";
t02="#!Temperature of the sample       T=%7.2f Kelvin                \n";
/*3="#-------------------------------------------------------------- \n";
t04="#\n";*/
    temperatur  = TEMPERATUR( iteration );
if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ){
    fprintf(fp,"%s",t01);fprintf(fp,t02,temperatur);
/* fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);*/
}
 
t01="#-------------------------------------------------------------- \n";
t02="#!Ion                    IONTYPE= %4s                          \n";
t03="#--------------------------------------------------------------\n";
t04="#                                                              \n";
t05="#!Lande factor of the ion   gJ= %8.6f                       \n";
t06="#                                                              \n";
t07="# Total angular momentum J of the                               \n";
t08="#!Spin - orbit - level     J= %4.1f                            \n";
t09="#                                                              \n";
t10="#!Electrons in 4f shell   Ne= %2d                              \n";
t11="#                                                              \n";
t12="#                                                              \n";
t13="#--------------------------------------------------------------\n";
t14="#                                                              \n";
t15="#  local point symmetry                                        \n";
t16="#! of the ion              PS= %10s                            \n";
t17="#                                                              \n";
t18="#! Symmetry_number            =%2d                              \n";
t19="#                                                              \n";
t20="#-------------------------------------------------------------- \n";
t21="#\n";
    ionname     = IONNAME( iteration );
    dimj        = DIMJ(    iteration );
    gj          = GJ(      iteration );
    ionennr     = isimplementiert( ionname );
    ionname     = IONENIMP[ ionennr ].ionname;
    elektronen4f= IONENIMP[ ionennr ].elektronen_in_vier_f;
    symmetrienr = SYMMETRIENR( kristallfeld );
/*  symname     = SYMLISTE[ symmetrienr ].symname; */
    einheitnr_out = EINHEITNROUT( iteration );
    einheit_out   = EINHEITIMP[ einheitnr_out ].einheit;
    entartung     = ewproblem->entartung;
    gi            = ewproblem->gi;
    ev            = ewproblem->eigenvektoren;
    ew            = ewproblem->eigenwerte;
    shift         = ewproblem->shift;
    anz_niveaus   = ANZ_ZE(entartung);
    macheps       = ewproblem->eps_machine;
    aJtb2         = aJtb_2(ewproblem,macheps);
 
    Bx =  B1(iteration);
    By =  B2(iteration);
    Bz =  B3(iteration);

    Bx2=  B1S(iteration);
    By2=  B2S(iteration);
    Bz2=  B3S(iteration);
 
   Bx += B1MOL(iteration);
    By += B2MOL(iteration);
    Bz += B3MOL(iteration);
    ps_null = is_parametersatz_null(iteration,symmetrienr,macheps);
    if(ps_null==JA && Bx==0.0 && By==0.0 && Bz==0.0)
      {symmetrienr = 10;/*symname = "sphaerisch";*/}
    if(ps_null==JA && (Bx!=0.0 || By!=0.0 || Bz!=0.0))
      {symmetrienr = 9;/*symname = "azimutal  ";*/}
 
 
 if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ){
   /*fprintf(fp,"%s",t01);*/fprintf(fp,t02,ionname);
   /* fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);*/
    fprintf(fp,t05,gj);fprintf(fp,"%s",t06);
    fprintf(fp,"%s",t07);fprintf(fp,t08,(DOUBLE)( (dimj-1)/2.0 ) );
  /*  fprintf(fp,"%s",t09);*/fprintf(fp,t10,elektronen4f);
  /*  fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
    fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);*/
  /*  fprintf(fp,"%s",t15);fprintf(fp,t16,symname);*/
  /*  fprintf(fp,"%s",t17);fprintf(fp,t18,symmetrienr);*/
  /*  fprintf(fp,"%s",t19);fprintf(fp,"%s",t20); fprintf(fp,"%s",t21);*/
 }
 
 
 
 
    symmetrienr  = SYMMETRIENR( kristallfeld ); /* stehen lassen !*/
    einheitnr_in = EINHEITNRIN( iteration );
    einheit_in   = EINHEITIMP[ einheitnr_in ].einheit;
    art          = EINGABEPARAMETERART(kristallfeld);
 if( *(FILENAME(kristallfeld)+16) == *(ORTHO+16) ){
  t01="#-------------------------------------------------------------- \n";
  t02="# Paramter          :     %2d von %2d                      \n";
  t03="#-------------------------------------------------------------- \n";
/*  fprintf(fp,"%s",t01);*/
  fprintf(fp,t02,ZEILE(iteration),MAX_ZEILE(iteration));
/*  fprintf(fp,"%s",t03);*/
 }
 
 if(ps_null == NEIN){
    switch(modus){
       case XW  :
       case AKQ :
       case BKQ :
       case LKQ :
                  if( (symmetrienr==8 &&
                       *(FILENAME(kristallfeld)+16) != *(ORTHO+16)) ||
                       *(FILENAME(kristallfeld)+16) == *(ORTHO+16)   )
                      parametersatz(fp,XW,kristallfeld,ionennr,einheit_in,art);
 
                  parametersatz(fp,AKQ,kristallfeld,ionennr,einheit_in,art);
                  parametersatz(fp,BKQ,kristallfeld,ionennr,einheit_in,art);
 
                if( *(FILENAME(kristallfeld)+16) == *(ORTHO+16) ) break;
 
       case DKQ :
       case WKQ :
       case VKQ :
                  if( (modus==VKQ || modus==WKQ || modus== DKQ || modus==LKQ)
                                  && isreell(symmetrienr,ionname) ){
                     if( symmetrienr == 8 )
                         parametersatz(fp,XW,kristallfeld,ionennr,einheit_in,art);
                     parametersatz(fp,AKQ,kristallfeld,ionennr,einheit_in,art);
                     parametersatz(fp,BKQ,kristallfeld,ionennr,einheit_in,art);
                  }
                  parametersatz(fp,DKQ,kristallfeld,ionennr,einheit_in,art);
                  parametersatz(fp,LKQ,kristallfeld,ionennr,einheit_in,art);
                  parametersatz(fp,VKQ,kristallfeld,ionennr,einheit_in,art);
                  parametersatz(fp,WKQ,kristallfeld,ionennr,einheit_in,art);
 
    }
}
 
 
if( Bx!=0.0 || By!=0.0 || Bz!=0.0 ){
t01="#\n";
t02="#-------------------------------------------------------------- \n";
t03="#                   Magnetic field B in Tesla .                    \n";
t04="#                                                              \n";
t05="#                   H = H   +  H       + H                     \n";
t06="#                        CF     mag_ex    mag_mol              \n";
t07="#                                                              \n";
t08="#--------------------------------------------------------------\n";
t09="#                                                              \n";
t10="#  H        = - g  my   J  B                                   \n";
t11="#   mag_ex       J   B  -  - ex                                \n";
t12="#                                                              \n";
t13="#--------------------------------------------------------------\n";
t14="#                                                              \n";
t15="#  H        = - 2 ( g - 1 ) my   J  B                          \n";
t16="#   mag_mol          J        B  -  - mol                      \n";
t17="#                                                              \n";
t18="#--------------------------------------------------------------\n";
 
t19="#                                                              \n";
t20="#! Bx = Bx_ex   =   %11.2f                                \n";
t21="#! By = By_ex   =   %11.2f                                \n";
t22="#! By = Bz_ex   =   %11.2f                                \n";
 
t23="#                                                              \n";
t24="#! Bx_mol  =   %11.2f                                     \n";
t25="#! By_mol  =   %11.2f                                     \n";
t26="#! Bz_mol  =   %11.2f                                     \n";
t27="#                                                              \n";
t28="#-------------------------------------------------------------- \n";
t29="#\n";
    if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) )  fprintf(fp,"%s",t01);
    fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
    fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
    fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
    fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,"%s",t15);fprintf(fp,"%s",t16);
    fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);
 
 if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ||
     (*(FILENAME(kristallfeld)+16) == *(ORTHO+16)&&!IS_MAGFIT(iteration)) ){
 
     if( *(FILENAME(kristallfeld)+16) == *(ORTHO+16) ){
         B1(iteration) = 0.0;
         B2(iteration) = 0.0;
         B3(iteration) = 0.0;
     }
 
    fprintf(fp,"%s",t19);
    fprintf(fp,t20,is_null(B1(iteration),0.001 ));
    fprintf(fp,t21,is_null(B2(iteration),0.001 ));
    fprintf(fp,t22,is_null(B3(iteration),0.001 ));
  }
 if( *(FILENAME(kristallfeld)+16) == *(ORTHO+16)&&IS_MAGFIT(iteration) ){
 
     richtung( iteration ); /* R1(),... setzen */
t19="#                                                              \n";
t20="# for the Approximation of the magnetic Momente thats in a      \n";
t21="# external Field B_ex in direction of  [ %5d, %5d, %5d].   \n";
    fprintf(fp,"%s",t19); fprintf(fp,"%s",t20);
    fprintf(fp,t21,R1(iteration),R2(iteration),R3(iteration));
  }
    fprintf(fp,"%s",t23);
    fprintf(fp,t24,is_null(B1MOL(iteration),0.001 ));
    fprintf(fp,t25,is_null(B2MOL(iteration),0.001 ));
    fprintf(fp,t26,is_null(B3MOL(iteration),0.001 ));fprintf(fp,"%s",t23);
    fprintf(fp,"%s",t27);fprintf(fp,"%s",t28);
    if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ) fprintf(fp,"%s",t29);
}
if( (Bx2!=0.0 || By2!=0.0 || Bz2!=0.0) && (symmetrienr == 10 || symmetrienr == 0) ){
t02="#-------------------------------------------------------------- \n";
t03="#                Anisotropy parameters in meV.                  \n";
t04="#-------------------------------------------------------------- \n";
t05="#         H= + Dx2 Jx ^ 2+ Dy2 Jy ^ 2+ Dz2 Jz ^ 2               \n";
t06="#! Dx2 =   %11.2f                                     \n";
t07="#! Dy2 =   %11.2f                                     \n";
t08="#! Dz2 =   %11.2f                                     \n";
t09="#                                                              \n";
t10="#-------------------------------------------------------------- \n";
    fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04); fprintf(fp,"%s",t05);
    fprintf(fp,t06,is_null(B1S(iteration),0.001 ));
    fprintf(fp,t07,is_null(B2S(iteration),0.001 ));
    fprintf(fp,t08,is_null(B3S(iteration),0.001 ));
    fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);fprintf(fp,"%s",t29);
}
 
t01="#-------------------------------------------------------------- \n";
t02="# Energy Eigenvalues are in  %6s.          	            \n";
t03="#--------------------------------------------------------------\n";
t04="#!Nr of different energy levels   noflevels=%2d                 \n";
t05="#!Energy shift   Eshift=             %17.4f \n";
t06="#                                                              \n";
t07="# Because of the calibration freedom the smallest energy       \n";
t08="# eigenvalue is shifted to zero. You can get the energy        \n";
t09="# eigen-value of the applied eigen-value problem               \n";
t10="# by shifting the energy by the added energy value above       \n";
t11="#                                                              \n";
 
t13 ="#! E(%2d)=   %11.4f        Degeneracy =  %2d -fold            \n";
t13s="#!*E(%2d)=   %11.4f        Degeneracy =  %2d -fold            \n";
 
t14="#                                                              \n";
t15="#--------------------------------------------------------------\n";
t15="# Them with  *  marked Energy eigenvalues have a non-          \n";
t16="# vanishing Matrix element of the Ground state E( 1).          \n";
t17="#-------------------------------------------------------------- \n";
t18="#\n";
 
    fprintf(fp,"%s",t01);fprintf(fp,t02,einheit_out);
    fprintf(fp,"%s",t03);fprintf(fp,t04,anz_niveaus);
    fprintf(fp,t05,shift);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);
    fprintf(fp,"%s",t08);fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
/*    fprintf(fp,"%s",t11);*/
 
 
    for( zeile=1 ; zeile<=anz_niveaus ; ++zeile){
       dummy = t13s;
       if( is_equal(R(aJtb2,zeile,1),0.0,0.005) ) dummy=t13;
       fprintf(fp,dummy,zeile,
       RV(ew,(INT)R(entartung,zeile,1))
       *EINHEITIMP[einheitnr_in].fek*EINHEITIMP[einheitnr_out].fke,
                  VALUE(gi,zeile) );
    }
  /*  fprintf(fp,"%s",t14);*/fprintf(fp,"%s",t15);
    fprintf(fp,"%s",t16);fprintf(fp,"%s",t17);
    if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ) fprintf(fp,"%s",t18);
 
 
 
t01="#-------------------------------------------------------------- \n";
t02="# The orthonormal characteristic Eigenvectors  |i,r>  with     \n";
t03="#                                                              \n";
t04="#             H |i,r>  = ( E  + E      ) |i,r>  , r=1 ... n    \n";
t05="#                           i    shift                     i   \n";
t06="#                                                              \n";
t07="#                            and                               \n";
t08="#                                                              \n";
t09="#                    <i,r|j,s> = D   D                         \n";
t10="#                                 ij  rs                       \n";
t11="#                                                              \n";
t12="#              D  = Kronecker- Delta function                  \n";
t13="#               ij                                             \n";
t14="#                                                              \n";
t15="# the |i,r> are also orthonormal .                             \n";
t16="# We consider in this program  Crystal Field spliting          \n";
t17="# at the lowest Spin-orbit-level (Ground state)                \n";
t18="#                                                              \n";
t19="# 2S+1                                                         \n";
t20="#     L  of the calculated Ion. The |i,r> are in               \n";
t21="#      J                                                       \n";
t22="#                                                              \n";
t23="# a more developed  Eigenfunction   			    \n";
t24="#                                                              \n";
t25="# system  [ | J  M   > ]               ,where                  \n";
t26="#                 J     M = -J,...,J                           \n";
t27="#                        J                                     \n";
t46="#                               ,                              \n";
t47="#              <  J  M   |  J  M  >  =  D    ,    .            \n";
t48="#                     J         J        M  M                  \n";
t49="#                                         J  J                 \n";
t28="#--------------------------------------------------------------\n";
t29="#                                                              \n";
 
t31p="# I%2d,%2d> =  ( %7.4f + %7.4f*i ) I%4.1f %5.1f>              \n";
t31m="# I%2d,%2d> =  ( %7.4f - %7.4f*i ) I%4.1f %5.1f>              \n";
t32p  ="#          + ( %7.4f + %7.4f*i ) I%4.1f %5.1f>              \n";
t32m  ="#          + ( %7.4f - %7.4f*i ) I%4.1f %5.1f>              \n";
t32mm ="#          - ( %7.4f + %7.4f*i ) I%4.1f %5.1f>              \n";
 
t36="# I%2d,%2d> =    %7.4f               I%4.1f %5.1f>              \n";
t37="# I%2d,%2d> =              %7.4f*i   I%4.1f %5.1f>              \n";
 
t38p="#          +   %7.4f               I%4.1f %5.1f>              \n";
t38m="#          -   %7.4f               I%4.1f %5.1f>              \n";
t39p="#          +             %7.4f*i   I%4.1f %5.1f>              \n";
t39m="#          -             %7.4f*i   I%4.1f %5.1f>              \n";
 
t43="#                                                              \n";
t44="#-------------------------------------------------------------- \n";
 
if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ){
fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);fprintf(fp,"%s",t05);
fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,"%s",t15);
fprintf(fp,"%s",t16);fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);fprintf(fp,"%s",t19);fprintf(fp,"%s",t20);
fprintf(fp,"%s",t21);fprintf(fp,"%s",t22);fprintf(fp,"%s",t23);fprintf(fp,"%s",t24);fprintf(fp,"%s",t25);
fprintf(fp,"%s",t26);fprintf(fp,"%s",t27);
fprintf(fp,"%s",t46);fprintf(fp,"%s",t47);fprintf(fp,"%s",t48);fprintf(fp,"%s",t49);
fprintf(fp,"%s",t28);fprintf(fp,"%s",t29);
}
 
if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ){
    j = (dimj -1)/2.0;
    for( zeile=1 ; zeile<=anz_niveaus ; ++zeile )
       for( spalte=1 ; spalte<= VALUE(gi,zeile) ; ++spalte ){
             v = MXSP( ev , (INT)R(entartung,zeile,spalte) );
             first_line = JA;
             for( i=1 ; i<=VRDIM(v) ; ++i ){
                re = RV(v,i);im = IV(v,i);
                re = is_null(re,1.0/1000);
                im = is_null(im,1.0/1000);
                mj = i - j - 1;
                if( re!=0.0 || im!=0.0 ) {
                   if( first_line ){
                       first_line = NEIN;
                       if( re==0.0 )
                          {fprintf(fp,t37,zeile,spalte,im,j,mj);
			  }
                       if( im==0.0 )
                         {fprintf(fp,t36,zeile,spalte,re,j,mj);
			  }
                       if( re!=0.0 && im!=0.0 ){
                          if( im>0 )
                            {fprintf(fp,t31p,zeile,spalte,re,im,j,mj);
			    }
                          else
                            {fprintf(fp,t31m,zeile,spalte,re,ABSD(im),j,mj);
			    }
                       }
                   }else{
                           if( re==0.0 ){
                             if( im>0 )
                                 {fprintf(fp,t39p,im,j,mj);}
                             else
                                 {fprintf(fp,t39m,ABSD(im),j,mj);}
                           }
                           if( im==0.0 ){
                              if( re>0 )
                                  {fprintf(fp,t38p,re,j,mj);}
                              else
                                 {fprintf(fp,t38m,ABSD(re),j,mj);}
                           }
                           if( re!=0.0 && im!=0.0 ){
                              if( im>0 )
                                 {fprintf(fp,t32p,re,im,j,mj);}
                              else if( re>0.0 )
                                 {fprintf(fp,t32m,re,ABSD(im),j,mj);}
                              else
                                 {fprintf(fp,t32mm,ABSD(re),ABSD(im),j,mj);}
                           }
                        }
                }
             }
             fprintf(fp,"%s",t43);
       }

    ewev=fopen_errchk("results/levels.cef","w"); /* output also  in uncommented format with states to be used in mcdiff.in */
  time_t curtime;
  struct tm *loctime;
   fprintf(ewev, "# output file of program so1ion/cfield ");
   curtime=time(NULL); loctime=localtime(&curtime); fprintf(ewev,"%s",asctime(loctime));
    fprintf(ewev,"#J=value {atom-file}  <Jx> <Jy>) <Jz>   \n");
    fprintf(ewev,"#!J= %4.1f sipffile=%s Jx=%6.3f Jy=%6.3f Jz=%6.3f\n",j,INFILE(kristallfeld),
  magnetm(mat_Jy,setup,ewproblem,kristallfeld,B1(iteration),B2(iteration),B3(iteration),temperatur)/gj,
  magnetm(mat_Jz,setup,ewproblem,kristallfeld,B1(iteration),B2(iteration),B3(iteration),temperatur)/gj,
  magnetm(mat_Jx,setup,ewproblem,kristallfeld,B1(iteration),B2(iteration),B3(iteration),temperatur)/gj);
    fprintf(ewev,"#! Eigenvalues= ");

    for( zeile=1 ; zeile<=anz_niveaus ; ++zeile )
       for( spalte=1 ; spalte<= VALUE(gi,zeile) ; ++spalte ){
       fprintf(ewev,"%g ",
       RV(ew,(INT)R(entartung,zeile,1))
       *EINHEITIMP[einheitnr_in].fek*EINHEITIMP[einheitnr_out].fke);
       }
    fprintf(ewev," %s\n",EINHEITIMP[einheitnr_out].einheit);
    fprintf(ewev,"# Eigenvectors [as columns]\n"/*,EINHEITIMP[einheitnr_out].einheit*/);
    fprintf(ewev,"#Real Part\n");
    for( i=1 ; i<=MXDIM(ev) ; ++i ){
     for( zeile=1 ; zeile<=anz_niveaus ; ++zeile )
       for( spalte=1 ; spalte<= VALUE(gi,zeile) ; ++spalte ){
             v = MXSP( ev , (INT)R(entartung,zeile,spalte) );
                re = RV(v,i);im = IV(v,i);
                re = is_null(re,1.0/1000);
                im = is_null(im,1.0/1000);
                             fprintf(ewev,"%+6.5f ",re);
                        }
	     fprintf(ewev,"\n");
             }
    fprintf(ewev,"#Imaginary Part\n");
    for( zeile=1 ; zeile<=anz_niveaus ; ++zeile )
       for( spalte=1 ; spalte<= VALUE(gi,zeile) ; ++spalte ){
//       fprintf(ewev,"0.000 ");
       }
//    fprintf(ewev,"\n");
    for( i=1 ; i<=MXDIM(ev) ; ++i ){
     for( zeile=1 ; zeile<=anz_niveaus ; ++zeile )
       for( spalte=1 ; spalte<= VALUE(gi,zeile) ; ++spalte ){
             v = MXSP( ev , (INT)R(entartung,zeile,spalte) );
                re = RV(v,i);im = IV(v,i);
                re = is_null(re,1.0/1000);
                im = is_null(im,1.0/1000);
                             fprintf(ewev,"%+6.5f ",im);
             }
	     fprintf(ewev,"\n");
       }


    fprintf(fp,"%s",t44);
    fclose(ewev); /*close short output file */  
    printf("Results results/levels.cef written ...\n");
}
 
 
 
if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ){
 
fprintf(fp,"#\n");
fprintf(fp,"#!magnetic moment(mb/f.u.): mx=%6.3f my=%6.3f mz=%6.3f\n", 
  magnetm(mat_Jx,setup,ewproblem,kristallfeld,B1(iteration),B2(iteration),B3(iteration),temperatur),
  magnetm(mat_Jy,setup,ewproblem,kristallfeld,B1(iteration),B2(iteration),B3(iteration),temperatur),
  magnetm(mat_Jz,setup,ewproblem,kristallfeld,B1(iteration),B2(iteration),B3(iteration),temperatur));
fprintf(fp,"#\n");
t01="#-------------------------------------------------------------- \n";
t02="#                                                              \n";
t03="#                  M A T R I X  E L E M E N T                  \n";
t04="#                  S I N G L E  C R Y S T A L                    \n";
t05="#                                                              \n";
t06="#--------------------------------------------------------------\n";
t07="#    Only the marix elements that are non-zero are written     \n";
t08="#                                                              \n";
t09="#--------------------------------------------------------------\n";
t10="#                              2 |            2 |            2 \n";
t11="#    a  <-->  b      |<a|J |b>|  |  |<a|J |b>|  |  |<a|J |b>|  \n";
t12="#                         x      |       y      |       z      \n";
t13="#--------------------------------+--------------+--------------\n";
t14="#   %2d  <--> %2d    %11.6f   |%11.6f   |%11.6f   \n";
t15="#-------------------------------------------------------------- \n";
fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
 
    for( i=1 ; i<= anz_niveaus ; ++i ){
        for( k=1 ; k<= i; ++k ){
            sumx=0.0;sumy=0.0;sumz=0.0;
            for( r=1 ; r<= VALUE(gi,i) ; ++r )
                for( s=1 ; s<= VALUE(gi,k) ; ++s ){
                     ev_ir = MXSP(ev, (INT)R(entartung,i,r) );
                     ev_ks = MXSP(ev, (INT)R(entartung,k,s) );
                     sumx += mat_Jx2( ev_ir , ev_ks , macheps);
                     sumy += mat_Jy2( ev_ir , ev_ks , macheps);
                     sumz += mat_Jz2( ev_ir , ev_ks , macheps);
                }
            sumx = is_null( sumx, 0.000001 );
            sumy = is_null( sumy, 0.000001 );
            sumz = is_null( sumz, 0.000001 );
            if( sumx!=0.0 || sumy!=0.0 || sumz!=0.0){
                fprintf(fp,"%s",t13);
                fprintf(fp,t14,i,k,sumx,sumy,sumz);
            }
        }
    }
    fprintf(fp,"%s",t15);
 
fprintf(fp,"#\n#\n");
t01="#-------------------------------------------------------------- \n";
t02="#                                                              \n";
t03="#                  M A T R I X  E L E M E N T S                \n";
t04="#                    P O L Y C R Y S T A L                     \n";
t05="#                                                              \n";
t06="#--------------------------------------------------------------\n";
t07="#                                                              \n";
t08="#                             ----                2            \n";
t09="# Matrix element       :      >     |<i,r|J |k,s>|             \n";
t10="#                             ----         T                   \n";
t11="#                              r,s                             \n";
t12="#                                                              \n";
t13="#                                                              \n";
t14="# for the transition :      E  ->  E                           \n";
t15="#                              i      k                        \n";
t16="#                                                              \n";
t17="#--------------------------------------------------------------\n";
t18="#                                                              \n";
t19="#                         SUm rule :                           \n";
t20="#                                                              \n";
t21="#             ----                2       2                    \n";
t22="#             >     |<i,r|J |k,s>|  = n  --- J(J+1)            \n";
t23="#             ----         T           i  3                    \n";
t24="#             k,r,s                                            \n";
t25="#-------------------------------------------------------------- \n";
fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);fprintf(fp,"%s",t05);
fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,"%s",t15);
fprintf(fp,"%s",t16);fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);fprintf(fp,"%s",t19);fprintf(fp,"%s",t20);
fprintf(fp,"%s",t21);fprintf(fp,"%s",t22);fprintf(fp,"%s",t23);fprintf(fp,"%s",t24);fprintf(fp,"%s",t25);
}/* end nicht ortho */
 
if( *(FILENAME(kristallfeld)+16) == *(ORTHO+16) ){
t01="#-------------------------------------------------------------- \n";
t02="# Polycrystal transition matrix elements                       \n";
t03="#-------------------------------------------------------------- \n";
fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);
}
    tabelle(fp,anz_niveaus,aJtb2);
 
if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ){
fprintf(fp,"#\n#\n");
t01="#-------------------------------------------------------------- \n";
t02="#               Transition intensities in barn.                \n";
t03="#                                                              \n";
t04="#                                                              \n";
t05="#                                                              \n";
t06="#                                                              \n";
t07="# =                    - E /T     -----                        \n";
t08="# |                   e   i       \\                   2        \n";
t09="# |     = const --------------     >    |<i,r|J |k,s>|         \n";
t10="# |             ----    - E /T    /            T               \n";
t11="# =             >    n e   i      -----                        \n";
t12="#  E -> E       ----  i            r,s                         \n";
t13="#   i    k       i                                             \n";
t14="#                                                              \n";
t15="#                            with                              \n";
t16="#                                                              \n";
t17="#                                                              \n";
t18="#                                -----                         \n";
t19="#                      2     2   \\                    2        \n";
t20="#        |<i,r|J |k,s>|   = ---   >     |<i,r|J |k,s>|         \n";
t21="#               T            3   /             u               \n";
t22="#                                -----                         \n";
t23="#                             u = x,y,z                        \n";
t24="#                                                              \n";
t25="#                                                              \n";
t26="#                             und                              \n";
t27="#                                                              \n";
t28="#                                   1          2               \n";
t29="#                  const  = 4*pi*( --- r   g  )                \n";
t30="#                                   2   0   J                  \n";
t31="#                                                              \n";
t32="#                                      -12                     \n";
t33="#                  r      = - 0.54 * 10    cm                  \n";
t34="#                   0                                          \n";
t35="#                                                              \n";
t36="#--------------------------------------------------------------\n";
t37="#                                                              \n";
t38="#                       1.Sum rule :                           \n";
t39="#                                                 - E /T       \n";
t40="#                                            n   e   i         \n";
t41="#  ----  =            2                       i                \n";
t42="#  >     |         = --- * const * J(J+1) * ----------------   \n";
t43="#  ----  =            3                     ----     - E /T    \n";
t44="#   k     E -> E                            >    n  e   i      \n";
t45="#          i    k                           ----  i            \n";
t46="#                                            i                 \n";
t47="#--------------------------------------------------------------\n";
t48="#                                                              \n";
t49="#                       2. sum rule :                          \n";
t50="#                                                              \n";
t51="#                                                              \n";
t52="#            ----  =            2                              \n";
t53="#            >     |         = --- * const * J(J+1)            \n";
t54="#            ----  =            3                              \n";
t55="#             k,i   E -> E                                     \n";
t56="#                    i    k                                    \n";
t57="#-------------------------------------------------------------- \n";
fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);fprintf(fp,"%s",t05);
fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);fprintf(fp,"%s",t15);
fprintf(fp,"%s",t16);fprintf(fp,"%s",t17);fprintf(fp,"%s",t18);fprintf(fp,"%s",t19);fprintf(fp,"%s",t20);
fprintf(fp,"%s",t21);fprintf(fp,"%s",t22);fprintf(fp,"%s",t23);fprintf(fp,"%s",t24);fprintf(fp,"%s",t25);
fprintf(fp,"%s",t26);fprintf(fp,"%s",t27);fprintf(fp,"%s",t28);fprintf(fp,"%s",t29);fprintf(fp,"%s",t30);
fprintf(fp,"%s",t31);fprintf(fp,"%s",t32);fprintf(fp,"%s",t33);fprintf(fp,"%s",t34);fprintf(fp,"%s",t35);
fprintf(fp,"%s",t36);fprintf(fp,"%s",t37);fprintf(fp,"%s",t38);fprintf(fp,"%s",t39);fprintf(fp,"%s",t40);
fprintf(fp,"%s",t41);fprintf(fp,"%s",t42);fprintf(fp,"%s",t43);fprintf(fp,"%s",t44);fprintf(fp,"%s",t45);
fprintf(fp,"%s",t46);fprintf(fp,"%s",t47);fprintf(fp,"%s",t48);fprintf(fp,"%s",t49);fprintf(fp,"%s",t50);
fprintf(fp,"%s",t51);fprintf(fp,"%s",t52);fprintf(fp,"%s",t53);fprintf(fp,"%s",t54);fprintf(fp,"%s",t55);
fprintf(fp,"%s",t56);fprintf(fp,"%s",t57);
}/* endif nicht ortho */
 
t01="#-------------------------------------------------------------- \n";
t02="#!Temperature of the sample               T= %7.2f Kelvin         \n";
    fprintf(fp,"%s",t01);fprintf(fp,t02,temperatur);
 
 
 zu_summe    = zustandssumme( einheitnr_in , ew , temperatur );
 gj          = IONENIMP[ ionennr ].gj;
 gesamte_intensitaet = 2.0/3.0*const*gj*gj* ((DOUBLE)(dimj*dimj-1)/4);
 
t01="#-------------------------------------------------------------- \n";
t02="#!parition function            Z  = %6.2f                 \n";
t03="#-------------------------------------------------------------- \n";
t04="#!Total_magnetic_scattering_intensity = %6.2f barn            \n";
t05="#-------------------------------------------------------------- \n";
fprintf(fp,"%s",t01);fprintf(fp,t02,zu_summe);fprintf(fp,"%s",t03);
fprintf(fp,t04,gesamte_intensitaet);fprintf(fp,"%s",t05);
 
if( *(FILENAME(kristallfeld)+16) == *(ORTHO+16) ){
t01="#-------------------------------------------------------------- \n";
t02="# Transition intensities Poly crystal                         \n";
t03="#-------------------------------------------------------------- \n";
fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);
}
 
 /* Intensitaet noch mit gi_ze*const*gj**2 * exp( -Ei/T)/Z multipizieren */
    faktor      = EINHEITIMP[ einheitnr_in].fek;
    for( ze=1 ; ze<=anz_niveaus ; ++ze ){
      energie = RV( ew , (INT)R(entartung,ze,1) );
      for( sp=1 ; sp<=anz_niveaus ; ++sp )
        {R(aJtb2,ze,sp) *= const*gj*gj
                           *exp_(-energie*faktor/temperatur)/zu_summe;
        }
    }
    tabelle(fp,anz_niveaus,aJtb2);
    if( anz_niveaus > 1)
        mittlung(fp,anz_niveaus,aJtb2,ew,entartung,einheitnr_in,
                    einheitnr_out,einheit_out);
t01="#-------------------------------------------------------------- \n";
t02="# Transition Energy (%6s) vs Intensity (barn)                  \n";
t03="#-------------------------------------------------------------- }\n";
fprintf(fp,"%s",t01);fprintf(fp,t02,einheit_out);fprintf(fp,"%s",t03);

/* Intensitaet noch als tabelle ausgeben */
    for( ze=1 ; ze<=anz_niveaus ; ++ze ){
      for( sp=1 ; sp<=anz_niveaus ; ++sp )
        {
      energie = (RV( ew , (INT)R(entartung,sp,1) )-RV( ew , (INT)R(entartung,ze,1) ))
                 *EINHEITIMP[einheitnr_in].fek*EINHEITIMP[einheitnr_out].fke;
         fprintf(fp,"%6.2f %6.2f\n",energie,R(aJtb2,ze,sp));
        }
    }

    free_mx(aJtb2);


/*    for( zeile=1 ; zeile<=anz_niveaus ; ++zeile){
       fprintf(fp,dummy,zeile,
       RV(ew,(INT)R(entartung,zeile,1))
       *EINHEITIMP[einheitnr_in].fek*EINHEITIMP[einheitnr_out].fke,
                  VALUE(gi,zeile) );
    }*/
 
if( *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ){
    fclose(fp);
}/* endif nicht ortho */
 
 
 if(IS_SUSZEPT(kristallfeld)==JA)
       raus_suszept(setup,iteration,ewproblem,kristallfeld);
 
 
 if(IS_KPOLY(kristallfeld)==JA){
       anf_feld   = ANFANG_FELD(kristallfeld);
       end_feld   = END_FELD(   kristallfeld);
       temp       = TEMP(       kristallfeld);
    raus_kpoly(setup,ewproblem,kristallfeld,anf_feld,end_feld,
               temp,ionname);
 }
 
 if(IS_ORTHO(kristallfeld)==JA){
    fit_ortho(setup,ewproblem,kristallfeld);
 }

  if( IS_MAGNETM(kristallfeld)==JA && IS_KPOLY(kristallfeld)==NEIN&&
     IS_ORTHO(  kristallfeld)==NEIN ){
       anf_feld   = ANFANG_FELD(kristallfeld);
       end_feld   = END_FELD(   kristallfeld);
       temp       = TEMP(       kristallfeld);
    raus_magnetm(setup,ewproblem,kristallfeld,anf_feld,end_feld,
                 temp,ionname);
 }
}

/*----------------------------------------------------------------------------
                                   raus_suszept()
-----------------------------------------------------------------------------*/
void raus_suszept(setup,iteration,ewproblem,kristallfeld)
    SETUP     *setup;
    ITERATION *iteration;
    EWPROBLEM *ewproblem;
    KRISTALLFELD *kristallfeld;
{
  FILE   *fp=0;
  INT    einheitnr_in,datensatz_nr=-1,anz_daten=0/*,dummy*/,i,anz_temp;
  INT    anzdatnr,datnr;
  INT    lesethetafile,anz,is_pbekannt();
  DOUBLE gj/*,t*/,macheps;
  DOUBLE suszept(),mat_Jx2(),mat_Jy2(),mat_Jz2(),zwischen;
  DOUBLE *x_s,*y_s,*z_s,*temp,temp_step;
  DOUBLE anf_temp,end_temp,lambda,theta,*xx,*ff;
  DOUBLE ttheta;
  CHAR   *t01,*t02,*t03;
  CHAR   *namethetafile;
 
  anf_temp   = ANFANG_TEMPERATUR(kristallfeld);
  end_temp   = END_TEMPERATUR(   kristallfeld);
  lambda     = LAMBDA(           kristallfeld);
  theta      = THETA(            kristallfeld);
  namethetafile = NAMETHETAFILE( kristallfeld);
  lesethetafile = LESETHETAFILE( kristallfeld);

  if(NUMMERIERUNG(setup)==JA)  datensatz_nr = 0;
 
  macheps      = ewproblem->eps_machine;
  einheitnr_in = EINHEITNRIN( iteration );
  gj           = GJ(iteration );
 
  if(lesethetafile==JA){
    NRPOSA(iteration) = i= 1;
    NRPOSE(iteration) = 1;
    FORMATPA(iteration)=7;   /* 7 mal 11 */
    FORMATPB(iteration)=11;
    FORMATNA(iteration)=2;   /* 2 mal 12 */
    FORMATNB(iteration)=12;
    anz = 200;
    xx  = DOUBLE_ALLOC(anz);
    ff  = DOUBLE_ALLOC(anz);
    lesedaten(fp,xx,ff,&anz,namethetafile,NUMMERIERUNG(setup),
              &i,1,is_pbekannt,&datnr,iteration,&anzdatnr);
    fclose(fp);
  }
 
  anz_temp     = ANZ_DATENPUNKTE(setup);  /*sollte >=2 sein */
  temp_step    = (end_temp - anf_temp)/(DOUBLE)(anz_temp-1);
 
 
  temp     = DOUBLE_ALLOC(anz_temp);
  for( i=1; i<= anz_temp ; ++i )
       VALUE(temp,i) = temp_step*(DOUBLE)(i-1)+anf_temp;
 
 
  t01="#%5d%5d %s\n";
  t02="%10.5e %10.5e\n";
 
  printf("calculating inverse Susceptibility in [100] ... \n");
  x_s = DOUBLE_ALLOC( anz_temp );
  for( i=1; i<= anz_temp ; ++i )
   VALUE(x_s,i)=suszept(mat_Jx2,ewproblem,einheitnr_in,VALUE(temp,i),gj);
 
  printf("calculating inverse Susceptibility %s ... \n",SUSZEPT);
  fp           = fopen_errchk( SUSZEPT, "w");
  anz_daten    = 0;
  for( i=1; i<= anz_temp ; ++i )
     if( !is_equal(VALUE(x_s,i),0.0,macheps) && VALUE(x_s,i) >0.0 )
                 ++anz_daten;
  if( anz_daten != 0 ){
    if(NUMMERIERUNG(setup)==JA) --datensatz_nr;
    t03="Paramagnetic inverse Susceptibility in [100] (in mol/emu)";
    fprintf(fp,t01,datensatz_nr,anz_daten,t03);
    for( i=1; i<= anz_temp; ++i )
       if( !is_equal(VALUE(x_s,i),0.0,macheps) && VALUE(x_s,i) >0.0 ){
          ttheta = theta;
          if(lesethetafile==JA) ttheta=spline(anz,xx,ff,VALUE(temp,i),macheps);
          fprintf(fp,t02,VALUE(temp,i)+ttheta,1.0/VALUE(x_s,i)-lambda);
       }
  }
  fclose(fp);
 
 
  printf("calculating inverse Susceptibility in [010] ...\n");
  y_s = DOUBLE_ALLOC( anz_temp );
  for( i=1; i<= anz_temp ; ++i )
   VALUE(y_s,i)=suszept(mat_Jy2,ewproblem,einheitnr_in,VALUE(temp,i),gj);
 
  printf("calculating inverse Susceptibility  %s ...... \n",SUSZEPT);
  fp           = fopen_errchk( SUSZEPT, "a");
  anz_daten    = 0;
  for( i=1; i<= anz_temp ; ++i )
     if( !is_equal(VALUE(y_s,i),0.0,macheps) && VALUE(y_s,i) >0.0 )
                 ++anz_daten;
  if( anz_daten != 0 ){
    if(NUMMERIERUNG(setup)==JA) --datensatz_nr;
    t03="Paramagnetic inverse Susceptibility in [010] (in mol/emu)";
    fprintf(fp,t01,datensatz_nr,anz_daten,t03);
    for( i=1; i<= anz_temp; ++i )
       if( !is_equal(VALUE(y_s,i),0.0,macheps) && VALUE(y_s,i) >0.0 ){
          ttheta = theta;
          if(lesethetafile==JA) ttheta=spline(anz,xx,ff,VALUE(temp,i),macheps);
          fprintf(fp,t02,VALUE(temp,i)+ttheta,1.0/VALUE(y_s,i)-lambda);
       }
  }
  fclose(fp);
 
 
 
  printf("calculating inverse Susceptibility in [001] ... \n");
  z_s = DOUBLE_ALLOC( anz_temp );
  for( i=1; i<= anz_temp ; ++i )
   VALUE(z_s,i)=suszept(mat_Jz2,ewproblem,einheitnr_in,VALUE(temp,i),gj);
 
  printf("calculating inverse Susceptibility  %s ... \n",SUSZEPT);
  fp           = fopen_errchk( SUSZEPT, "a");
  anz_daten    = 0;
  for( i=1; i<= anz_temp ; ++i )
     if( !is_equal(VALUE(z_s,i),0.0,macheps) && VALUE(z_s,i) >0.0 )
                 ++anz_daten;
  if( anz_daten != 0 ){
    if(NUMMERIERUNG(setup)==JA) --datensatz_nr;
    t03="Paramagnetic inverse Susceptibility in [001] (in mol/emu)";
    fprintf(fp,t01,datensatz_nr,anz_daten,t03);
    for( i=1; i<= anz_temp; ++i )
       if( !is_equal(VALUE(z_s,i),0.0,macheps) && VALUE(z_s,i) >0.0 ){
          ttheta = theta;
          if(lesethetafile==JA) ttheta=spline(anz,xx,ff,VALUE(temp,i),macheps);
          fprintf(fp,t02,VALUE(temp,i)+ttheta,1.0/VALUE(z_s,i)-lambda);
       }
  }
  fclose(fp);
 
 
 
printf("calculating inverse Susceptibility of a Polycrystal ... \n");
anz_daten    = 0;
for( i=1; i<= anz_temp ; ++i ){
    if( !is_equal(VALUE(x_s,i),0.0,macheps) && VALUE(x_s,i) >0.0 )
      if( !is_equal(VALUE(y_s,i),0.0,macheps) && VALUE(y_s,i) >0.0 )
        if( !is_equal(VALUE(z_s,i),0.0,macheps) && VALUE(z_s,i) >0.0 )
           ++anz_daten;
}
printf("calculating inverse Susceptibility  %s ... \n",SUSZEPT);
fp           = fopen_errchk( SUSZEPT, "a");
if( anz_daten != 0 ){
 if(NUMMERIERUNG(setup)==JA) --datensatz_nr;
 t03="Param. inv. Suszeptib. fuer ein Polykristall(mol/emu)";
 fprintf(fp,t01,datensatz_nr,anz_daten,t03);
 for( i=1; i<= anz_temp; ++i ){
    if( !is_equal(VALUE(x_s,i),0.0,macheps) && VALUE(x_s,i) >0.0 )
      if( !is_equal(VALUE(y_s,i),0.0,macheps) && VALUE(y_s,i) >0.0 )
        if( !is_equal(VALUE(z_s,i),0.0,macheps) && VALUE(z_s,i) >0.0 ){
          zwischen  = 1.0/( 1.0/VALUE(x_s,i)-lambda);
          zwischen += 1.0/( 1.0/VALUE(y_s,i)-lambda);
          zwischen += 1.0/( 1.0/VALUE(z_s,i)-lambda);
          fprintf(fp,t02,VALUE(temp,i),3.0/zwischen);
        }
  }
}
 
 
 t01="#\n#Alle Suszeptibilitaeten sind Null. Daher wurde %s nicht angelegt.\n";
 if( datensatz_nr<= 0 && anz_daten == 0)
     printf(t01,SUSZEPT);
 
  raus_kommentar(fp,anf_temp,end_temp,lambda,theta,lesethetafile,
                 namethetafile);
 
  fclose(fp);

  free_(x_s);
  free_(y_s);
  free_(z_s);
  free_(temp);
  if(lesethetafile==JA){
     free_(xx);
     free_(ff);
  }
 
}
/*----------------------------------------------------------------------------
                                   raus_kpoly()
-----------------------------------------------------------------------------*/
void raus_kpoly( setup,ewproblem,kristallfeld,anf_feld,end_feld,
                  temp,ionname)
    SETUP        *setup;
    EWPROBLEM    *ewproblem;
    KRISTALLFELD *kristallfeld;
    DOUBLE       anf_feld,end_feld,temp;
    CHAR         *ionname;
{
  FILE   *fp;
  INT    datensatz_nr=-1,anz_daten=0/*,dummy*/,i,anz_feld/*,r1,r2,r3*/;
  DOUBLE /*t,*/macheps;
  DOUBLE magnetm()/*,b_norm*/,sqrt();
  DOUBLE *x_s/*,*y_s,*z_s*/,*feld,feld_step,b1,b2,b3,mx,my,mz;
  CHAR   *t01,*t02/*,*t03,*t04*/,*t05;
  ITERATION *iteration;
  COMHES *comhes;
  MATRIX *matrix;
  BRUCH  *is_rational()/*,*z1,*z2,*z3*/;
  KOMPLEX *mat_Jx(),*mat_Jy(),*mat_Jz();
  DOUBLE sqrt();
  LONG /*hauptnenner,*/ggt_l();
 
  if(NUMMERIERUNG(setup)==JA)  datensatz_nr = 0;
 
  macheps   = ewproblem->eps_machine;
  comhes    = ewproblem->comhes;
  matrix    = comhes   ->matrix;
  ewproblem = set_ewproblem( setup,ewproblem,matrix,macheps,NOSPACE);
 
  macheps      = ewproblem->eps_machine;
  anz_feld     = ANZ_DATENPUNKTE(setup);  /*sollte >=2 sein */
  feld_step    = (end_feld - anf_feld)/(DOUBLE)(anz_feld-1);
  iteration    = ITERATION(kristallfeld);
  feld         = DOUBLE_ALLOC(anz_feld);
  t01          = "%5d%5d %s %4s %s \n";
  t05          = "%5d%5d Magnetic Moment von %4s in the polycrystal average\n";
  t02          = "%10.5e %10.5e\n";
 
  for( i=1; i<= anz_feld ; ++i )
       VALUE(feld,i) = feld_step*(DOUBLE)(i-1)+anf_feld;
 
  printf("calculating polycrystalline average (cubic)... \n"/*,MAGNETM*/);
  x_s = DOUBLE_ALLOC( anz_feld );
  for( i=1; i<= anz_feld ; ++i ){
       ++anz_daten;
       b1  = 1.0;
       b2  = 0.0;
       b3  = 0.0;
       mx  = magnetm(mat_Jx,setup,ewproblem,kristallfeld,
          (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
       my  = magnetm(mat_Jy,setup,ewproblem,kristallfeld,
          (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
       mz  = magnetm(mat_Jz,setup,ewproblem,kristallfeld,
          (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
       VALUE(x_s,i) = 6.0*sqrt(mx*mx + my*my + mz*mz);
       b1  = 1.0;
       b2  = 1.0;
       b3  = 0.0;
       mx  = magnetm(mat_Jx,setup,ewproblem,kristallfeld,
          (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
       my  = magnetm(mat_Jy,setup,ewproblem,kristallfeld,
          (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
       mz  = magnetm(mat_Jz,setup,ewproblem,kristallfeld,
          (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
       VALUE(x_s,i) += 12.0*sqrt(mx*mx + my*my + mz*mz);
       b1  = 1.0;
       b2  = 1.0;
       b3  = 1.0;
       mx  = magnetm(mat_Jx,setup,ewproblem,kristallfeld,
          (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
       my  = magnetm(mat_Jy,setup,ewproblem,kristallfeld,
          (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
       mz  = magnetm(mat_Jz,setup,ewproblem,kristallfeld,
          (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
      VALUE(x_s,i) +=  8.0*sqrt(mx*mx + my*my + mz*mz);
      VALUE(x_s,i) /=  26.0;
}
 
 
      printf("Berechnetes Moment auf %s rausschreiben ... \n",MAGNETM);
      fp = fopen_errchk( KPOLY, "w");
      if(NUMMERIERUNG(setup)==JA) --datensatz_nr;
      fprintf(fp,t05,datensatz_nr,anz_feld,ionname);
      for( i=1; i<= anz_feld; ++i )
           fprintf(fp,t02,VALUE(feld,i),VALUE(x_s,i) );
 
     t01="\n%s wurde nicht angelegt. Magnetnetisches Moment ist Null.\n";
      if( datensatz_nr<=0 && anz_daten == 0)
           printf(t01,KPOLY);
 
     free_(x_s);
     free_(feld);
 
     raus_kkommentar(fp,iteration,anf_feld,end_feld,temp,ionname);
 
     fclose(fp);
}
/*----------------------------------------------------------------------------
                                   richtung()
-----------------------------------------------------------------------------*/
void richtung( iteration )
  ITERATION *iteration;
{
  BRUCH  *is_rational(),*z1,*z2,*z3;
  LONG   hauptnenner,ggt_l();
  DOUBLE b1,b2,b3,b_norm,sqrt();
 
  b1        = B1(iteration);
  b2        = B2(iteration);
  b3        = B3(iteration);
  b_norm = sqrt( b1*b1 + b2*b2 + b3*b3 );
  if( !null(b_norm,MACHEPSFACT*accuracy())  ){
        b1 /= b_norm; b2 /= b_norm; b3 /= b_norm; }
  else{ b1  =  0.0;   b2  =  0.0;   b3  =  1.0;   }
 
 /* Richtung von B ausrechnen  */
    z1 = is_rational(b1*b1);  hauptnenner = ggt_l(ZAEHLER(z1),NENNER(z1));
    ZAEHLER(z1) /= hauptnenner; NENNER( z1) /= hauptnenner;
    z2 = is_rational(b2*b2);  hauptnenner = ggt_l(ZAEHLER(z2),NENNER(z2));
    ZAEHLER(z2) /= hauptnenner; NENNER( z2) /= hauptnenner;
    z3 = is_rational(b3*b3);  hauptnenner = ggt_l(ZAEHLER(z3),NENNER(z3));
    ZAEHLER(z3) /= hauptnenner; NENNER( z3) /= hauptnenner;
 
    hauptnenner = 1;
    if( NENNER(z1) > hauptnenner )  hauptnenner = NENNER(z1);
    if( NENNER(z2) > hauptnenner )  hauptnenner = NENNER(z2);
    if( NENNER(z3) > hauptnenner )  hauptnenner = NENNER(z3);
 
    ZAEHLER(z1) *= (hauptnenner/NENNER(z1));
    ZAEHLER(z2) *= (hauptnenner/NENNER(z2));
    ZAEHLER(z2) *= (hauptnenner/NENNER(z2));
 
    R1(iteration) = DSIGN(b1)*(LONG)sqrt( (DOUBLE)ZAEHLER(z1) );
    R2(iteration) = DSIGN(b2)*(LONG)sqrt( (DOUBLE)ZAEHLER(z2) );
    R3(iteration) = DSIGN(b3)*(LONG)sqrt( (DOUBLE)ZAEHLER(z3) );
    free_(z1);free_(z2);free_(z3);
}
/*----------------------------------------------------------------------------
                                   raus_magnetm()
-----------------------------------------------------------------------------*/
void raus_magnetm( setup,ewproblem,kristallfeld,anf_feld,end_feld,
                  temp,ionname)
    SETUP        *setup;
    EWPROBLEM    *ewproblem;
    KRISTALLFELD *kristallfeld;
    DOUBLE       anf_feld,end_feld,temp;
    CHAR         *ionname;
{
  FILE   *fp;
  INT    datensatz_nr=-1,anz_daten=0/*,dummy*/,i,anz_feld,r1,r2,r3;
  DOUBLE /*t,*/macheps;
  KOMPLEX *mat_Jx(),*mat_Jy(),*mat_Jz();
  DOUBLE magnetm(),b_norm,sqrt();
  DOUBLE *x_s,*y_s,*z_s,*feld,feld_step,b1,b2,b3;
  DOUBLE b1mol,b2mol,b3mol,bmol_norm;
  CHAR   *t01,*t02/*,*t03*/,*t04,*t05;
  ITERATION *iteration;
  COMHES *comhes;
  MATRIX *matrix;
  DOUBLE sqrt();
 
  if(NUMMERIERUNG(setup)==JA)  datensatz_nr = 0;
 
  macheps   = ewproblem->eps_machine;
  comhes    = ewproblem->comhes;
  matrix    = comhes   ->matrix;
  ewproblem = set_ewproblem( setup,ewproblem,matrix,macheps,NOSPACE);
 
  macheps   = ewproblem->eps_machine;
  anz_feld  = ANZ_DATENPUNKTE(setup);  /*sollte >=2 sein */
  feld_step = (end_feld - anf_feld)/(DOUBLE)(anz_feld-1);
  iteration = ITERATION(kristallfeld);
  feld      = DOUBLE_ALLOC(anz_feld);
  t01       = "%5d%5d %s %4s %s \n";
  t04       = "#%5d%5d absolute value of magnetic moment of %4s  for B in [%2g,%2g,%2g]\n";
  t05       = "#%5d%5d magnetic moment of %4s in [%2d,%2d,%2d] for B in [%2g,%2g,%2g]\n";
  t02       = "%10.5e %10.5e\n";
  b1        = B1(iteration);
  b2        = B2(iteration);
  b3        = B3(iteration);
  b1mol     = B1MOL(iteration);
  b2mol     = B2MOL(iteration);
  b3mol     = B3MOL(iteration);
 
  b_norm    = sqrt( b1   *b1    + b2   *b2    + b3   *b3 );
  bmol_norm = sqrt( b1mol*b1mol + b2mol*b2mol + b3mol*b3mol );
  if( !null(b_norm,macheps)  ){
        b1 /= b_norm; b2 /= b_norm; b3 /= b_norm; }
  else{
           b1=0.0;
           b2=0.0;
           b3=1.0;
  }
 
  richtung( iteration );  /* r1,r2,r3 setzen */
  r1 = R1(iteration);
  r2 = R2(iteration);
  r3 = R3(iteration);
 
    if( null(b_norm,macheps) && null(bmol_norm,macheps) ){
        b1 = 1.0; b2 = 0.0; b3 = 0.0;
        r1 = 1.0  ; r2 = 0.0  ; r3 = 0.0  ;
    }
    for( i=1; i<= anz_feld ; ++i )
          VALUE(feld,i) = feld_step*(DOUBLE)(i-1)+anf_feld;
    printf("Calculating magnetic moment in [100] for B in [%2g,%2g,%2g] ... \n", b1,b2,b3);
    x_s = DOUBLE_ALLOC( anz_feld );
    for( i=1; i<= anz_feld ; ++i ){
         ++anz_daten;
         VALUE(x_s,i)=magnetm(mat_Jx,setup,ewproblem,kristallfeld,
                              (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
    }
    printf("Berechnetes Moment auf %s rausschreiben ... \n",MAGNETM);
    fp = fopen_errchk( MAGNETM, "w");
    if(NUMMERIERUNG(setup)==JA) --datensatz_nr;
    fprintf(fp,t05,datensatz_nr,anz_feld,ionname,1,0,0,b1,b2,b3);
    for( i=1; i<= anz_feld; ++i )
         fprintf(fp,t02,VALUE(feld,i),VALUE(x_s,i) );
 
 
    if( null(b_norm,macheps) &&  null(bmol_norm,macheps) ){
        b1 = 0.0; b2 = 1.0; b3 = 0.0;
        r1 = 0  ; r2 = 1  ; r3 = 0  ;
    }
    printf("calculating magnetic moment in [010] for B in [%2g,%2g,%2g] ... \n", b1,b2,b3);
    y_s = DOUBLE_ALLOC( anz_feld );
    for( i=1; i<= anz_feld ; ++i ){
         ++anz_daten;
         VALUE(y_s,i)=magnetm(mat_Jy,setup,ewproblem,kristallfeld,
                             (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
    }
    printf("Writing results to %s ... \n",MAGNETM);
    if(NUMMERIERUNG(setup)==JA) --datensatz_nr;
    fprintf(fp,t05,datensatz_nr,anz_feld,ionname,0,1,0,b1,b2,b3);
    for( i=1; i<= anz_feld; ++i )
         fprintf(fp,t02,VALUE(feld,i),VALUE(y_s,i) );
 
 
    if( null(b_norm,macheps) &&  null(bmol_norm,macheps) ){
        b1 = 0.0; b2 = 0.0; b3 = 1.0;
        r1 = 0  ; r2 = 0  ; r3 = 1  ;
    }
    printf("calculating magnetic moment in [001] for B in [%2g,%2g,%2g] ... \n", b1,b2,b3);
    z_s = DOUBLE_ALLOC( anz_feld );
    for( i=1; i<= anz_feld ; ++i ){
         ++anz_daten;
         VALUE(z_s,i)=magnetm(mat_Jz,setup,ewproblem,kristallfeld,
                              (b1*VALUE(feld,i)),(b2*VALUE(feld,i)),(b3*VALUE(feld,i)),temp);
    }
    printf("Berechnetes Moment auf %s rausschreiben ... \n",MAGNETM);
    if(NUMMERIERUNG(setup)==JA) --datensatz_nr;
    fprintf(fp,t05,datensatz_nr,anz_feld,ionname,0,0,1,b1,b2,b3);
    for( i=1; i<= anz_feld; ++i )
         fprintf(fp,t02,VALUE(feld,i),VALUE(z_s,i) );

    if( !null(b_norm,macheps) ||  !null(bmol_norm,macheps) )
    {    
     printf("Berechne Betrag des magnetisches Moments fuer B in [%2g,%2g,%2g] ... \n", b1,b2,b3);
     printf("writing results to %s ... \n",MAGNETM);
     if(NUMMERIERUNG(setup)==JA) --datensatz_nr;
     fprintf(fp,t04,datensatz_nr,anz_feld,ionname,b1,b2,b3);
     for( i=1; i<= anz_feld; ++i )
         fprintf(fp,t02,VALUE(feld,i),sqrt( VALUE(x_s,i)*VALUE(x_s,i)+VALUE(y_s,i)*VALUE(y_s,i)+VALUE(z_s,i)*VALUE(z_s,i) ) );
   }
 
    t01="\n%s was not stored. Magnetic Moment is zero.\n";
    if( datensatz_nr <= 0 && anz_daten == 0)  printf(t01,MAGNETM);
 
    free_(x_s);
    free_(y_s);
    free_(z_s);
    free_(feld);
 
    raus_mkommentar(fp,iteration,anf_feld,end_feld,temp,r1,r2,r3,ionname);
 
    fclose(fp);
}
/*----------------------------------------------------------------------------
                                magnetm()
-----------------------------------------------------------------------------*/
/*          ---                                         */
/* (-) removed from this formula MR 28.9.08 */
/* m  =  g  >    w   <ir|J |ir>     c=x,y,z             */
/*  c     j ---   i       c                             */
/*         i,r                                          */
/*                                                         */
/*                                                         */
/*  w  =  exp( -E / T )  /  Z(T)                           */
/*   i           i                                         */
/*                                                         */
 
 
DOUBLE magnetm(mat_Ji,setup,ewproblem,kristallfeld,Bx,By,Bz,t)
    KOMPLEX   *(*mat_Ji)();
    SETUP     *setup;
    EWPROBLEM *ewproblem;
    KRISTALLFELD *kristallfeld;
    DOUBLE    Bx,By,Bz; /* angelegetes aeusseres feld*/
    DOUBLE    t; /* angelegete  Temperatur    */
{
    DOUBLE    Bxmol,Bymol,Bzmol; /* Molekularfeld   */
    ITERATION *iteration;
    INT    i/*,k*/,r/*,s*/,anz_niveaus,*gi;
    VEKTOR *ev_ir;
    VEKTOR *ew;
    MATRIX *ev,*bmag;
    MATRIX *entartung;
    DOUBLE faktor,exp_(),macheps,sumr=0.0,sumi=0.0,zusumme=0.0;
    DOUBLE ew_i,wi,gj,myB;
    KOMPLEX *mat;
    INT     einheitnr_in;
    CHAR  art;
 
 
    iteration    = ITERATION( kristallfeld );
    einheitnr_in = EINHEITNRIN( iteration );
    bmag         = HMAG( iteration );
    gj           = GJ(iteration );
    myB          = EINHEITIMP[ einheitnr_in ].myB;
    art          = EINGABEPARAMETERART(kristallfeld);
    Bxmol        = B1MOL(iteration);
    Bymol        = B2MOL(iteration);
    Bzmol        = B3MOL(iteration);
 
  
    HMAG(iteration) = calc_iBmag( bmag,gj,myB,Bx,By,Bz,Bxmol,Bymol,Bzmol);

    ewproblem       = solve(setup,ewproblem,NEIN,kristallfeld,art);
/* here changed to NOSPACE to NEIN MR okt 2002 - because NOSPACE leads
to overwriting the eigenvectors and eigenvalue (intended to savememory) -
this has the consequence that the entartung gets different in different
fields and the diagonalize routine cannot deal (in the NOSPACE=overwrite
modus) with repeated calls leading to different entartungen */
    macheps   = ewproblem->eps_machine;
    entartung = ewproblem->entartung;
    gi        = ewproblem->gi;
    ev        = ewproblem->eigenvektoren;
    ew        = ewproblem->eigenwerte;
 
 
    faktor        = EINHEITIMP[ einheitnr_in].fek;
    anz_niveaus   = ANZ_ZE(entartung);
 
    for( i=1 ; i<=VRDIM(ew) ; ++i )
        zusumme += exp_( - faktor*RV(ew,i) / t );
 
    for( i=1 ; i<= anz_niveaus ; ++i )
                for( r=1 ; r<= VALUE(gi,i); ++r ){
                   ew_i  = RV(ew,(INT)R(entartung,i,r))*faktor / t;
                   wi    = exp_( - ew_i )/zusumme;
                   ev_ir = MXSP(ev, (INT)R(entartung,i,r) );
                   mat   = (*mat_Ji)(ev_ir,ev_ir);
                     sumr += is_null(wi*RT(mat),macheps);
                     sumi += is_null(wi*IT(mat),macheps);
                   free_(mat);
                }
 
 
    if( sumi!= 0.0 )
      printf("\n imaginary part m_c = %e is not zero!\n",sumi);
 
 /* (-) removed from this formula MR 28.9.08 */
    return(gj*sumr);
 
}
/*----------------------------------------------------------------------------
                                   raus_mkommentar()
-----------------------------------------------------------------------------*/
void raus_mkommentar(fp,iteration,anf_feld,end_feld,temp,b1,b2,b3,ionname)
  FILE       *fp;
  ITERATION  *iteration;
  DOUBLE     anf_feld,end_feld,temp;
  INT        b1,b2,b3;
  CHAR       *ionname;
{
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11,*t12,*t13,*t14,*t15,*t16,*t17,*t18,*t19,*t20;
    CHAR *t21,*t22,*t23,*t24,*t25,*t26,*t27,*t28,*t29,*t30;
    CHAR *t31,*t32,*t33,*t34,*t35,*t36,*t37,*t38/*,*t39,*t40;
    CHAR *t41,*t42,*t43,*t44,*t45*/,*t46,*t47,*t48,*t49,*t50;
    CHAR *t51,*t52,*t53,*t54,*t55,*t56,*t57,*t58,*t59,*t60;
    CHAR *t61/*,*t62,*t63,*t64,*t65,*t66,*t67,*t68,*t69,*t70;
    CHAR *t71,*t72,*t73,*t74,*t75,*t76,*t77,*t78,*t79,*t80*/;
 
 
 
t01="# -------------------------------------------------------------- \n";
t02="#|    comment about the calculation of the magnetic moment      \n";
t03="#|    of the Rare Earth Ion at the Temperature temp and in      \n";
t04="#|    Field range : anf_feld until end_feld                     \n";
t05="# -------------------------------------------------------------- \n";
t06="#|                                                              \n";
t07="#|  The magnetic Moment  m    is                                \n";
t08="#|                       i                                      \n";
t09="#|                                                              \n";
t10="#|                  ----                                        \n";
/* minus removed from this formula MR 24.9.08 */
t11="#|  m   =    g  *   >     w   <n| J |n>     with   i = a,b,c    \n";
t12="#|   i        J     ----   n       i                            \n";
t13="#|                    n                                         \n";
t14="#|                                                              \n";
t15="#|                                                              \n";
t16="#|                             and                              \n";
t17="#|                                                              \n";
t18="#|              -E /T                                           \n";
t19="#|             e  n                                             \n";
t20="#|  w   =  ------------      ,   Z(T) = partition sum           \n";
t21="#|   n         Z(T)                                             \n";
t22="#|                                                              \n";
t23="#|                                                              \n";
t24="#|  calculated.                                                  \n";
t25="#|                                                              \n";
t26="#|                                                              \n";
t27="#|  where  :                                                    \n";
t28="#|                                                              \n";
t29="#|  a,b,c   :   the thress crystal directions                   \n";
t30="#|                                                              \n";
t31="#|  E       :   crystal Field energy to the state  |n>          \n";
t32="#|   n                                                          \n";
t33="#|                                                              \n";
t34="#|  g       :   Lande factor                                     \n";
t35="#|   J                                                          \n";
t36="#|                                                              \n";
t37="#|                                                              \n";
t38="# -------------------------------------------------------------- \n";
t46="#! start_field       = %6.2f Tesla                             \n"  ;
t47="#! End_field           = %6.2f Tesla                             \n"  ;
t48="#|                                                              \n";
t49="#! temp              = %6.2f Kelvin                            \n";
t50="#|                                                              \n";
t51="#|                                                              \n";
t52="#| The magnetic %4s - moment in der Field direction             \n"  ;
t53="#|                                                              \n";
t54="#|         B  =  (%2d,%2d,%2d)                                     \n";
t55="#|         -                                                    \n";
t56="#| calculated.                                                   \n";
t57="#|                                                              \n";
t58="#! Bx_mol            = %6.2f Tesla                             \n"  ;
t59="#! By_mol            = %6.2f Tesla                             \n"  ;
t60="#! Bz_mol            = %6.2f Tesla                             \n"  ;
t61="# -------------------------------------------------------------- \n";
 
fprintf(fp,"%s",t01);
fprintf(fp,"%s",t02);
fprintf(fp,"%s",t03);
fprintf(fp,"%s",t04);
fprintf(fp,"%s",t05);
fprintf(fp,"%s",t06);
fprintf(fp,"%s",t07);
fprintf(fp,"%s",t08);
fprintf(fp,"%s",t09);
fprintf(fp,"%s",t10);
 
fprintf(fp,"%s",t11);
fprintf(fp,"%s",t12);
fprintf(fp,"%s",t13);
fprintf(fp,"%s",t14);
fprintf(fp,"%s",t15);
fprintf(fp,"%s",t16);
fprintf(fp,"%s",t17);
fprintf(fp,"%s",t18);
fprintf(fp,"%s",t19);
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
fprintf(fp,"%s",t32);
fprintf(fp,"%s",t33);
fprintf(fp,"%s",t34);
fprintf(fp,"%s",t35);
fprintf(fp,"%s",t36);
fprintf(fp,"%s",t37);
fprintf(fp,"%s",t38);
fprintf(fp,t46,is_null(anf_feld,1/1000) );
fprintf(fp,t47,is_null(end_feld,1/1000) );
fprintf(fp,"%s",t48);
fprintf(fp,t49,temp);
fprintf(fp,"%s",t50);
fprintf(fp,"%s",t51);
fprintf(fp,t52,ionname);
fprintf(fp,"%s",t53);
fprintf(fp,t54,b1,b2,b3);
fprintf(fp,"%s",t55);
fprintf(fp,"%s",t56);
fprintf(fp,"%s",t57);
fprintf(fp,t58,is_null(B1MOL(iteration),1/1000) );
fprintf(fp,t59,is_null(B2MOL(iteration),1/1000) );
fprintf(fp,t60,is_null(B3MOL(iteration),1/1000) );
fprintf(fp,"%s",t61);
}
/*----------------------------------------------------------------------------
                                   raus_kkommentar()
-----------------------------------------------------------------------------*/
void raus_kkommentar(fp,iteration,anf_feld,end_feld,temp,ionname)
  FILE      *fp;
  ITERATION *iteration;
  DOUBLE    anf_feld,end_feld,temp;
  CHAR      *ionname;
{
    UNUSED_PARAMETER(ionname);

    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11,*t12,*t13,*t14,*t15,*t16,*t17,*t18,*t19,*t20;
    CHAR *t21,*t22,*t23,*t24,*t25,*t26,*t27/*,*t28,*t29,*t30;
    CHAR *t31,*t32,*t33,*t34,*t35,*t36,*t37,*t38,*t39,*t40;
    CHAR *t41,*t42,*t43,*t44,*t45,*t46,*t47,*t48,*t49,*t50;
    CHAR *t51,*t52,*t53,*t54,*t55,*t56,*t57,*t58,*t59,*t60;
    CHAR *t61,*t62,*t63,*t64,*t65,*t66,*t67,*t68,*t69,*t70;
    CHAR *t71,*t72,*t73,*t74,*t75,*t76,*t77,*t78,*t79,*t80*/;
 
 
 
t01="# -------------------------------------------------------------- \n";
t02="#|    comment about the calculation of the magnetic moment      \n";
t03="#|    of the Rare Earth Ion in a polycrystalline                \n";
t04="#|    average to the Temperature temp and in the Field range:   \n";
t05="#|    anf_feld bis end_feld                                     \n";
t06="# -------------------------------------------------------------- \n";
t07="#|                                                              \n";
t08="#|  The magnetic moment  m      is                              \n";
t09="#|                        poly                               \n";
t10="#|                                                              \n";
t11="#|             1                                                \n";
t12="#|  m     =  ---- ( 6 m      +  12 m       + 8 m     )          \n";
t13="#|   poly     26       [100]        [110]       [111]           \n";
t14="#|                                                              \n";
t15="#|                                                              \n";
t16="#|  calculated.                                                  \n";
t17="#|                                                              \n";
t18="# -------------------------------------------------------------- \n";
t19="#! start_field       = %6.2f Tesla                             \n"  ;
t20="#! End_field         = %6.2f Tesla                             \n"  ;
t21="#|                                                              \n";
t22="#! temp              = %6.2f Kelvin                            \n";
t23="#|                                                              \n";
t24="#! Bx_mol            = %6.2f Tesla                             \n"  ;
t25="#! By_mol            = %6.2f Tesla                             \n"  ;
t26="#! Bz_mol            = %6.2f Tesla                             \n"  ;
t27="# -------------------------------------------------------------- \n";
 
fprintf(fp,"%s",t01);
fprintf(fp,"%s",t02);
fprintf(fp,"%s",t03);
fprintf(fp,"%s",t04);
fprintf(fp,"%s",t05);
fprintf(fp,"%s",t06);
fprintf(fp,"%s",t07);
fprintf(fp,"%s",t08);
fprintf(fp,"%s",t09);
fprintf(fp,"%s",t10);
 
fprintf(fp,"%s",t11);
fprintf(fp,"%s",t12);
fprintf(fp,"%s",t13);
fprintf(fp,"%s",t14);
fprintf(fp,"%s",t15);
fprintf(fp,"%s",t16);
fprintf(fp,"%s",t17);
fprintf(fp,"%s",t18);
 
 
fprintf(fp,t19,anf_feld);
fprintf(fp,t20,end_feld);
fprintf(fp,"%s",t21);
fprintf(fp,t22,temp);
fprintf(fp,"%s",t23);
fprintf(fp,t24,B1MOL(iteration) );
fprintf(fp,t25,B2MOL(iteration) );
fprintf(fp,t26,B3MOL(iteration) );
fprintf(fp,"%s",t27);
}
/*----------------------------------------------------------------------------
                                   raus_kommentar()
-----------------------------------------------------------------------------*/
void raus_kommentar(fp,anf_temp,end_temp,lambda,theta,lesethetafile,
                 namethetafile)
  FILE   *fp;
  DOUBLE anf_temp,end_temp,lambda,theta;
  INT    lesethetafile;
  CHAR  *namethetafile;
{
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11,*t12,*t13,*t14,*t15,*t16,*t17,*t18,*t19,*t20;
    CHAR *t21,*t22,*t23,*t24,*t25,*t26,*t27,*t28,*t29,*t30;
    CHAR *t31,*t32,*t33,*t34,*t35,*t36,*t37,*t38,*t39,*t40;
    CHAR *t41,*t42,*t43,*t44,*t45,*t46,*t47,*t48,*t49,*t50;
    CHAR *t51,*t52,*t53,*t54,*t55,*t56,*t57,*t58,*t59,*t60;
    CHAR *t61,*t62,*t63,*t64,*t65,*t66,*t67,*t68,*t69/*,*t70*/;
    CHAR *t71,*t72,*t73,*t74,*t75,*t76,*t77,*t78,*t79,*t80;
    CHAR *t81,*t82,*t83,*t84/*,*t85,*t86,*t87,*t88,*t89,*t90;
    CHAR *t91,*t92,*t93,*t94,*t95,*t96,*t97,*t98,*t99*/;
 
t01="# -------------------------------------------------------------- \n";
t02="#|                                                              \n";
t03="#|Comment about the calculation of the inverse Susceptibility   \n";
t04="#|                                                              \n";
t05="# -------------------------------------------------------------- \n";
t06="#|                         CF                                   \n";
t07="#|  The susceptability  X     is                                \n";
t08="#|                         i                                    \n";
t09="#|                                                              \n";
t10="#|   CF 1        2     ----                2                    \n";
t11="#|  X = - *X *  g  *   >    w   |<m| J |n>|   with  i = a,b,c   \n";
t12="#|   i  T   o    J     ----  mn       i                         \n";
t13="#|                    n,m                                       \n";
t14="#|                                                              \n";
t15="#|                                                              \n";
t16="#|                             und                              \n";
t17="#|                                                              \n";
t18="#|                                                              \n";
t19="#|                (E  - E )/T                                   \n";
t20="#|               e  m    n      -    1                          \n";
t21="#|  w   =  w   ---------------------------  with   w   =  w     \n";
t22="#|   mn     m         (E  - E )/T                   nn     n    \n";
t23="#|                      m    n                                  \n";
t24="#|                                                              \n";
t25="#|                                                              \n";
t26="#|                             and                              \n";
t27="#|                                                              \n";
t28="#|              -E /T                                           \n";
t29="#|             e  m                                             \n";
t30="#|  w   =  ------------      ,   Z(T) = partition sum           \n";
t31="#|   m         Z(T)                                             \n";
t32="#|                                                              \n";
t33="#|                                                              \n";
t34="#|  calculates.                                                 \n";
t35="#|                                                              \n";
t36="#|                                                              \n";
t37="#|  with  :                                                     \n";
t38="#|                                                              \n";
t39="#|  a,b,c   :   the three crystal directions                    \n";
t40="#|                                                              \n";
t41="#|  E       :   Crystal Field energy to the state |m> in Kelvin \n";
t42="#|   m                                                          \n";
t43="#|                                                              \n";
t44="#|  g       :   Lande factor                                    \n";
t45="#|   J                                                          \n";
t46="#|                                                              \n";
t47="#|  T       :   Temperature in Kelvin                           \n";
t48="#|                                                              \n";
t49="#|              my0    2                                        \n";
t50="#|      X   =   ---  my  N                                      \n";
t51="#|       o      4pi    B  A                                     \n";
t52="#|                                                              \n";
t53="#|                                                              \n";
t54="#|                                                              \n";
t55="#|          =  0.375150...  emu / mol                           \n";
t56="#|                                                              \n";
t57="# -------------------------------------------------------------- \n";
t58="#|    The Molecular field constant lambda is by                 \n";
t59="#|                                                              \n";
t60="#|                CF                                            \n";
t61="#|    1/X  =  1/X    -  lambda   as defined.                    \n";
t62="#|       i       i                                              \n";
t63="#|                                                              \n";
t64="#|    for T > theta  is  1/X   proportional to  T + theta .     \n";
t65="#|                            i                                 \n";
t66="#|                                                              \n";
t67="# -------------------------------------------------------------- \n";
t68="#|    The susceptability  of the polycrsyalline sample comes    \n";
t69="#|    from                                                      \n"; /*
t70="#|                                                              \n"; */
t71="#|            1                                                 \n";
t72="#|      X  =  - ( X  + X  +  X  ) .                             \n";
t73="#|       p    3    a    b     c                                 \n";
t74="#|                                                              \n";
t75="#|                                                              \n";
t76="#|                                                              \n";
t77="# -------------------------------------------------------------- \n";
t78="#! start_temperature = %4.0f                                    \n";
t79="#! end_temperature     = %4.0f                                  \n";
t80="#|                                                              \n";
t81="#! lambda            = %7.2f mol/emu                          \n";
t82="#! theta             = %7.2f Kelvin                           \n";
t83="#| theta(T) is read from file %16s .            \n";
t84="# -------------------------------------------------------------- \n";
 
fprintf(fp,"%s",t01);
fprintf(fp,"%s",t02);
fprintf(fp,"%s",t03);
fprintf(fp,"%s",t04);
fprintf(fp,"%s",t05);
fprintf(fp,"%s",t06);
fprintf(fp,"%s",t07);
fprintf(fp,"%s",t08);
fprintf(fp,"%s",t09);
fprintf(fp,"%s",t10);
 
fprintf(fp,"%s",t11);
fprintf(fp,"%s",t12);
fprintf(fp,"%s",t13);
fprintf(fp,"%s",t14);
fprintf(fp,"%s",t15);
fprintf(fp,"%s",t16);
fprintf(fp,"%s",t17);
fprintf(fp,"%s",t18);
fprintf(fp,"%s",t19);
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
fprintf(fp,"%s",t32);
fprintf(fp,"%s",t33);
fprintf(fp,"%s",t34);
fprintf(fp,"%s",t35);
fprintf(fp,"%s",t36);
fprintf(fp,"%s",t37);
fprintf(fp,"%s",t38);
fprintf(fp,"%s",t39);
fprintf(fp,"%s",t40);
 
fprintf(fp,"%s",t41);
fprintf(fp,"%s",t42);
fprintf(fp,"%s",t43);
fprintf(fp,"%s",t44);
fprintf(fp,"%s",t45);
fprintf(fp,"%s",t46);
fprintf(fp,"%s",t47);
fprintf(fp,"%s",t48);
fprintf(fp,"%s",t49);
fprintf(fp,"%s",t50);
 
fprintf(fp,"%s",t51);
fprintf(fp,"%s",t52);
fprintf(fp,"%s",t53);
fprintf(fp,"%s",t54);
fprintf(fp,"%s",t55);
fprintf(fp,"%s",t56);
fprintf(fp,"%s",t57);
fprintf(fp,"%s",t58);
fprintf(fp,"%s",t59);
fprintf(fp,"%s",t60);
 
fprintf(fp,"%s",t61);
fprintf(fp,"%s",t62);
fprintf(fp,"%s",t63);
fprintf(fp,"%s",t64);
fprintf(fp,"%s",t65);
fprintf(fp,"%s",t66);
fprintf(fp,"%s",t67);
fprintf(fp,"%s",t68);
fprintf(fp,"%s",t69);
fprintf(fp,"%s",t60);
 
fprintf(fp,"%s",t71);
fprintf(fp,"%s",t72);
fprintf(fp,"%s",t73);
fprintf(fp,"%s",t74);
fprintf(fp,"%s",t75);
fprintf(fp,"%s",t76);
fprintf(fp,"%s",t77);
fprintf(fp,t78,anf_temp);
fprintf(fp,t79,end_temp);
fprintf(fp,"%s",t80);
fprintf(fp,t81,lambda);
if(lesethetafile == NEIN ) fprintf(fp,t82,theta);
else                       fprintf(fp,t83,namethetafile);
fprintf(fp,"%s",t84);

}
/*----------------------------------------------------------------------------
                                   printspalte()
-----------------------------------------------------------------------------*/
INT printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im)
  INT  modus;
  FILE *fp;
  CHAR *t05,*t06,*t07,*t46,*t47;
  INT k,q;
  DOUBLE re,im;
{CHAR *tc,*ts;
 tc="#!%c%1d%1d  =  %16.6f                                     \n";
 ts="#!%c%1d%1ds =              %16.6f                         \n" ;
 if (modus==BKQ || modus=='B' || modus==AKQ || modus=='A')
  {
        if( re!=0.0 ){fprintf(fp,tc,modus,k,q,re);}
        if( im!=0.0 ){fprintf(fp,ts,modus,k,q,im);}
  }
  else
  {
   if(modus!=LKQ && modus!='L'){
      if( re!=0.0 && im!=0.0 ){fprintf(fp,t05,modus,k,q,re,im);}
      if( re!=0.0 && im==0.0 ){fprintf(fp,t06,modus,k,q,re);   }
      if( re==0.0 && im!=0.0 ){fprintf(fp,t07,modus,k,q,im);   }
   }
   else{
         if( re!=0.0 ){fprintf(fp,t46,modus,k,q,re);}
         if( im!=0.0 ){fprintf(fp,t47,modus,k,q,im);}
   }
  }
   return(JA);
}
/*----------------------------------------------------------------------------
                                   parametersatz()
-----------------------------------------------------------------------------*/
void parametersatz(fp,modus,kristallfeld,ionennr,einheit,eingabeart)
    FILE *fp;
    CHAR modus;
    KRISTALLFELD *kristallfeld;
    INT  ionennr;
    CHAR *einheit;
    CHAR eingabeart;
{
    CHAR    *t01,*t02,*t02a,*t02b,*t02l,*t03,*t04,*t05,*t06,*t07/*,*t08,*t09,*t10*/,*t12,*t12a;
    CHAR    *t16,*t22,*t26/*,*t33,*t34*/,*t36/*,*t45*/,*t46,*t47;
    DOUBLE  v20r,v21r,v22r,v40r,v41r,v42r,v43r,v44r;
    DOUBLE  v20i,v21i,v22i,v40i,v41i,v42i,v43i,v44i;
    DOUBLE  v60r,v61r,v62r,v63r,v64r,v65r,v66r;
    DOUBLE  v60i,v61i,v62i,v63i,v64i,v65i,v66i;
    DOUBLE    re,im,x,w,b40,b60,f2,f4,f6;
    ITERATION *iteration;
    INT       e_4f;
    INT       k,q,f;
    INT       zwei_j;
 
    iteration = ITERATION(   kristallfeld );
    zwei_j    = DIMJ(iteration) - 1;
 
if(eingabeart==modus &&  *(FILENAME(kristallfeld)+16) != *(ORTHO+16) ){
  t01="#-------------------------------------------------------------- \n";
  t02="#                     Parameters                     \n";
 /* fprintf(fp,"%s",t01);
  fprintf(fp,"%s",t02);*/
}
if( modus == XW      &&  *(FILENAME(kristallfeld)+16) == *(ORTHO+16) ){
  t01="#-------------------------------------------------------------- \n";
  t02="#                  Next cubic point                   \n";
  fprintf(fp,"%s",t01);
  fprintf(fp,"%s",t02);
}
 
t01="#-------------------------------------------------------------- \n";
t02="# Parameter           :  %ckq   in  %6s                        \n";
t02a="# (NOT the same Vlm as in Hutchings p255 or Elliot and Stevens)\n";
t02b="#                        Bkq are the Stevens Parameters  - see Hutchings Solid State Physics 16 (1964) 227\n";
t02l="#                        Lkq are  Wybourne normalised Parameters (see McPhase Manual)\n";
t12="# Parameter           :  %ckq   in  %6s/a0**k                  \n";
t12a="# (compare Hutchings Solid State Physics 16 (1964) 227, p 255 eq 5.5)                             \n";
t22="# Parameter           :  x,W    in  %6s                        \n";
t03="#--------------------------------------------------------------\n";
t04="#                                                              \n";
t05="#!%c%1d%1d =   %16.6f  +  %16.6f*i              \n" ;
t06="#!%c%1d%1d =   %16.6f                                     \n";
t07="#!%c%1d%1d =              %16.6f*i                        \n" ;
t46="#!%c%1d%1d  =  %16.6f                                     \n";
t47="#!%c%1d-%1d =              %16.6f                         \n" ;
t16="#!   x     =   %16.6f                                \n";
t36="#-1<=x<= 1                                           \n";
t26="#!   W     =   %16.6f                                \n";
    if( modus == SIN )  modus = VKQ;
    modus = up(modus);
 
 
 if( modus == up(XW) ){
 
    b40 = RT( V40(iteration) );
    b60 = RT( V60(iteration) );
    f4  = F4( ionennr        );
    f6  = F6( ionennr        );
 
    if( b40==0.0 && b60==0.0 ){
       fprintf(fp,"%s",t01);
       fprintf(fp,t22,einheit);
 /*      fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);*/
       fprintf(fp,"%s",t36);fprintf(fp,t26,0.0);
/*       fprintf(fp,"%s",t04);*/
       fprintf(fp,"%s",t01);
 /*      if(ZEILE(iteration) ==0 )fprintf(fp,"\n\n");*/
       return;
    }
    else {
 
    if( b40!=0.0 && b60==0.0){
             x = 1; w = b40*f4/x;
             fprintf(fp,"%s",t01);
             fprintf(fp,t22,einheit);
/*             fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);*/
             fprintf(fp,t16,x);fprintf(fp,t26,w);
/*             fprintf(fp,"%s",t04);*/
/*             fprintf(fp,"%s",t01);*/
/*             if(ZEILE(iteration) ==0 )fprintf(fp,"\n\n");*/
             x =-1; w = b40*f4/x;
 /*            fprintf(fp,"%s",t01);*/
             fprintf(fp,t22,einheit);
             fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
             fprintf(fp,t16,x);fprintf(fp,t26,w);
/*             fprintf(fp,"%s",t04);*/
             fprintf(fp,"%s",t01);
/*             if(ZEILE(iteration) ==0 )fprintf(fp,"\n\n");*/
             return;
    }
    else
 
    if( b40==0.0 && b60!=0.0 ){
        x = 0.0;
        w = b60*f6;
    }
    else {

//  if( b40!=0.0 && b60!=0.0 ){
         if( b40*f4==b60*f6 )
              x = 0.5;
         else{  if( b40*f4/b60/f6 < 0.0 )    /* dann x<0 */
                      x = b40*f4/( b60*f6 - b40*f4 );
                else  x = b40*f4/( b60*f6 + b40*f4 );
             }
         w = b40*f4/x;
//  }
    } 
    }
 
    fprintf(fp,"%s",t01);
    fprintf(fp,t22,einheit);
/*    fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);*/
    fprintf(fp,t16,x);fprintf(fp,t26,w);
/*    fprintf(fp,"%s",t04);*/
    fprintf(fp,"%s",t01);
/*    if(ZEILE(iteration) ==0 )fprintf(fp,"\n\n");*/
 
    return;
 }
 else{
    fprintf(fp,"%s",t01);
    if( modus==up(AKQ) || modus == up(WKQ) ) fprintf(fp,t12,modus,einheit);
    else fprintf(fp,t02,modus,einheit);

    if( modus==up(VKQ)) fprintf(fp,"%s",t02a);
    if( modus==up(BKQ)) fprintf(fp,"%s",t02b);
    if( modus==up(LKQ)) fprintf(fp,"%s",t02l);
    if( modus==up(AKQ)) fprintf(fp,"%s",t12a);
/*    fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);*/
    e_4f = E4f( ionennr );
    if(alpha_J[e_4f]==0.0 || beta_J[e_4f]==0.0 || gamma_J[e_4f]==0.0){
       if(alpha_J[e_4f]==0.0) alpha_J[e_4f]=1.0;
       if(beta_J[e_4f]==0.0)  beta_J[e_4f]=1.0;
       if(gamma_J[e_4f]==0.0) gamma_J[e_4f]=1.0;
    }
 }
    v20r = RT( V20(iteration) );v20i = IT( V20(iteration) );
    v21r = RT( V21(iteration) );v21i = IT( V21(iteration) );
    v22r = RT( V22(iteration) );v22i = IT( V22(iteration) );
 
    v40r = RT( V40(iteration) );v40i = IT( V40(iteration) );
    v41r = RT( V41(iteration) );v41i = IT( V41(iteration) );
    v42r = RT( V42(iteration) );v42i = IT( V42(iteration) );
    v43r = RT( V43(iteration) );v43i = IT( V43(iteration) );
    v44r = RT( V44(iteration) );v44i = IT( V44(iteration) );
 
    v60r = RT( V60(iteration) );v60i = IT( V60(iteration) );
    v61r = RT( V61(iteration) );v61i = IT( V61(iteration) );
    v62r = RT( V62(iteration) );v62i = IT( V62(iteration) );
    v63r = RT( V63(iteration) );v63i = IT( V63(iteration) );
    v64r = RT( V64(iteration) );v64i = IT( V64(iteration) );
    v65r = RT( V65(iteration) );v65i = IT( V65(iteration) );
    v66r = RT( V66(iteration) );v66i = IT( V66(iteration) );
 
    e_4f = E4f( ionennr );
    switch(modus){
 
       case AKQ :
       case WKQ :
       case 'A' :
       case 'W' :
                  f2 = A0_BOHR * A0_BOHR;
                  if( zwei_j >= 2 ){
                     v20r /= alpha_J[ e_4f ] * r2(ionennr);
                     v21r /= alpha_J[ e_4f ] * r2(ionennr);
                     v22r /= alpha_J[ e_4f ] * r2(ionennr);
 
                     v20i /= alpha_J[ e_4f ] * r2(ionennr);
                     v21i /= alpha_J[ e_4f ] * r2(ionennr);
                     v22i /= alpha_J[ e_4f ] * r2(ionennr);
                  }
 
                  f4 = f2 * f2;
                  if( zwei_j >= 4 ){
                     v40r /=  beta_J[ e_4f ] * r4(ionennr);
                     v41r /=  beta_J[ e_4f ] * r4(ionennr);
                     v42r /=  beta_J[ e_4f ] * r4(ionennr);
                     v43r /=  beta_J[ e_4f ] * r4(ionennr);
                     v44r /=  beta_J[ e_4f ] * r4(ionennr);
 
                     v40i /=  beta_J[ e_4f ] * r4(ionennr);
                     v41i /=  beta_J[ e_4f ] * r4(ionennr);
                     v42i /=  beta_J[ e_4f ] * r4(ionennr);
                     v43i /=  beta_J[ e_4f ] * r4(ionennr);
                     v44i /=  beta_J[ e_4f ] * r4(ionennr);
                  }
 
                  f6 = f4 * f2;
                  if( zwei_j >= 6 ){
                     v60r /= gamma_J[ e_4f ] * r6(ionennr);
                     v61r /= gamma_J[ e_4f ] * r6(ionennr);
                     v62r /= gamma_J[ e_4f ] * r6(ionennr);
                     v63r /= gamma_J[ e_4f ] * r6(ionennr);
                     v64r /= gamma_J[ e_4f ] * r6(ionennr);
                     v65r /= gamma_J[ e_4f ] * r6(ionennr);
                     v66r /= gamma_J[ e_4f ] * r6(ionennr);
 
                     v60i /= gamma_J[ e_4f ] * r6(ionennr);
                     v61i /= gamma_J[ e_4f ] * r6(ionennr);
                     v62i /= gamma_J[ e_4f ] * r6(ionennr);
                     v63i /= gamma_J[ e_4f ] * r6(ionennr);
                     v64i /= gamma_J[ e_4f ] * r6(ionennr);
                     v65i /= gamma_J[ e_4f ] * r6(ionennr);
                     v66i /= gamma_J[ e_4f ] * r6(ionennr);
                  }
 
       case BKQ :
       case 'B' :
                  if( modus == BKQ || modus == 'B' || modus=='A' || modus==AKQ){
                     v21r *=  2 * omegan1n(1);
                     v22r *=  2 * omegan0n(2);
 
                     v41r *=  2 * omegan3n(1);
                     v42r *=  2 * omegan2n(2);
                     v43r *=  2 * omegan1n(3);
                     v44r *=  2 * omegan0n(4);
 
                     v61r *=  2 * omegan5n(1);
                     v62r *=  2 * omegan4n(2);
                     v63r *=  2 * omegan3n(3);
                     v64r *=  2 * omegan2n(4);
                     v65r *=  2 * omegan1n(5);
                     v66r *=  2 * omegan0n(6);

                     v21i *=  -2 * omegan1n(1);
                     v22i *=  -2 * omegan0n(2);
 
                     v41i *=  -2 * omegan3n(1);
                     v42i *=  -2 * omegan2n(2);
                     v43i *=  -2 * omegan1n(3);
                     v44i *=  -2 * omegan0n(4);
 
                     v61i *=  -2 * omegan5n(1);
                     v62i *=  -2 * omegan4n(2);
                     v63i *=  -2 * omegan3n(3);
                     v64i *=  -2 * omegan2n(4);
                     v65i *=  -2 * omegan1n(5);
                     v66i *=  -2 * omegan0n(6);
                  }
 
 
       case LKQ :
       case 'L' : if(modus == LKQ || modus == 'L' ){
 
                     /*  die Vkq auf die Lkq umrechnen        */
                     /*                                       */
                     /*             q                         */
                     /* D     = (-1)  (  L      + i L      )  */
                     /*  k,|q|            k,|q|      k,-|q|   */
                     /*                                       */
 
                     v21r *=  -1;
                     v41r *=  -1;
                     v43r *=  -1;
                     v61r *=  -1;
                     v63r *=  -1;
                     v65r *=  -1;
 
                     v21i *=  -1;
                     v41i *=  -1;
                     v43i *=  -1;
                     v61i *=  -1;
                     v63i *=  -1;
                     v65i *=  -1;
 
                  }
 
       case DKQ :
       case 'D' :
                 if(modus== LKQ || modus== DKQ || modus== 'L' || modus== 'D'){
                  /* Vkq auf Dkq umrechnen            */
                  /*                                  */
                  /* Dkq = Vkq / epsilon_kq / theta_k */
                  /*                                  */
                  /*                                  */
 
                   if( zwei_j >= 2 ){
                     v20r /=  epn2n(0) * alpha_J[ e_4f ];
                     v21r /=  epn1n(1) * alpha_J[ e_4f ];
                     v22r /=  epn0n(2) * alpha_J[ e_4f ];
                   }
                   if( zwei_j >= 4 ){
                     v40r /=  epn4n(0) * beta_J[ e_4f ];
                     v41r /=  epn3n(1) * beta_J[ e_4f ];
                     v42r /=  epn2n(2) * beta_J[ e_4f ];
                     v43r /=  epn1n(3) * beta_J[ e_4f ];
                     v44r /=  epn0n(4) * beta_J[ e_4f ];
                   }
                   if( zwei_j >= 6 ){
                     v60r /=  epn6n(0) * gamma_J[ e_4f ];
                     v61r /=  epn5n(1) * gamma_J[ e_4f ];
                     v62r /=  epn4n(2) * gamma_J[ e_4f ];
                     v63r /=  epn3n(3) * gamma_J[ e_4f ];
                     v64r /=  epn2n(4) * gamma_J[ e_4f ];
                     v65r /=  epn1n(5) * gamma_J[ e_4f ];
                     v66r /=  epn0n(6) * gamma_J[ e_4f ];
                   }
 
 
                   if( zwei_j >= 2 ){
                     v21i /=  epn1n(1) * alpha_J[ e_4f ];
                     v22i /=  epn0n(2) * alpha_J[ e_4f ];
                   }
                   if( zwei_j >= 4 ){
                     v41i /=  epn3n(1) * beta_J[ e_4f ];
                     v42i /=  epn2n(2) * beta_J[ e_4f ];
                     v43i /=  epn1n(3) * beta_J[ e_4f ];
                     v44i /=  epn0n(4) * beta_J[ e_4f ];
                   }
                   if( zwei_j >= 6 ){
                     v61i /=  epn5n(1) * gamma_J[ e_4f ];
                     v62i /=  epn4n(2) * gamma_J[ e_4f ];
                     v63i /=  epn3n(3) * gamma_J[ e_4f ];
                     v64i /=  epn2n(4) * gamma_J[ e_4f ];
                     v65i /=  epn1n(5) * gamma_J[ e_4f ];
                     v66i /=  epn0n(6) * gamma_J[ e_4f ];
                   }
 
                }
 
       case VKQ :
       case 'V' :
                 f=NEIN;
                 re = v20r;
                 im = v20i;
                 k  = 2;
                 q  = 0;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v21r;
                 im = v21i;
                 if( VOR21(iteration)==-1.0 ){
                     im=re; re=0.0;}
                 k  = 2;
                 q  = 1;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v22r;
                 im = v22i;
                 if( VOR22(iteration)==-1.0 ){
                     im=re; re=0.0;}
                 k  = 2;
                 q  = 2;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 if( f==JA ) { f=NEIN; fprintf(fp,"%s",t04); }
 
                 re = v40r;
                 im = v40i;
                 k  = 4;
                 q  = 0;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v41r;
                 im = v41i;
                 k  = 4;
                 q  = 1;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v42r;
                 im = v42i;
                 if( VOR42(iteration)==-1.0 ){
                     im=re; re=0.0;}
                 k  = 4;
                 q  = 2;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v43r;
                 im = v43i;
                 if( VOR43(iteration)==-1.0 ){
                     im=re; re=0.0;}
                 k  = 4;
                 q  = 3;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
 
                 re = v44r;
                 im = v44i;
                 if( VOR44(iteration)==-1.0 ){
                     im=re; re=0.0;}
                 k  = 4;
                 q  = 4;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 if( f==JA ) { f=NEIN; fprintf(fp,"%s",t04); }
 
                 re = v60r;
                 im = v60i;
                 k  = 6;
                 q  = 0;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v61r;
                 im = v61i;
                 k  = 6;
                 q  = 1;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v62r;
                 im = v62i;
                 if( VOR62(iteration)==-1.0 ){
                     im=re; re=0.0;}
                 k  = 6;
                 q  = 2;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v63r;
                 im = v63i;
                 if( VOR63(iteration)==-1.0 ){
                     im=re; re=0.0;}
                 k  = 6;
                 q  = 3;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v64r;
                 im = v64i;
                 if( VOR64(iteration)==-1.0 ){
                     im=re; re=0.0;}
                 k  = 6;
                 q  = 4;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v65r;
                 im = v65i;
                 k  = 6;
                 q  = 5;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 re = v66r;
                 im = v66i;
                 if( VOR66(iteration)==-1.0 ){
                     im=re; re=0.0;}
                 k  = 6;
                 q  = 6;
                 f=printspalte(modus,fp,t05,t06,t07,t46,t47,k,q,re,im);
 
                 if( f==JA ) { f=NEIN; /*fprintf(fp,"%s",t04); */}
    }
 
    fprintf(fp,"%s",t01);
/*    if(ZEILE(iteration) ==0 )fprintf(fp,"\n\n");*/
 
 
}
/*----------------------------------------------------------------------------
                                   tabelle()
-----------------------------------------------------------------------------*/
void tabelle(fp,anz_niveaus,inten)
  FILE *fp;
  INT  anz_niveaus;
  MATRIX *inten;
{
    CHAR /* *t01,*/*t02,*t03,*t04,*t05,*t06/*,*t07,*t08,*t09,*t10;
    CHAR *t11*/,*t12,*t13/*,*t14*/,*t15,*t16/*,*t17,*t18,*t19,*t20;
    CHAR *t21,*t22,*t23,*t24*/,*t25/*,*t26,*t27,*t28,*t29,*t00*/;
    CHAR *nf_3s(),*f_3s();
    INT /*i,*/test;
    INT zeile,spalte;
    DOUBLE  sum;
 
    kopf(fp,anz_niveaus,JA);/*  wurde Kopf zum ersten Mal ausgegeben : JA */
 
/*
t01="#|    |" ;t11="   |";*/
t02="#| E  |" ;t12="%3s|";
t03="#|  %d |";t13="%3s|";
t04="#|  %2d|";/*t14="   |";*/
t05="#|----|" ;t15="---|"; t25="--- ";
t06="# -----" ;t16="----";
 
    for( zeile=1; zeile<=anz_niveaus ;++zeile ){
       if( zeile==(test=(INT)(anz_niveaus/2)+1) )
          if(test>7)
             kopf(fp,anz_niveaus,NEIN);
        fprintf(fp,"%s",t02);
        sum = 0.0;
        for( spalte=1; spalte<=anz_niveaus ;++spalte ){
             sum += R(inten,zeile,spalte);
             fprintf(fp,t12,f_3s( R(inten,zeile,spalte) ));
        }
        fprintf(fp,t12,f_3s( sum ) );
        fprintf(fp,"\n");
 
        if(zeile<=9) fprintf(fp,t03,zeile);  else fprintf(fp,t04,zeile);
        for( spalte=1; spalte<=anz_niveaus ;++spalte )
           fprintf(fp,t13,nf_3s(R(inten,zeile,spalte)));
        fprintf(fp,t13,nf_3s(sum) );
        fprintf(fp,"\n");
 
        if(zeile< anz_niveaus) fprintf(fp,"%s",t05);  else fprintf(fp,"%s",t06);
        for( spalte=1; spalte<=anz_niveaus+1 ;++spalte )
           if(zeile< anz_niveaus) fprintf(fp,"%s",t15);
           else if(spalte<=anz_niveaus) fprintf(fp,"%s",t16); else fprintf(fp,"%s",t25);
        fprintf(fp,"\n");
 
    }
    fprintf(fp,"#\n");
 
}
/*----------------------------------------------------------------------------
                                   mittlung()
-----------------------------------------------------------------------------*/
void mittlung(fp,anz_niveaus,inten,ew,entartung,
            einheitnr_in,einheitnr_out,einheit_out)
  FILE *fp;
  INT  anz_niveaus;
  MATRIX *inten;
  VEKTOR    *ew;
  MATRIX    *entartung;
  INT       einheitnr_in,einheitnr_out;
  CHAR      *einheit_out;
{
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10,*t11,*t12;
/*  INT i; */
    INT zeile/*,spalte,gi_sum,gi_zeile*/;
    DOUBLE  ew_zeile,ew_zeilem1,int_qe,faktor;
    DOUBLE  int_sum_v,int_sum_g;
    DOUBLE  rel_err_v,rel_err_g;
    DOUBLE  ew_mittel_v,ew_mittel_g;
    DOUBLE  delta_ew_v[17],int_ew_v[17];
    DOUBLE                 int_ew_g[17];
 
    t01="#-----------------------------------------------------------\n" ;
    t02="#!Total_quasielastic_intensity =      %11.2f barn           \n";
    t03="#-----------------------------------------------------------\n";
    t04="#                  Neutron-Energy-loss                      \n";
    t05="#!middle_position_of_the_energy        = %11.2f %6s \n";
    t06="#!relative_error_in_the_middl_Position  = %11.2f %%      \n"    ;
    t07="#!Intensity_of_the_middle_position   = %11.2f barn   \n";
    t08="#                  Neutron-Energy-Gain                 \n";
    t09="#!middle position of the energy        = %11.2f %6s \n";
    t10="#!relative_error_in_the_middl_Position  = %11.2f %%      \n"    ;
    t11="#!Intensity_of_the_middle_position   = %11.2f barn   \n";
    t12="#-----------------------------------------------------------\n";
 
/* gesamte quasielastische Intensitaet berechnen*/
    int_qe = 0.0;
    for( zeile=1 ; zeile<=anz_niveaus ; ++zeile)
       int_qe += R(inten,zeile,zeile);
 
/* energiedifferenzen der leiter berechnen */
    faktor = EINHEITIMP[einheitnr_in].fek*EINHEITIMP[einheitnr_out].fke;
    for( zeile=2 ; zeile<=anz_niveaus ; ++zeile){
       ew_zeile   = RV(ew,(INT)R(entartung,zeile  ,1))*faktor;
       ew_zeilem1 = RV(ew,(INT)R(entartung,zeile-1,1))*faktor;
       delta_ew_v[zeile] = ew_zeile - ew_zeilem1;
    }
 
/* Uebergangsintensitaeten berechnen */
    int_sum_v = 0.0;
    int_sum_g = 0.0;
    for( zeile=2 ; zeile<=anz_niveaus ; ++zeile){
       int_ew_v[zeile] = R(inten,zeile-1,zeile);
       int_sum_v      += int_ew_v[zeile];
       int_ew_g[zeile] = R(inten,zeile,zeile-1);
       int_sum_g      += int_ew_g[zeile];
    }
 
/* gewichtete mittlere Energieposition berechnen */
    ew_mittel_v = 0.0;
    ew_mittel_g = 0.0;
    for( zeile=2 ; zeile<=anz_niveaus ; ++zeile){
         ew_mittel_v += int_ew_v[zeile] * delta_ew_v[zeile];
         ew_mittel_g -= int_ew_g[zeile] * delta_ew_v[zeile];
    }
    if( int_sum_v != 0.0 ) ew_mittel_v /= int_sum_v;
    if( int_sum_g != 0.0 ) ew_mittel_g /= int_sum_g;
 
 
/* gewichteten Fehler in der mittlere Energieposition berechnen */
    rel_err_v = 0.0;
    rel_err_g = 0.0;
    for( zeile=2 ; zeile<=anz_niveaus ; ++zeile){
         rel_err_v += int_ew_v[zeile]
                    * ABSD( delta_ew_v[zeile]-ew_mittel_v);
         rel_err_g += int_ew_g[zeile]
                    * ABSD(-delta_ew_v[zeile]-ew_mittel_g);
    }
    if( int_sum_v   != 0.0 ) rel_err_v /= int_sum_v;
    if( ew_mittel_v != 0.0 ) rel_err_v /= ABSD(ew_mittel_v);
    rel_err_v *= 100.0;
    if( int_sum_g   != 0.0 ) rel_err_g /= int_sum_g;
    if( ew_mittel_g != 0.0 ) rel_err_g /= ABSD(ew_mittel_g);
    rel_err_g *= 100.0;
 
 
    fprintf(fp,"%s",t01);
    fprintf(fp,t02,is_null(int_qe,0.001));
    fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
    fprintf(fp,t05,is_null(ew_mittel_v,0.001),einheit_out);
    fprintf(fp,t06,is_null(rel_err_v,0.001));
    fprintf(fp,t07,is_null(int_sum_v,0.001));
    fprintf(fp,"%s",t08);
    fprintf(fp,t09,is_null(ew_mittel_g,0.001),einheit_out);
    fprintf(fp,t10,is_null(rel_err_g,0.001));
    fprintf(fp,t11,is_null(int_sum_g,0.001));
    fprintf(fp,"%s",t12);
    fprintf(fp,"\n");
 
 
}
/*----------------------------------------------------------------------------
                                   kopf()
-----------------------------------------------------------------------------*/
void kopf(fp,anz_niveaus,flag)  /* tabellenkopf drucken */
  FILE *fp;
  INT  anz_niveaus;
  INT  flag;        /* JA oder NEIN */
{
    CHAR *t01,*t02,*t03,*t04,*t05/*,*t06,*t07,*t08,*t09,*t10*/;
    CHAR *t11,*t12,*t13,*t14,*t15/*,*t16,*t17,*t18,*t19,*t20*/;
    CHAR *t21/*,*t22,*t23*/,*t24/*,*t25,*t26,*t27,*t28,*t29,*t00;
    CHAR *t31*/,*t32,*t33,*t34,*t35/*,*t36,*t37,*t38,*t39,*t40*/;
    INT i;
 
t01="#------" ;t11="----" ;t21="--- ";
t02="# \\ E |";t12="   |" ;
t03="#E \\ k|";t13="E  |" ;
t04="# i \\ |";t14=" %d |"; t24=" %2d|";
t05="#----|" ;t15="---|" ;
 
t21="--- ";
t32="Zei|";
t33="len|";
t34="sum|";
t35="---|";
 
       if( flag==JA ){
         fprintf(fp,"%s",t01);
         for( i=1 ; i<=anz_niveaus ;++i )
            fprintf(fp,"%s",t11);
         fprintf(fp,"%s",t21);
         fprintf(fp,"\n");
       }
 
    fprintf(fp,"%s",t02);
    for( i=1 ; i<=anz_niveaus ;++i )
         fprintf(fp,"%s",t12);
    fprintf(fp,"%s",t32);
    fprintf(fp,"\n");
 
    fprintf(fp,"%s",t03);
    for( i=1 ; i<=anz_niveaus ;++i )
         fprintf(fp,"%s",t13);
    fprintf(fp,"%s",t33);
    fprintf(fp,"\n");
 
    fprintf(fp,"%s",t04);
    for( i=1 ; i<=anz_niveaus ;++i )
         if( i<=9 ) fprintf(fp,t14,i); else fprintf(fp,t24,i);
    fprintf(fp,"%s",t34);
    fprintf(fp,"\n");
 
    fprintf(fp,"%s",t05);
    for( i=1 ; i<=anz_niveaus ;++i )
         fprintf(fp,"%s",t15);
    fprintf(fp,"%s",t35);
    fprintf(fp,"\n");
 
 
}
/*----------------------------------------------------------------------------
                                   f_3s()
-----------------------------------------------------------------------------*/
CHAR *f_3s( zahl )  /*   123.23  ->  "123"    */
  DOUBLE zahl;
 
{
    CHAR *s;
    CHAR i_toc();
    INT  z1,z2,z3/*,z99*/;
 
    s = STRING_ALLOC(4);
 
    z1    = (INT)(zahl+0.005);  /* zahl = 123.9918   ->  z1 = 123  */
    zahl -= z1;
    zahl *= 100;
/*  z99   = (INT)(zahl); */
 
    z2  = z1 - (z1/100)*100;    /* z2 = 23 */
    z3  = z2 - (z2/10 )* 10;    /* z3 = 3  */
    z2 /= 10;
    z1 /=100;
 
    if( z1 > 9 ){
       *s         = '*';
       VALUE(s,1) = '*';
       VALUE(s,2) = '*';
       VALUE(s,3) = '\0';
       return(s);
    }
 
 
       *s         = i_toc(z1);
       VALUE(s,1) = i_toc(z2);
       VALUE(s,2) = i_toc(z3);
       VALUE(s,3) = '\0';
 
       if( z1!=0 )
          return(s);
 
       if( z2!=0 ){
           *s = ' ';
           return(s);
       }
 
       if( z3!=0 ){
          *s         = ' ';
          VALUE(s,1) = ' ';
          return(s);
       }
 
       *s         = ' ';
       VALUE(s,1) = ' ';
       VALUE(s,2) = ' ';
       return(s);
 
}
/*----------------------------------------------------------------------------
                                  nf_3s()
-----------------------------------------------------------------------------*/
CHAR *nf_3s( zahl )  /*  1000.23   ->  ".23"    */
   DOUBLE zahl;
 
{
    CHAR *s;
    CHAR i_toc();
    INT  z1=0,z2=0,z99;
 
    s = STRING_ALLOC(4);
 
    zahl += 0.005;
    zahl -= (INT)zahl;         /* zahl = 123.9918   ->  .9918  */
    zahl *= 100;
    z99   = (INT)(zahl);
    z1  = z99/10;
    z2  = z99 - z1*10;
 
    if( z1!=0 || z2!=0 ){
 
       *s   = '.';
       VALUE(s,1) = i_toc(z1);
       VALUE(s,2) = i_toc(z2);
    }
 
    VALUE(s,3) = '\0';
 
    return(s);
}
/*----------------------------------------------------------------------------
                              i_toc()
-----------------------------------------------------------------------------*/
CHAR i_toc( i )     /*  3 -> '3'   : integer to character */
  INT i;
{
    switch(i){
        case 1 : return( '1' );
        case 2 : return( '2' );
        case 3 : return( '3' );
        case 4 : return( '4' );
        case 5 : return( '5' );
        case 6 : return( '6' );
        case 7 : return( '7' );
        case 8 : return( '8' );
        case 9 : return( '9' );
       default : return( '0' );
    }
}
/*----------------------------------------------------------------------------
                              aJtb_2()
-----------------------------------------------------------------------------*/
     /*                  ---                      2             */
     /* Matrixelememte : >   |<E ;tauº J ºE ;mue>|              */
     /*                  ---    i       T  k                    */
     /*                tau,mue                                  */
MATRIX *aJtb_2( ewproblem,macheps)
      EWPROBLEM *ewproblem;
      DOUBLE    macheps;
{
     MATRIX *aJtb_2,*mx_alloc();
     MATRIX *entartung,*ev;
     INT    groesse,*gi;
     INT    sp,ze;
     DOUBLE sum_mat_Jt2();
 
     entartung = ewproblem->entartung;
     groesse   = ANZ_ZE( entartung );  /* =Anzahl der verschiedenen Niveaus*/
     aJtb_2    = mx_alloc( groesse,groesse );
 
     /*obere Dreiecksmatrix in Matrix aJtb_2 berechnen */
     ev = ewproblem->eigenvektoren;
     gi = ewproblem->gi;
     for( ze=1 ; ze<=groesse ; ++ze )
        for( sp=ze ; sp<=groesse ; ++sp )
            R(aJtb_2,ze,sp) = sum_mat_Jt2(ev,entartung,ze,VALUE(gi,ze)
                                                      ,sp,VALUE(gi,sp)
                                                      ,macheps );
 
 
     /* untere Dreiecksmatrix in Matrix aJtb_2 durch */
     /* Spieglung an der Diagonalen berechnen        */
     for( ze=2 ; ze<=groesse ; ++ze )
       for( sp=1  ; sp<=ze-1 ; ++sp )
          R(aJtb_2,ze,sp) = R(aJtb_2,sp,ze);
 
   return( aJtb_2 );
}
/*----------------------------------------------------------------------------
                              zustandssumme()
-----------------------------------------------------------------------------*/
DOUBLE zustandssumme( einheitnr_in , ew , temperatur )
     INT    einheitnr_in;
     VEKTOR *ew;        /* eigenwerte,kleinster auf 0 geshiftet */
     DOUBLE temperatur; /* in Kelvin */
{
   DOUBLE zustandssumme = 0.0;
   DOUBLE faktor,exp_();
   INT    i;
 
   faktor = EINHEITIMP[ einheitnr_in].fek;
 
   for ( i=1 ; i<=VRDIM(ew) ; ++i )
       zustandssumme += exp_( - faktor*RV(ew,i) / temperatur );
   return( zustandssumme );
}
/*----------------------------------------------------------------------------
                              sum_mat_Jt2()
-----------------------------------------------------------------------------*/
/*  ---                     2*/
/*  >   |<E ,tauºJ ºE ,mue>| */
/*  ---    i      T  k       */
 
DOUBLE sum_mat_Jt2(ev,entartung,zeile,gi_ze,spalte,gi_sp,macheps)
    MATRIX *ev;
    MATRIX *entartung;
    INT    zeile ,gi_ze;
    INT    spalte,gi_sp;
    DOUBLE macheps;
{
    INT    tau,mue;
    VEKTOR *ev_tau,*ev_mue;
    DOUBLE sum=0.0;
    DOUBLE mat_Jt2();
 
    for( tau=1 ; tau<= gi_ze ; ++tau )
       for( mue=1 ; mue<= gi_sp ; ++mue ){
            ev_tau = MXSP(ev, (INT)R(entartung,zeile,tau) );
            ev_mue = MXSP(ev, (INT)R(entartung,spalte,mue) );
            sum   += mat_Jt2( ev_tau,ev_mue,macheps );
       }
 
    return(sum);
}
/*----------------------------------------------------------------------------
                                suszept()
-----------------------------------------------------------------------------*/
/*  ---                  2                                 */
/*  >    w   |<ir|J |ks>|          c=x,y,z                 */
/*  ---   ik       c                                       */
/* i,r,k,s                                                 */
/*                                                         */
/*                                                         */
/*  w  =  exp( -E / T )  /  Z(T)                           */
/*   i           i                                         */
/*                                                         */
/*             exp( (Ei-Ek)/T )  -  1                      */
/*  w   =  w   -----------------------                     */
/*   ik     i       (Ei-Ek)/T                              */
/*                                                         */
 
 
DOUBLE suszept(mat_Ji2,ewproblem,einheitnr_in,t,gj)
    DOUBLE   (*mat_Ji2)();
    EWPROBLEM *ewproblem;
    INT       einheitnr_in;
    DOUBLE    t; /* temperatur in kelvin */
    DOUBLE    gj;
{
    INT    i,k,r,s,anz_niveaus,*gi;
    VEKTOR *ev_ir,*ev_ks;
    VEKTOR *ew;
    MATRIX *ev;
    MATRIX *entartung;
    DOUBLE faktor,exp_(),macheps,sum=0.0,zusumme=0.0;
    DOUBLE ew_i,ew_k,wik;
 
 
    macheps   = ewproblem->eps_machine;
    entartung = ewproblem->entartung;
    gi        = ewproblem->gi;
    ev        = ewproblem->eigenvektoren;
    ew        = ewproblem->eigenwerte;
    faktor        = EINHEITIMP[ einheitnr_in].fek;
    anz_niveaus   = ANZ_ZE(entartung);
 
    for( i=1 ; i<=VRDIM(ew) ; ++i )
        zusumme += exp_( - faktor*RV(ew,i) / t );
 
    for( i=1 ; i<= anz_niveaus ; ++i )
        for( k=1 ; k<= i; ++k )
            for( r=1 ; r<= VALUE(gi,i) ; ++r )
                for( s=1 ; s<= VALUE(gi,k); ++s ){
                   ew_i  = RV(ew,(INT)R(entartung,i,r))*faktor / t;
                   ew_k  = RV(ew,(INT)R(entartung,k,s))*faktor / t;
                   wik = 1.0/zusumme;
                   if( i==k  ) wik *= exp_( - ew_i );
                   else        wik *= 2.0*(exp_(-ew_k)-exp_(-ew_i))/(ew_i-ew_k);
                   ev_ks = MXSP(ev, (INT)R(entartung,k,s) );
                   ev_ir = MXSP(ev, (INT)R(entartung,i,r) );
                   sum  += wik*(*mat_Ji2)(ev_ir,ev_ks,macheps);
                }
 
 
    return(sum*gj*gj*XHI_0/t);
 
}
/*----------------------------------------------------------------------------
                              mat_Jt2()
-----------------------------------------------------------------------------*/
                              /*                          2  */
DOUBLE mat_Jt2(a,b,macheps)   /* Matrixelement  |<aºJ ºb>|   */
    VEKTOR *a,*b;             /*                     T       */
    DOUBLE macheps;           /*  a,b  = ºa>,ºb> Spaltenvektoren */
 
{
   KOMPLEX *mat_Jx();
   KOMPLEX *mat_Jy();
   KOMPLEX *mat_Jz();
   KOMPLEX *ckon();
   KOMPLEX *cmult();
   KOMPLEX *cadd();
   KOMPLEX *aJxb,*aJyb,*aJzb,*f;
   KOMPLEX *aJxb_star,*aJyb_star,*aJzb_star;
   KOMPLEX *aJxb_norm2,*aJyb_norm2,*aJzb_norm2;
   KOMPLEX *aJxyb,*aJxyzb,*aJtb;
   DOUBLE  erg;
 
   aJxb       = mat_Jx(a,b);
   aJxb_star  = ckon( aJxb);
   aJxb_norm2 = cmult(aJxb,aJxb_star);
   free_(aJxb);
   free_(aJxb_star);
 
   aJyb       = mat_Jy(a,b);
   aJyb_star  = ckon( aJyb);
   aJyb_norm2 = cmult(aJyb,aJyb_star);
   free_(aJyb);
   free_(aJyb_star);
 
   aJzb       = mat_Jz(a,b);
   aJzb_star  = ckon( aJzb);
   aJzb_norm2 = cmult(aJzb,aJzb_star);
   free_(aJzb);
   free_(aJzb_star);
 
   aJxyb   = cadd(  aJxb_norm2 , aJyb_norm2 );
   aJxyzb  = cadd(  aJxyb      , aJzb_norm2 );
   free_(aJxb_norm2);
   free_(aJyb_norm2);
   free_(aJzb_norm2);
   free_(aJxyb);
 
   f = KX_ALLOC(1);
   RT(f) = 2.0/3.0;
   IT(f) = 0.0;
 
   aJtb = cmult( f , aJxyzb );
   free_( aJxyzb );
   free_( f );
 
   if( !is_equal(IT(aJtb),0.0,macheps) ){
       printf("Unexpected error in mat_Jt2() in Intensity.c\n");
       printf("              2             \n");
       printf("IT( |<aºJ ºb>|  ) = %f != 0 \n",IT(aJtb) );
       printf("         T                  \n\n");
       exit(0);
   }
   erg = RT(aJtb);
   free_(aJtb);
                   /*           2                  2                2   */
 return(erg);      /* |<a|J |b>| = 2/3*( |<a|J |b>| +... +|<a|J |b>|  ) */
                   /*      T                  x                z        */
}
/*----------------------------------------------------------------------------
                              mat_Jx2()
-----------------------------------------------------------------------------*/
                              /*                          2  */
DOUBLE mat_Jx2(a,b,macheps)   /* Matrixelement  |<aºJ ºb>|   */
    VEKTOR *a,*b;             /*                     x       */
    DOUBLE macheps;           /*  a,b  = ºa>,ºb> Spaltenvektoren */
 
{
   KOMPLEX *mat_Jx();
   KOMPLEX *ckon();
   KOMPLEX *cmult();
   KOMPLEX *cadd();
   KOMPLEX *aJxb;
   KOMPLEX *aJxb_star;
   KOMPLEX *aJxb_norm2;
   DOUBLE  erg;
 
   aJxb       = mat_Jx(a,b);
   aJxb_star  = ckon( aJxb);
   aJxb_norm2 = cmult(aJxb,aJxb_star);
 
   if( !is_equal(IT(aJxb_norm2),0.0,macheps) ){
       printf("Unexpected error in mat_Jx2() in Intensity.c\n");
       printf("              2             \n");
       printf("IT( |<aºJ ºb>|  ) = %f != 0 \n",IT(aJxb_norm2) );
       printf("         x                  \n\n");
       exit(0);
   }
   erg        = RT(aJxb_norm2);
   free_(aJxb);
   free_(aJxb_star);
   free_(aJxb_norm2);
 
   return(erg);
 
}
/*----------------------------------------------------------------------------
                              mat_Jy2()
-----------------------------------------------------------------------------*/
                              /*                          2  */
DOUBLE mat_Jy2(a,b,macheps)   /* Matrixelement  |<aºJ ºb>|   */
    VEKTOR *a,*b;             /*                     y       */
    DOUBLE macheps;           /*  a,b  = ºa>,ºb> Spaltenvektoren */
 
{
   KOMPLEX *mat_Jy();
   KOMPLEX *ckon();
   KOMPLEX *cmult();
   KOMPLEX *cadd();
   KOMPLEX *aJyb;
   KOMPLEX *aJyb_star;
   KOMPLEX *aJyb_norm2;
   DOUBLE  erg;
 
   aJyb       = mat_Jy(a,b);
   aJyb_star  = ckon( aJyb);
   aJyb_norm2 = cmult(aJyb,aJyb_star);
   if( !is_equal(IT(aJyb_norm2),0.0,macheps) ){
       printf("Unexpected error in mat_Jy2() in Intensity.c\n");
       printf("              2             \n");
       printf("IT( |<aºJ ºb>|  ) = %f != 0 \n",IT(aJyb_norm2) );
       printf("         y                  \n\n");
       exit(0);
   }
   erg        = RT(aJyb_norm2);
   free_(aJyb);
   free_(aJyb_star);
   free_(aJyb_norm2);
 
   return(erg);
 
}
/*----------------------------------------------------------------------------
                              mat_Jz2()
-----------------------------------------------------------------------------*/
                              /*                          2  */
DOUBLE mat_Jz2(a,b,macheps)   /* Matrixelement  |<aºJ ºb>|   */
    VEKTOR *a,*b;             /*                     z       */
    DOUBLE macheps;           /*  a,b  = ºa>,ºb> Spaltenvektoren */
 
{
   KOMPLEX *mat_Jz();
   KOMPLEX *ckon();
   KOMPLEX *cmult();
   KOMPLEX *cadd();
   KOMPLEX *aJzb;
   KOMPLEX *aJzb_star;
   KOMPLEX *aJzb_norm2;
   DOUBLE  erg;
 
   aJzb       = mat_Jz(a,b);
   aJzb_star  = ckon( aJzb);
   aJzb_norm2 = cmult(aJzb,aJzb_star);
   if( !is_equal(IT(aJzb_norm2),0.0,macheps) ){
       printf("Unexpected error in mat_Jz2() in Intensity.c\n");
       printf("              2             \n");
       printf("IT( |<aºJ ºb>|  ) = %f != 0 \n",IT(aJzb_norm2) );
       printf("         z                  \n\n");
       exit(0);
   }
   erg        = RT(aJzb_norm2);
   free_(aJzb);
   free_(aJzb_star);
   free_(aJzb_norm2);
 
   return(erg);
 
}
/*----------------------------------------------------------------------------
                              mat_Jx()
-----------------------------------------------------------------------------*/
KOMPLEX *mat_Jx(a,b)  /* Matrixelement  <aºJ ºb>  */
    VEKTOR *a,*b;     /*                    x     */
                      /*  a,b  = ºa>,ºb> Spaltenvektoren */
{
   KOMPLEX *mat_Jp();
   KOMPLEX *mat_Jm();
   KOMPLEX *cadd();
   KOMPLEX *cmult();
   KOMPLEX *aJpb,*aJmb,*aJxb,*aJxb05,*f;
 
   f = KX_ALLOC(1);
   RT(f) = 0.5;
   IT(f) = 0.0;
 
   aJpb = mat_Jp(a,b);
   aJmb = mat_Jm(a,b);
 
   aJxb   = cadd(  aJpb , aJmb );
   aJxb05 = cmult( f    , aJxb );
 
   free_( f  );
   free_(aJpb);
   free_(aJmb);
   free_(aJxb);
 
   return( aJxb05 );  /* <a|Jx|b> = 1/2*( <a|J-|b> + <a|J+|b> ) */
}
/*----------------------------------------------------------------------------
                              mat_Jy()
-----------------------------------------------------------------------------*/
KOMPLEX *mat_Jy(a,b)  /* Matrixelement  <aºJ ºb>  */
    VEKTOR *a,*b;     /*                    y     */
                      /*  a,b  = ºa>,ºb> Spaltenvektoren */
{
   KOMPLEX *mat_Jp();
   KOMPLEX *mat_Jm();
   KOMPLEX *csub();
   KOMPLEX *cmult();
   KOMPLEX *aJpb,*aJmb,*aJyb,*aJybi05,*f;
 
   f = KX_ALLOC(1);
   RT(f) = 0.0;
   IT(f) = 0.5;
 
   aJmb = mat_Jm(a,b);
   aJpb = mat_Jp(a,b);
 
   aJyb    = csub(  aJmb , aJpb );
   aJybi05 = cmult( f    , aJyb );
 
   free_( f  );
   free_(aJpb);
   free_(aJmb);
   free_(aJyb);
 
   return( aJybi05 );  /* <a|Jy|b> = i/2*( <a|J-|b> - <a|J+|b> ) */
}
/*----------------------------------------------------------------------------
                              mat_Jp()
-----------------------------------------------------------------------------*/
KOMPLEX *mat_Jp(a,b)  /* Matrixelement  <aºJ ºb>  */
    VEKTOR *a,*b;     /*                    +     */
                      /*  a,b  = ºa>,ºb> Spaltenvektoren */
 
{
    INT dimj,n;
    KOMPLEX *vskalar();
    VEKTOR  *Jpb,*vr_alloc();
    KOMPLEX *aJpb;
 
    #include "define_j.c"
    dimj = VRDIM(b);
 
    Jpb = vr_alloc( dimj );
    for( n=dimj ; n>=2 ; --n  ){      /*   J+ºb> */
       RV(Jpb,n) = RV(b,n-1)*JP(nj-1);
       IV(Jpb,n) = IV(b,n-1)*JP(nj-1);
    }
    RV(Jpb,1)=IV(Jpb,1)=0.0;
 
    aJpb = vskalar(a,Jpb);
    free_vr(Jpb);
 
    return( aJpb );
}
/*----------------------------------------------------------------------------
                              mat_Jm()
-----------------------------------------------------------------------------*/
KOMPLEX *mat_Jm(a,b)  /* Matrixelement  <aºJ ºb>  */
    VEKTOR *a,*b;     /*                    -     */
                      /*  a,b  = ºa>,ºb> Spaltenvektoren */
 
{
    INT dimj,n;
    KOMPLEX *vskalar();
    VEKTOR  *Jmb,*vr_alloc();
    KOMPLEX *aJmb;
 
    #include "define_j.c"
    dimj = VRDIM(b);
 
    Jmb = vr_alloc( dimj );
    for( n=1 ; n<dimj ; ++n  ){      /*   J-ºb> */
       RV(Jmb,n) = RV(b,n+1)*JM(nj+1);
       IV(Jmb,n) = IV(b,n+1)*JM(nj+1);
    }
    RV(Jmb,dimj)=IV(Jmb,dimj)=0.0;
 
    aJmb = vskalar(a,Jmb);
    free_vr(Jmb);
 
    return( aJmb );
}
/*----------------------------------------------------------------------------
                              mat_Jz()
-----------------------------------------------------------------------------*/
KOMPLEX *mat_Jz(a,b)  /* Matrixelement  <aºJ ºb>  */
    VEKTOR *a,*b;     /*                    z     */
                      /*  a,b  = ºa>,ºb> Spaltenvektoren */
 
{
    INT dimj,n;
    KOMPLEX *vskalar();
    VEKTOR  *Jzb,*vr_alloc();
    KOMPLEX *aJzb;
 
    #include "define_j.c"
    dimj = VRDIM(b);
 
    Jzb = vr_alloc( dimj );
 
    for( n=1 ; n<=dimj ; ++n  ){      /*   Jzºb> */
       RV(Jzb,n) = RV(b,n)*nj;
       IV(Jzb,n) = IV(b,n)*nj;
    }
 
    aJzb  = vskalar(a,Jzb);
    free_vr(Jzb);
 
    return( aJzb );
}
 
/*----------------------------------------------------------------------------
                         is_parametersatz_null()
-----------------------------------------------------------------------------*/
INT is_parametersatz_null(iter,symmetrienr,macheps)
        ITERATION *iter;
        INT       symmetrienr;
        DOUBLE    macheps;
{
    INT zwei_j;
    INT flag = 0;
 
    zwei_j = DIMJ(iter) - 1;
 
 
    if(  zwei_j >2 ){
       switch( symmetrienr ){
 
               case 0 : if( null(RT( V21(iter) ),macheps)) ++flag;
                        if( null(IT( V22(iter) ),macheps)) ++flag;
               case 1 :
               case 2 : if( null(RT( V22(iter) ),macheps)) ++flag;
               case 3 :
               case 4 :
               case 5 :
               case 6 :
               case 7 : if( null(RT( V20(iter) ),macheps)) ++flag;
               case 8 :
               case 9 :
               case 10: break;
       }
    }
 
 
 
    if(  zwei_j >4 ){
        switch( symmetrienr ){
            case  9:
            case 10: break;
 
            case 0 : if( null(RT( V41(iter) ),macheps)) ++flag;
                     if( null(IT( V41(iter) ),macheps)) ++flag;
                     if( null(RT( V43(iter) ),macheps)) ++flag;
                     if( null(IT( V43(iter) ),macheps)) ++flag;
 
            case 1 : if( null(IT( V42(iter) ),macheps)) ++flag;
                     if( null(IT( V44(iter) ),macheps)) ++flag;
 
            case 2 : if( null(RT( V42(iter) ),macheps)) ++flag;
            case 3 :
            case 4 : if( null(RT( V44(iter) ),macheps)) ++flag;
 
            case 5 :
            case 6 : if( symmetrienr == 5 || symmetrienr == 6)
                     if( null(RT( V43(iter) ),macheps)) ++flag;
 
            case 7 :
            case 8 : if( null(RT( V40(iter) ),macheps)) ++flag;
                     if( symmetrienr == 8)
                        if( null(RT( V44(iter) ),macheps)) ++flag;
                     break;
      }
    }
 
 
    if(  zwei_j >6 ){
        switch( symmetrienr ){
            case  9:
            case 10: break;
 
            case 0 : if( null(RT( V61(iter) ),macheps)) ++flag;
                     if( null(IT( V61(iter) ),macheps)) ++flag;
                     if( null(RT( V63(iter) ),macheps)) ++flag;
                     if( null(IT( V63(iter) ),macheps)) ++flag;
                     if( null(RT( V65(iter) ),macheps)) ++flag;
                     if( null(IT( V65(iter) ),macheps)) ++flag;
 
            case 1 : if( null(IT( V62(iter) ),macheps)) ++flag;
                     if( null(IT( V64(iter) ),macheps)) ++flag;
                     if( null(IT( V66(iter) ),macheps)) ++flag;
 
            case 2 : if( null(RT( V62(iter) ),macheps)) ++flag;
                     if( null(RT( V66(iter) ),macheps)) ++flag;
 
            case 3 : if( symmetrienr == 3)
                        if( null(IT( V64(iter) ),macheps)) ++flag;
 
            case 4 : if( null(RT( V64(iter) ),macheps)) ++flag;
 
            case 5 : if( symmetrienr == 5 ){
                        if( null(IT( V63(iter) ),macheps)) ++flag;
                        if( null(IT( V66(iter) ),macheps)) ++flag;
                     }
 
            case 6 : if( symmetrienr == 5 || symmetrienr == 6)
                     if( null(RT( V63(iter) ),macheps)) ++flag;
 
            case 7 : if(symmetrienr==5||symmetrienr==6||symmetrienr==7)
                        if( null(RT( V66(iter) ),macheps)) ++flag;
 
 
            case 8 : if( null(RT( V60(iter) ),macheps)) ++flag;
                     if( symmetrienr == 8)
                        if( null(RT( V64(iter) ),macheps)) ++flag;
                     break;
       }
     }
 
 
    if(  zwei_j >6 ){
        switch( symmetrienr ){
            case 0 : if( flag == 26 )  return(JA);
                     return(NEIN);
            case 1 : if( flag == 14 )  return(JA);
                     return(NEIN);
            case 2 : if( flag ==  9 )  return(JA);
                     return(NEIN);
            case 3 : if( flag ==  6 )  return(JA);
                     return(NEIN);
            case 4 : if( flag ==  5 )  return(JA);
                     return(NEIN);
            case 5 : if( flag ==  8 )  return(JA);
                     return(NEIN);
            case 6 : if( flag ==  6 )  return(JA);
                     return(NEIN);
            case 7 : if( flag ==  4 )  return(JA);
                     return(NEIN);
            case 8 : if( flag ==  4 )  return(JA);
                     return(NEIN);
            case 9 :
            case 10: return(JA);
        }
    }
    if(  zwei_j >4 ){
        switch( symmetrienr ){
            case 0 : if( flag == 13 )  return(JA);
                     return(NEIN);
            case 1 : if( flag ==  7 )  return(JA);
                     return(NEIN);
            case 2 : if( flag ==  5 )  return(JA);
                     return(NEIN);
            case 3 : if( flag ==  3 )  return(JA);
                     return(NEIN);
            case 4 : if( flag ==  3 )  return(JA);
                     return(NEIN);
            case 5 : if( flag ==  3 )  return(JA);
                     return(NEIN);
            case 6 : if( flag ==  3 )  return(JA);
                     return(NEIN);
            case 7 : if( flag ==  2 )  return(JA);
                     return(NEIN);
            case 8 : if( flag ==  2 )  return(JA);
                     return(NEIN);
            case 9 :
            case 10: return(JA);
        }
    }
    if(  zwei_j >=2 ){ /* changed from >2 to >=2 on 28.7.2010 MR to allow S=1 */
        switch( symmetrienr ){
            case 0 : if( flag ==  4 )  return(JA);
                     return(NEIN);
            case 1 : if( flag ==  2 )  return(JA);
                     return(NEIN);
            case 2 : if( flag ==  2 )  return(JA);
                     return(NEIN);
            case 3 : if( flag ==  1 )  return(JA);
                     return(NEIN);
            case 4 : if( flag ==  1 )  return(JA);
                     return(NEIN);
            case 5 : if( flag ==  1 )  return(JA);
                     return(NEIN);
            case 6 : if( flag ==  1 )  return(JA);
                     return(NEIN);
            case 7 : if( flag ==  1 )  return(JA);
                     return(NEIN);
            case 8 : if( flag ==  0 )  return(JA);
                     return(NEIN);
            case 9 :
            case 10: return(JA);
        }
    }
    return(NEIN);
}
/*------------------------------------------------------------------------------
ENDEMODUL    I N T E N S I T    C
------------------------------------------------------------------------------*/
