 
/*-----------------------------------------------------------------------------
 
                                M   O  D  U  L
 
                              ORTHORHOMBISCH C
 
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
extern DOUBLE    exp();              /* definiert in math.c */
extern VEKTOR    *vr_normalisieren();/* definiert in KOMPLEX.C */
 
extern MINIMUM   *amoeba();          /* definiert in  MINIMA.C */
extern MINIMUM   *powell();          /* definiert in  MINIMA.C */
extern MINIMUM   *frprmn();          /* definiert in  MINIMA.C */
extern MINIMUM   *va05a_();          /* definiert in  MINIMA.C */
extern MINIMUM   *fitnr5();          /* definiert in  MINIMA.C */
 
extern EWPROBLEM *diagonalisiere();  /* definiert in DIAHERMX.C*/
extern INT       is_equal();         /* definiert in DIAHERMX.C */
extern DOUBLE    accuracy();         /* definiert in DIAHERMX.C */
extern DOUBLE    is_null();          /* definiert in DIAHERMX.C */
extern INT       null();             /* definiert in DIAHERMX.C */
extern INT       write_title();      /* definiert in DIAHERMX.C */
extern INT       *sort();            /* definiert in DIAHERMX.C */
 
extern INT       output();           /* definiert in INTENSIT.C */
extern MATRIX    *aJtb_2();          /* definiert in INTENSIT.C */
extern DOUBLE    zustandssumme();    /* definiert in INTENSIT.C */
extern DOUBLE    suszept();          /* definiert in INTENSIT.C */
extern DOUBLE    mat_Jx2();          /* definiert in INTENSIT.C */
extern DOUBLE    mat_Jy2();          /* definiert in INTENSIT.C */
extern DOUBLE    mat_Jz2();          /* definiert in INTENSIT.C */
extern DOUBLE    magnetm();          /* definiert in INTENSIT.C */
extern KOMPLEX   *mat_Jx();          /* definiert in INTENSIT.C */
extern KOMPLEX   *mat_Jy();          /* definiert in INTENSIT.C */
extern KOMPLEX   *mat_Jz();          /* definiert in INTENSIT.C */
 
extern VEKTOR    *vr_alloc();        /* definiert in MATRIX.C */
extern INT       free_vr();          /* definiert in MATRIX.C */
extern INT       free_mx();          /* definiert in MATRIX.C */
extern MATRIX    *mx_alloc();        /* definiert in MATRIX.C  */
 
extern DOUBLE    norm_operator();    /* definiert in STEVENS.C */
extern DOUBLE    pow_();             /* definiert in STEVENS.C */
extern MATRIX    *stevkq();          /* definiert in STEVENS.C */
extern MATRIX    *stevks();          /* definiert in STEVENS.C */
extern INT       free_Pkq();         /* definiert in STEVENS.C */
extern INT       free_Okq();         /* definiert in STEVENS.C */
extern STEVENS   *calc_Pkq();        /* definiert in STEVENS.C */
 
extern EINHEIT   EINHEITIMP[];       /* definiert in CFIELD.C*/
extern IONEN     IONENIMP[];         /* definiert in CFIELD.C*/
extern ITERATION *auswahlregel();    /* definiert in CFIELD.C */
extern INT       isimplementiert();  /* definiert in CFIELD.C.*/
extern IONEN     IONENIMP[];         /* definiert in CFIELD.C*/
extern DOUBLE    omegan0n();         /* definiert in CFIELD.C  */
extern DOUBLE    omegan1n();         /* definiert in CFIELD.C  */
extern DOUBLE    omegan2n();         /* definiert in CFIELD.C  */
extern DOUBLE    omegan3n();         /* definiert in CFIELD.C  */
extern DOUBLE    omegan4n();         /* definiert in CFIELD.C  */
extern DOUBLE    omegan5n();         /* definiert in CFIELD.C  */
extern DOUBLE    omegan6n();         /* definiert in CFIELD.C  */
extern DOUBLE    spline();           /* definiert in SPLINE.C*/
extern MATRIX    *calcBmol();        /* definiert in CFIELD.C  */
 
extern FIT       FITIMP[];           /* definiert in Minima.c  */
 
extern NEBENBEDINGUNG *neben_read(); /* definiert in EINGABE.C */
extern FILE *fopen_errchk();         /* definiert in EINGABE.C*/ 
 
#define r0    (-1.91*_e*_e/_m )/*  * 10**(-12) cm  = -0.54*10(-12)cm */
#define const (r0*r0*pi)       /* in barn            2   */
                               /* spaeter noch mit gj  multiplizieren*/
 
/*----------------------------------------------------------------------------
   Internal function declarations
-----------------------------------------------------------------------------*/
void op_norm(STEVENS *stev, DOUBLE macheps);
void par_uebersicht(EWPROBLEM *ewproblem, FILE *fp, VEKTOR *y, MATRIX *p, ITERATION *iter, INT ionennr);
void xi_uebersicht(EWPROBLEM *ewproblem, FILE *fp, MATRIX *xi, ITERATION *iter, INT ionennr);
void grad_uebersicht(EWPROBLEM *ewproblem, FILE *fp, VEKTOR *grad, ITERATION *iter, INT ionennr);
void mpar_kopf(FILE *fp, EXPONENT *exp, INT anz_saetze, CHAR *einheit, CHAR *ionname);
void par_kopf(FILE *fp, EXPONENT *exp, INT anz_saetze, CHAR *einheit, CHAR *ionname);
void mpar_tab(FILE *fp, EXPONENT *exp, VEKTOR *p, INT zeile, INT max_zeile);
void par_tab(INT is_feld, FILE *fp, EXPONENT *exp, VEKTOR *Bkq, DOUBLE chi2, VEKTOR *xW, INT zeile, INT max_zeile);
INT exp_par(DOUBLE z);
INT optimal(INT max_stellen, INT exponent);
DOUBLE mat_Chi2(EWPROBLEM *ewproblem, ITERATION *iteration, MATRIX *aJtb2, MATRIX *intensit_exp, INT i_t, INT i_e);

/*----------------------------------------------------------------------------
                                 exp_()
-----------------------------------------------------------------------------*/
DOUBLE exp_(z)
    DOUBLE z;
{
  #define expppp 100.0
  DOUBLE exp();
  return(z>expppp? exp(expppp): exp(z) );
}
/*----------------------------------------------------------------------------
                                 fit_ortho()
-----------------------------------------------------------------------------*/
void fit_ortho(setup,ewproblem,kristallfeld)
    SETUP        *setup;
    EWPROBLEM    *ewproblem;
    KRISTALLFELD *kristallfeld;
{
    #define DIMFIX 12  /* 9 + 3 */
 
    INT       symmetrienr=2;
 
    NEBENBEDINGUNG *neben_read(),*neben;
    ITERATION *iter,*hamilton();
    STEVENS   *stevens,*calc_Okq(),*calc_Pkq();
    MATRIX    *xi;
    INT       ionennr, dimj,einheitnr_out,i/*,k*/,zeile;
    CHAR      cs,*ionname;
    CHAR   /* *text,*/*t01,*t02,*t03,*t04,*t05,*t06,*t07;
    DOUBLE    macheps/*,ftol*/;
    DOUBLE    Chi2();
    DOUBLE    no2p0,no2p2,no4p2,no6p2,no6p6;
    DOUBLE    no4p05,no4m05,no6p21,no6m21,fret;
    DOUBLE /* bmolx,bmoly,bmolz,*/b_r,b_theta,b_phi;
 
    VEKTOR    *vr_alloc(), *Vkq,*Bkq/*,*xW*/,*v0,*q0,*p0/*,*fix*/;
    VEKTOR    *Vkq_Bkq(),*Bkq_v(),*v_xW();
    VEKTOR    *Bkq_Vkq(),*v_Bkq(),*rphi_v();
    VEKTOR    *y,*v1_v2(),*v2_v1(),*pv, *grad;
 
    MINIMUM   *amoeba(),*fit,*powell(),*va05a_(),*fitnr5();
 
    MATRIX    *p;
    FILE      *fp, *fopen();
 
 
    iter          =  ITERATION(kristallfeld);
    ionname       =  IONNAME( iter );
    ionennr       =  isimplementiert( ionname );
    einheitnr_out =  EINHEITNROUT(  iter );
    dimj          =  DIMJ(   iter);
    macheps       =  ewproblem->eps_machine;
 
    free_Pkq( PKQ(iter) );     /* Matrizen (Jn|Pkq(J)|mJ) holen und */
                               /* freigeben, da nicht gebraucht     */
                           /* Matrizen (Jn|Okq(J)|mJ) und           */
                           /* (Jn|Ok0(J)+s*Ok4(J)|mJ) initialisieren*/
    PKQ(iter) = calc_Okq(dimj,iter);
 
    stevens   = PKQ(iter);
    op_norm(stevens,macheps);          /* Operatornormen berechnen */
    iter  =  auswahlregel(iter,symmetrienr); /* Den Vkq-Parameter die */
                                             /* Symmetrienummer 2     */
                                             /* aufzwingen            */
 
    NEBENBEDINGUNG(iter) = neben_read(setup,CHI2,kristallfeld,einheitnr_out,macheps);
 
    Vkq = vr_alloc( DIMFIX );
 
    RV(Vkq, 1) =  RT( V20(iter) );
    RV(Vkq, 2) =  RT( V22(iter) );
    RV(Vkq, 3) =  RT( V40(iter) );
    RV(Vkq, 4) =  RT( V42(iter) );
    RV(Vkq, 5) =  RT( V44(iter) );
    RV(Vkq, 6) =  RT( V60(iter) );
    RV(Vkq, 7) =  RT( V62(iter) );
    RV(Vkq, 8) =  RT( V64(iter) );
    RV(Vkq, 9) =  RT( V66(iter) );
    RV(Vkq,10) =  BMOL(  iter );
    RV(Vkq,11) =  PHI(   iter );  /* sind in  */
    RV(Vkq,12) =  THETA( iter );  /* radiants */
 
 
    Bkq      = Vkq_Bkq(Vkq,stevens,ionennr,macheps);free_vr(Vkq);
 
    v0       = Bkq_v(  Bkq,stevens,ionennr,macheps);free_vr(Bkq);
    V0(iter) = v0;
 
/********************************************************************/
 neben = NEBENBEDINGUNG( iter ); /*
 fix   = FIX( neben ); */
 
/********************************************************************/
 
 
    cs = EINGABEPARAMETERART(kristallfeld);
    IS_SUSZEPT( kristallfeld ) = NEIN;
    IS_MAGNETM( kristallfeld ) = NEIN;
    IS_KPOLY(   kristallfeld ) = NEIN;
    IS_ORTHO(   kristallfeld ) = NEIN;
    FILENAME(   kristallfeld ) = ORTHO;
    printf("Ergebnisse auf %s herausschreiben...\n",ORTHO);
 
 
    fp = fopen_errchk(FILENAME(kristallfeld),"w");
    write_title(fp);
 
t01=" -------------------------------------------------------------- \n";
t02="|             CRYSTAL FIELD FIT -------- OUTPUT               |\n";
t03=" -------------------------------------------------------------- \n";
t04="\n";
 
    fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
 
/********************************************************************/
 
if( ANZ_VAR(neben) != 0  && ANZAHL(neben) != 0 ){
  printf("Iterating ...\n");
  switch( FITROUTINENNR(iter) ){
    case 1:
               fit = MINIMUM_ALLOC(1);
               fit -> anz_wiederholung = 0;
               fit -> max_wiederholung = 3;
               p0  = v2_v1( iter,v0,macheps);
               fit = amoeba(fit,setup,ewproblem,iter,p0,Chi2);
               free_vr(p0);
            break;
    case 2:
               fit = MINIMUM_ALLOC(1);
               fit -> anz_wiederholung = 0;
               fit -> max_wiederholung = 3;
               p0  = v2_v1( iter,v0,macheps);
               fit = powell(fit,setup,ewproblem,iter,p0,Chi2);
               free_vr(p0);
            break;
    case 3:
               fit = MINIMUM_ALLOC(1);
               fit -> anz_wiederholung = 0;
               fit -> max_wiederholung = 3;
               p0  = v2_v1( iter,v0,macheps);
               fit = frprmn(fit,setup,ewproblem,iter,p0,Chi2);
               free_vr(p0);
            break;
    case 4:
               fit = MINIMUM_ALLOC(1);
               fit -> anz_wiederholung = 0;
               fit -> max_wiederholung = 3;
               p0  = v2_v1( iter,v0,macheps);
 printf("am in 4  \n");
               fit = va05a_(fit,setup,ewproblem,iter,p0,Chi2);
 printf("am out of 4  \n");
               free_vr(p0);
            break;
//  case 5:
    default:
 printf("am in 5  \n");
               fit = MINIMUM_ALLOC(1);
               fit -> anz_wiederholung = 0;
               fit -> max_wiederholung = 3;
               p0  = v2_v1( iter,v0,macheps);
               fit = fitnr5(fit,setup,ewproblem,iter,p0,Chi2);
               free_vr(p0);
            break;
 
  }
//}
//if( ANZ_VAR(neben) != 0  && ANZAHL(neben) != 0 ){
t01=" --------------------------------------------------------------------- \n";
t02="| Fit routine         : %14s                                |\n";
t03="| Iteration           : %14s                                |\n";
t04="| Iterationsverhalten : %14s within %3d                     |\n";
t05="| Iteration step  : %4d <= %7d                               |\n";
t06=" --------------------------------------------------------------------- \n";
t07="\n";
    fprintf(fp,"%s",t01);
    fprintf(fp,t02,FITIMP[FITROUTINENNR(iter)].fitname);
    fprintf(fp,t03,TEXT1(fit));
    fprintf(fp,t04,TEXT2(fit),MAX_WIEDER(fit));
    fprintf(fp,t05,ITER_STEPS(fit),FITMAX(iter));
    fprintf(fp,"%s",t06);
    fprintf(fp,"%s",t07);
 
 
    if( FITROUTINENNR(iter) == 1  ){
 
        p  = MATRIX(fit);
        y  = VEKTOR(fit);
 
        par_uebersicht(ewproblem,fp,y,p,iter,ionennr);
 
        p0 = vr_alloc( ANZ_VAR(neben) );
        for(zeile=1; zeile<=ANZ_ZE(p); ++zeile){
 
             for(i=1; i<=ANZ_SP(p); ++i)
                 RV(p0,i) = R(p,zeile,i);
 
             q0        = v1_v2(iter,p0);
             iter      = hamilton( iter, q0 );
             ewproblem = diagonalisiere( ewproblem,HAMILTONIAN(iter),
                                         NOSPACE,setup);
 
             Bkq  = v_Bkq(    q0,stevens,ionennr,macheps);
             Vkq  = Bkq_Vkq( Bkq,stevens,ionennr,macheps);free_vr(Bkq);
 
             RT( V20(iter) )  =  RV(Vkq,1);
             RT( V22(iter) )  =  RV(Vkq,2);
             RT( V40(iter) )  =  RV(Vkq,3);
             RT( V42(iter) )  =  RV(Vkq,4);
             RT( V44(iter) )  =  RV(Vkq,5);
             RT( V60(iter) )  =  RV(Vkq,6);
             RT( V62(iter) )  =  RV(Vkq,7);
             RT( V64(iter) )  =  RV(Vkq,8);
             RT( V66(iter) )  =  RV(Vkq,9);
             BMOL(  iter   )  =  b_r     = RV(Vkq,10);
             PHI(   iter   )  =  b_phi   = RV(Vkq,11);
             THETA( iter   )  =  b_theta = RV(Vkq,12);
             B1MOL(iter)      =  b_r * sin(b_theta) * cos(b_phi);
             B2MOL(iter)      =  b_r * sin(b_theta) * sin(b_phi);
             B3MOL(iter)      =  b_r * cos(b_theta);
             free_vr(Vkq);
 
             ZEILE(iter)     = zeile;
             MAX_ZEILE(iter) = ANZ_ZE(p);
             fclose(fp);
             output(setup,ewproblem,kristallfeld,cs,ORTHO);
 
             free_vr(q0);
        }/* end zeile */
        ZEILE(iter)     = 0;
        MAX_ZEILE(iter) = 0;
 
        free_vr(p0);    /* == free_vr(VEKTOR(fit)) */
        free_mx(MATRIX(fit));
        free_(fit);
 
 
    }/* Ende Output Fitroutine 1 */
 
 
    if( FITROUTINENNR(iter) == 2 || FITROUTINENNR(iter)==3 ||
        FITROUTINENNR(iter) == 4 || FITROUTINENNR(iter)==5   ){
 
//      if( FITROUTINENNR(iter)==2 ) xi = MATRIX(fit);
        pv   = P_VEKTOR( fit);
//      if( FITROUTINENNR(iter)==3  || FITROUTINENNR(iter)==4 )
//         grad = XI_VEKTOR(fit);
        fret = FRET(     fit);
 
        p  = mx_alloc(1,VRDIM(pv) );
        for(i=1; i<=VRDIM(pv); ++i)
            R(p,1,i) = RV(pv,i);
 
        y  = vr_alloc(1);
        RV(y,1) = fret;
 
        par_uebersicht(ewproblem,fp,y,p,iter,ionennr);
        if( FITROUTINENNR(iter)==2 ) {
            xi = MATRIX(fit);
            xi_uebersicht(ewproblem,fp,xi,iter,ionennr); }
        if( FITROUTINENNR(iter)==3 || FITROUTINENNR(iter)==4 ) {
            grad = XI_VEKTOR(fit);
            grad_uebersicht(ewproblem,fp,grad,iter,ionennr);
            free_vr(grad); }
 
        free_mx(p);
        free_vr(y);
//      if( FITROUTINENNR(iter)==3 || FITROUTINENNR(iter)==4 )
//          free_vr(grad);
 
        q0        = v1_v2(iter,pv    );
        iter      = hamilton( iter, q0 );
        ewproblem = diagonalisiere( ewproblem,HAMILTONIAN(iter),
                                         NOSPACE,setup);
 
        Bkq  = v_Bkq(    q0,stevens,ionennr,macheps);
        Vkq  = Bkq_Vkq( Bkq,stevens,ionennr,macheps);free_vr(Bkq);
 
        RT( V20(iter) )  =  RV(Vkq,1);
        RT( V22(iter) )  =  RV(Vkq,2);
        RT( V40(iter) )  =  RV(Vkq,3);
        RT( V42(iter) )  =  RV(Vkq,4);
        RT( V44(iter) )  =  RV(Vkq,5);
        RT( V60(iter) )  =  RV(Vkq,6);
        RT( V62(iter) )  =  RV(Vkq,7);
        RT( V64(iter) )  =  RV(Vkq,8);
        RT( V66(iter) )  =  RV(Vkq,9);
        BMOL(  iter   )  =  b_r     = RV(Vkq,10);
        PHI(   iter   )  =  b_phi   = RV(Vkq,11);
        THETA( iter   )  =  b_theta = RV(Vkq,12);
        B1MOL(iter)      =  b_r * sin(b_theta) * cos(b_phi);
        B2MOL(iter)      =  b_r * sin(b_theta) * sin(b_phi);
        B3MOL(iter)      =  b_r * cos(b_theta);
        free_vr(Vkq);
 
        ZEILE(iter)     = 1;
        MAX_ZEILE(iter) = 1;
        fclose(fp);
        output(setup,ewproblem,kristallfeld,cs,ORTHO);
        free_vr(q0);
        ZEILE(iter)     = 0;
        MAX_ZEILE(iter) = 0;
 
        free_vr(VEKTOR( fit)); /* == free_vr(pv) == free_vr(p0) */
        if( FITROUTINENNR(iter)==2 ) free_mx(MATRIX( fit));
        free_(fit);
 
 
    }/* Output Fitroutine 2  bis 5 */
}/* end if( ANZ_VAR(neben) != 0  && ANZAHL(neben) != 0 ) */
else{
 
t01=" -------------------------------------------------------------- \n";
t02="| No satisfactory iteration has been found, because:           |\n";
t03="| The number of operatirs to vary  is %3d , or        |\n";
t04="| The number of conditions is           %3d .           |\n";
t05=" -------------------------------------------------------------- \n";
t06="\n";
    fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);
    fprintf(fp,t03,ANZ_VAR(neben));fprintf(fp,t04,ANZAHL(neben));
    fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);
}
 
    no2p0  = N_O2P0(stevens);
    no2p2  = N_O2P2(stevens);
    no4p2  = N_O4P2(stevens);
    no6p2  = N_O6P2(stevens);
    no6p6  = N_O6P6(stevens);
    no4p05 = N_O4P05(stevens);
    no4m05 = N_O4M05(stevens);
    no6m21 = N_O6M21(stevens);
    no6p21 = N_O6P21(stevens);
    free_Okq( PKQ(iter) );/* Matrizen (Jn|Okq(J)|mJ) holen     */
                          /* und freigeben, da nicht gebraucht */
    printf("calculating Pkq(J) ...\n");
    PKQ(iter) = calc_Pkq(dimj);/* Matrizen (Jn|Pkq(J)|mJ) berechnen*/
    stevens = PKQ(iter);
    N_O2P0(stevens) = no2p0;
    N_O2P2(stevens) = no2p2;
    N_O4P2(stevens) = no4p2;
    N_O6P2(stevens) = no6p2;
    N_O6P6(stevens) = no6p6;
    N_O4P05(stevens)= no4p05;
    N_O4M05(stevens)= no4m05;
    N_O6M21(stevens)= no6m21;
    N_O6P21(stevens)= no6p21;
 
    free_vr(v0);
}
/*----------------------------------------------------------------------------
                        v1_v2()
----------------------------------------------------------------------------*/
VEKTOR *v1_v2(iter,v)  /* den Vektor v auf Dimension vom Vektor v0   */
    ITERATION *iter;     /* erweitern, mit Hilfe vom Vektor fix        */
    VEKTOR    *v;
{
    NEBENBEDINGUNG *neben;
    VEKTOR         *fix,*vr_alloc(),*p0,*v0;
    INT            i,k;
 
    neben = NEBENBEDINGUNG( iter );
    fix   = FIX( neben );
    v0    = V0(  iter  );
    p0    = vr_alloc( VRDIM(fix)  );
 
 
    k = 1;
    for(i=1; i<=VRDIM(fix); ++i)
        if( k<= ANZ_VAR(neben) ){
            if( RV(fix,i) == 0.0 ){ /* i-te komponente ist variabel */
                  RV(p0,i) = RV(v,k);
                  if( RV(v0,i) != 0.0 ) RV(p0,i) *= ABSD(RV(v0,i));
                  ++k;
            }
            else  RV(p0,i) = RV(v0,i);
        }
        else  RV(p0,i) = RV(v0,i);
 
    return(p0);
}
/*----------------------------------------------------------------------------
                        v2_v1()
----------------------------------------------------------------------------*/
VEKTOR *v2_v1(iter,v,macheps)
    ITERATION *iter;
    VEKTOR    *v;
    DOUBLE    macheps;
{
    NEBENBEDINGUNG *neben;
    VEKTOR         *fix,*vr_alloc(),*p0/*,*v0*/;
    INT            i,k;
 
    neben = NEBENBEDINGUNG( iter );
    fix   = FIX( neben );
    p0    = vr_alloc( ANZ_VAR(neben) );
 
    k = 1;
    for(i=1; i<=VRDIM(fix); ++i)
      if( RV(fix,i) == 0.0  )  /* i-te komponente ist variabel */
         if( k<= ANZ_VAR(neben) ){
            if(!null(RV(v,i),macheps) ) RV(p0,k) = RV(v,i)/ABSD(RV(v,i));
            else                        RV(p0,k) = 1.0;
            ++k;
         }
 
    return(p0);
}
/*----------------------------------------------------------------------------
                             par_uebersicht()
----------------------------------------------------------------------------*/
void par_uebersicht(ewproblem,fp,y,p,iter,ionennr)
   EWPROBLEM *ewproblem;
   FILE   *fp;
   VEKTOR *y;  /* chi2 */
   MATRIX *p;
   ITERATION *iter;
   INT       ionennr;
{
 
/* CHAR      *t01,*t02,*t03,*t04,*t05,*t06,*t07;
   CHAR      *t08,*t09,*t10,*t11,*t12; */
 
   DOUBLE    macheps;
   INT            zeile,spalte,anz_par_saetze,einheitnr_in;
   INT            exp_B20,exp_B22;
   INT            exp_B40,exp_B42,exp_B44;
   INT            exp_B60,exp_B62,exp_B64,exp_B66;
   INT            exp_Bx,exp_By,exp_Bz;
   INT         /* exp_Bmol,*/is_feld;
   INT            exp_x,exp_W,exp_chi2;
   VEKTOR         *p0,*ps,*Bkq,*v_Bkq(),*xW,*v_xW(),*vr_alloc();
   STEVENS        *stevens;
   CHAR           *einheit_in,*ionname;
   DOUBLE         chi2/*,x,W*/,pow__(),c;
   EXPONENT       *exp;
 
   c             = 180.0/pi;
   ionname       = IONNAME( iter );
   macheps       = ewproblem->eps_machine;
   stevens       = PKQ( iter );
   ps            = vr_alloc(ANZ_SP(p));
   einheitnr_in  = EINHEITNRIN( iter );
   einheit_in    = EINHEITIMP[  einheitnr_in ].einheit;
 
 
   exp = EXP_ALLOC(1);
   exp -> B2i = -1000;
   exp -> B4i = -1000;
   exp -> B6i = -1000;
   exp -> Bmol= -1000;
   exp -> chi2= -1000;
   exp -> x   = -1000;
   exp -> W   = -1000;
   exp_B20    = -1000;
   exp_B22    = -1000;
   exp_B40    = -1000;
   exp_B42    = -1000;
   exp_B44    = -1000;
   exp_B60    = -1000;
   exp_B62    = -1000;
   exp_B64    = -1000;
   exp_B66    = -1000;
   exp_Bx     = -1000;
   exp_By     = -1000;
   exp_Bz     = -1000;
   exp_chi2   = -1000;
   exp_x      = -1000;
   exp_W      = -1000;
 
   is_feld = NEIN;
   for( zeile=1; zeile<=ANZ_ZE(p); ++zeile){
 
        for( spalte=1; spalte<=ANZ_SP(p); ++spalte)
             RV(ps,spalte) = R(p,zeile,spalte);
        p0   = v1_v2( iter,ps );
        Bkq  = v_Bkq(p0,stevens,ionennr,macheps);
        xW   = v_xW( p0,stevens,ionennr,macheps);
        chi2 = RV(y,zeile);
 
if( ! is_equal(RV(Bkq,1),0.0,macheps ) ) exp_B20 = exp_par( RV(Bkq,1) );
if( ! is_equal(RV(Bkq,2),0.0,macheps ) ) exp_B22 = exp_par( RV(Bkq,2) );
if( ! is_equal(RV(Bkq,3),0.0,macheps ) ) exp_B40 = exp_par( RV(Bkq,3) );
if( ! is_equal(RV(Bkq,4),0.0,macheps ) ) exp_B42 = exp_par( RV(Bkq,4) );
if( ! is_equal(RV(Bkq,5),0.0,macheps ) ) exp_B44 = exp_par( RV(Bkq,5) );
if( ! is_equal(RV(Bkq,6),0.0,macheps ) ) exp_B60 = exp_par( RV(Bkq,6) );
if( ! is_equal(RV(Bkq,7),0.0,macheps ) ) exp_B62 = exp_par( RV(Bkq,7) );
if( ! is_equal(RV(Bkq,8),0.0,macheps ) ) exp_B64 = exp_par( RV(Bkq,8) );
if( ! is_equal(RV(Bkq,9),0.0,macheps ) ) exp_B66 = exp_par( RV(Bkq,9) );
if( ! is_equal(RV(Bkq,10),0.0,macheps) ) exp_Bx  = exp_par( RV(Bkq,10));
if( ! is_equal(RV(Bkq,11)*c,0.0,macheps)) exp_By  = exp_par( RV(Bkq,11)*c);
if( ! is_equal(RV(Bkq,12)*c,0.0,macheps)) exp_Bz  = exp_par( RV(Bkq,12)*c);
if( ! is_equal( chi2    ,0.0,macheps ) ) exp_chi2= exp_par(   chi2    );
if( ! is_equal(RV(xW ,1),0.0,macheps ) ) exp_x   = exp_par(  RV(xW,1) );
if( ! is_equal(RV(xW ,2),0.0,macheps ) ) exp_W   = exp_par(  RV(xW,2) );
 
        if( exp_B20 > exp->B2i ) exp->B2i = exp_B20;
        if( exp_B22 > exp->B2i ) exp->B2i = exp_B22;
        if( exp_B40 > exp->B4i ) exp->B4i = exp_B40;
        if( exp_B42 > exp->B4i ) exp->B4i = exp_B42;
        if( exp_B44 > exp->B4i ) exp->B4i = exp_B44;
        if( exp_B60 > exp->B6i ) exp->B6i = exp_B60;
        if( exp_B62 > exp->B6i ) exp->B6i = exp_B62;
        if( exp_B64 > exp->B6i ) exp->B6i = exp_B64;
        if( exp_B66 > exp->B6i ) exp->B6i = exp_B66;
        if( exp_Bx  > exp->Bmol) exp->Bmol= exp_Bx;
        if( exp_By  > exp->Bmol) exp->Bmol= exp_By;
        if( exp_Bz  > exp->Bmol) exp->Bmol= exp_Bz;
        if( exp_chi2> exp->chi2) exp->chi2= exp_chi2;
        if( exp_x   > exp->x   ) exp->x   = exp_x;
        if( exp_W   > exp->W   ) exp->W   = exp_W;
 
        if( !is_equal( RV(Bkq,10),0.0,macheps) ) is_feld=JA;
        if( !is_equal( RV(Bkq,11),0.0,macheps) ) is_feld=JA;
        if( !is_equal( RV(Bkq,12),0.0,macheps) ) is_feld=JA;
 
        free_vr(Bkq);
        free_vr(xW);
        free_vr(p0);
   }
 
if( exp->B2i == -1000 ) exp->B2i = 0;
if( exp->B4i == -1000 ) exp->B4i = 0;
if( exp->B6i == -1000 ) exp->B6i = 0;
if( exp->Bmol== -1000 ) exp->Bmol= 0;
if( exp->chi2== -1000 ) exp->chi2= 0;
if( exp->x   == -1000 ) exp->x   = 0;
if( exp->W   == -1000 ) exp->W   = 0;
 
   exp->B2i =  optimal(5,exp->B2i);
   exp->B4i =  optimal(5,exp->B4i);
   exp->B6i =  optimal(5,exp->B6i);
   exp->Bmol=  optimal(8,exp->Bmol);
   exp->chi2=  optimal(5,exp->chi2);
   exp->x   =  optimal(5,exp->x  );
   exp->W   =  optimal(6,exp->W  );
 
   exp->facB2i = 1.0/pow__(10.0,exp->B2i);
   exp->facB4i = 1.0/pow__(10.0,exp->B4i);
   exp->facB6i = 1.0/pow__(10.0,exp->B6i);
   exp->facBmol= 1.0/pow__(10.0,exp->Bmol);
   exp->facchi2= 1.0/pow__(10.0,exp->chi2);
   exp->facx   = 1.0/pow__(10.0,exp->x   );
   exp->facW   = 1.0/pow__(10.0,exp->W   );
 
 
   anz_par_saetze = ANZ_ZE(p);
   par_kopf(fp,exp,anz_par_saetze,einheit_in,ionname);
   for( zeile=1; zeile<=ANZ_ZE(p); ++zeile){
        for( spalte=1; spalte<=ANZ_SP(p); ++spalte)
             RV(ps,spalte) = R(p,zeile,spalte);
        p0   = v1_v2( iter,ps );
        Bkq  = v_Bkq(p0,stevens,ionennr,macheps);
        xW   = v_xW( p0,stevens,ionennr,macheps);
        chi2 = RV(y,zeile);
 
        par_tab( is_feld,fp,exp,Bkq,chi2,xW,zeile,anz_par_saetze );
 
        free_vr(Bkq);
        free_vr(xW);
        free_vr(p0);
   }
 
 if( is_feld==JA ){
   anz_par_saetze = ANZ_ZE(p);
   mpar_kopf(fp,exp,anz_par_saetze,einheit_in,ionname);
   for( zeile=1; zeile<=ANZ_ZE(p); ++zeile){
        for( spalte=1; spalte<=ANZ_SP(p); ++spalte)
             RV(ps,spalte) = R(p,zeile,spalte);
        p0   = v1_v2( iter,ps );
 
        mpar_tab( fp,exp,p0,zeile,anz_par_saetze );
 
        free_vr(p0);
   }
 }
 
   free_vr(ps);
   free_(exp);
 
 
 
}
/*----------------------------------------------------------------------------
                              xi_uebersicht()
----------------------------------------------------------------------------*/
void xi_uebersicht(ewproblem,fp,xi,iter,ionennr)
   EWPROBLEM *ewproblem;
   FILE      *fp;
   MATRIX    *xi;
   ITERATION *iter;
   INT       ionennr;
{
   UNUSED_PARAMETER(ewproblem);
   UNUSED_PARAMETER(ionennr);

   CHAR      *t01,*t02,*t03,*t04,*t05,*t06,*t07;
   CHAR      *t08,*t09,*t10,*t11/*,*t12*/;
   VEKTOR    *v,*vr_alloc(),*v0_old,*v0/*,*fix*/,*v1_v2(),*w;
   INT       ze,sp;
 
t01=" ----------------------------------------------------------------------------\n";
t02="| Found conjugated directions                                                |\n";
t03="|                                                                            |\n";
t04="| started with   e  =  ( 0,...,0, 1 , 0,...,0 )  i=1,..., n                  |\n";
t05="|                       i                i                                   |\n";
t06="|----------------------------------------------------------------------------|\n";
t07="|vi |             Crystal field parameter                 || molecular field |\n";
t08="|   |                                                     ||                 |\n";
t09="|   | B20 | B22 | B40 | B42 | B44 | B60 | B62 | B64 | B66 || Bx  | By  | Bz  |\n";
t10="|e%2d|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f||%5.2f|%5.2f|%5.2f|\n";
t11=" ----------------------------------------------------------------------------\n";
 
 fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
 fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
 fprintf(fp,"%s",t09);
 
 v0_old   = V0(  iter  );
 v0 = vr_alloc( VRDIM(v0_old) );
 
 V0(iter) = v0;
 
 v = vr_alloc( MXDIM(xi) );
 for( sp=1; sp<=MXDIM(xi); ++sp){
 
     for( ze=1; ze<=MXDIM(xi); ++ze)
          RV(v,ze) = R(xi,ze,sp);
 
     v = vr_normalisieren(v);
     w = v1_v2( iter,v );
     fprintf(fp,t10,sp, RV(w,1), RV(w, 2), RV(w, 3), RV(w, 4),
                        RV(w,5), RV(w, 6), RV(w, 7), RV(w, 8),
                        RV(w,9), RV(w,10), RV(w,11), RV(w,12)  );
 
     free_vr(w);
 }
 free_vr(v);
 
 fprintf(fp,"%s",t11);
 free_vr(v0);
 V0(iter) = v0_old;
 
}
/*----------------------------------------------------------------------------
                              grad_uebersicht()
----------------------------------------------------------------------------*/
void grad_uebersicht(ewproblem,fp,grad,iter,ionennr)
   EWPROBLEM *ewproblem;
   FILE      *fp;
   VEKTOR    *grad;
   ITERATION *iter;
   INT       ionennr;
{
   UNUSED_PARAMETER(ionennr);

   CHAR      *t01,*t02,*t03,*t04,*t05,*t06,*t07;
   CHAR      *t08,*t09,*t10,*t11,*t12;
   VEKTOR /* *v,*/*vr_alloc(),*v0_old,*v0/*,*fix*/,*v1_v2(),*w;
/* INT       ze,sp; */
   EXPONENT  *exp;
   INT            exp_B20,exp_B22;
   INT            exp_B40,exp_B42,exp_B44;
   INT            exp_B60,exp_B62,exp_B64,exp_B66;
   INT            exp_Bx,exp_By,exp_Bz;
/* INT            exp_Bmol; */
   DOUBLE         macheps,pow__();
   STEVENS        *stev;
   VEKTOR         *Bkq;
   macheps = ewproblem->eps_machine;
   stev    = PKQ(iter);
 
 
t01=" ----------------------------------------------------------------------------\n";
t02="| Gradient v                                                                 |\n";
t03="|          -                                                                 |\n";
t04="|----------------------------------------------------------------------------|\n";
t05="|   |              Crsyatal field parameter               || molecular field |\n";
t06="|   |                                                     ||                 |\n";
t07="| v | B20 | B22 | B40 | B42 | B44 | B60 | B62 | B64 | B66 || Bx  | By  | Bz  |\n";
t08="| - |           |                 |                       ||                 |\n";
t09="|   |     (%3d) |        (%3d)    |            (%3d)      ||       (%3d)     |\n";
t10="|   |  /10      |     /10         |         /10           ||    /10          |\n";
t11="|   |%5.0f|%5.0f|%5.0f|%5.0f|%5.0f|%5.0f|%5.0f|%5.0f|%5.0f||%5.0f|%5.0f|%5.0f|\n";
t12=" ----------------------------------------------------------------------------\n";
 
 fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
 fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);
 
 v0_old   = V0(  iter  );
 v0 = vr_alloc( VRDIM(v0_old) );
 
 V0(iter) = v0;
 
 w    = v1_v2( iter,grad );
 Bkq  = vr_alloc( VRDIM(w) );
 
 /*  d               d                                             */
 /*  ---   =  ||O  ||---                                           */
 /*  dB          20  dv                                            */
 /*    20              5                                           */
 /*                                                                */
 RV(Bkq,1) = RV(w,5)*N_O2P0(stev);
 
 /*  d               d                                             */
 /*  ---   =  ||O  ||---                                           */
 /*  dB          22  dw                                            */
 /*    22              6                                           */
 /*                                                                */
 RV(Bkq,2) = RV(w,6)*N_O2P2(stev);
 
 /*  d       1            d        1            d                  */
 /*  ---   = -||O  +5O  ||---   +  -||O  -5O  ||---                */
 /*  dB      2   40   44  dw       2   40   44  dw                 */
 /*    40                   1                     3                */
 /*                                                                */
 RV(Bkq,3) = 0.5*( RV(w,1)*N_O4P05(stev) + RV(w,3)*N_O4M05(stev) );
 
 /*  d               d                                             */
 /*  ---   =  ||O  ||---                                           */
 /*  dB          42  dw                                            */
 /*    42              7                                           */
 /*                                                                */
 RV(Bkq,4) = RV(w,7)*N_O4P2(stev);
 
 /*  d        1               d                   d                */
 /*  ---   = ---( ||O  +5O  ||---  -  ||O  -5O  ||--- )            */
 /*  dB      10      40   44  dv         40   44  dv               */
 /*    44                       1                   3              */
 /*                                                                */
 RV(Bkq,5) = 0.1*( RV(w,1)*N_O4P05(stev) - RV(w,3)*N_O4M05(stev) );
 
 /*  d       1            d        1            d                  */
 /*  ---   = -||O -21O  ||---   +  -||O +21O  ||---                */
 /*  dB      2   60   64  dw       2   60   64  dw                 */
 /*    60                   2                     4                */
 /*                                                                */
 RV(Bkq,6) = 0.5*( RV(w,2)*N_O6M21(stev) + RV(w,4)*N_O6P21(stev) );
 
 /*  d               d                                             */
 /*  ---   =  ||O  ||---                                           */
 /*  dB          62  dw                                            */
 /*    62              8                                           */
 /*                                                                */
 RV(Bkq,7) = RV(w,8)*N_O6P2(stev);
 
 /*  d        1               d                   d                */
 /*  ---   = ---( ||O +21O  ||---  -  ||O -21O  ||--- )            */
 /*  dB      42      60   64  dv         60   64  dv               */
 /*    64                       4                   2              */
 /*                                                                */
 RV(Bkq,8)=1.0/42.0*(RV(w,4)*N_O6P21(stev)-RV(w,2)*N_O6M21(stev));
 
 /*  d               d                                             */
 /*  ---   =  ||O  ||---                                           */
 /*  dB          66  dw                                            */
 /*    66              9                                           */
 /*                                                                */
 RV(Bkq,9) = RV(w,9)*N_O6P6(stev);
 RV(Bkq,10)= RV(w,10);
 RV(Bkq,11)= RV(w,11);
 RV(Bkq,12)= RV(w,12);
 
 
 exp = EXP_ALLOC(1);
 exp -> B2i = -1000;
 exp -> B4i = -1000;
 exp -> B6i = -1000;
 exp -> Bmol= -1000;
 
 exp_B20    = -1000;
 exp_B22    = -1000;
 exp_B40    = -1000;
 exp_B42    = -1000;
 exp_B44    = -1000;
 exp_B60    = -1000;
 exp_B62    = -1000;
 exp_B64    = -1000;
 exp_B66    = -1000;
 exp_Bx     = -1000;
 exp_By     = -1000;
 exp_Bz     = -1000;
 
if( ! is_equal(RV(Bkq,1),0.0,macheps ) ) exp_B20 = exp_par( RV(Bkq,1) );
if( ! is_equal(RV(Bkq,2),0.0,macheps ) ) exp_B22 = exp_par( RV(Bkq,2) );
if( ! is_equal(RV(Bkq,3),0.0,macheps ) ) exp_B40 = exp_par( RV(Bkq,3) );
if( ! is_equal(RV(Bkq,4),0.0,macheps ) ) exp_B42 = exp_par( RV(Bkq,4) );
if( ! is_equal(RV(Bkq,5),0.0,macheps ) ) exp_B44 = exp_par( RV(Bkq,5) );
if( ! is_equal(RV(Bkq,6),0.0,macheps ) ) exp_B60 = exp_par( RV(Bkq,6) );
if( ! is_equal(RV(Bkq,7),0.0,macheps ) ) exp_B62 = exp_par( RV(Bkq,7) );
if( ! is_equal(RV(Bkq,8),0.0,macheps ) ) exp_B64 = exp_par( RV(Bkq,8) );
if( ! is_equal(RV(Bkq,9),0.0,macheps ) ) exp_B66 = exp_par( RV(Bkq,9) );
if( ! is_equal(RV(Bkq,10),0.0,macheps) ) exp_Bx  = exp_par( RV(Bkq,10));
if( ! is_equal(RV(Bkq,11),0.0,macheps) ) exp_By  = exp_par( RV(Bkq,11));
if( ! is_equal(RV(Bkq,12),0.0,macheps) ) exp_Bz  = exp_par( RV(Bkq,12));
 
if( exp_B20 > exp->B2i ) exp->B2i = exp_B20;
if( exp_B22 > exp->B2i ) exp->B2i = exp_B22;
if( exp_B40 > exp->B4i ) exp->B4i = exp_B40;
if( exp_B42 > exp->B4i ) exp->B4i = exp_B42;
if( exp_B44 > exp->B4i ) exp->B4i = exp_B44;
if( exp_B60 > exp->B6i ) exp->B6i = exp_B60;
if( exp_B62 > exp->B6i ) exp->B6i = exp_B62;
if( exp_B64 > exp->B6i ) exp->B6i = exp_B64;
if( exp_B66 > exp->B6i ) exp->B6i = exp_B66;
if( exp_Bx  > exp->Bmol) exp->Bmol= exp_Bx;
if( exp_By  > exp->Bmol) exp->Bmol= exp_By;
if( exp_Bz  > exp->Bmol) exp->Bmol= exp_Bz;
 
 
if( exp->B2i == -1000 ) exp->B2i = 0;
if( exp->B4i == -1000 ) exp->B4i = 0;
if( exp->B6i == -1000 ) exp->B6i = 0;
if( exp->Bmol== -1000 ) exp->Bmol= 0;
 
exp->B2i =  optimal(5,exp->B2i);
exp->B4i =  optimal(5,exp->B4i);
exp->B6i =  optimal(5,exp->B6i);
exp->Bmol=  optimal(5,exp->Bmol);
 
exp->facB2i = 1.0/pow__(10.0,exp->B2i);
exp->facB4i = 1.0/pow__(10.0,exp->B4i);
exp->facB6i = 1.0/pow__(10.0,exp->B6i);
exp->facBmol= 1.0/pow__(10.0,exp->Bmol);
 
fprintf(fp,t09,exp->B2i,exp->B4i,exp->B6i,exp->Bmol);
fprintf(fp,"%s",t10);
 
fprintf(fp,t11,      is_null((RV(Bkq,1)*exp->facB2i ),0.5),
                     is_null((RV(Bkq,2)*exp->facB2i ),0.5),
                     is_null((RV(Bkq,3)*exp->facB4i ),0.5),
                     is_null((RV(Bkq,4)*exp->facB4i ),0.5),
                     is_null((RV(Bkq,5)*exp->facB4i ),0.5),
                     is_null((RV(Bkq,6)*exp->facB6i ),0.5),
                     is_null((RV(Bkq,7)*exp->facB6i ),0.5),
                     is_null((RV(Bkq,8)*exp->facB6i ),0.5),
                     is_null((RV(Bkq,9)*exp->facB6i ),0.5),
                     is_null((RV(Bkq,10)*exp->facBmol ),0.5),
                     is_null((RV(Bkq,11)*exp->facBmol ),0.5),
                     is_null((RV(Bkq,12)*exp->facBmol ),0.5) );
 
 
 free_vr(w);
 
 fprintf(fp,"%s",t12);
 free_vr(v0);
 free_vr(Bkq);
 free_(exp);
 V0(iter) = v0_old;
 
}
/*----------------------------------------------------------------------------
                              optimal()
----------------------------------------------------------------------------*/
INT optimal(max_stellen,exponent)
   INT max_stellen,exponent;
{
 
 INT exp_kritisch;
 exp_kritisch = max_stellen-2;
 
 if(exponent >=0 ){
    if( exponent != exp_kritisch  ) return ( exponent-exp_kritisch );
    if( exponent == exp_kritisch  ) return ( exp_kritisch );
 
 }
 return ( exponent - exp_kritisch );
 
}
/*----------------------------------------------------------------------------
                              pow__()
----------------------------------------------------------------------------*/
DOUBLE pow__(z,n)
   DOUBLE z;
   INT    n;
{
  DOUBLE pow_();
 
  if( n>=0 ) return( pow_(z,n) );
 
  return( 1.0/pow_(z,-n) );
 
}
/*----------------------------------------------------------------------------
                              exp_par()
----------------------------------------------------------------------------*/
INT exp_par(z)
   DOUBLE z;
{
   DOUBLE log(),modf()/*,dummy*/;
   DOUBLE    *iptr,i;
 
    iptr  = (DOUBLE*)c_alloc( 1, sizeof( DOUBLE ) );
   *iptr  = 0.0;
 
   if( z!=0.0) /*dummy  =*/ modf(log(ABSD(z))/log(10.0),iptr );
 
   i = *iptr;
   free_(iptr);
 
   return( (INT)i );
}
/*----------------------------------------------------------------------------
                              par_tab()
----------------------------------------------------------------------------*/
void par_tab( is_feld,fp,exp,Bkq,chi2,xW,zeile,max_zeile)
  INT    is_feld;
  FILE   *fp;
  EXPONENT *exp;
  VEKTOR *Bkq;
  DOUBLE chi2;
  VEKTOR *xW;
  INT    zeile,max_zeile;
{
  CHAR *t01,*t02;
  DOUBLE is_null();
 
t02="|%2d|%5.0f|%5.0f|%5.0f|%5.0f|%5.0f|%5.0f|%5.0f|%5.0f|%5.0f|%5.0f|%6.0f|%6.0f|\n";
fprintf(fp,t02,zeile,is_null((RV(Bkq,1)*exp->facB2i ),0.5),
                     is_null((RV(Bkq,2)*exp->facB2i ),0.5),
                     is_null((RV(Bkq,3)*exp->facB4i ),0.5),
                     is_null((RV(Bkq,4)*exp->facB4i ),0.5),
                     is_null((RV(Bkq,5)*exp->facB4i ),0.5),
                     is_null((RV(Bkq,6)*exp->facB6i ),0.5),
                     is_null((RV(Bkq,7)*exp->facB6i ),0.5),
                     is_null((RV(Bkq,8)*exp->facB6i ),0.5),
                     is_null((RV(Bkq,9)*exp->facB6i ),0.5),
                     is_null((chi2     *exp->facchi2),0.5),
                     is_null((RV(xW,1) *exp->facx   ),0.5),
                     is_null((RV(xW,2) *exp->facW   ),0.5)     );
 
if( zeile == max_zeile && is_feld==NEIN){
t01=" ----------------------------------------------------------------------------\n";
  fprintf(fp,"%s",t01);
}
else if( zeile == max_zeile && is_feld==JA){
t01="|--------------------------------------------------------|-------------------\n";
  fprintf(fp,"%s",t01);
}
 
}
/*----------------------------------------------------------------------------
                             mpar_tab()
----------------------------------------------------------------------------*/
void mpar_tab( fp,exp,p,zeile,max_zeile)
  FILE   *fp;
  EXPONENT *exp;
  VEKTOR *p;
  INT    zeile,max_zeile;
{
  CHAR *t01,*t02;
  DOUBLE is_null(),c;
 
c = 180.0/pi;
t02="|%2d|%8.0f|%8.0f|%8.0f|                          |\n";
fprintf(fp,t02,zeile,is_null((RV(p,10)  *exp->facBmol),0.5),
                     is_null((RV(p,11)*c*exp->facBmol),0.5),
                     is_null((RV(p,12)*c*exp->facBmol),0.5) );
 
 
if( zeile == max_zeile ){
t01=" --------------------------------------------------------";
  fprintf(fp,"%s",t01);
  fprintf(fp,"\n");
}
 
 
}
/*----------------------------------------------------------------------------
                              par_kopf()
----------------------------------------------------------------------------*/
void par_kopf(fp,exp,anz_saetze,einheit,ionname)
  FILE       *fp;
  EXPONENT   *exp;
  INT        anz_saetze;
  CHAR       *einheit,*ionname;
{
  CHAR *t01,*t02,*t03,*t04,*t05,*t06;
  CHAR *t07,*t08,*t09,*t10,*t11/*,*t12*/;
 
t01=" -------------------------------------------------------------- \n";
t02="| It was found %2d Bkq-Parameters  which fitted         |\n";
t03="| the constraints.                                   |\n";
t04=" -------------------------------------------------------------- \n";
t05="\n";
    fprintf(fp,"%s",t01);fprintf(fp,t02,anz_saetze);
    fprintf(fp,"%s",t03);
    fprintf(fp,"%s",t04);
    fprintf(fp,"%s",t05);
t01=" ----------------------------------------------------------------------------\n";
t02="|                                                        |     | next        |\n";
t03="|     Bkq - Parameters in %6s.                   |             | cubic       |\n";
t04="|                                                        |     | Point:      |\n";
t05="|     ion : %4s                                          |     | W in %6s |\n";
t06="|--------------------------------------------------------|-----|-------------|\n";
t07="|  | B20 | B22 | B40 | B42 | B44 | B60 | B62 | B64 | B66 | Chi2|  x   | W    |\n";
t08="|  |           |                 |                       |     |      |      |\n";
t09="|  |    (%3d)  |       (%3d)     |         (%3d)         |(%3d)| (%3d)| (%3d)|\n";
t10="|  | /10       |    /10          |      /10              |/10  |/10   |/10   |\n";
t11="|--|-----------|-----------------|-----------------------|-----|------|------|\n";
 
 
 fprintf(fp,"%s",t01);
 fprintf(fp,"%s",t02);
 fprintf(fp,t03,einheit);
 fprintf(fp,"%s",t04);
 fprintf(fp,t05,ionname,einheit);
 fprintf(fp,"%s",t06);
 fprintf(fp,"%s",t07);
 fprintf(fp,"%s",t08);
 fprintf(fp,t09,exp->B2i,exp->B4i,exp->B6i,exp->chi2,exp->x,exp->W);
 fprintf(fp,"%s",t10);
 fprintf(fp,"%s",t11);
 
}
/*----------------------------------------------------------------------------
                             mpar_kopf()
----------------------------------------------------------------------------*/
void mpar_kopf(fp,exp,anz_saetze,einheit,ionname)
  FILE       *fp;
  EXPONENT   *exp;
  INT        anz_saetze;
  CHAR       *einheit,*ionname;
{
  UNUSED_PARAMETER(anz_saetze);
  UNUSED_PARAMETER(einheit);
  UNUSED_PARAMETER(ionname);

  CHAR *t01,*t02,*t03,*t04,*t05/*,*t06;
  CHAR *t07,*t08,*t09,*t10,*t11,*t12*/;
t01="|  | B      | Phi    | theta  | molecular field  B        |\n";
t02="|  |  mol   |        |        |                  -mol     |\n";
t03="|  |             (%3d)        |  in Tesla.               |\n";
t04="|  |          /10             |                          |\n";
t05="|--|--------|--------|--------|  Phi und theta in Grad.  |\n";
 
 fprintf(fp,"%s",t01);
 fprintf(fp,"%s",t02);
 fprintf(fp,t03,exp->Bmol);
 fprintf(fp,"%s",t04);
 fprintf(fp,"%s",t05);
 
}
/*----------------------------------------------------------------------------
                                  Chi2()
----------------------------------------------------------------------------*/
DOUBLE   Chi2(setup,ewproblem,iteration,v)/*ist funk in amoeba() in Minima.c */
  SETUP     *setup;
  EWPROBLEM *ewproblem;
  ITERATION *iteration;
  VEKTOR    *v;
{
   INT  /* *gi,*/anz_niveaus,einheitnr_in,einheitnr_out;
   MATRIX  *entartung/*,*ev*/,*aJtb_2(),*aJtb2,*intensit_exp,*d_intensit_exp;
   MATRIX  *ueber, *mx_alloc();
   DOUBLE/*shift,*/macheps,temperatur,e_calc,e_exp,gj,deltae;
   DOUBLE  zu_summe,zustandssumme(),faktor,energie,i_calc,i_exp;
   INT     anz_k,anz_n;
 
   CHAR           *ionname;
   VEKTOR      /* *fix,*v0,*/*w,*vr_alloc(),*ew,*ew_exp,*d_ew_exp;
   VEKTOR         *v1_v2();
   KOMPLEX        *mat_Jx(),*mat_Jy(),*mat_Jz();
   NEBENBEDINGUNG *neben;
   INT            free_vr(),i,k,ii,ionennr/*,diff*/;
   INT            datnr,anz_dat,ipos,*nummer,posanzahl,*sort(),*numcomp;
   DOUBLE         pos_icin,pos_icqe,*pos_e,*pos_i,dummy,*werte,*werti;
   DOUBLE         mat_Jx2(),mat_Jy2(),mat_Jz2();
   DOUBLE         chi2,chi2_z/*,chi2_n*/,norm,normi,dum;
   DOUBLE         chi2_1/*,chi2_1z,chi2_1n*/,chi2_g,chi2_p;
   DOUBLE         chi2_a,chi2_b,chi2_c;
/* DOUBLE         chi2_az,chi2_bz,chi2_cz,chi2_pz;
   DOUBLE         chi2_an,chi2_bn,chi2_cn,chi2_pn; */
   DOUBLE         dumx,dumy,dumz,dump,b1,b2,b3,b_norm,sqrt();
   DOUBLE         mag(),e_dummy,i_dummy;
   EWPROBLEM      *diagonalisiere();
   ITERATION      *hamilton();
   DOUBLE         ttheta;
 
/**********************************************************************/
/* angenommener Fehler in der Suszeptibilitaet und der Magnetisierung */
/* und Position                                                       */
/**********************************************************************/
   #define SUSERROR 0.10
   #define MAGERROR 0.10
   #define POSERROR 0.10  /*  10% in Position              bei Posfit*/
   #define INTERROR 0.10  /*  10% in Intensitaet           bei Posfit*/
   #define INTPOS   0.10  /* 0.1*chi2_int + 0.9*chi2_Pos   bei Posfit*/
   #define POSMAX   200
/**********************************************************************/
 
 
 
   neben     = NEBENBEDINGUNG( iteration );
   w         = v1_v2(    iteration, v );
   iteration = hamilton( iteration, w );
   ewproblem = diagonalisiere( ewproblem,HAMILTONIAN(iteration),
                                NOSPACE,setup);
 
   entartung     = ewproblem->entartung;
/* gi            = ewproblem->gi;
   ev            = ewproblem->eigenvektoren; */
   ew            = ewproblem->eigenwerte;
/* shift         = ewproblem->shift; */
   macheps       = ewproblem->eps_machine;
   anz_niveaus   = ANZ_ZE(entartung);
   aJtb2         = aJtb_2(ewproblem,macheps );
   einheitnr_out = EINHEITNROUT(  iteration );
   einheitnr_in  = EINHEITNRIN(   iteration );
   ionname       = IONNAME( iteration );
   ionennr       = isimplementiert( ionname );
 
     anz_n         = ANZ_VAR(         neben );
     ew_exp        =   EIGENWERTE(    neben );
   d_ew_exp        = D_EIGENWERTE(    neben );
     intensit_exp  =   INTENSITAETEN( neben );
   d_intensit_exp  = D_INTENSITAETEN( neben );
 
   faktor          = EINHEITIMP[ einheitnr_in].fek;
 
chi2 = 0.0;
 
if(IS_EIGENWERT(iteration)){
/****************************************/
/*    Chi2   = Chi2_Energieeigenwerte   */
/****************************************/
   chi2_z = 0.0;
   anz_k  = 0;
   norm   = 0.0;
   for( i=1; i<= VRDIM(ew_exp); ++i )
     if( i<= anz_niveaus ){
        e_exp  = RV(ew_exp,i);
        if( !is_equal(e_exp,0.0,macheps) ){
            e_calc = RV(ew,(INT)R(entartung,i,1))
                     *EINHEITIMP[einheitnr_in].fek
                     *EINHEITIMP[einheitnr_out].fke;
            normi   = 1.0/RV(d_ew_exp,i)/RV(d_ew_exp,i);
            chi2_z += (e_exp-e_calc)*(e_exp-e_calc)*normi;
            norm  += normi;
            ++anz_k;
        }
     }
   if( !is_equal(norm,0.0,macheps) )  chi2_z /= norm;
   if( anz_k-anz_n-1 > 0           )  chi2_z /= anz_k-anz_n-1;
 
   chi2 += WEW(iteration)*chi2_z;
}
 
 
/****************************************/
/*    Chi2   = Chi2_Intensitaeten       */
/****************************************/
if( IS_INTENSITAET(iteration) ){
   /*    Chi2  += Chi2_Intensitaeten   */
   temperatur= INTTEMP(iteration);
   norm      = 0.0;
   zu_summe  = zustandssumme( einheitnr_in , ew , temperatur );
   gj        = IONENIMP[ ionennr ].gj;
   chi2_1    = 0.0;
   anz_k  = 0;
   for( i=1; i<= MXDIM(intensit_exp); ++i )     /* i->k */
       for( k=1; k<= MXDIM(intensit_exp); ++k )
           if( i<= anz_niveaus && k<= anz_niveaus ){
               i_exp     = R(intensit_exp,i,k);
               if( !is_equal(i_exp,0.0,macheps) ){
                     energie   = RV( ew , (INT)R(entartung,i,1) );
                     i_calc    = R(aJtb2,i,k) * const*gj*gj
                               *exp_(-energie*faktor/temperatur)/zu_summe;
                     normi     = 1.0/R(d_intensit_exp,i,k)
                                    /R(d_intensit_exp,i,k);
                     chi2_1   += (i_exp-i_calc)*(i_exp-i_calc)*normi;
                     norm     += normi;
                     ++anz_k;
               }
            }
   if( !is_equal(norm  ,0.0,macheps) )  chi2_1  /= norm;
   if( anz_k-anz_n-1 > 0           )    chi2_1  /= anz_k-anz_n-1;
 
   chi2  += WINT(iteration)*chi2_1;
}
 
 
/****************************************/
/*    Chi2   = Chi2_Matrixelemente      */
/****************************************/
else if( IS_MATRIXELEMENT(iteration) ){
   /*    Chi2  += Chi2_Matrixelemente  */
   norm      = 0.0;
   chi2_1    = 0.0;
   anz_k  = 0;
   for( i=1; i<= MXDIM(intensit_exp); ++i )     /* i->k */
       for( k=1; k<= MXDIM(intensit_exp); ++k )
           if( i<= anz_niveaus && k<= anz_niveaus ){
               i_exp     = R(intensit_exp,i,k);
               if( !is_equal(i_exp,0.0,macheps) ){
                     i_calc    = R(aJtb2,i,k);
                     normi     = 1.0/R(d_intensit_exp,i,k)
                                    /R(d_intensit_exp,i,k);
                     chi2_1   += (i_exp-i_calc)*(i_exp-i_calc)*normi;
                     norm     += normi;
                     ++anz_k;
               }
            }
   if( !is_equal(norm  ,0.0,macheps) )  chi2_1 /= norm;
   if( anz_k-anz_n-1 > 0           )    chi2_1 /= anz_k-anz_n-1;
 
   chi2  += WMAT(iteration)*chi2_1;
}
 
 
/********************************************/
/*    Chi2   = Chi2_Positionen+Intensitaeten*/
/********************************************/
if( IS_POSFIT(iteration) ){
   werte     = DOUBLE_ALLOC( POSMAX );
   werti     = DOUBLE_ALLOC( POSMAX );
   ueber     = mx_alloc(anz_niveaus,anz_niveaus);
   gj        = IONENIMP[ ionennr ].gj;
   chi2_g    = 0.0;
   for( datnr=1; datnr<=POSDATANZ(iteration); ++datnr ){
 
      /* experimentelle Daten holen */
      anz_dat    = VALUE(POSANZ( iteration),datnr);
      temperatur = VALUE(POST(   iteration),datnr);
      pos_icin   = VALUE(POSICIN(iteration),datnr);
      pos_icqe   = VALUE(POSICQE(iteration),datnr);
      pos_e      = VALUE(POSE(   iteration),datnr);
      pos_i      = VALUE(POSI(   iteration),datnr);
 
      /* Zustandssumme mittels des theoretischen Spektrums berechnen */
      zu_summe   = zustandssumme( einheitnr_in ,ew,temperatur );
 
      /* Uebergangsmatrixelemente auf Matrix ueber  uebertragen */
      for(i=1; i<=anz_niveaus; ++i)
         for(k=1; k<=anz_niveaus; ++k)
             R(ueber,i,k) = R(aJtb2,i,k);
 
         /* Uebergangsintensitaeten berechnen */
         for(i=1; i<=anz_niveaus; ++i){
            energie  = RV( ew , (INT)R(entartung,i,1) );
            for(k=1; k<=anz_niveaus; ++k)
                R(ueber,i,k) *=  const*gj*gj
                             *exp_(-energie*faktor/temperatur)/zu_summe;
         }
 
      /* Linien berechnen, welche man bei der Temperatur sieht */
      for(i=1; i<=POSMAX; ++i){
          VALUE(werte,i) = 0.0;
          VALUE(werti,i) = 0.0;
      }
      ipos  = 1;
      dummy = 0.1; /* alle linien notieren, deren Intensitaet >= 0.1 b */
      for(i=1; i<=anz_niveaus; ++i)
         for(k=i; k<=anz_niveaus ; ++k)
             if( R(ueber,i,k) >=dummy && ipos<= POSMAX){
                 deltae=RV(ew,(INT)R(entartung,k,1))-RV(ew,(INT)R(entartung,i,1));
                 deltae *= faktor*FKelvinmeV;
                 for(ii=1; ii<=ipos-1; ++ii)
                    if( is_equal(VALUE(werte,ii),deltae,0.001) ){
                        VALUE(werti,ii) += R(ueber,i,k);
                        goto weiter;
                    }
                 VALUE(werte,ipos)  = deltae;
                 VALUE(werti,ipos)  = R(ueber,i,k);
                 ++ipos;
                 weiter:  ;
             }
      posanzahl = ipos-1;
      /* nach aufsteigender energie ordnen */
      nummer = INT_ALLOC(posanzahl);
      nummer = sort(werte,nummer,posanzahl);
/*
for(i=1; i<=anz_dat; ++i)
    printf("e_exp=%12.6f\n",VALUE(pos_e,i) );
printf("\n");
for(i=1; i<=posanzahl; ++i)
    printf("e_calc=%12.6f\n",VALUE(werte,VALUE(nummer,i)));
*/
      /* bestimme welcher theoretischer Uebergang */
      /* am Besten  zum experimentellen passt     */
      numcomp= INT_ALLOC(anz_dat);
      for(i=1;i<=anz_dat;++i){
         chi2_p   = 1000000.0;
         e_exp  = VALUE( pos_e,i); /* stets in meV */
         i_exp  = VALUE( pos_i,i); /* stets in barn*/
         for(k=1;k<=posanzahl;++k){
            e_calc = VALUE( werte, VALUE(nummer,k) );
            i_calc = VALUE( werti, VALUE(nummer,k) );
            e_dummy = (e_calc-e_exp)*(e_calc-e_exp);
            if( !is_equal(e_exp,0.0,macheps) ) e_dummy /= POSERROR*e_exp*POSERROR*e_exp;
            i_dummy = (i_calc-i_exp)*(i_calc-i_exp);
            if( !is_equal(i_exp,0.0,macheps) ) i_dummy /= INTERROR*i_exp*INTERROR*i_exp;
            dummy = e_dummy*(1.0-INTPOS) + i_dummy*INTPOS;
            if( dummy < chi2_p ){
                chi2_p           = dummy;
                VALUE(numcomp,i) = k;
            }
         }
      }
/*
for(i=1; i<=anz_dat; ++i)
    printf("numc(%2d)=%2d  ",i,VALUE(numcomp,i));
printf("\n");
printf("\n");
*/
      /* chi2 fuer die Linienenergien berechnen */
      chi2_p = 0.0;
      anz_k  = 0;
      norm   = 0.0;
      for(i=1;i<=anz_dat;++i){
         e_exp   = VALUE( pos_e,i);
         e_calc  = VALUE( werte, VALUE(numcomp,i) );
         if( !is_equal(e_exp,0.0,macheps) ){
            normi   = POSERROR*e_exp;
            if( !is_equal(normi,0.0,macheps) ) normi = 1.0/normi;
            else                               normi = 1.0;
            normi  *= normi;
            chi2_p += (e_calc-e_exp)*(e_calc-e_exp)*normi;
            norm   += normi;
            ++anz_k;
         }
      }
      if( !is_equal(norm ,0.0,macheps) )  chi2_p  /= norm;
      if( anz_k-anz_n-1 > 0           )   chi2_p /= anz_k-anz_n-1;
      chi2_g += chi2_p*(1.0 - INTPOS);
 
      /* chi2 fuer die Linienintensitaeten  berechnen */
      chi2_p = 0.0;
      anz_k  = 0;
      norm   = 0.0;
      for(i=1;i<=anz_dat;++i){
         e_exp   = VALUE( pos_i,i);
         e_calc  = VALUE( werti, VALUE(numcomp,i) );
         dummy = pos_icin;
         if( is_equal(e_calc,0.0,macheps) )  dummy=pos_icqe;
         e_calc *= dummy;
         if( !is_equal(e_exp,0.0,macheps) ){
            normi   = INTERROR*e_exp;
            if( !is_equal(normi,0.0,macheps) ) normi = 1.0/normi;
            else                               normi = 1.0;
            normi  *= normi;
            chi2_p += (e_calc-e_exp)*(e_calc-e_exp)*normi;
            norm   += normi;
            ++anz_k;
         }
      }
      if( !is_equal(norm ,0.0,macheps) )  chi2_p  /= norm;
      if( anz_k-anz_n-1 > 0           )   chi2_p  /= anz_k-anz_n-1;
      chi2_g += chi2_p*INTPOS;
 
      free_(numcomp);
      free_(nummer);
   }
   chi2  += WPOS(iteration)*chi2_g;
   free_mx(ueber);
   free_(werti);
   free_(werte);
}
 
 
/***********************************************/
/*    Chi2   = Chi2_inverse_Suszeptibilitaeten */
/***********************************************/
if( IS_SUSFIT(iteration) ){
 
 
   chi2_a = 0.0;
   anz_k  = 0;
   if( IS_SUSA(iteration) ){
       norm   = 0.0;
       for( i=1; i<= ANZSA(iteration); ++i){
            ttheta = 0.0;
            if(LESETHETAFILE(iteration)==JA)
               ttheta=spline(THETAANZ(iteration),THETAX(iteration),
                             THETAF(iteration),VALUE(SUSAT(iteration),i),
                             macheps);
            if( VALUE(SUSAT(iteration),i)-ttheta > 0.0 )
                dum = suszept(mat_Jx2,ewproblem,einheitnr_in,
                              VALUE( SUSAT(iteration),i )-ttheta,
                              GJ(          iteration) );
            else dum = 0.0;
 
 
          if( !is_equal(dum,0.0,macheps) && dum>0.0 ){
            if(VALUE(SUSASUS(iteration),i) >  0.0){
              normi   = 1.0/(SUSERROR*VALUE(SUSASUS(iteration),i));
              dum     = 1.0/dum - LAMBDA(iteration);
              dum    -=     VALUE( SUSASUS(iteration),i );
              dum    *= dum*normi*normi;
              chi2_a += dum;
              norm   += normi*normi;
            ++anz_k;
            }
          }
       }
      if( !is_equal(norm,0.0,macheps) )  chi2_a /= norm;
      if( anz_k-anz_n-1 > 0           )  chi2_a  /= anz_k-anz_n-1;
   }
 
   chi2_b = 0.0;
   anz_k  = 0;
   if( IS_SUSB(iteration) ){
       norm   = 0.0;
       for( i=1; i<= ANZSB(iteration); ++i){
            ttheta = 0.0;
            if(LESETHETAFILE(iteration)==JA)
               ttheta=spline(THETAANZ(iteration),THETAX(iteration),
                             THETAF(iteration),VALUE(SUSBT(iteration),i),
                             macheps);
            if( VALUE(SUSBT(iteration),i)-ttheta > 0.0 )
                dum = suszept(mat_Jy2,ewproblem,einheitnr_in,
                              VALUE( SUSBT(iteration),i )-ttheta,
                              GJ(          iteration) );
            else dum = 0.0;
          if( !is_equal(dum,0.0,macheps) && dum>0.0){
            if(VALUE(SUSBSUS(iteration),i) >  0.0){
              normi   = 1.0/(SUSERROR*VALUE(SUSBSUS(iteration),i));
              dum     = 1.0/dum - LAMBDA(iteration);
              dum    -=     VALUE( SUSBSUS(iteration),i );
              dum    *= dum*normi*normi;
              chi2_b += dum;
              norm   += normi*normi;
            ++anz_k;
            }
          }
       }
      if( !is_equal(norm,0.0,macheps) )  chi2_b /= norm;
      if( anz_k-anz_n-1 > 0           )  chi2_b  /= anz_k-anz_n-1;
   }
 
   chi2_c = 0.0;
   anz_k  = 0;
   if( IS_SUSC(iteration) ){
       norm   = 0.0;
       for( i=1; i<= ANZSC(iteration); ++i){
            ttheta = 0.0;
            if(LESETHETAFILE(iteration)==JA)
               ttheta=spline(THETAANZ(iteration),THETAX(iteration),
                             THETAF(iteration),VALUE(SUSCT(iteration),i),
                             macheps);
            if( VALUE(SUSCT(iteration),i)-ttheta > 0.0 )
                dum = suszept(mat_Jz2,ewproblem,einheitnr_in,
                              VALUE( SUSCT(iteration),i )-ttheta,
                              GJ(          iteration) );
            else dum = 0.0;
          if( !is_equal(dum,0.0,macheps) && dum>0.0 ){
            if(VALUE(SUSCSUS(iteration),i) >  0.0){
              normi   = 1.0/(SUSERROR*VALUE(SUSCSUS(iteration),i));
              dum     = 1.0/dum - LAMBDA(iteration);
              dum    -=     VALUE( SUSCSUS(iteration),i );
              dum    *= dum*normi*normi;
              chi2_c += dum;
              norm   += normi*normi;
            ++anz_k;
            }
          }
       }
      if( !is_equal(norm,0.0,macheps) )  chi2_c /= norm;
      if( anz_k-anz_n-1 > 0           )  chi2_c  /= anz_k-anz_n-1;
   }
 
   chi2_p = 0.0;
   anz_k  = 0;
   if( IS_SUSP(iteration) ){
       norm   = 0.0;
       for( i=1; i<= ANZSP(iteration); ++i){
            ttheta = 0.0;
            if(LESETHETAFILE(iteration)==JA)
               ttheta=spline(THETAANZ(iteration),THETAX(iteration),
                             THETAF(iteration),VALUE(SUSPT(iteration),i),
                             macheps);
            if( VALUE(SUSPT(iteration),i)-ttheta > 0.0 )
            dumx = suszept(mat_Jx2,ewproblem,einheitnr_in,
                           VALUE( SUSPT(iteration),i)-ttheta,
                           GJ(iteration) );
            else dumx = 0.0;
            if( !is_equal(dumx,0.0,macheps) ) dumx = 1.0/dumx - LAMBDA(iteration);
            if( VALUE(SUSPT(iteration),i)-ttheta > 0.0 )
            dumy = suszept(mat_Jy2,ewproblem,einheitnr_in,
                           VALUE( SUSPT(iteration),i)-ttheta,
                           GJ(iteration) );
            else dumy = 0.0;
            if( !is_equal(dumy,0.0,macheps) ) dumy = 1.0/dumy - LAMBDA(iteration);
            if( VALUE(SUSPT(iteration),i)-ttheta > 0.0 )
            dumz = suszept(mat_Jz2,ewproblem,einheitnr_in,
                           VALUE( SUSPT(iteration),i)-ttheta,
                           GJ(iteration) );
            else dumz = 0.0;
            if( !is_equal(dumz,0.0,macheps) ) dumz = 1.0/dumz - LAMBDA(iteration);
 
            if( !is_equal(dumx,0.0,macheps) && !is_equal(dumy,0.0,macheps) &&
                !is_equal(dumz,0.0,macheps) && dumx>0.0 && dumy>0.0 && dumz>0.0){
                dump = 1.0/3.0*(1.0/dumx + 1.0/dumy + 1.0/dumz);
                if( !is_equal(dump,0.0,macheps) ){
                  if(VALUE(SUSPSUS(iteration),i) >  0.0){
                    normi = 1.0/(SUSERROR*VALUE(SUSPSUS(iteration),i));
                    dump  = 1.0/dump;
                    dump -=     VALUE( SUSPSUS(iteration),i);
                    dump *= dump*normi*normi;
                    chi2_p += dump;
                    norm   += normi*normi;
                    ++anz_k;
                  }
                }
            }
       }
      if( !is_equal(norm,0.0,macheps) )  chi2_p /= norm;
      if( anz_k-anz_n-1 > 0           )  chi2_p  /= anz_k-anz_n-1;
   }
 
  chi2   += (chi2_a + chi2_b + chi2_c + chi2_p)*WSUS(iteration);
}
 
 
/***********************************************/
/*    Chi2   = Chi2_Magnetisierungskurven      */
/***********************************************/
if( IS_MAGFIT(iteration) ){
 
   chi2_a    = 0.0;
   anz_k  = 0;
   if( IS_MAGA(iteration) ){
       b1        = 1.0;
       b2        = 0.0;
       b3        = 0.0;
       norm      = 0.0;
       temperatur= MAGATEMP(iteration);
       for( i=1; i<= ANZMA(iteration) ; ++i ){
         if(VALUE(MAGAMAG(iteration),i) != 0.0){
             normi= 1.0/(MAGERROR*VALUE(MAGAMAG(iteration),i));
             dum  = mag(mat_Jx,setup,ewproblem,iteration,v,
                        (b1*VALUE(MAGAB(iteration),i)),
                        (b2*VALUE(MAGAB(iteration),i)),
                        (b3*VALUE(MAGAB(iteration),i)),temperatur);
              dum    -=     VALUE( MAGAMAG(iteration),i );
              dum    *= dum*normi*normi;
              chi2_a += dum;
              norm   += normi*normi;
              ++anz_k;
         }
       }
      if( !is_equal(norm,0.0,macheps) )  chi2_a /= norm;
      if( anz_k-anz_n-1 > 0           )  chi2_a  /= anz_k-anz_n-1;
   }
 
 
   chi2_b    = 0.0;
   anz_k  = 0;
   if( IS_MAGB(iteration) ){
       b1        = 0.0;
       b2        = 1.0;
       b3        = 0.0;
       temperatur= MAGBTEMP(iteration);
       norm      = 0.0;
       for( i=1; i<= ANZMB(iteration) ; ++i ){
         if(VALUE(MAGBMAG(iteration),i) != 0.0){
             normi= 1.0/(MAGERROR*VALUE(MAGBMAG(iteration),i));
             dum  = mag(mat_Jy,setup,ewproblem,iteration,v,
                        (b1*VALUE(MAGBB(iteration),i)),
                        (b2*VALUE(MAGBB(iteration),i)),
                        (b3*VALUE(MAGBB(iteration),i)),temperatur);
              dum    -=     VALUE( MAGBMAG(iteration),i );
              dum    *= dum*normi*normi;
              chi2_b += dum;
              norm   += normi*normi;
              ++anz_k;
         }
       }
      if( !is_equal(norm,0.0,macheps) )  chi2_b /= norm;
      if( anz_k-anz_n-1 > 0           )  chi2_b  /= anz_k-anz_n-1;
   }
 
   chi2_c    = 0.0;
   anz_k  = 0;
   if( IS_MAGC(iteration) ){
       b1        = 0.0;
       b2        = 0.0;
       b3        = 1.0;
       norm      = 0.0;
       temperatur= MAGCTEMP(iteration);
       for( i=1; i<= ANZMC(iteration) ; ++i ){
         if(VALUE(MAGCMAG(iteration),i) != 0.0){
             normi= 1.0/(MAGERROR*VALUE(MAGCMAG(iteration),i));
             dum  = mag(mat_Jz,setup,ewproblem,iteration,v,
                        (b1*VALUE(MAGCB(iteration),i)),
                        (b2*VALUE(MAGCB(iteration),i)),
                        (b3*VALUE(MAGCB(iteration),i)),temperatur);
              dum    -=     VALUE( MAGCMAG(iteration),i );
              dum    *= dum*normi*normi;
              chi2_c += dum;
              norm   += normi*normi;
              ++anz_k;
         }
       }
      if( !is_equal(norm,0.0,macheps) )  chi2_c /= norm;
      if( anz_k-anz_n-1 > 0           )  chi2_c  /= anz_k-anz_n-1;
   }
 
   chi2_p    = 0.0;
   anz_k  = 0;
   if( IS_MAGP(iteration)  ){
     b1        = B1(iteration);
     b2        = B2(iteration);
     b3        = B3(iteration);
     temperatur= MAGPTEMP(iteration);
     norm      = 0.0;
     b_norm = sqrt( b1*b1 + b2*b2 + b3*b3 );
     if( !is_equal(b_norm,0.0,macheps)  ){
           b1/=b_norm; b2/=b_norm; b3/=b_norm; }
     else{ b1 = 0.0; b2 = 0.0; b3 = 1.0;     }
     if( !(b1==0.0 && b2==0.0 && b3==1.0 && IS_MAGC(iteration)) ){
       for( i=1; i<= ANZMP(iteration) ; ++i ){
         if(VALUE(MAGPMAG(iteration),i) != 0.0){
             dumx = mag(mat_Jx,setup,ewproblem,iteration,v,
                        (b1*VALUE(MAGPB(iteration),i)),
                        (b2*VALUE(MAGPB(iteration),i)),
                        (b3*VALUE(MAGPB(iteration),i)),temperatur);
             dumy = mag(mat_Jy,setup,ewproblem,iteration,v,
                        (b1*VALUE(MAGPB(iteration),i)),
                        (b2*VALUE(MAGPB(iteration),i)),
                        (b3*VALUE(MAGPB(iteration),i)),temperatur);
             dumz = mag(mat_Jz,setup,ewproblem,iteration,v,
                        (b1*VALUE(MAGPB(iteration),i)),
                        (b2*VALUE(MAGPB(iteration),i)),
                        (b3*VALUE(MAGPB(iteration),i)),temperatur);
 
              normi   = 1.0/(MAGERROR*VALUE(MAGPMAG(iteration),i));
              dum     = sqrt( dumx*dumx + dumy*dumy + dumz*dumz );
              dum    -=     VALUE( MAGPMAG(iteration),i );
              dum    *= dum*normi*normi;
              chi2_p += dum;
              norm   += normi*normi;
              ++anz_k;
         }
       }
     }
      if( !is_equal(norm,0.0,macheps) )  chi2_p /= norm;
      if( anz_k-anz_n-1 > 0           )  chi2_p  /= anz_k-anz_n-1;
   }
 
   chi2   += (chi2_a + chi2_b + chi2_c + chi2_p)*WMAG(iteration);
}
 
 
   free_vr(w);
   free_mx(aJtb2);
 
 
   return(chi2);
}
 
/*----------------------------------------------------------------------------
                                mag()
-----------------------------------------------------------------------------*/
/*          ---                                         */
/* m  = -g  >    w   <ir|J |ir>     c=x,y,z             */
/*  c     j ---   i       c                             */
/*         i,r                                          */
/*                                                         */
/*                                                         */
/*  w  =  exp( -E / T )  /  Z(T)                           */
/*   i           i                                         */
/*                                                         */
 
 
DOUBLE mag(mat_Ji,setup,ewproblem,iteration,v,Bx,By,Bz,t)
    KOMPLEX   *(*mat_Ji)();
    SETUP     *setup;
    EWPROBLEM *ewproblem;
    ITERATION *iteration;
    VEKTOR    *v;
    DOUBLE    Bx,By,Bz; /* angelegetes aeusseres feld*/
    DOUBLE    t; /* angelegete  Temperatur    */
{
    INT     i/*,k*/,r/*,s*/,anz_niveaus,*gi;
    ITERATION  *hamilton();
    VEKTOR  *ev_ir,*w,*v1_v2();
    VEKTOR  *ew;
    MATRIX  *ev/*,*bmag*/;
    MATRIX  *entartung;
    DOUBLE  faktor,exp_(),macheps,sumr=0.0,sumi=0.0,zusumme=0.0;
    DOUBLE  ew_i,wi,gj/*,myB*/;
    KOMPLEX *mat;
    INT     einheitnr_in;
 
 
    einheitnr_in = EINHEITNRIN( iteration );
/*  bmag         = HMAG( iteration ); */
    gj           = GJ(iteration );
/*  myB          = EINHEITIMP[ einheitnr_in ].myB; */
 
 
    w          = v1_v2( iteration,v );
    RV(w,10)  += Bx*gj/(2.0*(gj-1));
    RV(w,11)  += By*gj/(2.0*(gj-1));
    RV(w,12)  += Bz*gj/(2.0*(gj-1));
 
    iteration = hamilton( iteration,w );
    free_vr(w);
 
    ewproblem = diagonalisiere( ewproblem,HAMILTONIAN(iteration), NOSPACE,setup);
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
      printf("\n Imaginary part of m_c = %e is not zero!\n",sumi);
 
 
    return(-gj*sumr);
 
}
 
/*------------------------------------------------------------------------------
                               hamilton()
--------------------------------------------------------------------------------
 
 H  =  H    +  H
        CF      mol
 
 
                ~              ~              ~              ~
 H    =  v1 * O4P05  +  v2 * O6M21  +  v3 * O4M05  +  v4 * O6P21
  CF
                              ~              ~              ~
                     +  v5 * O20    +  v6 * O22    +  v7 * O42
 
                                             ~              ~
                                    +  v8 * O62    +  v9 * O66
 
 
 H    = myB 2(g -1) J B
  mol          J    - - mol
 
------------------------------------------------------------------------------*/
ITERATION *hamilton(i,v)
    ITERATION *i;
    VEKTOR    *v;
{
    STEVENS *s;
    MATRIX  *h,*mag,*calcBmol();
    INT     n,m;
    DOUBLE  myB;
    DOUBLE  bmolx,bmoly,bmolz;
 
    h    = HAMILTONIAN(i); /* Speicher fuer GesamtHamiltonian holen */
    s    = PKQ( i);        /* (Jn|Pkq(J)|mJ) holen */
    mag  = HMAG(i);        /* (Jn|Hmag|mJ holen  */
 
    myB  = EINHEITIMP[ EINHEITNRIN(i) ].myB;
 
    /* Speicher des Gesamthamiltonians h loeschen */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) = 0.0;
              I(h,n,m) = 0.0;
         }
 
 
    /* h += h             */
    /*       kristallfeld */
    for( n=DIMJ(i) ; n>=1 ; --n )
       for( m=DIMJ(i) ; m>=1 ; --m ){
         if(VOR44(i)==1.0) R(h,n,m)+=RV(v,1)*R(O4P05(s),n,m)/N_O4P05(s);
         else              I(h,n,m)+=RV(v,1)*R(O4P05(s),n,m)/N_O4P05(s);
         if(VOR64(i)==1.0) R(h,n,m)+=RV(v,2)*R(O6M21(s),n,m)/N_O6M21(s);
         else              I(h,n,m)+=RV(v,2)*R(O6M21(s),n,m)/N_O6M21(s);
         if(VOR44(i)==1.0) R(h,n,m)+=RV(v,3)*R(O4M05(s),n,m)/N_O4M05(s);
         else              I(h,n,m)+=RV(v,3)*R(O4M05(s),n,m)/N_O4M05(s);
         if(VOR64(i)==1.0) R(h,n,m)+=RV(v,4)*R(O6P21(s),n,m)/N_O6P21(s);
         else              I(h,n,m)+=RV(v,4)*R(O6P21(s),n,m)/N_O6P21(s);
 
         R(h,n,m) += RV(v,5)*R( O2P0(s),n,m )/N_O2P0(s);
 
         if(VOR22(i)==1.0) R(h,n,m)+=RV(v,6)*R(O2P2(s),n,m)/N_O2P2(s);
         else              I(h,n,m)+=RV(v,6)*R(O2P2(s),n,m)/N_O2P2(s);
         if(VOR42(i)==1.0) R(h,n,m)+=RV(v,7)*R(O4P2(s),n,m)/N_O4P2(s);
         else              I(h,n,m)+=RV(v,7)*R(O4P2(s),n,m)/N_O4P2(s);
         if(VOR62(i)==1.0) R(h,n,m)+=RV(v,8)*R(O6P2(s),n,m)/N_O6P2(s);
         else              I(h,n,m)+=RV(v,8)*R(O6P2(s),n,m)/N_O6P2(s);
         if(VOR66(i)==1.0) R(h,n,m)+=RV(v,9)*R(O6P6(s),n,m)/N_O6P6(s);
         else              I(h,n,m)+=RV(v,9)*R(O6P6(s),n,m)/N_O6P6(s);
       }
 
    /* Speicher von Hmag loeschen */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(mag,n,m) = 0.0;
              I(mag,n,m) = 0.0;
         }
 
    bmolx   = RV(v,10) * sin( RV(v,12) ) * cos( RV(v,11) );
    bmoly   = RV(v,10) * sin( RV(v,12) ) * sin( RV(v,11) );
    bmolz   = RV(v,10) * cos( RV(v,12) );
 
    mag  = calcBmol( DIMJ(i) , mag     , 2.0*(GJ(i)-1.0), myB,
                     bmolx   , bmoly   , bmolz    );
 
    /* h += mag           */
    /*       Magnetfeld   */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) += R(mag,n,m);
              I(h,n,m) += I(mag,n,m);
         }
 
 
    return( i );
}
/*----------------------------------------------------------------------------
                                 op_norm()
----------------------------------------------------------------------------*/
void op_norm( stev,macheps )
    STEVENS  *stev;
    DOUBLE   macheps;
{
    DOUBLE norm_operator();
 
    N_O2P0( stev )  = norm_operator( O2P0(stev) ,macheps );
    N_O2P2( stev )  = norm_operator( O2P2(stev) ,macheps );
    N_O4P2( stev )  = norm_operator( O4P2(stev) ,macheps );
    N_O6P2( stev )  = norm_operator( O6P2(stev) ,macheps );
    N_O6P6( stev )  = norm_operator( O6P6(stev) ,macheps );
 
    N_O4P05( stev ) = norm_operator( O4P05(stev) ,macheps );
    N_O4M05( stev ) = norm_operator( O4M05(stev) ,macheps );
 
    N_O6P21( stev ) = norm_operator( O6P21(stev) ,macheps );
    N_O6M21( stev ) = norm_operator( O6M21(stev) ,macheps );
 
 
}
/*----------------------------------------------------------------------------
                                 calc_Okq()
----------------------------------------------------------------------------*/
STEVENS *calc_Okq( dimj,i )
    INT dimj;
    ITERATION *i;
{
    STEVENS *stevens;
    MATRIX  *stevks(),*stevkq();
 
    stevens = STEVENS_ALLOC(1);
    DIMJ(stevens) = dimj;
 
    O2P0(stevens)  = stevkq(2,0,dimj);
    if( VOR22(i)==1.0 ) O2P2(stevens)  =  stevkq(2, 2,dimj);
    else                O2P2(stevens)  =  stevkq(2,-2,dimj);
    if( VOR42(i)==1.0 ) O4P2(stevens)  =  stevkq(4, 2,dimj);
    else                O4P2(stevens)  =  stevkq(4,-2,dimj);
    if( VOR62(i)==1.0 ) O6P2(stevens)  =  stevkq(6, 2,dimj);
    else                O6P2(stevens)  =  stevkq(6,-2,dimj);
    if( VOR66(i)==1.0 ) O6P6(stevens)  =  stevkq(6, 6,dimj);
    else                O6P6(stevens)  =  stevkq(6,-6,dimj);
    O4P05(stevens) = stevks(4,  5,dimj,i);
    O4M05(stevens) = stevks(4, -5,dimj,i);
    O6P21(stevens) = stevks(6, 21,dimj,i);
    O6M21(stevens) = stevks(6,-21,dimj,i);
 
 
    return( stevens );
}
/*----------------------------------------------------------------------------
                                 Bkq_v()
------------------------------------------------------------------------------
*******************        1   2   3   4   5   6   7   8   9  t
* Transformation: *  Bkq=(B20,B22,B40,B42,B44,B60,B62,B64,B66)
*******************                      t
                        -> (v_1,....,v_9)
-----------------------------------------------------------------------------*/
VEKTOR *Bkq_v(Bkq,stev,ionennr,macheps)
    VEKTOR    *Bkq;
    STEVENS   *stev;
    INT       ionennr;
    DOUBLE    macheps;
{
   UNUSED_PARAMETER(ionennr);
   UNUSED_PARAMETER(macheps);

   VEKTOR *v,*vr_alloc();
   DOUBLE B20,B22,B40,B42,B44,B60,B62,B64,B66;
   DOUBLE B1mol,B2mol,B3mol;
 
   B20   =  RV(Bkq,1);
   B22   =  RV(Bkq,2);
   B40   =  RV(Bkq,3);
   B42   =  RV(Bkq,4);
   B44   =  RV(Bkq,5);
   B60   =  RV(Bkq,6);
   B62   =  RV(Bkq,7);
   B64   =  RV(Bkq,8);
   B66   =  RV(Bkq,9);
   B1mol =  RV(Bkq,10);
   B2mol =  RV(Bkq,11);
   B3mol =  RV(Bkq,12);
 
   v   =  vr_alloc( VRDIM(Bkq) );
 
   RV(v, 1) = 1.0/10.0*( 5.0*B40+B44 )*N_O4P05(stev);
   RV(v, 2) = 1.0/42.0*(21.0*B60-B64 )*N_O6M21(stev);
   RV(v, 3) = 1.0/10.0*( 5.0*B40-B44 )*N_O4M05(stev);
   RV(v, 4) = 1.0/42.0*(21.0*B60+B64 )*N_O6P21(stev);
 
   RV(v, 5) = B20*N_O2P0(stev);
   RV(v, 6) = B22*N_O2P2(stev);
   RV(v, 7) = B42*N_O4P2(stev);
   RV(v, 8) = B62*N_O6P2(stev);
   RV(v, 9) = B66*N_O6P6(stev);
 
   RV(v,10) = B1mol;
   RV(v,11) = B2mol;
   RV(v,12) = B3mol;
 
   return( v );
}
/*----------------------------------------------------------------------------
                                 Bkq_Vkq()
------------------------------------------------------------------------------
*******************        1   2   3   4   5   6   7   8   9  t
* Transformation: *  Bkq=(B20,B22,B40,B42,B44,B60,B62,B64,B66)
*******************
                        -> Vkq
-----------------------------------------------------------------------------*/
VEKTOR *Bkq_Vkq(Bkq,stev,ionennr,macheps)
    VEKTOR    *Bkq;
    STEVENS   *stev;
    INT       ionennr;
    DOUBLE    macheps;
{
   UNUSED_PARAMETER(stev);
   UNUSED_PARAMETER(ionennr);
   UNUSED_PARAMETER(macheps);

   VEKTOR *Vkq,*vr_alloc();
 
   Vkq =  vr_alloc( VRDIM(Bkq) );
 
     /* Bkq auf Vkq umrechnen       */
     /*                             */
     /*       | Bk0  ,         q =0 */
     /* Vkq = |                     */
     /*       | Bkq/2/omegakq, q >0 */
     /*                             */
 
 
     RV( Vkq,1 )  = RV(Bkq,1);
     RV( Vkq,2 )  = RV(Bkq,2) / 2 / omegan0n(2);
 
     RV( Vkq,3 )  = RV(Bkq,3);
     RV( Vkq,4 )  = RV(Bkq,4) / 2 / omegan2n(2);
     RV( Vkq,5 )  = RV(Bkq,5) / 2 / omegan0n(4);
 
     RV( Vkq,6 )  = RV(Bkq,6);
     RV( Vkq,7 )  = RV(Bkq,7) / 2 / omegan4n(2);
     RV( Vkq,8 )  = RV(Bkq,8) / 2 / omegan2n(4);
     RV( Vkq,9 )  = RV(Bkq,9) / 2 / omegan0n(6);
 
     RV( Vkq,10)  = RV(Bkq,10);
     RV( Vkq,11)  = RV(Bkq,11);
     RV( Vkq,12)  = RV(Bkq,12);
 
     return( Vkq );
}
/*----------------------------------------------------------------------------
                                 Vkq_Bkq()
------------------------------------------------------------------------------
*******************        1   2   3   4   5   6   7   8   9  t
* Transformation: *  Vkq=(V20,V22,V40,V42,V44,V60,V62,V64,V66)
*******************
                        -> Bkq
-----------------------------------------------------------------------------*/
VEKTOR *Vkq_Bkq(Vkq,stev,ionennr,macheps)
    VEKTOR    *Vkq;
    STEVENS   *stev;
    INT       ionennr;
    DOUBLE    macheps;
{
   UNUSED_PARAMETER(stev);
   UNUSED_PARAMETER(ionennr);
   UNUSED_PARAMETER(macheps);

   VEKTOR *Bkq,*vr_alloc();
 
   Bkq =  vr_alloc( VRDIM(Vkq) );
 
     /* Bkq auf Vkq umrechnen       */
     /*                             */
     /*       | Vk0  ,         q =0 */
     /* Bkq = |                     */
     /*       | Vkq*2*omegakq, q >0 */
     /*                             */
 
 
     RV( Bkq,1 )  = RV(Vkq,1);
     RV( Bkq,2 )  = RV(Vkq,2) * 2 * omegan0n(2);
 
     RV( Bkq,3 )  = RV(Vkq,3);
     RV( Bkq,4 )  = RV(Vkq,4) * 2 * omegan2n(2);
     RV( Bkq,5 )  = RV(Vkq,5) * 2 * omegan0n(4);
 
     RV( Bkq,6 )  = RV(Vkq,6);
     RV( Bkq,7 )  = RV(Vkq,7) * 2 * omegan4n(2);
     RV( Bkq,8 )  = RV(Vkq,8) * 2 * omegan2n(4);
     RV( Bkq,9 )  = RV(Vkq,9) * 2 * omegan0n(6);
 
     RV( Bkq,10)  = RV(Vkq,10);
     RV( Bkq,11)  = RV(Vkq,11);
     RV( Bkq,12)  = RV(Vkq,12);
 
     return( Bkq );
}
/*----------------------------------------------------------------------------
                                 v_Bkq()
------------------------------------------------------------------------------
*******************                   t
* Transformation: *  v = (v_1,...,v_9) -> Bkq
*******************
-----------------------------------------------------------------------------*/
VEKTOR *v_Bkq(v,stev,ionennr,macheps)
    VEKTOR    *v;
    STEVENS   *stev;
    INT       ionennr;
    DOUBLE    macheps;
{
   UNUSED_PARAMETER(ionennr);
   UNUSED_PARAMETER(macheps);

   VEKTOR *Bkq,*vr_alloc();
   DOUBLE B20,B22,B40,B42,B44,B60,B62,B64,B66;
 
 
   B20 = RV(v,5)/N_O2P0(stev);
   B22 = RV(v,6)/N_O2P2(stev);
   B42 = RV(v,7)/N_O4P2(stev);
   B62 = RV(v,8)/N_O6P2(stev);
   B66 = RV(v,9)/N_O6P6(stev);
 
   B40 =        RV(v,1)/N_O4P05(stev) + RV(v,3)/N_O4M05(stev);
   B44 =  5.0*( RV(v,1)/N_O4P05(stev) - RV(v,3)/N_O4M05(stev) );
 
   B60 =        RV(v,2)/N_O6M21(stev) + RV(v,4)/N_O6P21(stev);
   B64 =-21.0*( RV(v,2)/N_O6M21(stev) - RV(v,4)/N_O6P21(stev) );
 
   Bkq       = vr_alloc( VRDIM(v) );
   RV(Bkq,1) = B20;
   RV(Bkq,2) = B22;
   RV(Bkq,3) = B40;
   RV(Bkq,4) = B42;
   RV(Bkq,5) = B44;
   RV(Bkq,6) = B60;
   RV(Bkq,7) = B62;
   RV(Bkq,8) = B64;
   RV(Bkq,9) = B66;
 
   RV(Bkq,10)= RV(v,10);
   RV(Bkq,11)= RV(v,11);
   RV(Bkq,12)= RV(v,12);
 
   return(Bkq);
}
/*----------------------------------------------------------------------------
                                 rphi_v()
------------------------------------------------------------------------------
*******************                               t                    t
* Transformation: *  rphi=( r,phi_1,...,phi_(n-1) ) -> v=( v_1,...,v_n )
*******************
-----------------------------------------------------------------------------*/
VEKTOR *rphi_v(rphi,stev,ionennr,macheps)
    VEKTOR    *rphi;
    STEVENS   *stev;
    INT       ionennr;
    DOUBLE    macheps;
{
    UNUSED_PARAMETER(stev);
    UNUSED_PARAMETER(ionennr);
    UNUSED_PARAMETER(macheps);

    VEKTOR *v,*vr_alloc();
    INT    zeile,j,n;
 
    n      = VRDIM(rphi);
    v      = vr_alloc( n );
 
 
    for( zeile=1; zeile<=n-3; ++zeile){
         RV(v,zeile) = RV(rphi,1);
 
         for( j=n; j>=zeile+1; --j)
            RV(v,zeile) *= cos( RV(rphi,j) );
 
         if( zeile > 1)
            RV(v,zeile) *= sin( RV(rphi,zeile) );
    }
    RV(v,n-2) = RV(rphi,n-2);
    RV(v,n-1) = RV(rphi,n-1);
    RV(v,n-0) = RV(rphi,n-0);
 
 
    return(v);
}
 
 
 
/*----------------------------------------------------------------------------
                                 abfrage()
------------------------------------------------------------------------------
*********************
* Quadrantenabfrage * fuer  r aus IR  und  phi_i aus [0,pi[
*********************
-----------------------------------------------------------------------------*/
QUADRANT *abfrage( rsin, rcos, macheps  )
    DOUBLE rsin,rcos,macheps;
{
    QUADRANT *q;
 
    q = Q_ALLOC(1);
 
    if( is_equal( rsin, 0.0, macheps) ){
        WINKEL(q) = 0.0;
        RADIUS(q) = rcos;
    }
    else{
          WINKEL(q) = acos( DSIGN(rsin)*rcos/sqrt(rsin*rsin+rcos*rcos) );
          RADIUS(q) = rsin / sin( WINKEL(q) );
        }
 
    return(q);
}
 
 
/*----------------------------------------------------------------------------
                                 v_rphi()
------------------------------------------------------------------------------
*******************                t                                  t
* Transformation: *  v=(v_1,...,v_n) -> rphi=( r,phi_1,...,phi_(n-1) )
*******************
-----------------------------------------------------------------------------*/
VEKTOR *v_rphi(v,stev,ionennr,macheps)
    VEKTOR    *v;
    STEVENS   *stev;
    INT       ionennr;
    DOUBLE    macheps;
{
    UNUSED_PARAMETER(stev);
    UNUSED_PARAMETER(ionennr);
 
    QUADRANT *q,*abfrage();
    VEKTOR   *rphi,*vr_alloc();
    DOUBLE   r;
    INT      zeile,n;
 
    n      = VRDIM( v );
    rphi   = vr_alloc( n );
 
    r      = RV(v,1);
    for( zeile=2; zeile<=n-3; ++zeile){
         q              = abfrage( RV(v,zeile), r, macheps );
         RV(rphi,zeile) = WINKEL(q);
         r              = RADIUS(q);
         free_(q);
    }
    RV(rphi,1) = r;
 
    RV(rphi,n-2) = RV(v,n-2);
    RV(rphi,n-1) = RV(v,n-1);
    RV(rphi,n-0) = RV(v,n-0);
 
    return(rphi);
}
 
 
/*----------------------------------------------------------------------------
                                 xW_v()
------------------------------------------------------------------------------
*******************            t                    t
* Transformation: *  xW=( x,W )   ->   v=( v_1,v_2 )
*******************
-----------------------------------------------------------------------------*/
VEKTOR *xW_v(xW,stev,ionennr,macheps)
    VEKTOR    *xW;
    STEVENS   *stev;
    INT       ionennr;
    DOUBLE    macheps;
{
    UNUSED_PARAMETER(macheps);

    VEKTOR  *v,*vr_alloc();
    DOUBLE  x,W;
 
    v = vr_alloc( VRDIM(xW) );
    x = RV(xW,1);
    W = RV(xW,2);
 
    RV(v,1) = W *        x  *N_O4P05( stev )/F4(ionennr);
    RV(v,2) = W *(1-ABSD(x))*N_O6M21( stev )/F6(ionennr);
 
    return( v );
}
 
/*----------------------------------------------------------------------------
                                 v_xW()
------------------------------------------------------------------------------
*******************                t             t
* Transformation: *  v =( v_1,v_2 ) -> xW = (x,W)
*******************
-----------------------------------------------------------------------------*/
VEKTOR *v_xW(v,stev,ionennr,macheps)
    VEKTOR    *v;
    STEVENS   *stev;
    INT       ionennr;
    DOUBLE    macheps;
{
    VEKTOR  *xW,*vr_alloc();
    DOUBLE  x,W,v1,v2,f4,f6,d;
    INT     isv1null,isv2null;
 
    xW = vr_alloc( VRDIM(v) );
    v1 = RV(v,1);
    v2 = RV(v,2);
    f4 = F4( ionennr );
    f6 = F6( ionennr );
    isv1null = is_equal( v1, 0.0, macheps );
    isv2null = is_equal( v2, 0.0, macheps );
 
    if( isv1null==JA && isv2null==JA  ){
        W = 0.0;
        x = 1.0;  /* x ist eigentlich beliebig */
    }
    if( isv2null==JA && isv1null==NEIN){
        x = 1.0;                   /* auch x = -1 moeglich */
        W = x * v1 * f4 / N_O4P05( stev );
    }
    if( isv2null==NEIN ){
        d = v1/v2 * f4/f6 * N_O6M21( stev )/N_O4P05( stev );
 
        if( is_equal(d, 0.0, macheps) ){
            x = 0.0;
            W = v2 / N_O6M21( stev ) * f6;
        }
        else{
              x = d  / ( 1.0 + ABSD(d) );
              W = v1 / N_O4P05( stev ) * f4 / x;
            }
    }
 
    RV(xW,1) = x;
    RV(xW,2) = W;
    return( xW );
}
 
/*------------------------------------------------------------------------------
                               strategie()
------------------------------------------------------------------------------*/
void strategie(ewproblem,iteration,p0)
    EWPROBLEM *ewproblem;
    ITERATION *iteration;
    VEKTOR    *p0;
{
    UNUSED_PARAMETER(ewproblem);
    UNUSED_PARAMETER(iteration);
    UNUSED_PARAMETER(p0);
/*
   if( IS_MATRIXELEMENT(iteration) || IS_INTENSITAET(iteration)
                                   || IS_POSFIT(iteration)       )
       return(vertausche_niveaus(ewproblem,iteration,p0));
 */
}
/*------------------------------------------------------------------------------
                               vertausche_niveaus()
------------------------------------------------------------------------------*/
INT vertausche_niveaus(ewproblem,iteration,p0)
    EWPROBLEM *ewproblem;
    ITERATION *iteration;
    VEKTOR    *p0;
{
  UNUSED_PARAMETER(p0);

  NEBENBEDINGUNG *neben;
  MATRIX  *entartung/*, *ev*/,*aJtb_2(), *aJtb2, *exp;
  MATRIX  *mx_alloc(), *mx_Chi2;
  VEKTOR/**ew, *fix, *v0, *neu,*/ *vr_alloc()/*,*v*/, *v1_v2();
  DOUBLE  macheps/*, shift, temperatur,chi2_e,chi2_t*/,mat_Chi2();
/*CHAR    *ionname;
  INT     ionennr, einheitnr_in, einheitnr_out, dum; */
  INT /*  *gi,*/ anz_niveaus, i_e, i_t, flag/*, i*/, ndim;
/*STEVENS *stevens; */
 
 
 
  entartung     = ewproblem->entartung;
/*gi            = ewproblem->gi;
  ev            = ewproblem->eigenvektoren;
  ew            = ewproblem->eigenwerte;
  shift         = ewproblem->shift; */
  macheps       = ewproblem->eps_machine;
  anz_niveaus   = ANZ_ZE(entartung);
  aJtb2         = aJtb_2(ewproblem,macheps );
/*einheitnr_out = EINHEITNROUT(  iteration );
  einheitnr_in  = EINHEITNRIN(   iteration ); */
 
/*ionname       = IONNAME( iteration );
  ionennr       = isimplementiert( ionname );
  temperatur    = TEMPERATUR(    iteration ); */
  neben         = NEBENBEDINGUNG(iteration);
/*stevens       = PKQ(           iteration);
  v0            = V0(            iteration );
  fix           = FIX(    neben     );*//* z.B. (0,1,0,1,0,1,1,1,0) */
  exp           = INTENSITAETEN(neben);
 
  ndim = MIN( MXDIM(exp), anz_niveaus );
 
  mx_Chi2 = mx_alloc( ndim,ndim ); /*i_t-> 1 ... ndim (i_e) */
  for(i_t=1; i_t<=ndim; ++i_t )
      for(i_e=1; i_e<=ndim; ++i_e )
          R(mx_Chi2,i_t,i_e)=mat_Chi2(ewproblem,iteration,aJtb2,exp,i_t,i_e);
 
  for(i_t=1; i_t<=ndim; ++i_t )
      for(i_e=1; i_e<=ndim; ++i_e )
          printf("mx(%d,%d) = %f\n",i_t,i_e,R(mx_Chi2,i_t,i_e));
 
 
 
 
 
  free_mx(mx_Chi2);
//free_vr(neu);
  free_mx(aJtb2);
 
  flag=0;
  if( flag==0 ) return( NEIN );
  else          return( JA   );
}
/*----------------------------------------------------------------------------
                                  mat_Chi2()
----------------------------------------------------------------------------*/
DOUBLE   mat_Chi2(ewproblem,iteration,aJtb2,intensit_exp,i_t,i_e)
  EWPROBLEM *ewproblem;
  ITERATION *iteration;
  MATRIX    *aJtb2, *intensit_exp;
  INT       i_t,i_e;
{
   INT  /* *gi,*/anz_niveaus,einheitnr_in/*,einheitnr_out*/;
   MATRIX  *entartung/*,*ev*/,*aJtb_2()/*,*d_intensit_exp*/;
   DOUBLE/*shift,*/macheps,temperatur/*,e_calc,e_exp*/,gj;
   DOUBLE  zu_summe,zustandssumme(),faktor,energie,i_calc,i_exp/*,chi2_1*/;
   VEKTOR  *ew;
 
   CHAR           *ionname;
   INT         /* i,*/k,ionennr/*,diff*/;
   DOUBLE         chi2/*,norm,normi*/;
 
 
   entartung     = ewproblem->entartung;
/* gi            = ewproblem->gi;
   ev            = ewproblem->eigenvektoren; */
   ew            = ewproblem->eigenwerte;
/* shift         = ewproblem->shift; */
   macheps       = ewproblem->eps_machine;
   anz_niveaus   = ANZ_ZE(entartung);
/* einheitnr_out = EINHEITNROUT(  iteration ); */
   einheitnr_in  = EINHEITNRIN(   iteration );
   temperatur    = TEMPERATUR(    iteration );
   ionname       = IONNAME( iteration );
   ionennr       = isimplementiert( ionname );
   chi2          = 0.0;
 
 
if( IS_INTENSITAET(iteration) ){
   /*    Chi2  += Chi2_Intensitaeten   */
   zu_summe  = zustandssumme( einheitnr_in , ew , temperatur );
   gj        = IONENIMP[ ionennr ].gj;
 
       for( k=1; k<= MXDIM(intensit_exp); ++k )
           if( i_t<=anz_niveaus&& i_e<=anz_niveaus && k<=anz_niveaus){
               i_exp     = R(intensit_exp,i_e,k);
               if( !is_equal(i_exp,0.0,macheps) ){
                     faktor    = EINHEITIMP[ einheitnr_in].fek;
                     energie   = RV( ew , (INT)R(entartung,i_e,1) );
                     i_calc    = R(aJtb2,i_t,k) * const*gj*gj
                               *exp_(-energie*faktor/temperatur)/zu_summe;
                     chi2     += (i_exp-i_calc)*(i_exp-i_calc);
               }
            }
}
else if( IS_MATRIXELEMENT(iteration) ){
   /*    Chi2  += Chi2_Matrixelemente  */
       for( k=1; k<= MXDIM(intensit_exp); ++k )
           if( i_e<= anz_niveaus && i_t<=anz_niveaus && k<= anz_niveaus){
               i_exp     = R(intensit_exp,i_e,k);
               if( !is_equal(i_exp,0.0,macheps) ){
                     i_calc    = R(aJtb2,i_t,k);
                     chi2     += (i_exp-i_calc)*(i_exp-i_calc);
               }
            }
}
 
 
 
   return(chi2);
}
/*------------------------------------------------------------------------------
                               menue()
------------------------------------------------------------------------------*/
INT menue(amoeba,setup,ewproblem,iteration,w,iter_steps,chi2)
    MINIMUM   *amoeba;
    SETUP     *setup;
    EWPROBLEM *ewproblem;
    ITERATION *iteration;
    VEKTOR    *w;
    INT       iter_steps;
    DOUBLE    chi2;
{
  UNUSED_PARAMETER(amoeba); 
  #define DIMFIX 12   /* VRDIM(fix) */
 
  CHAR **bild/*, *bildzeile*/,*ionname;
  INT  ionennr;
  INT  einheitnr_in, einheitnr_out;
  INT  anz_zeilen  = 25,zeile;
/*INT  anz_spalten = 80,spalte; */
 
  DOUBLE  /*r,*/x,W/*,shift*/,macheps,gj,faktor,energie,zu_summe;
  DOUBLE  phi[DIMFIX+1];
  DOUBLE  e[18],pe[18],i01[18],i02[18],e_calc,e_exp;
  DOUBLE  p1[18],p2[18],i_calc,i_exp,temperatur;
  CHAR    c[DIMFIX+1];
  INT     g[18],i/*,ze*/,sp/*,i_v,i_v0,flag*/;
  INT     *gi,anz_niveaus;
  MATRIX  *entartung/*,*ev*/,*aJtb_2(),*aJtb2;
  VEKTOR  *ew/*,*v0*/,*fix,*v1_v2(),*v;
  CHAR    *einheit;
  NEBENBEDINGUNG *neben;
  VEKTOR  *ew_exp,*rphi,*xW,*v_rphi(),*v_xW();
  MATRIX  *intensit_exp;
  STEVENS *stevens;
  DOUBLE  b1,b2,b3;
 
  b1 = B1(iteration);
  b2 = B2(iteration);
  b3 = B3(iteration);
  B1(iteration)=0.0;
  B2(iteration)=0.0;
  B3(iteration)=0.0;
 
  v         = v1_v2(    iteration, w );
  iteration = hamilton( iteration, v );
  ewproblem = diagonalisiere( ewproblem,HAMILTONIAN(iteration),
                                NOSPACE,setup);
 
 
  B1(iteration)=b1;
  B2(iteration)=b2;
  B3(iteration)=b3;
 
  entartung     = ewproblem->entartung;
  gi            = ewproblem->gi;
/*ev            = ewproblem->eigenvektoren;*/
  ew            = ewproblem->eigenwerte;
/*shift         = ewproblem->shift; */
  macheps       = ewproblem->eps_machine;
  anz_niveaus   = ANZ_ZE(entartung);
  aJtb2         = aJtb_2(ewproblem,macheps );
  einheitnr_out = EINHEITNROUT(  iteration );
  einheitnr_in  = EINHEITNRIN(   iteration );
  einheit       = EINHEITIMP[ einheitnr_out ].einheit;
 
  ionname       = IONNAME( iteration );
  ionennr       = isimplementiert( ionname );
  temperatur    = TEMPERATUR(    iteration );
  neben         = NEBENBEDINGUNG(iteration);
  stevens       = PKQ(           iteration);
/*v0            = V0(            iteration ); */
  fix           = FIX(    neben     );  /* z.B. (0,1,0,1,0,1,1,1,0) */
 
 
  xW    =  v_xW(  v,stevens,ionennr,macheps);
  x     =  is_null(RV(xW,1),0.00005);
  W     =  is_null(RV(xW,2),0.00005);
  rphi  =  v_rphi(v,stevens,ionennr,macheps);
 
  RV(rphi,11) *= 180.0/pi;  /* phi auf   grad umrechnen */
  RV(rphi,12) *= 180.0/pi;  /* theta auf grad umrechnen */
 
  fix  = FIX( neben );
  for(i=1; i<=VRDIM(fix); ++i)
      if( RV(fix,i) == 0.0 ) c[i]   = ' ';
      else                   c[i]   = '*';
 
 
  phi[1] = is_null(RV(rphi,1),0.00005);
  for(i=2; i<=VRDIM(rphi)-3; ++i)
     phi[i] = is_null(RV(rphi,i)*180.0/pi,0.00005);
  for(i=VRDIM(rphi)-2; i<=VRDIM(rphi); ++i)
     phi[i] = is_null(RV(rphi,i),0.00005);
 
  for(i=1; i<=anz_niveaus; ++i)
     g[i] = VALUE(gi,i);
 
 
  for( i=1 ; i<=anz_niveaus ; ++i){
       e[i] = RV(ew,(INT)R(entartung,i,1))
       *EINHEITIMP[einheitnr_in].fek*EINHEITIMP[einheitnr_out].fke;
       e[i] = is_null(e[i],0.005);
 
   }
 
 
 
 zu_summe    = zustandssumme( einheitnr_in , ew , temperatur );
 gj          = IONENIMP[ ionennr ].gj;
if( IS_INTENSITAET(iteration) || IS_POSFIT(iteration)){
 /* Intensitaet noch mit gi_ze*const*gj**2 * exp( -Ei/T)/Z multipizieren */
    faktor      = EINHEITIMP[ einheitnr_in].fek;
    energie     = RV( ew , (INT)R(entartung,1,1) );
      for( sp=1 ; sp<=anz_niveaus ; ++sp ){
        i01[sp] = R(aJtb2,1,sp) * const*gj*gj
                           *exp_(-energie*faktor/temperatur)/zu_summe;
        i01[sp] = is_null(i01[sp],0.05);
 
        }
if( anz_niveaus >= 2){
 /* Intensitaet noch mit gi_ze*const*gj**2 * exp( -Ei/T)/Z multipizieren */
    faktor      = EINHEITIMP[ einheitnr_in].fek;
    energie     = RV( ew , (INT)R(entartung,2,1) );
      for( sp=1 ; sp<=anz_niveaus ; ++sp ){
        i02[sp] = R(aJtb2,2,sp) * const*gj*gj
                           *exp_(-energie*faktor/temperatur)/zu_summe;
        i02[sp] = is_null(i02[sp],0.05);
      }
}
}/* end if intensitaet */
else if( IS_MATRIXELEMENT(iteration) || IS_EIGENWERT(iteration)
         || IS_SUSFIT(iteration)     || IS_MAGFIT(iteration) ){
      for( sp=1 ; sp<=anz_niveaus ; ++sp ){
           i01[sp] = R(aJtb2,1,sp);
           i01[sp] = is_null(i01[sp],0.05);
      }
if( anz_niveaus >= 2){
      for( sp=1 ; sp<=anz_niveaus ; ++sp ){
           i02[sp] = R(aJtb2,2,sp);
           i02[sp] = is_null(i02[sp],0.05);
      }
}
}/* end if intensitaet */
 
 
  ew_exp  = EIGENWERTE( neben );
  for(i=1; i<=anz_niveaus; ++i){
   e_calc = RV(ew,(INT)R(entartung,i,1))
            *EINHEITIMP[einheitnr_in].fek*EINHEITIMP[einheitnr_out].fke;
   e_exp  = RV(ew_exp,i);
   if( is_equal(e_exp,0.0,macheps) )
      pe[i] = -1.0;
   else pe[i] = ABSD(100.0*(e_calc - e_exp)/e_exp);
 
   if(i==1) pe[i] = 0.0;
   pe[i] = is_null(pe[i],0.5);
   if( pe[i]>9999.0 ) pe[i] = 9999.0;
 
  }
 
 
 
  intensit_exp = INTENSITAETEN( neben );
  for(i=1; i<=anz_niveaus; ++i){
   i_calc = i01[i];
   i_exp  = R(intensit_exp,1,i);
   if( is_equal(i_exp,0.0,macheps) )
      p1[i] = -1.0;
   else p1[i] = ABSD(100.0*(i_calc - i_exp)/i_exp);
   p1[i] = is_null(p1[i],0.5);
   if( p1[i]>9999.0 ) p1[i] = 9999.0;
  }
 
 
if( anz_niveaus >= 2){
  for(i=1; i<=anz_niveaus; ++i){
   i_calc = i02[i];
   i_exp  = R(intensit_exp,2,i);
   if( is_equal(i_exp ,0.0,macheps) )
      p2[i] = -1.0;
   else p2[i] = ABSD(100.0*(i_calc - i_exp)/i_exp);
   p2[i] = is_null(p2[i],0.5);
   if( p2[i]>9999.0 ) p2[i] = 9999.0;
  }
 
}
 
 
bild = STRINGP_ALLOC( anz_zeilen );
 
VALUE(bild, 1) =
" -----------------------------------------------------------------------\n";
if( IS_INTENSITAET(iteration) || IS_POSFIT(iteration) )
VALUE(bild, 2) =
"|Temperatur: %9.2fK|  |  |Energie  | |dE||| 1->i | |dI|| 2->i | |dI||\n";
else if( IS_MATRIXELEMENT(iteration) || IS_EIGENWERT(iteration)
         || IS_SUSFIT(iteration)     || IS_MAGFIT(iteration)  )
VALUE(bild, 2) =
"|Temperatur: %9.2fK|  |  |Energie  | |dE||| 1->i | |dM|| 2->i | |dM||\n";
 
if( IS_INTENSITAET(iteration) || IS_POSFIT(iteration) )
VALUE(bild, 3) =
"|Iterationsschritt:%4d| i|gi|in %6s|in %% ||in b  |in %% |in b  |in %% |\n"
;
else if( IS_MATRIXELEMENT(iteration) || IS_EIGENWERT(iteration)
         || IS_SUSFIT(iteration)     || IS_MAGFIT(iteration)  )
VALUE(bild, 3) =
"|Iterationsschritt:%4d| i|gi|in %6s|in %% ||      |in %% |      |in %% |\n"
;
 
VALUE(bild, 4) =
"|----------------------|---------------------||------------|------------|\n";
 
if( anz_niveaus >= 2 )
VALUE(bild, 5) =
"|  r | %13.4f |%c| 1 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
else
VALUE(bild, 5) =
"|  r | %13.4f |%c| 1 %2d %9.2f %4.0f || %5.1f %4.0f |            |\n";
 
if( anz_niveaus >= 2 ){
VALUE(bild, 6) =
"|Phi1| %13.4f |%c| 2 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild, 6) =
"|Phi1| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 3 ){
VALUE(bild, 7) =
"|Phi2| %13.4f |%c| 3 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild, 7) =
"|Phi2| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 4 ){
VALUE(bild, 8) =
"|Phi3| %13.4f |%c| 4 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild, 8) =
"|Phi3| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 5 ){
VALUE(bild, 9) =
"|Phi4| %13.4f |%c| 5 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild, 9) =
"|Phi4| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 6 ){
VALUE(bild,10) =
"|Phi5| %13.4f |%c| 6 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,10) =
"|Phi5| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 7 ){
VALUE(bild,11) =
"|Phi6| %13.4f |%c| 7 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,11) =
"|Phi6| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 8 ){
VALUE(bild,12) =
"|Phi7| %13.4f |%c| 8 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,12) =
"|Phi7| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 9 ){
VALUE(bild,13) =
"|Phi8| %13.4f |%c| 9 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,13) =
"|Phi8| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 10){
VALUE(bild,14) =
"|B  m| %13.4f |%c|10 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,14) =
"|B  m| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 11){
VALUE(bild,15) =
"|Phio| %13.4f |%c|11 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,15) =
"|Phio| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 12){
VALUE(bild,16) =
"|Ttal| %13.4f |%c|12 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,16) =
"|Ttal| %13.4f |%c|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 13){
VALUE(bild,17) =
"|----------------------|13 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,17) =
"|----------------------|                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 14){
VALUE(bild,18) =
"|kubischer Startpunkt: |14 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,18) =
"|kubischer Startpunkt: |                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 15){
VALUE(bild,19) =
"|(W in %6s)         |15 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,19) =
"|(W in %6s)         |                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 16){
VALUE(bild,20) =
"|                      |16 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,20) =
"|                      |                     ||            |            |\n";
}
 
 
if( anz_niveaus >= 17){
VALUE(bild,21) =
"|  x = %13.4f   |17 %2d %9.2f %4.0f || %5.1f %4.0f | %5.1f %4.0f |\n";
}
else{
VALUE(bild,21) =
"|  x = %13.4f   |                     ||            |            |\n";
}
 
 
VALUE(bild,22) =
"|  W = %13.4f   |------------------------------------------------|\n";
 
 
VALUE(bild,23) =
"|----------------------|sus_a %c|sus_b %c|sus_c %c|sus_p %c|ew %c|int %c|mat %c|\n";
 
 
VALUE(bild,24) =
"|Chi2| %16.4e|mag_a %c|mag_b %c|mag_c %c|mag_B %c|posfit %c        |\n";
 
VALUE(bild,25) =
" ----------------------------------------------------------------------- \n";
 
 
clearscreen;
for(zeile=1; zeile<= 2; ++zeile)
  printf(VALUE(bild,zeile),temperatur);
 
  printf(VALUE(bild, 3),iter_steps,einheit);
  printf("%s",VALUE(bild, 4));
 
 
for(i=1; i<=12; ++i)
if( i <= anz_niveaus )
  printf(VALUE(bild,i+4),phi[i],c[i],g[i],e[i],pe[i],i01[i],p1[i],i02[i],p2[i]);
else
  printf(VALUE(bild,i+4),phi[i],c[i]);
 
 
for(i=13; i<=14; ++i)
if( i <= anz_niveaus )
  printf(VALUE(bild,i+4),g[i],e[i],pe[i],i01[i],p1[i],i02[i],p2[i]);
else
  printf("%s",VALUE(bild,i+4) );
 
 
for(i=15; i<=15; ++i)
if( i <= anz_niveaus )
  printf(VALUE(bild,i+4),einheit,g[i],e[i],pe[i],i01[i],p1[i],i02[i],p2[i]);
else
  printf(VALUE(bild,i+4),einheit );
 
 
for(i=16; i<=16; ++i)
if( i <= anz_niveaus )
  printf(VALUE(bild,i+4),        g[i],e[i],pe[i],i01[i],p1[i],i02[i],p2[i]);
else
  printf("%s",VALUE(bild,i+4) );
 
for(i=17; i<=17; ++i)
if( i <= anz_niveaus )
  printf(VALUE(bild,i+4),x,      g[i],e[i],pe[i],i01[i],p1[i],i02[i],p2[i]);
else
  printf(VALUE(bild,i+4),x );
 
  printf(VALUE(bild,22),W );
  printf(VALUE(bild,23),(IS_SUSA(iteration)? '*':' '),
                        (IS_SUSB(iteration)? '*':' '),
                        (IS_SUSC(iteration)? '*':' '),
                        (IS_SUSP(iteration)? '*':' '),
                        (IS_EIGENWERT(    iteration)? '*':' '),
                        (IS_INTENSITAET(  iteration)? '*':' '),
                        (IS_MATRIXELEMENT(iteration)? '*':' ') );
  printf(VALUE(bild,24),chi2,(IS_MAGA(iteration)? '*':' '),
                             (IS_MAGB(iteration)? '*':' '),
                             (IS_MAGC(iteration)? '*':' '),
                             (IS_MAGP(iteration)? '*':' '),
                             (IS_POSFIT(iteration)? '*':' ')  );
 
  printf("%s",VALUE(bild,25) );
 
 
 
  free_vr(v);
  free_vr(xW);
  free_vr(rphi);
  free_mx(aJtb2);
  return 0;
}
/*------------------------------------------------------------------------------
                               if_stabil()
------------------------------------------------------------------------------*/
void if_stabil(amoeba,ewproblem,iteration,w,iter_steps,chi2)
    MINIMUM   *amoeba;
    EWPROBLEM *ewproblem;
    ITERATION *iteration;
    VEKTOR    *w;
    INT       iter_steps;
    DOUBLE    chi2;
{
  UNUSED_PARAMETER(iter_steps); 
 
  DOUBLE  macheps;
  INT     i_v,i_v0,flag;
  VEKTOR  *fix;
  NEBENBEDINGUNG *neben;
  VEKTOR  *rphi,*v_rphi(),*v,*v0;
  INT  ionennr;
  STEVENS *stevens;
  CHAR *ionname;
 
  macheps       = ewproblem->eps_machine;
  v0            = V0( iteration );
  neben         = NEBENBEDINGUNG(iteration);
  stevens       = PKQ(           iteration);
  ionname       = IONNAME( iteration );
 
  ionennr       = isimplementiert( ionname );
  fix           = FIX(    neben     );  /* z.B. (0,1,0,1,0,1,1,1,0) */
 
 
  v   = vr_alloc( VRDIM(v0) );
  for( i_v0=1,i_v=1; i_v0<=VRDIM(v0); ++i_v0 )
       if( RV(fix,i_v0) == 1.0 )   RV(v,i_v0)=RV(v0,i_v0);
       else                      { RV(v,i_v0)=RV(w,i_v); ++i_v;}
  rphi  =  v_rphi(v,stevens,ionennr,macheps);
 
 
 
   flag                       = 0;
   if( is_equal( chi2       - amoeba->chi2,0.0,macheps)) ++flag;
       amoeba->chi2 = chi2;
 
   if( is_equal( RV(rphi,1) - amoeba->r,0.0,macheps)) ++flag;
       amoeba->r = RV(rphi,1);
 
   if( is_equal( RV(rphi,2) - amoeba->phi1,0.0,macheps)) ++flag;
       amoeba->phi1 = RV(rphi,2);
 
   if( is_equal( RV(rphi,3) - amoeba->phi2,0.0,macheps)) ++flag;
       amoeba->phi2 = RV(rphi,3);
 
   if( is_equal( RV(rphi,4) - amoeba->phi3,0.0,macheps)) ++flag;
       amoeba->phi3 = RV(rphi,4);
 
   if( is_equal( RV(rphi,5) - amoeba->phi4,0.0,macheps)) ++flag;
       amoeba->phi4 = RV(rphi,5);
 
   if( is_equal( RV(rphi,6) - amoeba->phi5,0.0,macheps)) ++flag;
       amoeba->phi5 = RV(rphi,6);
 
   if( is_equal( RV(rphi,7) - amoeba->phi6,0.0,macheps)) ++flag;
       amoeba->phi6 = RV(rphi,7);
 
   if( is_equal( RV(rphi,8) - amoeba->phi7,0.0,macheps)) ++flag;
       amoeba->phi7 = RV(rphi,8);
 
   if( is_equal( RV(rphi,9) - amoeba->phi8,0.0,macheps)) ++flag;
       amoeba->phi8 = RV(rphi,9);
 
   if( flag == 10){
      TEXT2(amoeba)                = STABLE;
      amoeba -> anz_wiederholung  += 1;
   }
   else{
      amoeba -> anz_wiederholung  = 0;
      TEXT2(amoeba)               = UNSTABLE;
      }
 
 
 
  free_vr(v);
  free_vr(rphi);
 
}
/*------------------------------------------------------------------------------
ENDEMODUL    O R T H O   C
------------------------------------------------------------------------------*/
