/*-----------------------------------------------------------------------------
 
                         H A U P T P R O G R A M M
 
                                 CFIELD C
 
-------------------------------------------------------------------------------
Literatur :
-----------
 
Welche der V  fuer die verschiedenen Symmetrien unabhaengig sind , wurde
            kq
der Tabelle 5 ( Seite 484 ) aus  "Handbook on Physics and Chemistry of
                                  RARE EARTHS ",Vol.5,North-Holland
                                  Publishing Company ,Amsterdam 1982
entnommen.
 
 
 
Die Konstanten der Elektronenmasse ,Elektronenladung,Lichtgeschwindigkeit
sowie das Plankschewirkungsquantum wurden dem HELLWEGE(1988):
"Einfuehrung in die Festkoerperphysik"  entnommen .
 
 
 Includedateien:  types     c
                  define_J  c
 
 
 
 Module        :  cfield    c
                  diahermx  c
                  eingabe   c
                  intensit  c
                  komplex   c
                  matrix    c
                  stevens   c
                  minima    c
                  theta     c
                  ortho     c
                  mb11a     c
                  va05a     c
                  spline    c
 
 
 
 
-----------------------------------------------------------------------------*/
 
 
/*----------------------------------------------------------------------------
Includedateien holen
-----------------------------------------------------------------------------*/
#include <stdlib.h> /* atoi()   definieren          */
#include <stdio.h>  /* printf() definieren          */
#include <ctype.h>  /* isupper(),isdigit() defin.   */
#include <math.h>   /* sqrt,acos,etc definieren     */
#define pi (4.0*atan(1.0))
#include "types.c"  /* selbstdefinierte Typen holen */
 
 
/*----------------------------------------------------------------------------
Extern definierte Funktionen
-----------------------------------------------------------------------------*/
extern EWPROBLEM *diagonalisiere(EWPROBLEM *ewproblem,MATRIX *matrix,INT overwrite,SETUP *setup);  /* definiert in DIAHERMX.C*/
extern INT       neben_create(CHAR modus,CHAR *name_chi2);     /* definiert in EINGABE.C */
 
extern STEVENS   *calc_Pkq(INT dimj);        /* definiert in STEVENS.C */
extern KOMPLEX   *Clm(INT l,INT m,DOUBLE theta,DOUBLE phi);             /* definiert in STEVENS.C */
extern MATRIX    *mx_alloc(INT anz_ze,INT anz_sp);        /* definiert in MATRIX.C  */
extern DOUBLE    pow_(DOUBLE x,INT n) ;             /* definiert in STEVENS.C */
extern SETUP     *cfield_setup();    /* definiert in DIAHERMX.C*/
extern INT info_ewproblem(EWPROBLEM * ewproblem);         /* definiert in DIAHERMX.C*/
extern INT info_Pkq(INT k,INT q,INT dimj,CHAR modus);               /* definiert in STEVENS.C */
extern INT info_Fkq(INT k,INT q,INT dimj);               /* definiert in STEVENS.C */
extern INT info_Gkq(INT k,INT q,INT dimj);               /* definiert in STEVENS.C */
extern INT info_STEVkq(INT k,INT q,INT dimj,CHAR modus);            /* definiert in STEVENS.C */
extern INT info_tensor_Clm(INT l,INT m,DOUBLE theta,DOUBLE phi);        /* definiert in STEVENS.C */
extern INT fkq_tabelle(INT anz_ionen);            /* definiert in STEVENS.C */
extern INT gkq_tabelle(INT anz_ionen);            /* definiert in STEVENS.C */
extern INT info_thetakq();           /* definiert in THETA.C   */
extern INT info_rn();                /* definiert in THETA.C   */
extern INT atoi(const char *);                   /* definiert in <stdlib.h>*/
extern INT output(SETUP *setup,EWPROBLEM *ewproblem,KRISTALLFELD *kristallfeld,CHAR  modus);                 /* definiert in INTENSITY.C*/
extern INT write_title(FILE *fp);            /* definiert in DIAHERMX.C*/
 
extern INT create_Akq(INT einheitnr_in,INT einheitnr_out,CHAR *ion,INT symmetrienr,CHAR modus,DOUBLE temp);             /* definiert in EINGABE.C*/
extern INT create_Bkq(INT einheitnr_in,INT einheitnr_out,CHAR *ion,INT symmetrienr,CHAR modus,DOUBLE temp);             /* definiert in EINGABE.C*/
extern INT create_Dkq(INT einheitnr_in,INT einheitnr_out,CHAR *ion,INT symmetrienr,CHAR modus,DOUBLE temp);             /* definiert in EINGABE.C*/
extern INT create_Lkq(INT einheitnr_in,INT einheitnr_out,CHAR *ion,INT symmetrienr,CHAR modus,DOUBLE temp);             /* definiert in EINGABE.C*/
extern INT create_Vkq(INT einheitnr_in,INT einheitnr_out,CHAR *ion,INT symmetrienr,CHAR modus,DOUBLE temp);             /* definiert in EINGABE.C*/
extern INT create_Wkq(INT einheitnr_in,INT einheitnr_out,CHAR *ion,INT symmetrienr,CHAR modus,DOUBLE temp);             /* definiert in EINGABE.C*/
extern INT create_xW(INT einheitnr_in,INT einheitnr_out,CHAR *ion,INT symmetrienr,CHAR modus,DOUBLE temp);              /* definiert in EINGABE.C*/
extern INT create_nn(CHAR   *name,CHAR modus,INT nn,CHAR   *ion,DOUBLE temp);              /* definiert in EINGABE.C*/
 
extern ITERATION *read_Akq(CHAR *name, INT vsymmetrienr_vor);        /* definiert in EINGABE.C*/
extern ITERATION *read_Bkq(CHAR *name, INT vsymmetrienr_vor);        /* definiert in EINGABE.C*/
extern ITERATION *read_Bkqnew(CHAR *ion);        /* definiert in EINGABE.C*/
extern ITERATION *read_Dkq(CHAR *name, INT vsymmetrienr_vor);        /* definiert in EINGABE.C*/
extern ITERATION *read_Lkq(CHAR *name, INT vsymmetrienr_vor);        /* definiert in EINGABE.C*/
extern ITERATION *read_Vkq(CHAR *name, INT vsymmetrienr_vor);        /* definiert in EINGABE.C*/
extern ITERATION *read_Wkq(CHAR *name, INT vsymmetrienr_vor);        /* definiert in EINGABE.C*/
extern ITERATION *read_xW(CHAR *name, INT vsymmetrienr_vor);         /* definiert in EINGABE.C*/
extern UMGEBUNG  *read_nn(CHAR *name);         /* definiert in EINGABE.C*/

extern FILE *fopen_errchk(const char * filename,const char * mode);         /* definiert in EINGABE.C*/
  
extern DOUBLE alpha_J[]; /* definiert in theta.c */
extern DOUBLE beta_J[];  /* definiert in theta.c */
extern DOUBLE gamma_J[]; /* definiert in theta.c */

extern IONEN IONENIMP[]; /* definiert in theta.c*/ 

INT isimplementiert(CHAR*);
void read_error(INT,FILE*,CHAR*);
INT coor_error(FILE*,CHAR*,INT,INT);
int strlen_own(CHAR*);
INT equal(CHAR*,CHAR*);
/*----------------------------------------------------------------------------
                                Programmname
-----------------------------------------------------------------------------*/
 
#define PROGRAMMNAME "Cfield"
 
/*----------------------------------------------------------------------------
Globale Arrays und deren defines
-----------------------------------------------------------------------------*/
_CALFUN CF;  /*globale struktur fuer routine CALFUN, */
             /* welche von der fortranroutine VA05A   */
             /* aufgerufen wird                       */
 

/*define ANZ_IONEN              (  sizeof(IONENIMP)/sizeof(IONEN)  )*/
#define ANZ_IONEN              30
 
SYMMETRIE SYMLISTE[]={
/* 0 */   {  0 , "triklin   " } ,
/* 1 */   {  1 , "monoklin  " } ,
/* 2 */   {  2 , "rhombisch " } ,
/* 3 */   {  3 , "tetragonal" } ,
/* 4 */   {  4 , "tetragonal" } ,
/* 5 */   {  5 , "trigonal  " } ,
/* 6 */   {  6 , "trigonal  " } ,
/* 7 */   {  7 , "hexagonal " } ,
/* 8 */   {  8 , "kubisch   " } ,
};
#define ANZ_SYMNR   ( sizeof(SYMLISTE) / sizeof(SYMMETRIE) )
 
INFO INFOLISTE[]={
/* 0 */   {"c[lm] l m theta phi"   },
/* 1 */   {"ep[silonkq]"           },
/* 2 */   {"ew[=eigenwertproblem] A[kq]|B[kq]|D[kq]|L[kq]|V[kq]|W[kq]|x[W]"},
/* 3 */   {"ew[=eigenwertproblem] s[i]ngleion] filename.type symmetrienr"},
/* 4 */   {"F[kq(J)]    k q J"      } ,
/* 5 */   {"F[kq-Tabelle]    "      } ,
/* 6 */   {"G[kq(J)]    k q J"      } ,
/* 7 */   {"G[kq-Tabelle]    "      } ,
/* 8 */   {"h[amilton]"             } ,
/* 9 */   {"k[onstanten]"           } ,
/*10 */   {"m[agnetfeld] J Bx By Bz"} ,
/*11 */   {"o[megakq]"              } ,
/*12 */   {"P[kq(J)]    G[GT]|n[orm] k q J"},
/*13 */   {"S[TEVkq(J)] G[GT]|n[orm] k q J"},
/*14 */   {"r[n]"                    },
/*15 */   {"sy[mmetrienummern]"     } ,
/*16 */   {"t[hetakq]"             },
/*17 */   {"V[ersionsnummer des Programmes]"},
/*18 */   {"v[lm] eV|cm-1 filename.filetype symmetrienummer"},
};
#define ANZ_INFOS   ( sizeof(INFOLISTE) / sizeof(INFO) )
 
 
BEFEHL BEFEHLLISTE[]={
/* 0 */{"Information (help)................................"  , "-i[nfo   ]"},
/* 1 */{"Create an input file.............................."  , "-c[reate ]"},
/* 2 */{"Read an input file................................"  , "-r[ead   ]"},
/* 3 */{"Calculate susceptability (in small applied fields)"  , "-s[uszept]"},
/* 4 */{"Calculate RE (Rare Earth) moment.................."  , "-m[oment ]"},
/* 5 */{"Polycrystalline average of the RE-Moments (cubic)"  , "-k[poly  ]"},
/* 6 */{"Rhombic crystal field fit (Symmetry number = 2)..."  , "-o[rtho  ]"},
/* 7 */{"Determine Chi-squared for crystal field fit......."  , "-b[eside ]"},
};
#define ANZ_BEFEHLE ( sizeof(BEFEHLLISTE) / sizeof(BEFEHL) )
 
 
EINHEIT EINHEITIMP[]={
/*  0 */   { "Kelvin"  , 'K' , 0.0 , 0.0 , 0.0},
/*  1 */   { "Kelvin"  , 'k' , 0.0 , 0.0 , 0.0},
/*  2 */   { "meV   "  , 'M' , 0.0 , 0.0 , 0.0},
/*  3 */   { "meV   "  , 'm' , 0.0 , 0.0 , 0.0},
/*  4 */   { "eV    "  , 'E' , 0.0 , 0.0 , 0.0},
/*  5 */   { "eV    "  , 'e' , 0.0 , 0.0 , 0.0},
/*  6 */   { "cm-1  "  , 'C' , 0.0 , 0.0 , 0.0},
/*  7 */   { "cm-1  "  , 'c' , 0.0 , 0.0 , 0.0},
/*  8 */   { "Thz   "  , 'T' , 0.0 , 0.0 , 0.0},
/*  9 */   { "Thz   "  , 't' , 0.0 , 0.0 , 0.0},
};
#define ANZ_EINHEIT            (  sizeof(EINHEITIMP)/sizeof(EINHEIT)  )
 
void init_einheit()
{
    INT i;
 
    for(i=0 ; i<(int)ANZ_EINHEIT ; ++i )
          switch(i){
             case 0:
             case 1: EINHEITIMP[i].myB = myBKelvin;
                     EINHEITIMP[i].fke = FKelvinKelvin;
                     EINHEITIMP[i].fek = FKelvinKelvin;
                     break;
             case 2:
             case 3: EINHEITIMP[i].myB = myBmeV;
                     EINHEITIMP[i].fke = FKelvinmeV;
                     EINHEITIMP[i].fek = FmeVKelvin;
                     break;
             case 4:
             case 5: EINHEITIMP[i].myB = myBeV;
                     EINHEITIMP[i].fke = FKelvineV;
                     EINHEITIMP[i].fek = FeVKelvin;
                     break;
 
             case 6:
             case 7: EINHEITIMP[i].myB = myBcm;
                     EINHEITIMP[i].fke = FKelvincm;
                     EINHEITIMP[i].fek = FcmKelvin;
                     break;
 
             case 8:
             case 9: EINHEITIMP[i].myB = myBThz;
                     EINHEITIMP[i].fke = FKelvinThz;
                     EINHEITIMP[i].fek = FThzKelvin;
                     break;
 
        default : clearscreen;
                  printf("mistake in init_einheit in %s.C !\n\n",PROGRAMMNAME);
                  printf("Unit %d not implemented.\n",i);
                  exit(0);
          }
}
 
INFO INFO_EINHEIT[]={
/*  0 */     {"K[elvin]"},
/*  1 */     {"m[eV]"   },
/*  2 */     {"e[V]"    },
/*  3 */     {"c[m-1]"  },
/*  4 */     {"T[hz]"   },
};
#define ANZ_INFO_EINHEIT       (  sizeof(INFO_EINHEIT)/sizeof(INFO)  )
 
/*-------------------------------*/
 
 
INT is_einheit_imp(CHAR c ) /* Bestimme Einheit anhand des Erkennugszeichens c */
{
   INT i;
 
   for(i=0; i<(int)ANZ_EINHEIT; ++i)
       if( c==EINHEITIMP[i].c )
             return( i );
 
   return( NICHTIMP ); /* Einheit nicht implementiert */
}
/*----------------------------------------------------------------------------
                                 Faktoren epsilon     definieren
                                                 kq
 
-----------------------------------------------------------------------------*/
DOUBLE fak(INT n)   /*   n!   */
{
   if( n==0 )
      return( 1.0 );
   return( (DOUBLE)n*fak(n-1) );
}
 
 
DOUBLE powi(DOUBLE x,INT n)  /* x**n */
{
   if( n==0 ) return( 1.0 );
   return( x*powi(x,n-1) );
}
 
#define gn0n(n) ( fak(2*(n))/fak(n) )
#define gn1n(n) ( gn0n(n+1) )
#define gn2n(n) ( gn0n(n+2)/(2.0*(n)+3) )
#define gn3n(n) ( gn0n(n+3)/(2.0*(n)+5) )
#define gn4n(n) ( gn0n(n+4)/(2.0*(n)+5)/(2.0*(n)+7) )
#define gn5n(n) ( gn0n(n+5)/(2.0*(n)+7)/(2.0*(n)+9) )
#define gn6n(n) ( gn0n(n+6)/(2.0*(n)+7)/(2.0*(n)+9)/(2.0*(n)+11)*3.0 )
 
#define dummy1(k,q) (  powi(-1.0,(q+ABS(q))/2)/powi(2.0,k)    )
#define dummy2(k,q) (1.0/sqrt((DOUBLE)(fak(k+ABS(q))*fak(k-ABS(q)))))
#define dummy(k,q)  dummy1(k,q)*dummy2(k,q)
 
DOUBLE  epn0n(INT n) {return( dummy(ABS(n)  ,n)*gn0n(ABS(n)) );}
DOUBLE  epn1n(INT n) {return( dummy(ABS(n)+1,n)*gn1n(ABS(n)) );}
DOUBLE  epn2n(INT n) {return( dummy(ABS(n)+2,n)*gn2n(ABS(n)) );}
DOUBLE  epn3n(INT n) {return( dummy(ABS(n)+3,n)*gn3n(ABS(n)) );}
DOUBLE  epn4n(INT n) {return( dummy(ABS(n)+4,n)*gn4n(ABS(n)) );}
DOUBLE  epn5n(INT n) {return( dummy(ABS(n)+5,n)*gn5n(ABS(n)) );}
DOUBLE  epn6n(INT n) {return( dummy(ABS(n)+6,n)*gn6n(ABS(n)) );}
 
 
/*------------------------*/
 
DOUBLE omegan0n(INT n)
{
   UNUSED_PARAMETER(n);
   return( 1.0 );
}
 
/*------------------------*/
DOUBLE omegan1n(INT n)
{
   UNUSED_PARAMETER(n);
   return( 1.0 );
}
 
/*------------------------*/
 
DOUBLE omegan2n(INT n)
{
   UNUSED_PARAMETER(n);
   return( 1.0 );
}
 
/*------------------------*/
 
DOUBLE omegan3n(INT n)
{
  INT k=0;
  INT z;
 
   if(n==0) return( 1.0 );
 
   do{
   z = 2 + 3*k;
   if( z==n )  return( 3.0 );
   if( z>n  )  return( 1.0 );
   k++;
   }while(1);
}
 
/*------------------------*/
 
DOUBLE omegan4n(INT n)
{
  INT k=0;
  INT z1;
  INT z2;
 
   if(n==0) return( 1.0 );
 
   do{
   z1 = 1 + 3*k;
   z2 = 2 + 3*k;
   if( z1==n || z2==n )  return( 3.0 );
   if( z1>n  || z2>n  )  return( 1.0 );
   k++;
   }while(1);
}
 
/*------------------------*/
 
DOUBLE omegan5n(INT n)
{
  INT k=0;
  INT z1;
  INT z2;
  INT z3;
  DOUBLE e1/*,e2*/;
 
   if(n==0) return( 1.0 );
 
   e1 = 1.0;
   do{
   z1 = 1 + 3*k;
   if( z1==n ) e1=3.0;
   k++;
   }while( z1<n );
 
   k=0;
/* e2 = 1.0; */
   do{
   z2 = 3 + 5*k;
   z3 = 4 + 5*k;
   if( z2==n || z3==n ) return( 5.0*e1 );
   if( z2> n || z3> n ) return( 1.0*e1 );
   k++;
   }while( 1 );
 
}
 
/*------------------------*/
 
DOUBLE omegan6n(INT n)
{
  INT k=0;
  INT z1;
  INT z2;
  INT z3;
 
   if(n==0) return( 1.0 );
 
   do{
   z1 = 2 + 5*k;
   z2 = 3 + 5*k;
   z3 = 4 + 5*k;
   if( z1==n || z2==n || z3==n )  return( 5.0 );
   if( z1> n || z2> n || z3> n )  return( 1.0 );
   k++;
   }while(1);
}
 
/*------------------------*/
 
/*------------------------------------------------------------------------------
 function for extracting some important parameters of ions - for INTERNAL 
 class ionpars
------------------------------------------------------------------------------*/
void getpar(char * iontype, int * dj, double * alpha, double * beta, double * gamma, double * lande,
       double * rh2, double * rh4,double * rh6, int * nof_electrons )
{
    KRISTALLFELD *kristallfeld,*init_iterationnew(CHAR*);
/*  SETUP        *setup; */
    ITERATION    *iteration;
    init_einheit();
/*  setup =  cfield_setup();*/ /* setup's aus SETUP-file holen */
    kristallfeld=init_iterationnew(iontype);/*sets up HMAG */
    iteration= ITERATION (kristallfeld);

    (*dj)=DIMJ(iteration);
    (*alpha)=alpha_J[E4f(IONENNUMMER(iteration))];
    (*beta)=beta_J[E4f(IONENNUMMER(iteration))];
    (*gamma)=gamma_J[E4f(IONENNUMMER(iteration))];
    (*lande)=GJ(iteration);
    (*rh2)=r2(IONENNUMMER(iteration));
    (*rh4)=r4(IONENNUMMER(iteration));
    (*rh6)=r6(IONENNUMMER(iteration));
    (*nof_electrons)= IONENIMP[ IONENNUMMER(iteration) ].elektronen_in_vier_f;
   
}

/*------------------------------------------------------------------------------
 mcphase function for calculating the crystal field matrix  - full diagonalization - for INTERNAL
 module cfield
------------------------------------------------------------------------------*/
/*cfield_mcphasnew is designed for operation without reading a Bkq.parameter file*/
void cfield_mcphasnew(char * iontype, double * Jxr,double * Jxi,  double * Jyr, double * Jyi, double * Jzr, double * Jzi,
                              double * mo22sr, double * mo22si,
                              double * mo21sr, double * mo21si,
                              double * mo20cr, double * mo20ci,
                              double * mo21cr, double * mo21ci,
                              double * mo22cr, double * mo22ci,

                              double * mo33sr, double * mo33si,
                              double * mo32sr, double * mo32si,
                              double * mo31sr, double * mo31si,
                              double * mo30cr, double * mo30ci,
                              double * mo31cr, double * mo31ci,
                              double * mo32cr, double * mo32ci,
                              double * mo33cr, double * mo33ci,

                              double * mo44sr, double * mo44si,
                              double * mo43sr, double * mo43si,
                              double * mo42sr, double * mo42si,
                              double * mo41sr, double * mo41si,
                              double * mo40cr, double * mo40ci,
                              double * mo41cr, double * mo41ci,
                              double * mo42cr, double * mo42ci,
                              double * mo43cr, double * mo43ci,
                              double * mo44cr, double * mo44ci,

                              double * mo55sr, double * mo55si,
                              double * mo54sr, double * mo54si,
                              double * mo53sr, double * mo53si,
                              double * mo52sr, double * mo52si,
                              double * mo51sr, double * mo51si,
                              double * mo50cr, double * mo50ci,
                              double * mo51cr, double * mo51ci,
                              double * mo52cr, double * mo52ci,
                              double * mo53cr, double * mo53ci,
                              double * mo54cr, double * mo54ci,
                              double * mo55cr, double * mo55ci,

                              double * mo66sr, double * mo66si,
                              double * mo65sr, double * mo65si,
                              double * mo64sr, double * mo64si,
                              double * mo63sr, double * mo63si,
                              double * mo62sr, double * mo62si,
                              double * mo61sr, double * mo61si,
                              double * mo60cr, double * mo60ci,
                              double * mo61cr, double * mo61ci,
                              double * mo62cr, double * mo62ci,
                              double * mo63cr, double * mo63ci,
                              double * mo64cr, double * mo64ci,
                              double * mo65cr, double * mo65ci,
                              double * mo66cr, double * mo66ci,

                              double * modxcr, double * modxci,
                              double * modycr, double * modyci,
                              double * modzcr, double * modzci,

                              double * modxdxcr, double * modxdxci,
                              double * modydycr, double * modydyci,
                              double * modzdzcr, double * modzdzci,
                              int * dj, 
                              double * alpha, double * beta, double * gamma, 
                              double * lande,
                              double * rh2,double * rh4,double * rh6, int * nof_electrons)
{
    KRISTALLFELD *kristallfeld,*init_iterationnew(CHAR*);
    EWPROBLEM /* *ewproblem,*/ *setuphcf( EWPROBLEM    *,INT, KRISTALLFELD *,CHAR );
    /*SETUP        *setup;*/
    INT          m,n;
    INT          i,j;
    DOUBLE       Bx,By,Bz,myB;
    ITERATION    *iteration;
    MATRIX       *calc_Bmag(INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

    init_einheit();


   /* setup = cfield_setup(); *//* setup's aus SETUP-file holen */
  
/* determine crystal field stevens operator matrices */

    kristallfeld=init_iterationnew(iontype);/*sets up HMAG */
							    
    iteration= ITERATION (kristallfeld);
    myB              = myBmeV;    

    (*dj)=DIMJ(iteration);
    (*alpha)=alpha_J[E4f(IONENNUMMER(iteration))];
    (*beta)=beta_J[E4f(IONENNUMMER(iteration))];
    (*gamma)=gamma_J[E4f(IONENNUMMER(iteration))];
    (*lande)=GJ(iteration);
    (*rh2)=r2(IONENNUMMER(iteration));
    (*rh4)=r4(IONENNUMMER(iteration));
    (*rh6)=r6(IONENNUMMER(iteration));
    (*nof_electrons)= IONENIMP[ IONENNUMMER(iteration) ].elektronen_in_vier_f;

    Bx=-1.0/GJ(iteration)/myB;By=0.0;Bz=0.0;
      HMAG(iteration)=calc_Bmag(DIMJ(iteration),GJ(iteration),myB,Bx,By,Bz);
    for (i=1;i<=DIMJ(iteration);++i)
     {for (j=1;j<=DIMJ(iteration);++j)
     {Jxr[30*(i-1)+j-1]=R(HMAG(iteration),i,j);Jxi[30*(i-1)+j-1]=I(HMAG(iteration),i,j);}
     }
     
     
    Bx=0.0;By=-1.0/GJ(iteration)/myB;Bz=0.0;
    HMAG(iteration)=calc_Bmag( DIMJ(iteration),GJ(iteration),myB,Bx,By,Bz);
    for (i=1;i<=DIMJ(iteration);++i)
     {for (j=1;j<=DIMJ(iteration);++j)
     {Jyr[30*(i-1)+j-1]=R(HMAG(iteration),i,j);Jyi[30*(i-1)+j-1]=I(HMAG(iteration),i,j);}
     }

    Bx=0.0;By=0.0; Bz=-1.0/GJ(iteration)/myB;
    HMAG(iteration)=calc_Bmag( DIMJ(iteration),GJ(iteration),myB,Bx,By,Bz);
    for (i=1;i<=DIMJ(iteration);++i)
     {for (j=1;j<=DIMJ(iteration);++j)
     {Jzr[30*(i-1)+j-1]=R(HMAG(iteration),i,j);Jzi[30*(i-1)+j-1]=I(HMAG(iteration),i,j);}
     }

    Bx=0.0;By=0.0; Bz=0.0;
    HMAG(iteration)=calc_Bmag( DIMJ(iteration),GJ(iteration),myB,Bx,By,Bz);
/*  ewproblem=*/setuphcf((EWPROBLEM*)0,NEIN,kristallfeld,BKQ);

    /* h += singleion anisotropy */
    for( n= DIMJ(iteration); n>=1 ; --n) for( m=DIMJ(iteration) ; m>=1 ; --m) {
      modxcr[30*(n-1)+m-1]=0.; modxci[30*(n-1)+m-1]=0.;
      modycr[30*(n-1)+m-1]=0.; modyci[30*(n-1)+m-1]=0.;
      modzcr[30*(n-1)+m-1]=0.; modzci[30*(n-1)+m-1]=0.;
      for( i=DIMJ(iteration) ; i>=1 ; --i){
              modxcr[30*(n-1)+m-1] += Jxr[30*(n-1)+i-1]*Jxr[30*(i-1)+m-1];
              modycr[30*(n-1)+m-1] -= Jyi[30*(n-1)+i-1]*Jyi[30*(i-1)+m-1];
              modzcr[30*(n-1)+m-1] += Jzr[30*(n-1)+i-1]*Jzr[30*(i-1)+m-1]; }
      }
    for( n= DIMJ(iteration); n>=1 ; --n) for( m=DIMJ(iteration) ; m>=1 ; --m) {
      modxdxcr[30*(n-1)+m-1]=0.; modxdxci[30*(n-1)+m-1]=0.;
      modydycr[30*(n-1)+m-1]=0.; modydyci[30*(n-1)+m-1]=0.;
      modzdzcr[30*(n-1)+m-1]=0.; modzdzci[30*(n-1)+m-1]=0.;
      for( i=DIMJ(iteration) ; i>=1 ; --i){
              modxdxcr[30*(n-1)+m-1] += modxcr[30*(n-1)+i-1]*modxcr[30*(i-1)+m-1];
              modydycr[30*(n-1)+m-1] += modycr[30*(n-1)+i-1]*modycr[30*(i-1)+m-1];
              modzdzcr[30*(n-1)+m-1] += modzcr[30*(n-1)+i-1]*modzcr[30*(i-1)+m-1]; }
      }

/*sets up hamiltonian using olms and solves ev problem (necessary step although CF is zero in this routine)*/

/* here we use the stevens operators for the quadrupolar hamiltonian  (NOT the Racah) */
/* ===============================================================
|                                                             |
|Kristallfeldparameter  B    in  meV     und  reell           |
|                        kq                                   |
|-----------------------------                                |
|        ---                 |                                |
| H   =  >     B   STEV (J)  |  mit  STEV (J)  hermitesch     |
|  KF    ---    kq     kq    |           kq                   |
|        k>=0                |                                |
|        q>=0                |                                |
|-----------------------------                                |
|                                                             |
| und    q=0 :  STEV (J) := P  (J)                            |
|                   k0       k0                               |
|                                                             |
|                                       +                     |
|        q>0 :  STEV (J)c := ( P  (J) + P  (J) )/2/omega      |
|                   kq         kq       kq             kq     |
|                                                             |
|                                       +                     |
|        q<0 :  STEV (J)s := ( P  (J) - P  (J) )/2/i/omega    |
|                   kq         kq       kq               kq   |
|                                                             |
*/

    for( n=DIMJ(iteration) ; n>=1 ; --n )
         for( m=DIMJ(iteration) ; m>=1 ; --m ){


             mo22sr[30*(n-1)+m-1]= ( I(P2P2(PKQ(iteration)),n,m )-I(P2M2(PKQ(iteration)),n,m))/2/omegan0n(2);
	     mo22si[30*(n-1)+m-1]=-( R(P2P2(PKQ(iteration)),n,m )-R(P2M2(PKQ(iteration)),n,m))/2/omegan0n(2);
             mo21sr[30*(n-1)+m-1]= ( I(P2P1(PKQ(iteration)),n,m )-I(P2M1(PKQ(iteration)),n,m))/2/omegan1n(1);
	     mo21si[30*(n-1)+m-1]=-( R(P2P1(PKQ(iteration)),n,m )-R(P2M1(PKQ(iteration)),n,m))/2/omegan1n(1);
	     mo20cr[30*(n-1)+m-1]=  R(P2P0(PKQ(iteration)),n,m );
	     mo20ci[30*(n-1)+m-1]=0;
             mo21cr[30*(n-1)+m-1]=( R(P2P1(PKQ(iteration)),n,m )+R(P2M1(PKQ(iteration)),n,m))/2/omegan1n(1);
	     mo21ci[30*(n-1)+m-1]=( I(P2P1(PKQ(iteration)),n,m )+I(P2M1(PKQ(iteration)),n,m))/2/omegan1n(1);
             mo22cr[30*(n-1)+m-1]=( R(P2P2(PKQ(iteration)),n,m )+R(P2M2(PKQ(iteration)),n,m))/2/omegan0n(2);
	     mo22ci[30*(n-1)+m-1]=( I(P2P2(PKQ(iteration)),n,m )+I(P2M2(PKQ(iteration)),n,m))/2/omegan0n(2);


             mo33sr[30*(n-1)+m-1]= ( I(P3P3(PKQ(iteration)),n,m )-I(P3M3(PKQ(iteration)),n,m))/2/omegan0n(3);
	     mo33si[30*(n-1)+m-1]=-( R(P3P3(PKQ(iteration)),n,m )-R(P3M3(PKQ(iteration)),n,m))/2/omegan0n(3);
             mo32sr[30*(n-1)+m-1]= ( I(P3P2(PKQ(iteration)),n,m )-I(P3M2(PKQ(iteration)),n,m))/2/omegan1n(2);
	     mo32si[30*(n-1)+m-1]=-( R(P3P2(PKQ(iteration)),n,m )-R(P3M2(PKQ(iteration)),n,m))/2/omegan1n(2);
             mo31sr[30*(n-1)+m-1]= ( I(P3P1(PKQ(iteration)),n,m )-I(P3M1(PKQ(iteration)),n,m))/2/omegan2n(1);
	     mo31si[30*(n-1)+m-1]=-( R(P3P1(PKQ(iteration)),n,m )-R(P3M1(PKQ(iteration)),n,m))/2/omegan2n(1);
	     mo30cr[30*(n-1)+m-1]=  R(P3P0(PKQ(iteration)),n,m );
	     mo30ci[30*(n-1)+m-1]=0;
             mo31cr[30*(n-1)+m-1]=( R(P3P1(PKQ(iteration)),n,m )+R(P3M1(PKQ(iteration)),n,m))/2/omegan2n(1);
	     mo31ci[30*(n-1)+m-1]=( I(P3P1(PKQ(iteration)),n,m )+I(P3M1(PKQ(iteration)),n,m))/2/omegan2n(1);
             mo32cr[30*(n-1)+m-1]=( R(P3P2(PKQ(iteration)),n,m )+R(P3M2(PKQ(iteration)),n,m))/2/omegan1n(2);
	     mo32ci[30*(n-1)+m-1]=( I(P3P2(PKQ(iteration)),n,m )+I(P3M2(PKQ(iteration)),n,m))/2/omegan1n(2);
             mo33cr[30*(n-1)+m-1]=( R(P3P3(PKQ(iteration)),n,m )+R(P3M3(PKQ(iteration)),n,m))/2/omegan0n(3);
	     mo33ci[30*(n-1)+m-1]=( I(P3P3(PKQ(iteration)),n,m )+I(P3M3(PKQ(iteration)),n,m))/2/omegan0n(3);

             mo44sr[30*(n-1)+m-1]= ( I(P4P4(PKQ(iteration)),n,m )-I(P4M4(PKQ(iteration)),n,m))/2/omegan0n(4);
	     mo44si[30*(n-1)+m-1]=-( R(P4P4(PKQ(iteration)),n,m )-R(P4M4(PKQ(iteration)),n,m))/2/omegan0n(4);
             mo43sr[30*(n-1)+m-1]= ( I(P4P3(PKQ(iteration)),n,m )-I(P4M3(PKQ(iteration)),n,m))/2/omegan1n(3);
	     mo43si[30*(n-1)+m-1]=-( R(P4P3(PKQ(iteration)),n,m )-R(P4M3(PKQ(iteration)),n,m))/2/omegan1n(3);
             mo42sr[30*(n-1)+m-1]= ( I(P4P2(PKQ(iteration)),n,m )-I(P4M2(PKQ(iteration)),n,m))/2/omegan2n(2);
	     mo42si[30*(n-1)+m-1]=-( R(P4P2(PKQ(iteration)),n,m )-R(P4M2(PKQ(iteration)),n,m))/2/omegan2n(2);
             mo41sr[30*(n-1)+m-1]= ( I(P4P1(PKQ(iteration)),n,m )-I(P4M1(PKQ(iteration)),n,m))/2/omegan3n(1);
	     mo41si[30*(n-1)+m-1]=-( R(P4P1(PKQ(iteration)),n,m )-R(P4M1(PKQ(iteration)),n,m))/2/omegan3n(1);
	     mo40cr[30*(n-1)+m-1]=  R(P4P0(PKQ(iteration)),n,m );
	     mo40ci[30*(n-1)+m-1]=0;
             mo41cr[30*(n-1)+m-1]=( R(P4P1(PKQ(iteration)),n,m )+R(P4M1(PKQ(iteration)),n,m))/2/omegan3n(1);
	     mo41ci[30*(n-1)+m-1]=( I(P4P1(PKQ(iteration)),n,m )+I(P4M1(PKQ(iteration)),n,m))/2/omegan3n(1);
             mo42cr[30*(n-1)+m-1]=( R(P4P2(PKQ(iteration)),n,m )+R(P4M2(PKQ(iteration)),n,m))/2/omegan2n(2);
	     mo42ci[30*(n-1)+m-1]=( I(P4P2(PKQ(iteration)),n,m )+I(P4M2(PKQ(iteration)),n,m))/2/omegan2n(2);
             mo43cr[30*(n-1)+m-1]=( R(P4P3(PKQ(iteration)),n,m )+R(P4M3(PKQ(iteration)),n,m))/2/omegan1n(3);
	     mo43ci[30*(n-1)+m-1]=( I(P4P3(PKQ(iteration)),n,m )+I(P4M3(PKQ(iteration)),n,m))/2/omegan1n(3);
             mo44cr[30*(n-1)+m-1]=( R(P4P4(PKQ(iteration)),n,m )+R(P4M4(PKQ(iteration)),n,m))/2/omegan0n(4);
	     mo44ci[30*(n-1)+m-1]=( I(P4P4(PKQ(iteration)),n,m )+I(P4M4(PKQ(iteration)),n,m))/2/omegan0n(4);

             mo55sr[30*(n-1)+m-1]= ( I(P5P5(PKQ(iteration)),n,m )-I(P5M5(PKQ(iteration)),n,m))/2/omegan0n(5);
	     mo55si[30*(n-1)+m-1]=-( R(P5P5(PKQ(iteration)),n,m )-R(P5M5(PKQ(iteration)),n,m))/2/omegan0n(5);
             mo54sr[30*(n-1)+m-1]= ( I(P5P4(PKQ(iteration)),n,m )-I(P5M4(PKQ(iteration)),n,m))/2/omegan1n(4);
	     mo54si[30*(n-1)+m-1]=-( R(P5P4(PKQ(iteration)),n,m )-R(P5M4(PKQ(iteration)),n,m))/2/omegan1n(4);
             mo53sr[30*(n-1)+m-1]= ( I(P5P3(PKQ(iteration)),n,m )-I(P5M3(PKQ(iteration)),n,m))/2/omegan2n(3);
	     mo53si[30*(n-1)+m-1]=-( R(P5P3(PKQ(iteration)),n,m )-R(P5M3(PKQ(iteration)),n,m))/2/omegan2n(3);
             mo52sr[30*(n-1)+m-1]= ( I(P5P2(PKQ(iteration)),n,m )-I(P5M2(PKQ(iteration)),n,m))/2/omegan3n(2);
	     mo52si[30*(n-1)+m-1]=-( R(P5P2(PKQ(iteration)),n,m )-R(P5M2(PKQ(iteration)),n,m))/2/omegan3n(2);
             mo51sr[30*(n-1)+m-1]= ( I(P5P1(PKQ(iteration)),n,m )-I(P5M1(PKQ(iteration)),n,m))/2/omegan4n(1);
	     mo51si[30*(n-1)+m-1]=-( R(P5P1(PKQ(iteration)),n,m )-R(P5M1(PKQ(iteration)),n,m))/2/omegan4n(1);
	     mo50cr[30*(n-1)+m-1]=  R(P5P0(PKQ(iteration)),n,m );
	     mo50ci[30*(n-1)+m-1]=0;
             mo51cr[30*(n-1)+m-1]=( R(P5P1(PKQ(iteration)),n,m )+R(P5M1(PKQ(iteration)),n,m))/2/omegan4n(1);
	     mo51ci[30*(n-1)+m-1]=( I(P5P1(PKQ(iteration)),n,m )+I(P5M1(PKQ(iteration)),n,m))/2/omegan4n(1);
             mo52cr[30*(n-1)+m-1]=( R(P5P2(PKQ(iteration)),n,m )+R(P5M2(PKQ(iteration)),n,m))/2/omegan3n(2);
	     mo52ci[30*(n-1)+m-1]=( I(P5P2(PKQ(iteration)),n,m )+I(P5M2(PKQ(iteration)),n,m))/2/omegan3n(2);
             mo53cr[30*(n-1)+m-1]=( R(P5P3(PKQ(iteration)),n,m )+R(P5M3(PKQ(iteration)),n,m))/2/omegan2n(3);
	     mo53ci[30*(n-1)+m-1]=( I(P5P3(PKQ(iteration)),n,m )+I(P5M3(PKQ(iteration)),n,m))/2/omegan2n(3);
             mo54cr[30*(n-1)+m-1]=( R(P5P4(PKQ(iteration)),n,m )+R(P5M4(PKQ(iteration)),n,m))/2/omegan1n(4);
	     mo54ci[30*(n-1)+m-1]=( I(P5P4(PKQ(iteration)),n,m )+I(P5M4(PKQ(iteration)),n,m))/2/omegan1n(4);
             mo55cr[30*(n-1)+m-1]=( R(P5P5(PKQ(iteration)),n,m )+R(P5M5(PKQ(iteration)),n,m))/2/omegan0n(5);
	     mo55ci[30*(n-1)+m-1]=( I(P5P5(PKQ(iteration)),n,m )+I(P5M5(PKQ(iteration)),n,m))/2/omegan0n(5);

             mo66sr[30*(n-1)+m-1]= ( I(P6P6(PKQ(iteration)),n,m )-I(P6M6(PKQ(iteration)),n,m))/2/omegan0n(6);
	     mo66si[30*(n-1)+m-1]=-( R(P6P6(PKQ(iteration)),n,m )-R(P6M6(PKQ(iteration)),n,m))/2/omegan0n(6);
             mo65sr[30*(n-1)+m-1]= ( I(P6P5(PKQ(iteration)),n,m )-I(P6M5(PKQ(iteration)),n,m))/2/omegan1n(5);
	     mo65si[30*(n-1)+m-1]=-( R(P6P5(PKQ(iteration)),n,m )-R(P6M5(PKQ(iteration)),n,m))/2/omegan1n(5);
             mo64sr[30*(n-1)+m-1]= ( I(P6P4(PKQ(iteration)),n,m )-I(P6M4(PKQ(iteration)),n,m))/2/omegan2n(4);
	     mo64si[30*(n-1)+m-1]=-( R(P6P4(PKQ(iteration)),n,m )-R(P6M4(PKQ(iteration)),n,m))/2/omegan2n(4);
             mo63sr[30*(n-1)+m-1]= ( I(P5P3(PKQ(iteration)),n,m )-I(P6M3(PKQ(iteration)),n,m))/2/omegan3n(3);
	     mo63si[30*(n-1)+m-1]=-( R(P6P3(PKQ(iteration)),n,m )-R(P6M3(PKQ(iteration)),n,m))/2/omegan3n(3);
             mo62sr[30*(n-1)+m-1]= ( I(P6P2(PKQ(iteration)),n,m )-I(P6M2(PKQ(iteration)),n,m))/2/omegan4n(2);
	     mo62si[30*(n-1)+m-1]=-( R(P6P2(PKQ(iteration)),n,m )-R(P6M2(PKQ(iteration)),n,m))/2/omegan4n(2);
             mo61sr[30*(n-1)+m-1]= ( I(P6P1(PKQ(iteration)),n,m )-I(P6M1(PKQ(iteration)),n,m))/2/omegan5n(1);
	     mo61si[30*(n-1)+m-1]=-( R(P6P1(PKQ(iteration)),n,m )-R(P6M1(PKQ(iteration)),n,m))/2/omegan5n(1);
	     mo60cr[30*(n-1)+m-1]=  R(P6P0(PKQ(iteration)),n,m );
	     mo60ci[30*(n-1)+m-1]=0;
             mo61cr[30*(n-1)+m-1]=( R(P6P1(PKQ(iteration)),n,m )+R(P6M1(PKQ(iteration)),n,m))/2/omegan5n(1);
	     mo61ci[30*(n-1)+m-1]=( I(P6P1(PKQ(iteration)),n,m )+I(P6M1(PKQ(iteration)),n,m))/2/omegan5n(1);
             mo62cr[30*(n-1)+m-1]=( R(P6P2(PKQ(iteration)),n,m )+R(P6M2(PKQ(iteration)),n,m))/2/omegan4n(2);
	     mo62ci[30*(n-1)+m-1]=( I(P6P2(PKQ(iteration)),n,m )+I(P6M2(PKQ(iteration)),n,m))/2/omegan4n(2);
             mo63cr[30*(n-1)+m-1]=( R(P6P3(PKQ(iteration)),n,m )+R(P6M3(PKQ(iteration)),n,m))/2/omegan3n(3);
	     mo63ci[30*(n-1)+m-1]=( I(P6P3(PKQ(iteration)),n,m )+I(P6M3(PKQ(iteration)),n,m))/2/omegan3n(3);
             mo64cr[30*(n-1)+m-1]=( R(P6P4(PKQ(iteration)),n,m )+R(P6M4(PKQ(iteration)),n,m))/2/omegan2n(4);
	     mo64ci[30*(n-1)+m-1]=( I(P6P4(PKQ(iteration)),n,m )+I(P6M4(PKQ(iteration)),n,m))/2/omegan2n(4);
             mo65cr[30*(n-1)+m-1]=( R(P6P5(PKQ(iteration)),n,m )+R(P6M5(PKQ(iteration)),n,m))/2/omegan1n(5);
	     mo65ci[30*(n-1)+m-1]=( I(P6P5(PKQ(iteration)),n,m )+I(P6M5(PKQ(iteration)),n,m))/2/omegan1n(5);
             mo66cr[30*(n-1)+m-1]=( R(P6P6(PKQ(iteration)),n,m )+R(P6M6(PKQ(iteration)),n,m))/2/omegan0n(6);
	     mo66ci[30*(n-1)+m-1]=( I(P6P6(PKQ(iteration)),n,m )+I(P6M6(PKQ(iteration)),n,m))/2/omegan0n(6);
       }
}


/*------------------------------------------------------------------------------

Bisherige Module  :
-------------------
                        main     c  : *

                        types    c  : Includefile
                        define_J c  : Includefile

                        matrix   c  : *
                        theta    c  : *
                        komplex  c  : *
                        diahermx c  : *
                        stevens  c  : *
			minima   c  : *
			spline   c  : *
			mb11a    c  : *
			ortho    c  : *
			va05a    c  : *
			intensit c  : *

 * : wird jeweils einzeln kompiliert und dann zusammengelinkt .
------------------------------------------------------------------------------*/
 
 
 
/*------------------------------------------------------------------------------
                                 isreell()
------------------------------------------------------------------------------*/
INT isreell(INT  symmetrienr ,CHAR * ion ) /* testet ob alle Kristallfeldparameter     */
                                 /* fuer das Ion ion bei der  Symmetrienummer*/
                                 /* symmetrienr i.a. reell sind              */
{
   INT ionennr,dimj,zwei_j;
 
   ionennr = isimplementiert( ion );
   dimj    = IONENIMP[ ionennr ].dimj;
   zwei_j  = dimj - 1;
 
 
   if( zwei_j >= 6 ){  switch( symmetrienr ){
                           case 8:
                           case 7:
                           case 6:
                           case 4:
                           case 2: return(JA);break;
 
                           case 5:
                           case 3:
                           case 1:
                           case 0: return(NEIN);break;
                       }
                    }
 
   if( zwei_j >= 4 ){  switch( symmetrienr ){
                           case 8:
                           case 7:
                           case 6:
                           case 5:
                           case 4:
                           case 3:
                           case 2: return(JA);break;
 
                           case 1:
                           case 0: return(NEIN);break;
                       }
                    }
 
   if( zwei_j >= 2 ){  switch( symmetrienr ){
                           case 8:
                           case 7:
                           case 6:
                           case 5:
                           case 4:
                           case 3:
                           case 2:
                           case 1: return(JA);break;
 
                           case 0: return(NEIN);break;
                       }
                    }
return(NEIN);
}
/*------------------------------------------------------------------------------
                                 solve()
------------------------------------------------------------------------------*/
/* Kristallfeldhamiltonian loesen */
EWPROBLEM *solve(SETUP        *setup,EWPROBLEM    *ewproblem,INT overwrite, KRISTALLFELD *kristallfeld,CHAR modus)
{
    EWPROBLEM *diagonalisiere(EWPROBLEM *ewproblem,MATRIX *matrix,INT overwrite,SETUP *setup);
 
    ITERATION *iteration;
    ITERATION *Vkq0( ITERATION *,INT);
    ITERATION *Vkq1( ITERATION *,INT);
    ITERATION *Vkq2( ITERATION *,INT);
    ITERATION *Vkq3( ITERATION *,INT);
    ITERATION *Vkq4( ITERATION *,INT);
    ITERATION *Vkq5( ITERATION *,INT);
    ITERATION *Vkq6( ITERATION *,INT);
    ITERATION *Vkq7( ITERATION *,INT);
    ITERATION *Vkq8( ITERATION *,INT);
 
    ITERATION *hamltn0(ITERATION *);
    ITERATION *hamltn1(ITERATION *);
    ITERATION *hamltn2(ITERATION *);
    ITERATION *hamltn3(ITERATION *);
    ITERATION *hamltn4(ITERATION *);
    ITERATION *hamltn5(ITERATION *);
    ITERATION *hamltn6(ITERATION *);
    ITERATION *hamltn7(ITERATION *);
    ITERATION *hamltn8(ITERATION *);
 
    int sym;
 
    iteration = ITERATION(kristallfeld);
    sym       =  SYMMETRIENR(kristallfeld);
    switch( sym  ){
 
         case 0 : switch(modus){
                     case SIN: iteration = Vkq0(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn0(iteration );
                  }
                  break;
 
         case 1 : switch(modus){
                     case SIN: iteration = Vkq1(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn1(iteration );
                  }
                  break;
         case 2 : switch(modus){
                     case SIN: iteration = Vkq2(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn2(iteration );
                  }
                  break;
 
         case 3 : switch(modus){
                     case SIN: iteration = Vkq3(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn3(iteration );
                  }
                  break;
 
         case 4 : switch(modus){
                     case SIN: iteration = Vkq4(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn4(iteration );
                  }
                  break;
 
         case 5 : switch(modus){
                     case SIN: iteration = Vkq5(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn5(iteration );
                  }
                  break;
 
         case 6 : switch(modus){
                     case SIN: iteration = Vkq6(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn6(iteration );
                  }
                  break;
 
         case 7 : switch(modus){
                     case SIN: iteration = Vkq7(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn7(iteration );
                  }
                  break;
 
         case 8 : switch(modus){
                     case SIN: iteration = Vkq8(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn8(iteration );
                  }
                  break;
    }
/* if( overwrite != NOSPACE )
    printf("diagonalising the hamiltonian...\n");*/
 
 if( overwrite == NOSPACE ){
    return(  diagonalisiere( ewproblem,HAMILTONIAN(iteration),NOSPACE,setup ));
 }
 else
    return(  diagonalisiere( (EWPROBLEM*)0,HAMILTONIAN(iteration),NEIN,setup ));
}

/*------------------------------------------------------------------------------
                                 setuphcf()
------------------------------------------------------------------------------*/
/* Kristallfeldhamiltonian loesen */
EWPROBLEM *setuphcf(   EWPROBLEM    *ewproblem,INT overwrite, KRISTALLFELD *kristallfeld,CHAR modus)
{
    UNUSED_PARAMETER(overwrite);

    EWPROBLEM *diagonalisiere(EWPROBLEM *ewproblem,MATRIX *matrix,INT overwrite,SETUP *setup);
 
    ITERATION *iteration;
    ITERATION *Vkq0( ITERATION *,INT);
    ITERATION *Vkq1( ITERATION *,INT);
    ITERATION *Vkq2( ITERATION *,INT);
    ITERATION *Vkq3( ITERATION *,INT);
    ITERATION *Vkq4( ITERATION *,INT);
    ITERATION *Vkq5( ITERATION *,INT);
    ITERATION *Vkq6( ITERATION *,INT);
    ITERATION *Vkq7( ITERATION *,INT);
    ITERATION *Vkq8( ITERATION *,INT);
 
    ITERATION *hamltn0(ITERATION *);
    ITERATION *hamltn1(ITERATION *);
    ITERATION *hamltn2(ITERATION *);
    ITERATION *hamltn3(ITERATION *);
    ITERATION *hamltn4(ITERATION *);
    ITERATION *hamltn5(ITERATION *);
    ITERATION *hamltn6(ITERATION *);
    ITERATION *hamltn7(ITERATION *);
    ITERATION *hamltn8(ITERATION *);
 
    int sym;
 
    iteration = ITERATION(kristallfeld);
    sym       =  SYMMETRIENR(kristallfeld);
    switch( sym  ){
 
         case 0 : switch(modus){
                     case SIN: iteration = Vkq0(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn0(iteration );
                  }
                  break;
 
         case 1 : switch(modus){
                     case SIN: iteration = Vkq1(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn1(iteration );
                  }
                  break;
         case 2 : switch(modus){
                     case SIN: iteration = Vkq2(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn2(iteration );
                  }
                  break;
 
         case 3 : switch(modus){
                     case SIN: iteration = Vkq3(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn3(iteration );
                  }
                  break;
 
         case 4 : switch(modus){
                     case SIN: iteration = Vkq4(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn4(iteration );
                  }
                  break;
 
         case 5 : switch(modus){
                     case SIN: iteration = Vkq5(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn5(iteration );
                  }
                  break;
 
         case 6 : switch(modus){
                     case SIN: iteration = Vkq6(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn6(iteration );
                  }
                  break;
 
         case 7 : switch(modus){
                     case SIN: iteration = Vkq7(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn7(iteration );
                  }
                  break;
 
         case 8 : switch(modus){
                     case SIN: iteration = Vkq8(iteration,sym);
                     case XW :
                     case AKQ:
                     case BKQ:
                     case DKQ:
                     case LKQ:
                     case VKQ:
                     case WKQ: iteration = hamltn8(iteration );
                  }
                  break;
    }
return ewproblem;
}

/*------------------------------------------------------------------------------
                                 show_()
------------------------------------------------------------------------------*/
void show_( EWPROBLEM *ew)
{
    info_ewproblem( ew );
}
 
/*------------------------------------------------------------------------------
                               hamltn0()
------------------------------------------------------------------------------*/
ITERATION *hamltn0( ITERATION *i)
{
    STEVENS *s;
    MATRIX  *h,*mag;
    INT     n,m;
 
    h    = HAMILTONIAN(i); /* Speicher fuer GesamtHamiltonian holen */
    s    = PKQ( i);        /* (Jn|Pkq(J)|mJ) holen */
    mag  = HMAG(i);        /* (Jn|Hmag|mJ holen  */
 
 
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
 
              R(h,n,m) += RT( V20(i) ) * R( P2P0(s),n,m );
              R(h,n,m) += RT( V40(i) ) * R( P4P0(s),n,m );
              R(h,n,m) += RT( V60(i) ) * R( P6P0(s),n,m );
 
              R(h,n,m) += RT(V21(i))*( R(P2P1(s),n,m)+R(P2M1(s),n,m) );
              I(h,n,m) += IT(V21(i))*( R(P2P1(s),n,m)-R(P2M1(s),n,m) );
 
              R(h,n,m) += RT(V22(i))*( R(P2P2(s),n,m)+R(P2M2(s),n,m) );
              I(h,n,m) += IT(V22(i))*( R(P2P2(s),n,m)-R(P2M2(s),n,m) );
 
 
 
              R(h,n,m) += RT(V41(i))*( R(P4P1(s),n,m)+R(P4M1(s),n,m) );
              I(h,n,m) += IT(V41(i))*( R(P4P1(s),n,m)-R(P4M1(s),n,m) );
 
              R(h,n,m) += RT(V42(i))*( R(P4P2(s),n,m)+R(P4M2(s),n,m) );
              I(h,n,m) += IT(V42(i))*( R(P4P2(s),n,m)-R(P4M2(s),n,m) );
 
              R(h,n,m) += RT(V43(i))*( R(P4P3(s),n,m)+R(P4M3(s),n,m) );
              I(h,n,m) += IT(V43(i))*( R(P4P3(s),n,m)-R(P4M3(s),n,m) );
 
              R(h,n,m) += RT(V44(i))*( R(P4P4(s),n,m)+R(P4M4(s),n,m) );
              I(h,n,m) += IT(V44(i))*( R(P4P4(s),n,m)-R(P4M4(s),n,m) );
 
 
 
 
              R(h,n,m) += RT(V61(i))*( R(P6P1(s),n,m)+R(P6M1(s),n,m) );
              I(h,n,m) += IT(V61(i))*( R(P6P1(s),n,m)-R(P6M1(s),n,m) );
 
              R(h,n,m) += RT(V62(i))*( R(P6P2(s),n,m)+R(P6M2(s),n,m) );
              I(h,n,m) += IT(V62(i))*( R(P6P2(s),n,m)-R(P6M2(s),n,m) );
 
              R(h,n,m) += RT(V63(i))*( R(P6P3(s),n,m)+R(P6M3(s),n,m) );
              I(h,n,m) += IT(V63(i))*( R(P6P3(s),n,m)-R(P6M3(s),n,m) );
 
              R(h,n,m) += RT(V64(i))*( R(P6P4(s),n,m)+R(P6M4(s),n,m) );
              I(h,n,m) += IT(V64(i))*( R(P6P4(s),n,m)-R(P6M4(s),n,m) );
 
              R(h,n,m) += RT(V65(i))*( R(P6P5(s),n,m)+R(P6M5(s),n,m) );
              I(h,n,m) += IT(V65(i))*( R(P6P5(s),n,m)-R(P6M5(s),n,m) );
 
              R(h,n,m) += RT(V66(i))*( R(P6P6(s),n,m)+R(P6M6(s),n,m) );
              I(h,n,m) += IT(V66(i))*( R(P6P6(s),n,m)-R(P6M6(s),n,m) );
 
 
         }
 
    /* h += h             */
    /*       Magnetfeld   */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) += R(mag,n,m);
              I(h,n,m) += I(mag,n,m);
         }
 
    /* h += singleion anisotropy */
    if( B1S(i)!=0.0 || B2S(i)!=0.0 || B3S(i)!=0.0 ) {
       #include "define_j.c"          /* mj,J2,J+,... definieren */
       MATRIX  *dx,*dy,*dz;
       DOUBLE d1=sqrt(fabs(B1S(i))),d2=sqrt(fabs(B2S(i))),d3=sqrt(fabs(B3S(i))),jm,jp,s1=1.,s2=1.,s3=1.;
       INT dimj=DIMJ(i),l; 
       if(B1S(i)<0) {s1=-1.;} if(B2S(i)<0) {s2=-1.;} if(B3S(i)<0) {s3=-1.;}
       dx = mx_alloc( dimj,dimj ); dy = mx_alloc( dimj,dimj ); dz = mx_alloc( dimj,dimj );
       for( n=DIMJ(i) ; n>=1 ; --n) for( m=DIMJ(i) ; m>=1 ; --m){
              jm=JM(mj)*D(nj,mj-1); jp=JP(mj)*D(nj,mj+1);
              R(dx,n,m) = d1*0.5*( jm+jp ); I(dy,n,m) = d2*0.5*( jm-jp ); R(dz,n,m) = d3*mj*D(nj,mj); }

       for( n=DIMJ(i) ; n>=1 ; --n ) for( m=DIMJ(i) ; m>=1 ; --m ) for( l=DIMJ(i) ; l>=1 ; --l){
              R(h,n,m) += ( s1*R(dx,n,l)*R(dx,l,m) - s2*I(dy,n,l)*I(dy,l,m) + s3*R(dz,n,l)*R(dz,l,m) ); }
       free(dx); free(dy); free(dz);
    }
if( B1SS(i)!=0.0 || B2SS(i)!=0.0 || B3SS(i)!=0.0 ) {
       #include "define_j.c"          /* mj,J2,J+,... definieren */
       MATRIX  *dx,*dy,*dz;
       MATRIX  *dxdx,*dydy,*dzdz;
       DOUBLE d1=sqrt(sqrt(fabs(B1SS(i)))),d2=sqrt(sqrt(fabs(B2SS(i)))),d3=sqrt(sqrt(fabs(B3SS(i)))),jm,jp,s1=1.,s2=1.,s3=1.;
       INT dimj=DIMJ(i),l; 
       if(B1SS(i)<0) {s1=-1.;} if(B2SS(i)<0) {s2=-1.;} if(B3SS(i)<0) {s3=-1.;}
       dx = mx_alloc( dimj,dimj ); dy = mx_alloc( dimj,dimj ); dz = mx_alloc( dimj,dimj );
       dxdx = mx_alloc( dimj,dimj ); dydy = mx_alloc( dimj,dimj ); dzdz = mx_alloc( dimj,dimj );
       for( n=DIMJ(i) ; n>=1 ; --n) for( m=DIMJ(i) ; m>=1 ; --m){
              jm=JM(mj)*D(nj,mj-1); jp=JP(mj)*D(nj,mj+1);
              R(dx,n,m) = d1*0.5*( jm+jp ); I(dy,n,m) = d2*0.5*( jm-jp ); R(dz,n,m) = d3*mj*D(nj,mj); }
       for( n=DIMJ(i) ; n>=1 ; --n ) for( m=DIMJ(i) ; m>=1 ; --m ) for( l=DIMJ(i) ; l>=1 ; --l){
              R(dxdx,n,m) += R(dx,n,l)*R(dx,l,m);
              R(dydy,n,m) -= I(dy,n,l)*I(dy,l,m); 
              R(dzdz,n,m) += R(dz,n,l)*R(dz,l,m); }

       for( n=DIMJ(i) ; n>=1 ; --n ) for( m=DIMJ(i) ; m>=1 ; --m ) for( l=DIMJ(i) ; l>=1 ; --l){
              R(h,n,m) += ( s1*R(dxdx,n,l)*R(dxdx,l,m) + s2*R(dydy,n,l)*R(dydy,l,m) + s3*R(dzdz,n,l)*R(dzdz,l,m) ); }
       free(dx); free(dy); free(dz);
       free(dxdx); free(dydy); free(dzdz);
    }
 
 
    return( i );
}
/*------------------------------------------------------------------------------
                               hamltn1()
------------------------------------------------------------------------------*/
ITERATION *hamltn1( ITERATION *i)
{
    STEVENS *s;
    MATRIX  *h,*mag;
    INT     n,m;
 
    h    = HAMILTONIAN(i); /* Speicher fuer GesamtHamiltonian holen */
    s    = PKQ( i);        /* (Jn|Pkq(J)|mJ) holen */
    mag  = HMAG(i);        /* (Jn|Hmag|mJ holen  */
 
 
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
 
              R(h,n,m) += RT( V20(i) ) * R( P2P0(s),n,m );
              R(h,n,m) += RT( V40(i) ) * R( P4P0(s),n,m );
              R(h,n,m) += RT( V60(i) ) * R( P6P0(s),n,m );
 
 
              R(h,n,m) += RT(V22(i))*( R(P2M2(s),n,m)+R(P2P2(s),n,m) );
              I(h,n,m) += IT(V22(i))*( R(P2M2(s),n,m)-R(P2P2(s),n,m) );
 
              R(h,n,m) += RT(V42(i))*( R(P4P2(s),n,m)+R(P4M2(s),n,m) );
              I(h,n,m) += IT(V42(i))*( R(P4P2(s),n,m)-R(P4M2(s),n,m) );
 
              R(h,n,m) += RT(V44(i))*( R(P4P4(s),n,m)+R(P4M4(s),n,m) );
              I(h,n,m) += IT(V44(i))*( R(P4P4(s),n,m)-R(P4M4(s),n,m) );
 
 
 
              R(h,n,m) += RT(V62(i))*( R(P6P2(s),n,m)+R(P6M2(s),n,m) );
              I(h,n,m) += IT(V62(i))*( R(P6P2(s),n,m)-R(P6M2(s),n,m) );
 
              R(h,n,m) += RT(V64(i))*( R(P6P4(s),n,m)+R(P6M4(s),n,m) );
              I(h,n,m) += IT(V64(i))*( R(P6P4(s),n,m)-R(P6M4(s),n,m) );
 
              R(h,n,m) += RT(V66(i))*( R(P6P6(s),n,m)+R(P6M6(s),n,m) );
              I(h,n,m) += IT(V66(i))*( R(P6P6(s),n,m)-R(P6M6(s),n,m) );
 
 
         }
 
    /* h += h             */
    /*       Magnetfeld   */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) += R(mag,n,m);
              I(h,n,m) += I(mag,n,m);
         }
 
 
    return( i );
}
/*------------------------------------------------------------------------------
                               hamltn2()
------------------------------------------------------------------------------*/
ITERATION *hamltn2( ITERATION *i)
{
    STEVENS *s;
    MATRIX  *h,*mag;
    INT     n,m;
 
    h    = HAMILTONIAN(i); /* Speicher fuer GesamtHamiltonian holen */
    s    = PKQ( i);        /* (Jn|Pkq(J)|mJ) holen */
    mag  = HMAG(i);        /* (Jn|Hmag|mJ holen  */
 
 
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
 
              R(h,n,m) += RT( V20(i) ) * R( P2P0(s),n,m );
              R(h,n,m) += RT( V40(i) ) * R( P4P0(s),n,m );
              R(h,n,m) += RT( V60(i) ) * R( P6P0(s),n,m );
 
              if(VOR22(i)==1.0){
                R(h,n,m) += RT(V22(i))*( R(P2M2(s),n,m)+R(P2P2(s),n,m) );
                I(h,n,m) += IT(V22(i))*( R(P2M2(s),n,m)-R(P2P2(s),n,m) );
              }else{
                I(h,n,m) += RT(V22(i))*( R(P2M2(s),n,m)-R(P2P2(s),n,m) );
              }
              if(VOR42(i)==1.0){
                R(h,n,m) += RT(V42(i))*( R(P4M2(s),n,m)+R(P4P2(s),n,m) );
                I(h,n,m) += IT(V42(i))*( R(P4M2(s),n,m)-R(P4P2(s),n,m) );
              }else{
                I(h,n,m) += RT(V42(i))*( R(P4M2(s),n,m)-R(P4P2(s),n,m) );
              }
              if(VOR44(i)==1.0){
                R(h,n,m) += RT(V44(i))*( R(P4M4(s),n,m)+R(P4P4(s),n,m) );
                I(h,n,m) += IT(V44(i))*( R(P4M4(s),n,m)-R(P4P4(s),n,m) );
              }else{
                I(h,n,m) += RT(V44(i))*( R(P4M4(s),n,m)-R(P4P4(s),n,m) );
              }
 
              if(VOR62(i)==1.0){
                R(h,n,m) += RT(V62(i))*( R(P6M2(s),n,m)+R(P6P2(s),n,m) );
                I(h,n,m) += IT(V62(i))*( R(P6M2(s),n,m)-R(P6P2(s),n,m) );
              }else{
                I(h,n,m) += RT(V62(i))*( R(P6M2(s),n,m)-R(P6P2(s),n,m) );
              }
              if(VOR64(i)==1.0){
                R(h,n,m) += RT(V64(i))*( R(P6M4(s),n,m)+R(P6P4(s),n,m) );
                I(h,n,m) += IT(V64(i))*( R(P6M4(s),n,m)-R(P6P4(s),n,m) );
              }else{
                I(h,n,m) += RT(V64(i))*( R(P6M4(s),n,m)-R(P6P4(s),n,m) );
              }
              if(VOR66(i)==1.0){
                R(h,n,m) += RT(V66(i))*( R(P6M6(s),n,m)+R(P6P6(s),n,m) );
                I(h,n,m) += IT(V66(i))*( R(P6M6(s),n,m)-R(P6P6(s),n,m) );
              }else{
                I(h,n,m) += RT(V66(i))*( R(P6M6(s),n,m)-R(P6P6(s),n,m) );
              }
 
 
         }
 
    /* h += h             */
    /*       Magnetfeld   */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) += R(mag,n,m);
              I(h,n,m) += I(mag,n,m);
         }
 
 
    return( i );
}
/*------------------------------------------------------------------------------
                               hamltn3()
------------------------------------------------------------------------------*/
ITERATION *hamltn3( ITERATION *i)
{
    STEVENS *s;
    MATRIX  *h,*mag;
    INT     n,m;
 
    h    = HAMILTONIAN(i); /* Speicher fuer GesamtHamiltonian holen */
    s    = PKQ( i);        /* (Jn|Pkq(J)|mJ) holen */
    mag  = HMAG(i);        /* (Jn|Hmag|mJ holen  */
 
 
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
 
              R(h,n,m) += RT( V20(i) ) * R( P2P0(s),n,m );
              R(h,n,m) += RT( V40(i) ) * R( P4P0(s),n,m );
              R(h,n,m) += RT( V60(i) ) * R( P6P0(s),n,m );
 
 
              if(VOR44(i)==1.0){
                R(h,n,m) += RT(V44(i))*( R(P4M4(s),n,m)+R(P4P4(s),n,m) );
                I(h,n,m) += IT(V44(i))*( R(P4M4(s),n,m)-R(P4P4(s),n,m) );
              }else{
                I(h,n,m) += RT(V44(i))*( R(P4M4(s),n,m)-R(P4P4(s),n,m) );
              }
 
 
 
              R(h,n,m) += RT(V64(i))*( R(P6P4(s),n,m)+R(P6M4(s),n,m) );
              I(h,n,m) += IT(V64(i))*( R(P6P4(s),n,m)-R(P6M4(s),n,m) );
 
 
         }
 
    /* h += h             */
    /*       Magnetfeld   */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) += R(mag,n,m);
              I(h,n,m) += I(mag,n,m);
         }
 
 
    return( i );
}
/*------------------------------------------------------------------------------
                               hamltn4()
------------------------------------------------------------------------------*/
ITERATION *hamltn4( ITERATION *i)
{
    STEVENS *s;
    MATRIX  *h,*mag;
    INT     n,m;
 
    h    = HAMILTONIAN(i); /* Speicher fuer GesamtHamiltonian holen */
    s    = PKQ( i);        /* (Jn|Pkq(J)|mJ) holen */
    mag  = HMAG(i);        /* (Jn|Hmag|mJ holen  */
 
 
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
 
              R(h,n,m) += RT( V20(i) ) * R( P2P0(s),n,m );
              R(h,n,m) += RT( V40(i) ) * R( P4P0(s),n,m );
              R(h,n,m) += RT( V60(i) ) * R( P6P0(s),n,m );
 
 
              if(VOR44(i)==1.0){
                R(h,n,m) += RT(V44(i))*( R(P4M4(s),n,m)+R(P4P4(s),n,m) );
                I(h,n,m) += IT(V44(i))*( R(P4M4(s),n,m)-R(P4P4(s),n,m) );
              }else{
                I(h,n,m) += RT(V44(i))*( R(P4M4(s),n,m)-R(P4P4(s),n,m) );
              }
 
              if(VOR64(i)==1.0){
                R(h,n,m) += RT(V64(i))*( R(P6M4(s),n,m)+R(P6P4(s),n,m) );
                I(h,n,m) += IT(V64(i))*( R(P6M4(s),n,m)-R(P6P4(s),n,m) );
              }else{
                I(h,n,m) += RT(V64(i))*( R(P6M4(s),n,m)-R(P6P4(s),n,m) );
              }
 
 
         }
 
    /* h += h             */
    /*       Magnetfeld   */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) += R(mag,n,m);
              I(h,n,m) += I(mag,n,m);
         }
 
 
    return( i );
}
/*------------------------------------------------------------------------------
                               hamltn5()
------------------------------------------------------------------------------*/
ITERATION *hamltn5( ITERATION *i)
{
    STEVENS *s;
    MATRIX  *h,*mag;
    INT     n,m;
 
    h    = HAMILTONIAN(i); /* Speicher fuer GesamtHamiltonian holen */
    s    = PKQ( i);        /* (Jn|Pkq(J)|mJ) holen */
    mag  = HMAG(i);        /* (Jn|Hmag|mJ holen  */
 
 
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
 
              R(h,n,m) += RT( V20(i) ) * R( P2P0(s),n,m );
              R(h,n,m) += RT( V40(i) ) * R( P4P0(s),n,m );
              R(h,n,m) += RT( V60(i) ) * R( P6P0(s),n,m );
 
              if(VOR43(i)==1.0){
                R(h,n,m) += RT(V43(i))*( R(P4M3(s),n,m)+R(P4P3(s),n,m) );
                I(h,n,m) += IT(V43(i))*( R(P4M3(s),n,m)-R(P4P3(s),n,m) );
              }else{
                I(h,n,m) += RT(V43(i))*( R(P4M3(s),n,m)-R(P4P3(s),n,m) );
              }
 
 
              R(h,n,m) += RT(V63(i))*( R(P6P3(s),n,m)+R(P6M3(s),n,m) );
              I(h,n,m) += IT(V63(i))*( R(P6P3(s),n,m)-R(P6M3(s),n,m) );
 
              R(h,n,m) += RT(V66(i))*( R(P6P6(s),n,m)+R(P6M6(s),n,m) );
              I(h,n,m) += IT(V66(i))*( R(P6P6(s),n,m)-R(P6M6(s),n,m) );
 
 
         }
 
    /* h += h             */
    /*       Magnetfeld   */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) += R(mag,n,m);
              I(h,n,m) += I(mag,n,m);
         }
 
 
    return( i );
}
/*------------------------------------------------------------------------------
                               hamltn6()
------------------------------------------------------------------------------*/
ITERATION *hamltn6( ITERATION *i)
{
    STEVENS *s;
    MATRIX  *h,*mag;
    INT     n,m;
 
    h    = HAMILTONIAN(i); /* Speicher fuer GesamtHamiltonian holen */
    s    = PKQ( i);        /* (Jn|Pkq(J)|mJ) holen */
    mag  = HMAG(i);        /* (Jn|Hmag|mJ holen  */
 
 
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
 
              R(h,n,m) += RT( V20(i) ) * R( P2P0(s),n,m );
              R(h,n,m) += RT( V40(i) ) * R( P4P0(s),n,m );
              R(h,n,m) += RT( V60(i) ) * R( P6P0(s),n,m );
 
 
              if(VOR43(i)==1.0){
                R(h,n,m) += RT(V43(i))*( R(P4M3(s),n,m)+R(P4P3(s),n,m) );
                I(h,n,m) += IT(V43(i))*( R(P4M3(s),n,m)-R(P4P3(s),n,m) );
              }else{
                I(h,n,m) += RT(V43(i))*( R(P4M3(s),n,m)-R(P4P3(s),n,m) );
              }
              if(VOR63(i)==1.0){
                R(h,n,m) += RT(V63(i))*( R(P6M3(s),n,m)+R(P6P3(s),n,m) );
                I(h,n,m) += IT(V63(i))*( R(P6M3(s),n,m)-R(P6P3(s),n,m) );
              }else{
                I(h,n,m) += RT(V63(i))*( R(P6M3(s),n,m)-R(P6P3(s),n,m) );
              }
              if(VOR66(i)==1.0){
                R(h,n,m) += RT(V66(i))*( R(P6M6(s),n,m)+R(P6P6(s),n,m) );
                I(h,n,m) += IT(V66(i))*( R(P6M6(s),n,m)-R(P6P6(s),n,m) );
              }else{
                I(h,n,m) += RT(V66(i))*( R(P6M6(s),n,m)-R(P6P6(s),n,m) );
              }
 
 
         }
 
    /* h += h             */
    /*       Magnetfeld   */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) += R(mag,n,m);
              I(h,n,m) += I(mag,n,m);
         }
 
 
    return( i );
}
/*------------------------------------------------------------------------------
                               hamltn7()
------------------------------------------------------------------------------*/
ITERATION *hamltn7( ITERATION *i)
{
    STEVENS *s;
    MATRIX  *h,*mag;
    INT     n,m;
 
    h    = HAMILTONIAN(i); /* Speicher fuer GesamtHamiltonian holen */
    s    = PKQ( i);        /* (Jn|Pkq(J)|mJ) holen */
    mag  = HMAG(i);        /* (Jn|Hmag|mJ holen  */
 
 
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
 
              R(h,n,m) += RT( V20(i) ) * R( P2P0(s),n,m );
              R(h,n,m) += RT( V40(i) ) * R( P4P0(s),n,m );
              R(h,n,m) += RT( V60(i) ) * R( P6P0(s),n,m );
 
              if(VOR66(i)==1.0){
                R(h,n,m) += RT(V66(i))*( R(P6M6(s),n,m)+R(P6P6(s),n,m) );
                I(h,n,m) += IT(V66(i))*( R(P6M6(s),n,m)-R(P6P6(s),n,m) );
              }else{
                I(h,n,m) += RT(V66(i))*( R(P6M6(s),n,m)-R(P6P6(s),n,m) );
              }
 
 
         }
 
    /* h += h             */
    /*       Magnetfeld   */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) += R(mag,n,m);
              I(h,n,m) += I(mag,n,m);
         }
 
 
    return( i );
}
/*------------------------------------------------------------------------------
                               hamltn8()
------------------------------------------------------------------------------*/
ITERATION *hamltn8( ITERATION *i)
{
    STEVENS *s;
    MATRIX  *h,*mag;
    INT     n,m;
 
    h    = HAMILTONIAN(i); /* Speicher fuer GesamtHamiltonian holen */
    s    = PKQ( i);        /* (Jn|Pkq(J)|mJ) holen */
    mag  = HMAG(i);        /* (Jn|Hmag|mJ holen  */
 
 
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
 
              R(h,n,m) += RT( V40(i) ) * R( P4P0(s),n,m );
              R(h,n,m) += RT( V60(i) ) * R( P6P0(s),n,m );
 
              R(h,n,m) += RT(V44(i))*( R(P4P4(s),n,m)+R(P4M4(s),n,m) );
 
              R(h,n,m) += RT(V64(i))*( R(P6P4(s),n,m)+R(P6M4(s),n,m) );
 
 
         }
 
    /* h += h             */
    /*       Magnetfeld   */
    for( n=DIMJ(i) ; n>=1 ; --n )
         for( m=DIMJ(i) ; m>=1 ; --m ){
              R(h,n,m) += R(mag,n,m);
              I(h,n,m) += I(mag,n,m);
         }
 
 
    return( i );
}
/*------------------------------------------------------------------------------
                                 Vkq0()
------------------------------------------------------------------------------*/
ITERATION *Vkq0(ITERATION *iter,INT sym)
{
    INT i ;
    ITERATION *auswahlregel(ITERATION *,INT);
 
    /*                                  */
    /*  Q_R   = Q(R)                    */
    /*                     *            */
    /*  CkqFS = beta  *   C             */
    /*              kq      kq          */
    /*                                  */
    /*   *                              */
    /*  C     = RT(C   ) -i* IT(C   )   */
    /*   kq         kq           kq     */
    /*                                  */
    /*                        k         */
    /*                      <r >        */
    /* Rk_k+1 = theta(J) * -------      */
    /*               k        k+1       */
    /*                       R          */
    /*                                  */
 
    RT( V20(iter) ) = 0.0;  IT( V20(iter) ) = 0.0;
    RT( V21(iter) ) = 0.0;  IT( V21(iter) ) = 0.0;
    RT( V22(iter) ) = 0.0;  IT( V22(iter) ) = 0.0;
 
    RT( V40(iter) ) = 0.0;  IT( V40(iter) ) = 0.0;
    RT( V41(iter) ) = 0.0;  IT( V41(iter) ) = 0.0;
    RT( V42(iter) ) = 0.0;  IT( V42(iter) ) = 0.0;
    RT( V43(iter) ) = 0.0;  IT( V43(iter) ) = 0.0;
    RT( V44(iter) ) = 0.0;  IT( V44(iter) ) = 0.0;
 
    RT( V60(iter) ) = 0.0;  IT( V60(iter) ) = 0.0;
    RT( V61(iter) ) = 0.0;  IT( V61(iter) ) = 0.0;
    RT( V62(iter) ) = 0.0;  IT( V62(iter) ) = 0.0;
    RT( V63(iter) ) = 0.0;  IT( V63(iter) ) = 0.0;
    RT( V64(iter) ) = 0.0;  IT( V64(iter) ) = 0.0;
    RT( V65(iter) ) = 0.0;  IT( V65(iter) ) = 0.0;
    RT( V66(iter) ) = 0.0;  IT( V66(iter) ) = 0.0;
 
    for( i=1 ; i<= ANZ_NN(iter) ; ++i ){
 
         RT(V20(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C20FS(iter,i));
         RT(V21(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C21FS(iter,i));
         RT(V22(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C22FS(iter,i));
 
         RT(V40(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C40FS(iter,i));
         RT(V41(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C41FS(iter,i));
         RT(V42(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C42FS(iter,i));
         RT(V43(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C43FS(iter,i));
         RT(V44(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C44FS(iter,i));
 
         RT(V60(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C60FS(iter,i));
         RT(V61(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C61FS(iter,i));
         RT(V62(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C62FS(iter,i));
         RT(V64(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C64FS(iter,i));
         RT(V65(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C65FS(iter,i));
         RT(V66(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C66FS(iter,i));
 
 
         IT(V22(iter)) += Q_R(iter,i) * R2_3(iter,i) * IT(C22FS(iter,i));
 
         IT(V41(iter)) += Q_R(iter,i) * R4_5(iter,i) * IT(C41FS(iter,i));
         IT(V42(iter)) += Q_R(iter,i) * R4_5(iter,i) * IT(C42FS(iter,i));
         IT(V43(iter)) += Q_R(iter,i) * R4_5(iter,i) * IT(C43FS(iter,i));
         IT(V44(iter)) += Q_R(iter,i) * R4_5(iter,i) * IT(C44FS(iter,i));
 
         IT(V61(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C61FS(iter,i));
         IT(V62(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C62FS(iter,i));
         IT(V63(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C63FS(iter,i));
         IT(V64(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C64FS(iter,i));
         IT(V65(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C65FS(iter,i));
         IT(V66(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C66FS(iter,i));
    }
 
 
    return( auswahlregel(iter,sym) );
}
 
/*------------------------------------------------------------------------------
                                 Vkq1()
------------------------------------------------------------------------------*/
ITERATION *Vkq1(ITERATION *iter,INT sym)
{
    INT i;
    ITERATION *auswahlregel(ITERATION *,INT);
 
    /*                                  */
    /*  Q_R   = Q(R)                    */
    /*                     *            */
    /*  CkqFS = beta  *   C             */
    /*              kq      kq          */
    /*                                  */
    /*   *                              */
    /*  C     = RT(C   ) -i* IT(C   )   */
    /*   kq         kq           kq     */
    /*                                  */
    /*                        k         */
    /*                      <r >        */
    /* Rk_k+1 = theta(J) * -------      */
    /*               k        k+1       */
    /*                       R          */
    /*                                  */
 
    RT( V20(iter) ) = 0.0;  IT( V20(iter) ) = 0.0;
    RT( V21(iter) ) = 0.0;  IT( V21(iter) ) = 0.0;
    RT( V22(iter) ) = 0.0;  IT( V22(iter) ) = 0.0;
 
    RT( V40(iter) ) = 0.0;  IT( V40(iter) ) = 0.0;
    RT( V41(iter) ) = 0.0;  IT( V41(iter) ) = 0.0;
    RT( V42(iter) ) = 0.0;  IT( V42(iter) ) = 0.0;
    RT( V43(iter) ) = 0.0;  IT( V43(iter) ) = 0.0;
    RT( V44(iter) ) = 0.0;  IT( V44(iter) ) = 0.0;
 
    RT( V60(iter) ) = 0.0;  IT( V60(iter) ) = 0.0;
    RT( V61(iter) ) = 0.0;  IT( V61(iter) ) = 0.0;
    RT( V62(iter) ) = 0.0;  IT( V62(iter) ) = 0.0;
    RT( V63(iter) ) = 0.0;  IT( V63(iter) ) = 0.0;
    RT( V64(iter) ) = 0.0;  IT( V64(iter) ) = 0.0;
    RT( V65(iter) ) = 0.0;  IT( V65(iter) ) = 0.0;
    RT( V66(iter) ) = 0.0;  IT( V66(iter) ) = 0.0;
 
    for( i=1 ; i<= ANZ_NN(iter) ; ++i ){
 
         RT(V20(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C20FS(iter,i));
         RT(V22(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C22FS(iter,i));
 
         RT(V40(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C40FS(iter,i));
         RT(V42(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C42FS(iter,i));
         RT(V44(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C44FS(iter,i));
 
         RT(V60(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C60FS(iter,i));
         RT(V62(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C62FS(iter,i));
         RT(V64(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C64FS(iter,i));
         RT(V66(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C66FS(iter,i));
 
 
         IT(V42(iter)) += Q_R(iter,i) * R4_5(iter,i) * IT(C42FS(iter,i));
         IT(V44(iter)) += Q_R(iter,i) * R4_5(iter,i) * IT(C44FS(iter,i));
 
         IT(V62(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C62FS(iter,i));
         IT(V64(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C64FS(iter,i));
         IT(V66(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C66FS(iter,i));
    }
    return( auswahlregel(iter,sym) );
}
/*------------------------------------------------------------------------------
                                 Vkq2()
------------------------------------------------------------------------------*/
ITERATION *Vkq2(ITERATION *iter,INT sym)
{
    INT i;
    ITERATION *auswahlregel(ITERATION *,INT);
 
    /*                                  */
    /*  Q_R   = Q(R)                    */
    /*                     *            */
    /*  CkqFS = beta  *   C             */
    /*              kq      kq          */
    /*                                  */
    /*   *                              */
    /*  C     = RT(C   ) -i* IT(C   )   */
    /*   kq         kq           kq     */
    /*                                  */
    /*                        k         */
    /*                      <r >        */
    /* Rk_k+1 = theta(J) * -------      */
    /*               k        k+1       */
    /*                       R          */
    /*                                  */
 
 
    RT( V20(iter) ) = 0.0;  IT( V20(iter) ) = 0.0;
    RT( V21(iter) ) = 0.0;  IT( V21(iter) ) = 0.0;
    RT( V22(iter) ) = 0.0;  IT( V22(iter) ) = 0.0;
 
    RT( V40(iter) ) = 0.0;  IT( V40(iter) ) = 0.0;
    RT( V41(iter) ) = 0.0;  IT( V41(iter) ) = 0.0;
    RT( V42(iter) ) = 0.0;  IT( V42(iter) ) = 0.0;
    RT( V43(iter) ) = 0.0;  IT( V43(iter) ) = 0.0;
    RT( V44(iter) ) = 0.0;  IT( V44(iter) ) = 0.0;
 
    RT( V60(iter) ) = 0.0;  IT( V60(iter) ) = 0.0;
    RT( V61(iter) ) = 0.0;  IT( V61(iter) ) = 0.0;
    RT( V62(iter) ) = 0.0;  IT( V62(iter) ) = 0.0;
    RT( V63(iter) ) = 0.0;  IT( V63(iter) ) = 0.0;
    RT( V64(iter) ) = 0.0;  IT( V64(iter) ) = 0.0;
    RT( V65(iter) ) = 0.0;  IT( V65(iter) ) = 0.0;
    RT( V66(iter) ) = 0.0;  IT( V66(iter) ) = 0.0;
 
 
    for( i=1 ; i<= ANZ_NN(iter) ; ++i ){
 
         RT(V20(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C20FS(iter,i));
         RT(V22(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C22FS(iter,i));
 
         RT(V40(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C40FS(iter,i));
         RT(V42(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C42FS(iter,i));
         RT(V44(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C44FS(iter,i));
 
         RT(V60(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C60FS(iter,i));
         RT(V62(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C62FS(iter,i));
         RT(V64(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C64FS(iter,i));
         RT(V66(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C66FS(iter,i));
 
    }
    return( auswahlregel(iter,sym) );
}
/*------------------------------------------------------------------------------
                                 Vkq3()
------------------------------------------------------------------------------*/
ITERATION *Vkq3(ITERATION *iter,INT sym)
{
    INT i;
    ITERATION *auswahlregel(ITERATION *,INT);
 
    /*                                  */
    /*  Q_R   = Q(R)                    */
    /*                     *            */
    /*  CkqFS = beta  *   C             */
    /*              kq      kq          */
    /*                                  */
    /*   *                              */
    /*  C     = RT(C   ) -i* IT(C   )   */
    /*   kq         kq           kq     */
    /*                                  */
    /*                        k         */
    /*                      <r >        */
    /* Rk_k+1 = theta(J) * -------      */
    /*               k        k+1       */
    /*                       R          */
    /*                                  */
 
 
    RT( V20(iter) ) = 0.0;  IT( V20(iter) ) = 0.0;
    RT( V21(iter) ) = 0.0;  IT( V21(iter) ) = 0.0;
    RT( V22(iter) ) = 0.0;  IT( V22(iter) ) = 0.0;
 
    RT( V40(iter) ) = 0.0;  IT( V40(iter) ) = 0.0;
    RT( V41(iter) ) = 0.0;  IT( V41(iter) ) = 0.0;
    RT( V42(iter) ) = 0.0;  IT( V42(iter) ) = 0.0;
    RT( V43(iter) ) = 0.0;  IT( V43(iter) ) = 0.0;
    RT( V44(iter) ) = 0.0;  IT( V44(iter) ) = 0.0;
 
    RT( V60(iter) ) = 0.0;  IT( V60(iter) ) = 0.0;
    RT( V61(iter) ) = 0.0;  IT( V61(iter) ) = 0.0;
    RT( V62(iter) ) = 0.0;  IT( V62(iter) ) = 0.0;
    RT( V63(iter) ) = 0.0;  IT( V63(iter) ) = 0.0;
    RT( V64(iter) ) = 0.0;  IT( V64(iter) ) = 0.0;
    RT( V65(iter) ) = 0.0;  IT( V65(iter) ) = 0.0;
    RT( V66(iter) ) = 0.0;  IT( V66(iter) ) = 0.0;
 
    for( i=1 ; i<= ANZ_NN(iter) ; ++i ){
 
         RT(V20(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C20FS(iter,i));
 
         RT(V40(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C40FS(iter,i));
         RT(V44(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C44FS(iter,i));
 
         RT(V60(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C60FS(iter,i));
         RT(V64(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C64FS(iter,i));
 
 
         IT(V64(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C64FS(iter,i));
    }
    return( auswahlregel(iter,sym) );
}
/*------------------------------------------------------------------------------
                                 Vkq4()
------------------------------------------------------------------------------*/
ITERATION *Vkq4(ITERATION *iter,INT sym)
{
    INT i;
    ITERATION *auswahlregel(ITERATION *,INT);
 
    /*                                  */
    /*  Q_R   = Q(R)                    */
    /*                     *            */
    /*  CkqFS = beta  *   C             */
    /*              kq      kq          */
    /*                                  */
    /*   *                              */
    /*  C     = RT(C   ) -i* IT(C   )   */
    /*   kq         kq           kq     */
    /*                                  */
    /*                        k         */
    /*                      <r >        */
    /* Rk_k+1 = theta(J) * -------      */
    /*               k        k+1       */
    /*                       R          */
    /*                                  */
 
 
    RT( V20(iter) ) = 0.0;  IT( V20(iter) ) = 0.0;
    RT( V21(iter) ) = 0.0;  IT( V21(iter) ) = 0.0;
    RT( V22(iter) ) = 0.0;  IT( V22(iter) ) = 0.0;
 
    RT( V40(iter) ) = 0.0;  IT( V40(iter) ) = 0.0;
    RT( V41(iter) ) = 0.0;  IT( V41(iter) ) = 0.0;
    RT( V42(iter) ) = 0.0;  IT( V42(iter) ) = 0.0;
    RT( V43(iter) ) = 0.0;  IT( V43(iter) ) = 0.0;
    RT( V44(iter) ) = 0.0;  IT( V44(iter) ) = 0.0;
 
    RT( V60(iter) ) = 0.0;  IT( V60(iter) ) = 0.0;
    RT( V61(iter) ) = 0.0;  IT( V61(iter) ) = 0.0;
    RT( V62(iter) ) = 0.0;  IT( V62(iter) ) = 0.0;
    RT( V63(iter) ) = 0.0;  IT( V63(iter) ) = 0.0;
    RT( V64(iter) ) = 0.0;  IT( V64(iter) ) = 0.0;
    RT( V65(iter) ) = 0.0;  IT( V65(iter) ) = 0.0;
    RT( V66(iter) ) = 0.0;  IT( V66(iter) ) = 0.0;
 
 
    for( i=1 ; i<= ANZ_NN(iter) ; ++i ){
 
         RT(V20(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C20FS(iter,i));
 
         RT(V40(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C40FS(iter,i));
         RT(V44(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C44FS(iter,i));
 
         RT(V60(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C60FS(iter,i));
         RT(V64(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C64FS(iter,i));
 
    }
    return( auswahlregel(iter,sym) );
}
/*------------------------------------------------------------------------------
                                 Vkq5()
------------------------------------------------------------------------------*/
ITERATION *Vkq5(ITERATION *iter,INT sym)
{
    INT i;
    ITERATION *auswahlregel(ITERATION *,INT);
 
    /*                                  */
    /*  Q_R   = Q(R)                    */
    /*                     *            */
    /*  CkqFS = beta  *   C             */
    /*              kq      kq          */
    /*                                  */
    /*   *                              */
    /*  C     = RT(C   ) -i* IT(C   )   */
    /*   kq         kq           kq     */
    /*                                  */
    /*                        k         */
    /*                      <r >        */
    /* Rk_k+1 = theta(J) * -------      */
    /*               k        k+1       */
    /*                       R          */
    /*                                  */
 
 
    RT( V20(iter) ) = 0.0;  IT( V20(iter) ) = 0.0;
    RT( V21(iter) ) = 0.0;  IT( V21(iter) ) = 0.0;
    RT( V22(iter) ) = 0.0;  IT( V22(iter) ) = 0.0;
 
    RT( V40(iter) ) = 0.0;  IT( V40(iter) ) = 0.0;
    RT( V41(iter) ) = 0.0;  IT( V41(iter) ) = 0.0;
    RT( V42(iter) ) = 0.0;  IT( V42(iter) ) = 0.0;
    RT( V43(iter) ) = 0.0;  IT( V43(iter) ) = 0.0;
    RT( V44(iter) ) = 0.0;  IT( V44(iter) ) = 0.0;
 
    RT( V60(iter) ) = 0.0;  IT( V60(iter) ) = 0.0;
    RT( V61(iter) ) = 0.0;  IT( V61(iter) ) = 0.0;
    RT( V62(iter) ) = 0.0;  IT( V62(iter) ) = 0.0;
    RT( V63(iter) ) = 0.0;  IT( V63(iter) ) = 0.0;
    RT( V64(iter) ) = 0.0;  IT( V64(iter) ) = 0.0;
    RT( V65(iter) ) = 0.0;  IT( V65(iter) ) = 0.0;
    RT( V66(iter) ) = 0.0;  IT( V66(iter) ) = 0.0;
 
 
    for( i=1 ; i<= ANZ_NN(iter) ; ++i ){
 
         RT(V20(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C20FS(iter,i));
 
         RT(V40(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C40FS(iter,i));
         RT(V43(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C43FS(iter,i));
 
         RT(V60(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C60FS(iter,i));
         RT(V63(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C63FS(iter,i));
         RT(V66(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C66FS(iter,i));
 
 
         IT(V63(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C63FS(iter,i));
         IT(V66(iter)) += Q_R(iter,i) * R6_7(iter,i) * IT(C66FS(iter,i));
    }
    return( auswahlregel(iter,sym) );
}
/*------------------------------------------------------------------------------
                                 Vkq6()
------------------------------------------------------------------------------*/
ITERATION *Vkq6(ITERATION *iter,INT sym)
{
    INT i;
    ITERATION *auswahlregel(ITERATION *,INT);
 
    /*                                  */
    /*  Q_R   = Q(R)                    */
    /*                     *            */
    /*  CkqFS = beta  *   C             */
    /*              kq      kq          */
    /*                                  */
    /*   *                              */
    /*  C     = RT(C   ) -i* IT(C   )   */
    /*   kq         kq           kq     */
    /*                                  */
    /*                        k         */
    /*                      <r >        */
    /* Rk_k+1 = theta(J) * -------      */
    /*               k        k+1       */
    /*                       R          */
    /*                                  */
 
 
    RT( V20(iter) ) = 0.0;  IT( V20(iter) ) = 0.0;
    RT( V21(iter) ) = 0.0;  IT( V21(iter) ) = 0.0;
    RT( V22(iter) ) = 0.0;  IT( V22(iter) ) = 0.0;
 
    RT( V40(iter) ) = 0.0;  IT( V40(iter) ) = 0.0;
    RT( V41(iter) ) = 0.0;  IT( V41(iter) ) = 0.0;
    RT( V42(iter) ) = 0.0;  IT( V42(iter) ) = 0.0;
    RT( V43(iter) ) = 0.0;  IT( V43(iter) ) = 0.0;
    RT( V44(iter) ) = 0.0;  IT( V44(iter) ) = 0.0;
 
    RT( V60(iter) ) = 0.0;  IT( V60(iter) ) = 0.0;
    RT( V61(iter) ) = 0.0;  IT( V61(iter) ) = 0.0;
    RT( V62(iter) ) = 0.0;  IT( V62(iter) ) = 0.0;
    RT( V63(iter) ) = 0.0;  IT( V63(iter) ) = 0.0;
    RT( V64(iter) ) = 0.0;  IT( V64(iter) ) = 0.0;
    RT( V65(iter) ) = 0.0;  IT( V65(iter) ) = 0.0;
    RT( V66(iter) ) = 0.0;  IT( V66(iter) ) = 0.0;
 
 
    for( i=1 ; i<= ANZ_NN(iter) ; ++i ){
 
         RT(V20(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C20FS(iter,i));
 
         RT(V40(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C40FS(iter,i));
         RT(V43(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C43FS(iter,i));
 
         RT(V60(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C60FS(iter,i));
         RT(V63(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C63FS(iter,i));
         RT(V66(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C66FS(iter,i));
    }
    return( auswahlregel(iter,sym) );
}
/*------------------------------------------------------------------------------
                                Vkq7()
------------------------------------------------------------------------------*/
ITERATION *Vkq7(ITERATION *iter,INT sym)
{
    INT i;
    ITERATION *auswahlregel(ITERATION *,INT);
 
    /*                                  */
    /*  Q_R   = Q(R)                    */
    /*                  *               */
    /*  CkqFS = f   *  C                */
    /*           kq     kq              */
    /*                                  */
    /*   *                              */
    /*  C     = RT(C   ) -i* IT(C   )   */
    /*   kq         kq           kq     */
    /*                                  */
    /*                        k         */
    /*                      <r >        */
    /* Rk_k+1 = theta(J) * -------      */
    /*               k        k+1       */
    /*                       R          */
    /*                                  */
 
 
    RT( V20(iter) ) = 0.0;  IT( V20(iter) ) = 0.0;
    RT( V21(iter) ) = 0.0;  IT( V21(iter) ) = 0.0;
    RT( V22(iter) ) = 0.0;  IT( V22(iter) ) = 0.0;
 
    RT( V40(iter) ) = 0.0;  IT( V40(iter) ) = 0.0;
    RT( V41(iter) ) = 0.0;  IT( V41(iter) ) = 0.0;
    RT( V42(iter) ) = 0.0;  IT( V42(iter) ) = 0.0;
    RT( V43(iter) ) = 0.0;  IT( V43(iter) ) = 0.0;
    RT( V44(iter) ) = 0.0;  IT( V44(iter) ) = 0.0;
 
    RT( V60(iter) ) = 0.0;  IT( V60(iter) ) = 0.0;
    RT( V61(iter) ) = 0.0;  IT( V61(iter) ) = 0.0;
    RT( V62(iter) ) = 0.0;  IT( V62(iter) ) = 0.0;
    RT( V63(iter) ) = 0.0;  IT( V63(iter) ) = 0.0;
    RT( V64(iter) ) = 0.0;  IT( V64(iter) ) = 0.0;
    RT( V65(iter) ) = 0.0;  IT( V65(iter) ) = 0.0;
    RT( V66(iter) ) = 0.0;  IT( V66(iter) ) = 0.0;
 
 
    for( i=1 ; i<= ANZ_NN(iter) ; ++i ){
 
         RT(V20(iter)) += Q_R(iter,i) * R2_3(iter,i) * RT(C20FS(iter,i));
 
         RT(V40(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C40FS(iter,i));
 
         RT(V60(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C60FS(iter,i));
         RT(V66(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C66FS(iter,i));
    }
    return( auswahlregel(iter,sym) );
}
/*------------------------------------------------------------------------------
                                 Vkq8()
------------------------------------------------------------------------------*/
ITERATION *Vkq8( ITERATION *iter,INT sym)
{
    INT i;
    ITERATION *auswahlregel(ITERATION *,INT);
 
    /*                                  */
    /*  Q_R   = Q(R)                    */
    /*                  *               */
    /*  CkqFS = f   *  C                */
    /*           kq     kq              */
    /*                                  */
    /*   *                              */
    /*  C     = RT(C   ) -i* IT(C   )   */
    /*   kq         kq           kq     */
    /*                                  */
    /*                        k         */
    /*                      <r >        */
    /* Rk_k+1 = theta(J) * -------      */
    /*               k        k+1       */
    /*                       R          */
    /*                                  */
 
 
 
    RT( V20(iter) ) = 0.0;  IT( V20(iter) ) = 0.0;
    RT( V21(iter) ) = 0.0;  IT( V21(iter) ) = 0.0;
    RT( V22(iter) ) = 0.0;  IT( V22(iter) ) = 0.0;
 
    RT( V40(iter) ) = 0.0;  IT( V40(iter) ) = 0.0;
    RT( V41(iter) ) = 0.0;  IT( V41(iter) ) = 0.0;
    RT( V42(iter) ) = 0.0;  IT( V42(iter) ) = 0.0;
    RT( V43(iter) ) = 0.0;  IT( V43(iter) ) = 0.0;
    RT( V44(iter) ) = 0.0;  IT( V44(iter) ) = 0.0;
 
    RT( V60(iter) ) = 0.0;  IT( V60(iter) ) = 0.0;
    RT( V61(iter) ) = 0.0;  IT( V61(iter) ) = 0.0;
    RT( V62(iter) ) = 0.0;  IT( V62(iter) ) = 0.0;
    RT( V63(iter) ) = 0.0;  IT( V63(iter) ) = 0.0;
    RT( V64(iter) ) = 0.0;  IT( V64(iter) ) = 0.0;
    RT( V65(iter) ) = 0.0;  IT( V65(iter) ) = 0.0;
    RT( V66(iter) ) = 0.0;  IT( V66(iter) ) = 0.0;
 
 
    for( i=1 ; i<= ANZ_NN(iter) ; ++i ){
 
         RT(V40(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C40FS(iter,i));
         RT(V44(iter)) += Q_R(iter,i) * R4_5(iter,i) * RT(C44FS(iter,i));
 
         RT(V60(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C60FS(iter,i));
         RT(V64(iter)) += Q_R(iter,i) * R6_7(iter,i) * RT(C64FS(iter,i));
    }
 
    if( ABSD( 2*RT(V44(iter))- 5*RT(V40(iter)) ) >EPS1 )
       read_error(30,(FILE*)0,(CHAR*)0);
    if( ABSD( 2*RT(V64(iter))+21*RT(V60(iter)) ) >EPS1 )
       read_error(31,(FILE*)0,(CHAR*)0);
 
    return( auswahlregel(iter,sym) );
}
/*------------------------------------------------------------------------------
                             init_iteration()
------------------------------------------------------------------------------*/
KRISTALLFELD *init_iterationnew(CHAR *iontype) /* [1] */
{
 
    ITERATION    *iteration;
    ITERATION    *read_Bkqnew(CHAR*);
    KRISTALLFELD *kristallfeld;
    iteration  = read_Bkqnew(iontype);
 
   kristallfeld              = KRISTALLFELD(1);
   ITERATION(kristallfeld)   = iteration;
   SYMMETRIENR(kristallfeld) = SYMMETRIENR( iteration );
   IONENNUMMER(kristallfeld) = IONENNUMMER( iteration );
   return( kristallfeld );
 
/* [1] : alles noetige fuer eine moegliche Iteration vorbereiten   */
/* [2] : Naechste Nachbarn und Magnetfeld von file filename lesen  */
/* [3] : Umgebungsionen fuer eine eventuelle Iteration vorbereiten */
 
}

KRISTALLFELD *init_iteration(CHAR *filename,INT symmetrienr,CHAR modus) /* [1] */
{
 
    UMGEBUNG     *umgebung ,   *read_nn(CHAR *name);
    ITERATION    *iteration,   *init_umgebung(CHAR*,UMGEBUNG*);
    ITERATION    *read_Vkq(CHAR*,INT),  *read_Bkq(CHAR*,INT) ,*read_Akq(CHAR*,INT);
    ITERATION    *read_Wkq(CHAR*,INT),  *read_xW(CHAR*,INT);
    ITERATION    *read_Dkq(CHAR*,INT),  *read_Lkq(CHAR*,INT);
    KRISTALLFELD *kristallfeld;
 
 
    switch(modus){
 
/* [2] */ case SIN:   umgebung  = read_nn( filename );
/* [3] */             iteration = init_umgebung( filename,umgebung );
                      SYMMETRIENR(iteration) = symmetrienr;
                      break;
 
          case XW :   iteration  = read_xW( filename ,symmetrienr);
                      break;
 
          case AKQ:   iteration  = read_Akq( filename ,symmetrienr);
                      break;
 
          case BKQ:   iteration  = read_Bkq( filename ,symmetrienr);
                      break;
 
          case DKQ:   iteration  = read_Dkq( filename ,symmetrienr);
                      break;
 
          case LKQ:   iteration  = read_Lkq( filename ,symmetrienr);
                      break;
 
          case VKQ:   iteration  = read_Vkq( filename ,symmetrienr);
                      break;
 
          case WKQ:   iteration  = read_Wkq( filename ,symmetrienr);
                      break;

          default:    fprintf(stderr,"init_iteration: Invalid parameter type %c\n",modus); exit(-1);
   }
 
   kristallfeld              = KRISTALLFELD(1);
   INFILE(kristallfeld)      = filename;
   ITERATION(kristallfeld)   = iteration;
   SYMMETRIENR(kristallfeld) = SYMMETRIENR( iteration );
   IONENNUMMER(kristallfeld) = IONENNUMMER( iteration );
   return( kristallfeld );
 
/* [1] : alles noetige fuer eine moegliche Iteration vorbereiten   */
/* [2] : Naechste Nachbarn und Magnetfeld von file filename lesen  */
/* [3] : Umgebungsionen fuer eine eventuelle Iteration vorbereiten */
 
}

/*------------------------------------------------------------------------------
                          auswahlregel()
------------------------------------------------------------------------------*/
ITERATION *auswahlregel(ITERATION *iter,INT symmetrienr)
{
    INT zwei_j;
    zwei_j = DIMJ(iter) - 1;
 
    switch((int)( zwei_j >2) ){
       case JA : switch( symmetrienr ){
                    case 8 : RT( V20(iter) ) = 0.0;
 
                    case 7 :
                    case 6 :
                    case 5 :
                    case 4 :
                    case 3 : RT( V22(iter) ) = 0.0;
                             IT( V22(iter) ) = 0.0;
 
                    case 2 :
                    case 1 : RT( V21(iter) ) = 0.0;
                             IT( V21(iter) ) = 0.0;
                             if(RT(V22(iter))!=0.0&&IT(V22(iter))!=0.0)
                                IT(V22(iter)) = 0.0;
 
                    case 0 : IT( V20(iter) ) = 0.0;
 /*                            if(RT(V21(iter))!=0.0&&IT(V21(iter))!=0.0)  // removed because it gave strange
                                IT(V21(iter)) = 0.0;                     // results in rotated systems
 */                }
                 break;
 
       default : RT( V20(iter) ) = 0.0; IT( V20(iter) ) = 0.0;
                 RT( V21(iter) ) = 0.0; IT( V21(iter) ) = 0.0;
                 RT( V22(iter) ) = 0.0; IT( V22(iter) ) = 0.0;
    }
 
 
    switch((int)( zwei_j >4) ){
       case JA : switch( symmetrienr ){
 
                    case 7 :
                    case 6 :
 
                    case 5 : RT( V44(iter) ) = 0.0;
                    case 8 :
                             IT( V44(iter) ) = 0.0;
 
                    case 4 :
                    case 3 :
                             RT( V42(iter) ) = 0.0;
                             IT( V42(iter) ) = 0.0;
 
                    case 2 :
                             if(RT(V44(iter))!=0.0&&IT(V44(iter))!=0.0)
                                IT(V44(iter)) = 0.0;
                             if(RT(V42(iter))!=0.0&&IT(V42(iter))!=0.0)
                                IT(V42(iter)) = 0.0;
 
                    case 1 :
                             if( symmetrienr==5||symmetrienr==6){
                               if(RT(V43(iter))!=0.0&&IT(V43(iter))!=0.0)
                                  IT(V43(iter)) = 0.0;
                             }
                             else{RT( V43(iter) ) = 0.0;
                                  IT( V43(iter) ) = 0.0;}
                             RT( V41(iter) ) = 0.0;
                             IT( V41(iter) ) = 0.0;
 
                    case 0 : IT( V40(iter) ) = 0.0;
                 }
                 break;
 
       default : RT( V40(iter) ) = 0.0; IT( V40(iter) ) = 0.0;
                 RT( V41(iter) ) = 0.0; IT( V41(iter) ) = 0.0;
                 RT( V42(iter) ) = 0.0; IT( V42(iter) ) = 0.0;
                 RT( V43(iter) ) = 0.0; IT( V43(iter) ) = 0.0;
                 RT( V44(iter) ) = 0.0; IT( V44(iter) ) = 0.0;
    }
 
 
    switch( (int)(zwei_j >6) ){
       case JA : switch( symmetrienr ){
 
                    case 7 :
                    case 6 :
                    case 5 :
                             RT( V64(iter) ) = 0.0;
                    case 8 :
                             IT( V64(iter) ) = 0.0;
                    case 4 :
                    case 3 :
                             RT( V62(iter) ) = 0.0;
                             IT( V62(iter) ) = 0.0;
                             if( symmetrienr!=5&&symmetrienr!=6&&symmetrienr!=7){
                                RT( V66(iter) ) = 0.0;
                                IT( V66(iter) ) = 0.0;
                             }else{
                              if( symmetrienr==5||symmetrienr==6)
                               if(RT(V66(iter))!=0.0&&IT(V66(iter))!=0.0)
                                  IT(V66(iter)) = 0.0;
                             }
                    case 2 :
                            if( symmetrienr==4||symmetrienr==2)
                              if(RT(V64(iter))!=0.0&&IT(V64(iter))!=0.0)
                                 IT(V64(iter)) = 0.0;
                            if(RT(V66(iter))!=0.0&&IT(V66(iter))!=0.0)
                               IT(V66(iter)) = 0.0;
                            if(RT(V62(iter))!=0.0&&IT(V62(iter))!=0.0)
                               IT(V62(iter)) = 0.0;
 
                    case 1 :
                             RT( V65(iter) ) = 0.0;
                             IT( V65(iter) ) = 0.0;
                             RT( V61(iter) ) = 0.0;
                             IT( V61(iter) ) = 0.0;
 
                            if( symmetrienr!=5&&symmetrienr!=6){
                               RT( V63(iter) ) = 0.0;
                               IT( V63(iter) ) = 0.0;
                            }else{
                              if( symmetrienr==6)
                               if(RT(V63(iter))!=0.0&&IT(V63(iter))!=0.0)
                                  IT(V63(iter)) = 0.0;
                             }
                    case 0 : IT( V60(iter) ) = 0.0;
                 }
                 break;
 
       default : RT( V60(iter) ) = 0.0; IT( V60(iter) ) = 0.0;
                 RT( V61(iter) ) = 0.0; IT( V61(iter) ) = 0.0;
                 RT( V62(iter) ) = 0.0; IT( V62(iter) ) = 0.0;
                 RT( V63(iter) ) = 0.0; IT( V63(iter) ) = 0.0;
                 RT( V64(iter) ) = 0.0; IT( V64(iter) ) = 0.0;
                 RT( V65(iter) ) = 0.0; IT( V65(iter) ) = 0.0;
                 RT( V66(iter) ) = 0.0; IT( V66(iter) ) = 0.0;
    }
 
 
    return( iter );
}
/*------------------------------------------------------------------------------
                          init_umgebung(CHAR*,UMGEBUNG*)
------------------------------------------------------------------------------*/
ITERATION *init_umgebung(  CHAR     *name,UMGEBUNG *umgebung )  /* Umgebungsionen fuer Iteration    */
                  /* vorbereiten                      */
                /* Stevensoperatoren initialisieren */
{
 
    ITERATION *iteration,*iter_alloc(INT,INT);
    UMGEBUNG  *trans_R_sphaerisch(CHAR*,UMGEBUNG*);
    UMGEBUNG  *trans_B_kartesisch(CHAR*,UMGEBUNG*);
    UMGEBUNG  *normiere_laengen(UMGEBUNG*);
    STEVENS   *calc_Pkq(INT);
    MATRIX    *calc_Bmag(INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);
    KOMPLEX   *Clm(INT,INT,DOUBLE,DOUBLE);
    DOUBLE    pow_(DOUBLE,INT),R,theta_R,phi_R,Bx,By,Bz,gj,myB;
    INT       anz_nn,i,_e4f,dimj;
 
 
 
    dimj      = IONENIMP[ IONENNR(umgebung) ].dimj;
    anz_nn    = ANZ_NN(umgebung);
    iteration = iter_alloc(dimj,anz_nn);
    ANZ_NN(iteration)     =  anz_nn;
    DIMJ(iteration)       =  dimj;
    GJ(iteration)         =  gj   = IONENIMP[ IONENNR(umgebung) ].gj;
    IONNAME(iteration)    =  IONENIMP[ IONENNR(umgebung) ].ionname;
    TEMPERATUR(iteration) =  TEMPERATUR(umgebung);
    PKQ(iteration)        =  calc_Pkq( dimj );
 
    umgebung =  trans_B_kartesisch( name,umgebung ); /* B = (Bx,By,Bz)   */
    Bx               = B1(umgebung);
    By               = B2(umgebung);
    Bz               = B3(umgebung);
 
    B1(iteration)    = B1(umgebung);
    B2(iteration)    = B2(umgebung);
    B3(iteration)    = B3(umgebung);
 
 
    myB              = myBmeV;
    HMAG(iteration)  = calc_Bmag( dimj,gj,myB,Bx,By,Bz );
 
    umgebung = trans_R_sphaerisch( name,umgebung );    /*  R = (R,theta,phi) */
    umgebung = normiere_laengen( umgebung );           /* |R|= |R|/a0        */
 
 
    /*                                  */
    /*  Q_R   = Q(R)                    */
    /*                        *         */
    /*  CkqFS = epsilon   * C           */
    /*                 kq    kq         */
    /*                                  */
    /*   *                              */
    /*  C     = RT(C   ) -i* IT(C   )   */
    /*   kq         kq           kq     */
    /*                                  */
    /*                           k      */
    /*                2S+1     <r >     */
    /* Rk_k+1 = theta(    L )*-------   */
    /*               k     J    k+1     */
    /*                         R        */
    /*                                  */
 
    _e4f     = e4f(umgebung);         /*  Elektronen des Aufions in 4f-Schale */
 
    for( i=1 ; i<=anz_nn ; ++i ){
         R       = X1(umgebung,i);
         theta_R = X3(umgebung,i);
         phi_R   = X2(umgebung,i);
 
         Q_R( iteration,i) = Q(umgebung,i);
         R2_3(iteration,i) = alpha_J[_e4f] * r2(IONENNR(umgebung)) / pow_( R , 3);
         R4_5(iteration,i) = beta_J[ _e4f] * r4(IONENNR(umgebung)) / pow_( R , 5);
         R6_7(iteration,i) = gamma_J[_e4f] * r6(IONENNR(umgebung)) / pow_( R , 7);
 
         R2_3(iteration,i) *= _E0programm*1000; /* Energieeinheit */
         R4_5(iteration,i) *= _E0programm*1000; /* meV            */
         R6_7(iteration,i) *= _E0programm*1000;
 
         EINHEITNROUT( iteration ) =  2; /* siehe EINHEITIMP[] */
         EINHEITNRIN(  iteration ) =  2; /* siehe EINHEITIMP[] */
 
            C20FS(iteration,i)   = Clm( 2,0,theta_R,phi_R );
         RT(C20FS(iteration,i)) *=      epn2n(0);
         IT(C20FS(iteration,i)) *=     -epn2n(0);
 
            C21FS(iteration,i)   = Clm( 2,1,theta_R,phi_R );
         RT(C21FS(iteration,i)) *=      epn1n(1);
         IT(C21FS(iteration,i)) *=     -epn1n(1);
 
            C22FS(iteration,i)   = Clm( 2,2,theta_R,phi_R );
         RT(C22FS(iteration,i)) *=      epn0n(2);
         IT(C22FS(iteration,i)) *=     -epn0n(2);
 
 
 
            C40FS(iteration,i)   = Clm( 4,0,theta_R,phi_R );
         RT(C40FS(iteration,i)) *=      epn4n(0);
         IT(C40FS(iteration,i)) *=     -epn4n(0);
 
            C41FS(iteration,i)   = Clm( 4,1,theta_R,phi_R );
         RT(C41FS(iteration,i)) *=      epn3n(1);
         IT(C41FS(iteration,i)) *=     -epn3n(1);
 
            C42FS(iteration,i)   = Clm( 4,2,theta_R,phi_R );
         RT(C42FS(iteration,i)) *=      epn2n(2);
         IT(C42FS(iteration,i)) *=     -epn2n(2);
 
            C43FS(iteration,i)   = Clm( 4,3,theta_R,phi_R );
         RT(C43FS(iteration,i)) *=      epn1n(3);
         IT(C43FS(iteration,i)) *=     -epn1n(3);
 
            C44FS(iteration,i)   = Clm( 4,3,theta_R,phi_R );
         RT(C44FS(iteration,i)) *=      epn0n(4);
         IT(C44FS(iteration,i)) *=     -epn0n(4);
 
 
 
            C60FS(iteration,i)   = Clm( 6,0,theta_R,phi_R );
         RT(C60FS(iteration,i)) *=      epn6n(0);
         IT(C60FS(iteration,i)) *=     -epn6n(0);
 
            C61FS(iteration,i)   = Clm( 6,1,theta_R,phi_R );
         RT(C61FS(iteration,i)) *=      epn5n(1);
         IT(C61FS(iteration,i)) *=     -epn5n(1);
 
            C62FS(iteration,i)   = Clm( 6,2,theta_R,phi_R );
         RT(C62FS(iteration,i)) *=      epn4n(2);
         IT(C62FS(iteration,i)) *=     -epn4n(2);
 
            C63FS(iteration,i)   = Clm( 6,3,theta_R,phi_R );
         RT(C63FS(iteration,i)) *=      epn3n(3);
         IT(C63FS(iteration,i)) *=     -epn3n(3);
 
            C64FS(iteration,i)   = Clm( 6,4,theta_R,phi_R );
         RT(C64FS(iteration,i)) *=      epn2n(4);
         IT(C64FS(iteration,i)) *=     -epn2n(4);
 
            C65FS(iteration,i)   = Clm( 6,5,theta_R,phi_R );
         RT(C65FS(iteration,i)) *=      epn1n(5);
         IT(C65FS(iteration,i)) *=     -epn1n(5);
 
            C66FS(iteration,i)   = Clm( 6,6,theta_R,phi_R );
         RT(C66FS(iteration,i)) *=      epn0n(6);
         IT(C66FS(iteration,i)) *=     -epn0n(6);
    }
 
    return( iteration );
}
/*------------------------------------------------------------------------------
                                iter_alloc(INT,INT)
------------------------------------------------------------------------------*/
ITERATION *iter_alloc(INT dimj,INT anz_nn) /* Speicher fuer Struktur  ITERATION*/
                  /* holen                            */
{
    ITERATION  *iteration_i;
    MATRIX     *mx_alloc(INT,INT);
    iteration_i = ITERATION_ALLOC(1);
 
 
    Q_R_P(  iteration_i ) = DOUBLE_ALLOC(anz_nn); /* Ladung der    */
                                                  /* Umgebungsionen*/
 
 
    R2_3_P( iteration_i ) = DOUBLE_ALLOC(anz_nn); /*  k                  k+1*/
    R4_5_P( iteration_i ) = DOUBLE_ALLOC(anz_nn); /*<r > / | R |            */
    R6_7_P( iteration_i ) = DOUBLE_ALLOC(anz_nn); /*                        */
 
    V20(    iteration_i ) = KX_ALLOC(1);
    V21(    iteration_i ) = KX_ALLOC(1);
    V22(    iteration_i ) = KX_ALLOC(1);
    V40(    iteration_i ) = KX_ALLOC(1);
    V41(    iteration_i ) = KX_ALLOC(1);
    V42(    iteration_i ) = KX_ALLOC(1);
    V43(    iteration_i ) = KX_ALLOC(1);
    V44(    iteration_i ) = KX_ALLOC(1);
    V60(    iteration_i ) = KX_ALLOC(1);
    V61(    iteration_i ) = KX_ALLOC(1);
    V62(    iteration_i ) = KX_ALLOC(1);
    V63(    iteration_i ) = KX_ALLOC(1);
    V64(    iteration_i ) = KX_ALLOC(1);
    V65(    iteration_i ) = KX_ALLOC(1);
    V66(    iteration_i ) = KX_ALLOC(1);
 
    RT( V20( iteration_i ) ) = 0.0;
    RT( V21( iteration_i ) ) = 0.0;
    RT( V22( iteration_i ) ) = 0.0;
    RT( V40( iteration_i ) ) = 0.0;
    RT( V41( iteration_i ) ) = 0.0;
    RT( V42( iteration_i ) ) = 0.0;
    RT( V43( iteration_i ) ) = 0.0;
    RT( V44( iteration_i ) ) = 0.0;
    RT( V60( iteration_i ) ) = 0.0;
    RT( V61( iteration_i ) ) = 0.0;
    RT( V62( iteration_i ) ) = 0.0;
    RT( V63( iteration_i ) ) = 0.0;
    RT( V64( iteration_i ) ) = 0.0;
    RT( V65( iteration_i ) ) = 0.0;
    RT( V66( iteration_i ) ) = 0.0;
 
    IT( V20( iteration_i ) ) = 0.0;
    IT( V21( iteration_i ) ) = 0.0;
    IT( V22( iteration_i ) ) = 0.0;
    IT( V40( iteration_i ) ) = 0.0;
    IT( V41( iteration_i ) ) = 0.0;
    IT( V42( iteration_i ) ) = 0.0;
    IT( V43( iteration_i ) ) = 0.0;
    IT( V44( iteration_i ) ) = 0.0;
    IT( V60( iteration_i ) ) = 0.0;
    IT( V61( iteration_i ) ) = 0.0;
    IT( V62( iteration_i ) ) = 0.0;
    IT( V63( iteration_i ) ) = 0.0;
    IT( V64( iteration_i ) ) = 0.0;
    IT( V65( iteration_i ) ) = 0.0;
    IT( V66( iteration_i ) ) = 0.0;
 
 
      C20_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1); /*  C    */
      C21_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1); /*   kq  */
      C22_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C40_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C41_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C42_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C43_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C44_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C60_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C61_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C62_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C63_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C64_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C65_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
      C66_P(  iteration_i ) = KX_P_ALLOC(anz_nn+1);
 
      VOR20(  iteration_i) = 1.0;
      VOR21(  iteration_i) = 1.0;
      VOR22(  iteration_i) = 1.0;
      VOR40(  iteration_i) = 1.0;
      VOR41(  iteration_i) = 1.0;
      VOR42(  iteration_i) = 1.0;
      VOR43(  iteration_i) = 1.0;
      VOR44(  iteration_i) = 1.0;
      VOR60(  iteration_i) = 1.0;
      VOR61(  iteration_i) = 1.0;
      VOR62(  iteration_i) = 1.0;
      VOR63(  iteration_i) = 1.0;
      VOR64(  iteration_i) = 1.0;
      VOR65(  iteration_i) = 1.0;
      VOR66(  iteration_i) = 1.0;
 
    HAMILTONIAN( iteration_i ) = mx_alloc( dimj,dimj );
 
    return( iteration_i );
}
/*------------------------------------------------------------------------------
                                    calc_Bmag()
------------------------------------------------------------------------------*/
MATRIX *calc_Bmag(INT dimj,DOUBLE gj,DOUBLE myB,DOUBLE Bx,DOUBLE By,DOUBLE Bz ) /*magnetischen Hamiltonian */
                         /* Hmag ausrechnen         */
{
    INT    m,n;
    MATRIX *bmag,*mx_alloc(INT,INT);
    DOUBLE jm,jp;
 
    #include "define_j.c"          /* mj,J2,J+,... definieren */
    bmag = mx_alloc( dimj,dimj );  /* Speicher fuer (J nj| Hmag |mj J)*/
 
    for( n=dimj ; n>=1 ; --n)
         for( m=dimj ; m>=1 ; --m){
              jm=JM(mj)*D(nj,mj-1);
              jp=JP(mj)*D(nj,mj+1);
/*       sign changed 24.9.08 because  Zeeman term has negative sign */
              R(bmag,n,m) = -gj*myB*(  0.5*Bx*( jm+jp ) + mj*Bz*D(nj,mj)  );
              I(bmag,n,m) = -gj*myB*   0.5*By*( jm-jp );
    }
    return( bmag );
}
/*------------------------------------------------------------------------------
                                    calc_Bmag_D()
------------------------------------------------------------------------------*/
MATRIX *calc_Bmag_D(INT dimj,DOUBLE gj,DOUBLE myB,DOUBLE Bx,DOUBLE By,DOUBLE Bz,DOUBLE Dx2,DOUBLE Dy2,DOUBLE Dz2,DOUBLE Dx4,DOUBLE Dy4,DOUBLE Dz4 ) /*magnetischen Hamiltonian */
                                   /* Hmag = - gJ muB J.B + simple anisotropy  H= + Dx2 Jx ^ 2+ Dy2 Jy ^ 2+ Dz2 Jz ^ 2 + Dx4 Jx ^ 4+ Dy4 Jy ^ 4+ Dz4 Jz ^ 4      */
{
    INT    m,n;
    MATRIX *bmag,*mx_alloc(INT,INT);
    DOUBLE jm,jp,jx2,jy2,jx4,jy4;

    #include "define_j.c"          /* mj,J2,J+,... definieren */
    bmag = mx_alloc( dimj,dimj );  /* Speicher fuer (J nj| Hmag |mj J)*/

 for( n=dimj ; n>=1 ; --n)
         for( m=dimj ; m>=1 ; --m){
              jm=JM(mj)*D(nj,mj-1);
              jp=JP(mj)*D(nj,mj+1);
              jx2=0.25*(JPMMP(mj)*D(nj,mj));
              jy2=jx2;
              jx4=jx2*jx2+0.0625*(JP2(mj-2)*JM2(mj)+JM2(mj+2)*JP2(mj));
              jy4=jx4;
              if (D(nj,mj-2)>0.5) {jx2+=0.25*JM2(mj);jy2-=0.25*JM2(mj);jx4+=0.0625*JM2(mj)*JPMMP(mj);jy4-=0.0625*JM2(mj)*JPMMP(mj);}
              if (D(nj,mj+2)>0.5) {jx2+=0.25*JP2(mj);jy2-=0.25*JP2(mj);jx4+=0.0625*JP2(mj)*JPMMP(mj);jy4-=0.0625*JP2(mj)*JPMMP(mj);}
              if (D(nj,mj-4)>0.5) {jx4+=0.0625*JM2(mj-2)*JM2(mj);jy4+=0.0625*JM2(mj-2)*JM2(mj);}
              if (D(nj,mj+4)>0.5) {jx4+=0.0625*JP2(mj+2)*JP2(mj);jy4+=0.0625*JP2(mj+2)*JP2(mj);}
/* sign changed 24.9.08 because zeeman term has negative sign */
              R(bmag,n,m) = -gj*myB*(  0.5*Bx*( jm+jp ) + mj*Bz*D(nj,mj)  )+Dz2*mj*mj*D(nj,mj)+Dx2*jx2+Dy2*jy2+Dz4*mj*mj*mj*mj*D(nj,mj)+Dx4*jx4+Dy4*jy4;
              I(bmag,n,m) = -gj*myB*   0.5*By*( jm-jp );
    }


    return( bmag );
}


/*------------------------------------------------------------------------------
                                    calcBmol()
------------------------------------------------------------------------------*/
MATRIX *calcBmol(INT dimj,MATRIX *bmag, DOUBLE gjs, DOUBLE myB, DOUBLE Bx, DOUBLE By, DOUBLE Bz )   /* gjs = 2(gJ-1) */
{
    INT    m,n;
    DOUBLE jm,jp;
 
    #include "define_j.c"          /* mj,J2,J+,... definieren */
 
    for( n=dimj ; n>=1 ; --n)
         for( m=dimj ; m>=1 ; --m){
              jm=JM(mj)*D(nj,mj-1);
              jp=JP(mj)*D(nj,mj+1);
 /*       sign changed 24.9.08 because  Zeeman term has negative sign */
             R(bmag,n,m) += -gjs*myB*(  0.5*Bx*( jm+jp ) + mj*Bz*D(nj,mj)  );
              I(bmag,n,m) += -gjs*myB*   0.5*By*( jm-jp );
    }
    return( bmag );
}
/*------------------------------------------------------------------------------
                                    calc_iBmag()
------------------------------------------------------------------------------*/
MATRIX *calc_iBmag(MATRIX * bmag, DOUBLE gj, DOUBLE myB, DOUBLE Bx, DOUBLE By, DOUBLE Bz, DOUBLE Bxmol, DOUBLE Bymol, DOUBLE Bzmol, DOUBLE Dx2, DOUBLE Dy2, DOUBLE Dz2, DOUBLE Dx4, DOUBLE Dy4, DOUBLE Dz4 )
{
    INT    m,n,dimj;
    DOUBLE jm,jp,jx2,jy2,jx4,jy4;
 
    #include "define_j.c"          /* mj,J2,J+,... definieren */
                                   /* <nj| A |mj>             */
    dimj = MXDIM(bmag);
/*    printf("dim=%i\n",dimj);*/
    for( n=dimj ; n>=1 ; --n)
         for( m=dimj ; m>=1 ; --m){
              jm=JM(mj)*D(nj,mj-1);
              jp=JP(mj)*D(nj,mj+1);
 /* sign changed 24.9.08 because zeeman term has negative sign */
             R(bmag,n,m) = -2.0*(gj-1.0)*myB*(0.5*Bxmol*( jm+jp )
                            + mj*Bzmol*D(nj,mj)  );
              I(bmag,n,m) = -2.0*(gj-1.0)*myB* 0.5*Bymol*( jm-jp );
         }
for( n=dimj ; n>=1 ; --n)
         for( m=dimj ; m>=1 ; --m){
              jm=JM(mj)*D(nj,mj-1);
              jp=JP(mj)*D(nj,mj+1);
              jx2=0.25*(JPMMP(mj)*D(nj,mj));
              jy2=jx2;
              jx4=jx2*jx2+0.0625*(JP2(mj-2)*JM2(mj)+JM2(mj+2)*JP2(mj));
              jy4=jx4;
              if (D(nj,mj-2)>0.5) {jx2+=0.25*JM2(mj);jy2-=0.25*JM2(mj);jx4+=0.0625*JM2(mj)*JPMMP(mj);jy4-=0.0625*JM2(mj)*JPMMP(mj);}
              if (D(nj,mj+2)>0.5) {jx2+=0.25*JP2(mj);jy2-=0.25*JP2(mj);jx4+=0.0625*JP2(mj)*JPMMP(mj);jy4-=0.0625*JP2(mj)*JPMMP(mj);}
              if (D(nj,mj-4)>0.5) {jx4+=0.0625*JM2(mj-2)*JM2(mj);jy4+=0.0625*JM2(mj-2)*JM2(mj);}
              if (D(nj,mj+4)>0.5) {jx4+=0.0625*JP2(mj+2)*JP2(mj);jy4+=0.0625*JP2(mj+2)*JP2(mj);}
/* sign changed 24.9.08 because zeeman term has negative sign */
              R(bmag,n,m) += -gj*myB*(  0.5*Bx*( jm+jp ) + mj*Bz*D(nj,mj)  )+Dz2*mj*mj*D(nj,mj)+Dx2*jx2+Dy2*jy2+Dz4*mj*mj*mj*mj*D(nj,mj)+Dx4*jx4+Dy4*jy4;
              I(bmag,n,m) += -gj*myB*   0.5*By*( jm-jp );
    }

    return( bmag );
}
/*------------------------------------------------------------------------------
                            normiere_laengen()
------------------------------------------------------------------------------*/
UMGEBUNG *normiere_laengen(UMGEBUNG *umgebung) /* Laengen auf Bohrschen Radius normieren*/
{
    INT i;
 
    for( i=1 ; i<=ANZ_NN(umgebung) ;++i )
         X1( umgebung,i ) /= A0_BOHR;
 
    return( umgebung );
}
/*------------------------------------------------------------------------------
                            trans_R_sphaerisch()
------------------------------------------------------------------------------*/
UMGEBUNG*trans_R_sphaerisch(CHAR     *name, UMGEBUNG *umgebung)/*Ortskoordinaten der Umgebungsionen*/
                     /* ins sphaerische Koordinatensystem */
                /* transformieren                    */
{
    INT anz_nn,i;
    SPHAERE *sphaere,*sphaerisch(CHAR*,INT,DOUBLE,DOUBLE,DOUBLE);
    DOUBLE  x,y,z;
 
    if( MODUS(umgebung)!= 'r' )   return( umgebung );
 
    anz_nn = ANZ_NN(umgebung);
    for( i=1 ; i<=anz_nn ; ++i ){
         x = X1(umgebung,i);
         y = X2(umgebung,i);
         z = X3(umgebung,i);
 
         sphaere=sphaerisch(name,i,x,y,z);/* ins sphaer. Koor.transformieren*/
 
         X1(umgebung,i) = SPH_R(    sphaere);
         X2(umgebung,i) = SPH_PHI(  sphaere);
         X3(umgebung,i) = SPH_THETA(sphaere);
 
         free_( sphaere );
    }
    return( umgebung );
}
/*------------------------------------------------------------------------------
                            trans_B_kartesisch()
------------------------------------------------------------------------------*/
UMGEBUNG*trans_B_kartesisch( CHAR     *name, UMGEBUNG *umgebung)/*Magnetfeld ins kartesische       */
                    /* Koordinatensystem transformieren */
{
    UNUSED_PARAMETER(name);

    DOUBLE  h,theta,phi;
   // DOUBLE  sin(),cos();
 
    if( MODUS(umgebung) == 'r' )   return( umgebung );
 
    h     = B1(umgebung);
    phi   = B2(umgebung)*pi/180.0;
    theta = B3(umgebung)*pi/180.0;
 
    if( MODUS(umgebung) == 'p' )  theta = pi/2;
 
    B1(umgebung) = h*sin(theta)*cos(phi);
    B2(umgebung) = h*sin(theta)*sin(phi);
    B3(umgebung) = h*cos(theta);
 
    return( umgebung );
}
/*------------------------------------------------------------------------------
                               sphaerisch()
------------------------------------------------------------------------------*/
SPHAERE *sphaerisch(CHAR   *name,INT ionnr, DOUBLE x, DOUBLE y, DOUBLE z)  /* ins sphaerische Koordinatensystem */
                      /* transformieren                    */
{
    SPHAERE *sphaere;
    DOUBLE  r,theta,phi;
    DOUBLE  sinphi,cosphi;
    //DOUBLE  sqrt(),acos(),atan();
 
 
    r     = sqrt( x*x + y*y + z*z );
 
    if( r==0.0 ){ coor_error( (FILE*)0,name,ionnr,4);
                  exit(1);
                }
 
    theta = acos( z/r )*180.0/pi;
    if( theta==0.0 || theta==180.0)  phi = 0.0;
    else{
         cosphi = x/r/sin(pi/180.0*theta);
         sinphi = y/r/sin(pi/180.0*theta);
 
         if( cosphi != 0.0 ){
              phi  = atan( ABSD(sinphi)/ABSD(cosphi) )*180.0/pi;
              if( sinphi >=0.0 && cosphi>0.0 )  ;
              else if( sinphi >=0.0 && cosphi<0.0 )  phi=180.0-phi;
              else if( sinphi <=0.0 && cosphi>0.0 )  phi=360.0-phi;
              else phi=180.0+phi;
         }
         else if( cosphi==0.0 && sinphi >= 0.0)  phi=90.0;  else  phi = 270.0;
    }
 
    sphaere            = SPHAERE_ALLOC(1);
    SPH_R(sphaere)     = r;
    SPH_PHI(sphaere)   = phi;
    SPH_THETA(sphaere) = theta;
 
    return( sphaere );
}
/*------------------------------------------------------------------------------
                               info_info()
------------------------------------------------------------------------------*/
void info_info()      /* Informationen ueber moegliche Infos ausgeben */
{
    INT i;
 
    clearscreen;
    printf("\nImplementierte Infos are :\n\n");
    for( i=0 ; i<  (int)ANZ_INFOS ; ++i )
         printf("%s -i[nfo] %s \n",PROGRAMMNAME,INFO(i) );
    printf("\n");
 
}
/*------------------------------------------------------------------------------
                               info_befehle()
------------------------------------------------------------------------------*/
void info_befehle()     /* informationen ueber moegliche Befehle  ausgeben */
{
 INT i;
 CHAR *t01,*t02,*t03;
 
    clearscreen;
 
t01=" ------------------------------------------------------------------ \n";
t02="|   %8s is a Crystal Field program.                         |\n";
t03=" ------------------------------------------------------------------ \n";
 
    printf("%s",t01);printf(t02,PROGRAMMNAME);printf("%s",t03);
    printf("\n");
 
 
    printf("\nImplemented commands are :\n\n");
    for( i=0 ; i<  (int)ANZ_BEFEHLE ; ++i )
         printf("%s : %s %s \n",BEFEHLKOM(i),PROGRAMMNAME,BEFEHL(i) );
    printf("\n");
 
}
/*------------------------------------------------------------------------------
                              info_magnetfeld()
------------------------------------------------------------------------------*/
void info_magnetfeld(INT dimj,DOUBLE Bx,DOUBLE By,DOUBLE Bz)  /* informiere ueber (Jn|Hmag|mJ)   */
                 /* Hmag in Einheiten von gJ*myB    */
{
    CHAR *name = "results/magnfeld.info";
 
    FILE   *fp,*fopen(const char*,const char*);
    MATRIX *bmag,*calc_Bmag(INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);
    DOUBLE gj,myB;
    INT   m,n;
    CHAR *t01,*t02,*t03,*t04,*t05,*t06;
    CHAR *t07,*t08,*t09,*t10,*t11,*t12;
    DOUBLE j;
    DOUBLE jm;
    DOUBLE jn;
 
    clearscreen;
    printf("Information given in the File %s .\n",name);
    fp = fopen_errchk(name,"w");
 
    gj   = 1.0;
    myB  = 1.0;
    bmag = calc_Bmag( dimj,gj,myB,Bx,By,Bz );
 
    t01 = "=========================================================\n";
    t02 = "| Table of (JM'|Hmag|MJ) in units of gJ*myB             |\n";
    t03 = "| d.h  Hmag =   (Jx,Jy,Jz)*(Bz,By,Bz)                   |\n";
    t04 = "---------------------------------------------------------\n";
    t05 = "| J  = %4.1f                                             |\n";
    t06 = "| Bx = %7.4f                                          |\n";
    t07 = "| By = %7.4f                                          |\n";
    t08 = "| Bz = %7.4f                                          |\n";
    t09 = "=========================================================\n";
    t10 = "| Realt (%3.1f %4.1f I Hmag I %4.1f %3.1f ) = %13.4f   |\n";
    t11 = "| Imagi (%3.1f %4.1f I Hmag I %4.1f %3.1f ) = %13.4f   |\n";
    t12 = "---------------------------------------------------------\n";
 
    j = ((DOUBLE)dimj-1)/2;
    fprintf(fp,"%s",t01);
    fprintf(fp,"%s",t02);
    fprintf(fp,"%s",t03);
    fprintf(fp,"%s",t04);
    fprintf(fp,t05,j);
    fprintf(fp,t06,Bx);
    fprintf(fp,t07,By);
    fprintf(fp,t08,Bz);
    fprintf(fp,"%s",t09);
    for( n=dimj ; n>=1 ; --n )
         for( m=dimj ; m>=1 ; --m ){
               jm = m-j-1;
               jn = n-j-1;
              fprintf(fp,t10,j,jn,jm,j,R(bmag,n,m));
              fprintf(fp,t11,j,jn,jm,j,I(bmag,n,m));
              fprintf(fp,"%s",t12);
         }
 
    fclose(fp);
}
/*------------------------------------------------------------------------------
                               info_symmetrien()
------------------------------------------------------------------------------*/
void info_symmetrien(CHAR*name)  /* Info ueber implementierte Symmetrien ausgeben */
{
 
    FILE *fp,*fopen(const char*,const char*);
 
    clearscreen;
    printf("Information given in the File %s.\n",name);
    fp = fopen_errchk(name,"w");
 
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"� Table of implemented Symmetries                         �\n");
    fprintf(fp,"===========================================================\n");
    fprintf(fp,"� +-Symmetrynr. � Point groups of the crystal field       �\n");
    fprintf(fp,"- V -------------+-----------------------------------------\n");
    fprintf(fp,"�   �            �  C    C                                �\n");
    fprintf(fp,"� 0 � triclinic    �   i    1                               �\n");
    fprintf(fp,"----+------------+-----------------------------------------\n");
    fprintf(fp,"�   �            �  C    C    C                           �\n");
    fprintf(fp,"� 1 � monoclinic   �   2    s    2h                         �\n");
    fprintf(fp,"----+------------+-----------------------------------------\n");
    fprintf(fp,"�   �            �  C    D    D                           �\n");
    fprintf(fp,"� 2 � rhombich  �   2v   2    2h                         �\n");
    fprintf(fp,"----+------------+-----------------------------------------\n");
    fprintf(fp,"�   �            �  C    S    C                           �\n");
    fprintf(fp,"� 3 �            �   4    4    4h                         �\n");
    fprintf(fp,"�---� tetragonal �-----------------------------------------\n");
    fprintf(fp,"�   �            �  D    C    D    D                      �\n");
    fprintf(fp,"� 4 �            �   4    4v   2d   4h                    �\n");
    fprintf(fp,"----+------------+-----------------------------------------\n");
    fprintf(fp,"�   �            �  C    S                                �\n");
    fprintf(fp,"� 5 �            �   3    6                               �\n");
    fprintf(fp,"----� trigonal   �-----------------------------------------\n");
    fprintf(fp,"�   �            �  D    C    D                           �\n");
    fprintf(fp,"� 6 �            �   3    3v   3d                         �\n");
    fprintf(fp,"----+------------+-----------------------------------------\n");
    fprintf(fp,"�   �            �  C    C    C    D    C    D    D       �\n");
    fprintf(fp,"� 7 � hexagonal  �   6    3h   6h   6    6v   3h   6h     �\n");
    fprintf(fp,"----+------------+-----------------------------------------\n");
    fprintf(fp,"�   �            �  T    T    T    O    O                 �\n");
    fprintf(fp,"� 8 � Cubic    �        d    h         h                �\n");
    fprintf(fp,"-----------------------------------------------------------\n");
    fprintf(fp,"\n");
 
    fclose(fp);
}
/*----------------------------------------------------------------------------
                           info_konstanten()
-----------------------------------------------------------------------------*/
void info_konstanten()   /* Liste der benutzten Naturkonstanten */
{
    CHAR *name = "results/konstntn.info";
    FILE *fp,*fopen(const char*,const char*);
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11,*t12,*t13,*t14,*t15;
    CHAR *s01,*s02,*s03,*s04,*s05,*s06,*s07;

 
    fp = fopen_errchk(name,"w");
    clearscreen;
    printf("Information given in the File %s.\n",name);
 
    t01 = "===========================================================\n"  ;
    t02 = "|Table of constants                                       |\n"  ;
    t03 = "===========================================================\n"  ;
    t04 = "| Planks constant             h  : %1.9f  e-34 Jsec |\n";
    t05 = "-----------------------------------------------------------\n";
    t06 = "| Elementary charge e              : %1.10f e-19 C    |\n";
    t07 = "-----------------------------------------------------------\n";
    t08 = "| Electron mass m              : %1.9f  e-31 Kg   |\n";
    t09 = "-----------------------------------------------------------\n";
    t10 = "| Velocity of light c         : %1.8f   e+08 m/sec|\n";
    t11 = "-----------------------------------------------------------\n";
    t12 = "| Bolzmann-constant k           : %1.5f      e-23 J/K  |\n";
    t13 = "-----------------------------------------------------------\n";
    t14 = "| Avogardo's-number Na          : %1.6f     e+23 1/mol|\n";
    t15 = "-----------------------------------------------------------\n";
 
    s01 = "=========================================================\n"  ;
    s02 = "|derived constants                                      |\n"  ;
    s03 = "=========================================================\n"  ;
    s04 = "| Bohr  Radius a0    :      %f Angstroem      |\n";
    s05 = "---------------------------------------------------------\n";
    s06 = "| Bohr Magneton myB :      %10.6e eV/Tesla   |\n";
    s07 = "---------------------------------------------------------\n";
 
    fprintf(fp,"%s",t01);
    fprintf(fp,"%s",t02);
    fprintf(fp,"%s",t03);
    fprintf(fp,t04,_h);
    fprintf(fp,"%s",t05);
    fprintf(fp,t06,_e);
    fprintf(fp,"%s",t07);
    fprintf(fp,t08,_m);
    fprintf(fp,"%s",t09);
    fprintf(fp,t10,_c);
    fprintf(fp,"%s",t11);
    fprintf(fp,t12,_k);
    fprintf(fp,"%s",t13);
    fprintf(fp,t14,_NA);
    fprintf(fp,"%s",t15);
 
    fprintf(fp,"%s",s01);
    fprintf(fp,"%s",s02);
    fprintf(fp,"%s",s03);
    fprintf(fp,s04,A0_BOHR);
    fprintf(fp,"%s",s05);
    fprintf(fp,s06,_myBplus);
    fprintf(fp,"%s",s07);
 
 
    fclose(fp);
}
/*----------------------------------------------------------------------------
                               info_Vlm()
-----------------------------------------------------------------------------*/
void info_Vlm(CHAR *filename,INT symmetrienr,CHAR *einheit)
{
    CHAR *name = "results/Vlm.info";
    FILE *fp,*fopen(const char*,const char*);
    KOMPLEX      *z;
    KRISTALLFELD *kf,*init_iteration(CHAR *,INT,CHAR);
/*  ITERATION    *iteration; */
    DOUBLE       renorm = E0_EINHEIT;
    INT          s;
    CHAR         *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    DOUBLE       rt,it;
 
 
    fp = fopen_errchk(name,"w");
    clearscreen;
    printf("Information given in the File %s.\n",name);
 
 
    if(VALUE(einheit,0)=='1') renorm *= 1e+26 / _h / _c ;
    kf = init_iteration( filename,symmetrienr,VKQ );
 
    switch( symmetrienr ){
         case 0 : /* iteration = */ Vkq0( ITERATION(kf),symmetrienr  );
                  break;
         case 1 : /* iteration = */ Vkq1( ITERATION(kf),symmetrienr  );
                  break;
         case 2 : /* iteration = */ Vkq2( ITERATION(kf),symmetrienr  );
                  break;
         case 3 : /* iteration = */ Vkq3( ITERATION(kf),symmetrienr  );
                  break;
         case 4 : /* iteration = */ Vkq4( ITERATION(kf),symmetrienr  );
                  break;
         case 5 : /* iteration = */ Vkq5( ITERATION(kf),symmetrienr  );
                  break;
         case 6 : /* iteration = */ Vkq6( ITERATION(kf),symmetrienr  );
                  break;
         case 7 : /* iteration = */ Vkq7( ITERATION(kf),symmetrienr  );
                  break;
         case 8 : /* iteration = */ Vkq8( ITERATION(kf),symmetrienr  );
                  break;
    }
 
    z=V20(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V20(ITERATION(kf))=z;
    z=V21(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V21(ITERATION(kf))=z;
    z=V22(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V22(ITERATION(kf))=z;
 
    z=V40(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V40(ITERATION(kf))=z;
    z=V41(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V41(ITERATION(kf))=z;
    z=V42(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V42(ITERATION(kf))=z;
    z=V43(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V43(ITERATION(kf))=z;
    z=V44(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V44(ITERATION(kf))=z;
 
    z=V60(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V60(ITERATION(kf))=z;
    z=V61(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V61(ITERATION(kf))=z;
    z=V62(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V62(ITERATION(kf))=z;
    z=V63(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V63(ITERATION(kf))=z;
    z=V64(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V64(ITERATION(kf))=z;
    z=V65(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V65(ITERATION(kf))=z;
    z=V66(ITERATION(kf));RT(z)*=renorm;IT(z)*=renorm;V66(ITERATION(kf))=z;
 
    t01 = "===============================================================";
    t02 = "|                                              ---            |";
    t03 = "|Crystal Field parameter  V [%4s] with H     = >   V  *  P (J)|\n";
    t04 = "|                         kq            krist  ---  kq   kq   |";
    t05 = "===============================================================";
    t06 = "|Symmetry : %s             Symmetry number : %1d   |\n";
    t07 = "|---------------------------------------------------------|\n";
    t08 = "|ReV%1d%1d = %13.6e            ImV%1d%1d = %13.6e   |\n"      ;
    t09 = "|  V%1d%1d = %13.6e                                    |\n"   ;
    t10 = "-----------------------------------------------------------\n";
 
    s=symmetrienr;
    fprintf(fp,"%s\n",t01);
    fprintf(fp,"%s\n",t02);
    fprintf(fp,t03,einheit);
    fprintf(fp,"%s\n",t04);
    fprintf(fp,"%s\n",t05);
    fprintf(fp,t06,SYMLISTE[s].symname,s);
    fprintf(fp,"%s",t07);
 
    if( s==1 || s==2 || s==3 || s==4 || s==5 || s==6 || s==7 ){
         rt = RT(  V20( ITERATION(kf) )  );
         fprintf(fp,t09,2,0,rt);}
    if( s==1 || s==2 ){
         rt = RT(  V22( ITERATION(kf) )  );
         fprintf(fp,t09,2,2,rt);}
 
 
         rt = RT(  V40( ITERATION(kf) )  );
         if( s!=8 )fprintf(fp,"%s",t07);
         fprintf(fp,t09,4,0,rt);
    if( s==1 ){
         rt = RT(  V42( ITERATION(kf) )  );
         it = IT(  V42( ITERATION(kf) )  );
         fprintf(fp,t08,4,2,rt,4,2,it);}
    if( s==2 ){
         rt = RT(  V42( ITERATION(kf) )  );
         fprintf(fp,t09,4,2,rt);}
    if( s==5 || s==6 ){
         rt = RT(  V43( ITERATION(kf) )  );
         fprintf(fp,t09,4,3,rt);}
    if( s==1 ){
         rt = RT(  V44( ITERATION(kf) )  );
         it = IT(  V44( ITERATION(kf) )  );
         fprintf(fp,t08,4,4,rt,4,4,it);}
    if( s==2 || s==3 || s==4 || s==8 ){
         rt = RT(  V44( ITERATION(kf) )  );
         fprintf(fp,t09,4,4,rt);}
 
         rt = RT(  V60( ITERATION(kf) )  );
         fprintf(fp,"%s",t07);
         fprintf(fp,t09,6,0,rt);
    if( s==1 ){
         rt = RT(  V62( ITERATION(kf) )  );
         it = IT(  V62( ITERATION(kf) )  );
         fprintf(fp,t08,6,2,rt,6,2,it);}
    if( s==2 ){
         rt = RT(  V62( ITERATION(kf) )  );
         fprintf(fp,t09,6,2,rt);}
    if( s==5 ){
         rt = RT(  V63( ITERATION(kf) )  );
         it = IT(  V63( ITERATION(kf) )  );
         fprintf(fp,t08,6,3,rt,6,3,it);}
    if( s==6 ){
         rt = RT(  V63( ITERATION(kf) )  );
         fprintf(fp,t09,6,3,rt);}
    if( s==1 || s==3 ){
         rt = RT(  V64( ITERATION(kf) )  );
         it = IT(  V64( ITERATION(kf) )  );
         fprintf(fp,t08,6,4,rt,6,4,it);}
    if( s==2 || s==4 || s==8 ){
         rt = RT(  V64( ITERATION(kf) )  );
         fprintf(fp,t09,6,4,rt);}
    if( s==1 || s==5 ){
         rt = RT(  V66( ITERATION(kf) )  );
         it = IT(  V66( ITERATION(kf) )  );
         fprintf(fp,t08,6,6,rt,6,6,it);}
    if( s==2 || s==6 || s==7 ) {
         rt = RT(  V66( ITERATION(kf) )  );
         fprintf(fp,t09,6,6,rt);}
 
         fprintf(fp,"%s",t10);
 
    fclose(fp);
}
/*----------------------------------------------------------------------------
                               info_epsilonkq()
-----------------------------------------------------------------------------*/
void info_epsilonkq()   /* Liste der Faktoren epsilonkq */
{
    CHAR *name = "results/epsilonkq.info";
    FILE *fp,*fopen(const char*,const char*);
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11,*t12,*t13,*t14,*t15,*t16,*t17,*t18,*t19,*t20;
    CHAR *t21,*t22,*t23,*t24,*t25,*t26,*t27,*t28,*t29,*t30;
    CHAR *t31,*t32,*t33,*t34,*t35,*t36,*t37,*t38,*t39,*t40;
    CHAR *t41,*t42,*t43,*t44,*t45,*t46,*t47,*t48,*t49,*t50;
    CHAR *s01,*s02,*s03,*s04,*s05,*s06,*s07;
 
    fp = fopen_errchk(name,"w");
    clearscreen;
    printf("Information given in the File %s .\n",name);
 
    t01 = "===========================================================";
    t02 = "|Table of epsilon Kq factors          |";
    t03 = "===========================================================";
    t04 = "|                                                         |";
    t05 = "| Given by :                                               |";
    t06 = "|                             2S+1                        |";
    t07 = "| C    = epsilon   *  theta (     L   ) * P  ( J )        |";
    t08 = "|  kq           kq         k        J      kq  -          |";
    t09 = "|---------------------------------------------------------|";
    t10 = "|                                                |q|+q    |";
    t11 = "|                                                -----    |";
    t12 = "|                                                  2      |";
    t13 = "|                  1          (2k)!          (-1)         |";
    t14 = "| epsilon   = ------------ * ------- * ------------------ |";
    t15 = "|        kq      h             k       sqrt([k-q]!*[k+q]!)|";
    t16 = "|                 k,|q|       2  k!                       |";
    t17 = "|                                                         |";
    t18 = "| h      =  h                                             |";
    t19 = "|  k,|q|     n+tau,n   Implemented are:                |";
    t20 = "|---------------------------------------------------------|";
    t21 = "|                                                         |";
    t22 = "|                                                         |";
    t23 = "| h      =  1                                             |";
    t24 = "|  n,n                                                    |";
    t25 = "|                                                         |";
    t26 = "|                                                         |";
    t27 = "| h      =  1                                             |";
    t28 = "|  n+1,n                                                  |";
    t29 = "|                                                         |";
    t30 = "|                                                         |";
    t31 = "| h      =  2n+3                                          |";
    t32 = "|  n+2,n                                                  |";
    t33 = "|                                                         |";
    t34 = "|                                                         |";
    t35 = "| h      =  2n+5                                          |";
    t36 = "|  n+3,n                                                  |";
    t37 = "|                                                         |";
    t38 = "|                                                         |";
    t39 = "| h      =  (2n+5)(2n+7)                                  |";
    t40 = "|  n+4,n                                                  |";
    t41 = "|                                                         |";
    t42 = "|                                                         |";
    t43 = "| h      =  (2n+7)(2n+9)                                  |";
    t44 = "|  n+5,n                                                  |";
    t45 = "|                                                         |";
    t46 = "|            (2n+7)(2n+9)(2n+11)                          |";
    t47 = "| h      =  ---------------------                         |";
    t48 = "|  n+6,n              3                                   |";
    t49 = "|                                                         |";
    t50 = "-----------------------------------------------------------";
 
 
    fprintf(fp,"%s\n",t01);
    fprintf(fp,"%s\n",t02);
    fprintf(fp,"%s\n",t03);
    fprintf(fp,"%s\n",t04);
    fprintf(fp,"%s\n",t05);
    fprintf(fp,"%s\n",t06);
    fprintf(fp,"%s\n",t07);
    fprintf(fp,"%s\n",t08);
    fprintf(fp,"%s\n",t09);
    fprintf(fp,"%s\n",t10);
 
    fprintf(fp,"%s\n",t11);
    fprintf(fp,"%s\n",t12);
    fprintf(fp,"%s\n",t13);
    fprintf(fp,"%s\n",t14);
    fprintf(fp,"%s\n",t15);
    fprintf(fp,"%s\n",t16);
    fprintf(fp,"%s\n",t17);
    fprintf(fp,"%s\n",t18);
    fprintf(fp,"%s\n",t19);
    fprintf(fp,"%s\n",t20);
 
    fprintf(fp,"%s\n",t21);
    fprintf(fp,"%s\n",t22);
    fprintf(fp,"%s\n",t23);
    fprintf(fp,"%s\n",t24);
    fprintf(fp,"%s\n",t25);
    fprintf(fp,"%s\n",t26);
    fprintf(fp,"%s\n",t27);
    fprintf(fp,"%s\n",t28);
    fprintf(fp,"%s\n",t29);
    fprintf(fp,"%s\n",t30);
 
    fprintf(fp,"%s\n",t31);
    fprintf(fp,"%s\n",t32);
    fprintf(fp,"%s\n",t33);
    fprintf(fp,"%s\n",t34);
    fprintf(fp,"%s\n",t35);
    fprintf(fp,"%s\n",t36);
    fprintf(fp,"%s\n",t37);
    fprintf(fp,"%s\n",t38);
    fprintf(fp,"%s\n",t39);
    fprintf(fp,"%s\n",t40);
 
    fprintf(fp,"%s\n",t41);
    fprintf(fp,"%s\n",t42);
    fprintf(fp,"%s\n",t43);
    fprintf(fp,"%s\n",t44);
    fprintf(fp,"%s\n",t45);
    fprintf(fp,"%s\n",t46);
    fprintf(fp,"%s\n",t47);
    fprintf(fp,"%s\n",t48);
    fprintf(fp,"%s\n",t49);
    fprintf(fp,"%s\n",t50);
 
 
     s01 = "======================== \n";
    s02 = "�        �             � \n";
    s03 = "�  k   q � epsilon     � \n";
    s04 = "�        �        kq   � \n";
    s05 = "======================== \n";
    s06 = "� %2d %2d � %12.9f �    \n";
    s07 = "------------------------ \n";


     fprintf(fp,"\n\n");
 
     fprintf(fp,"%s",s01);
     fprintf(fp,"%s",s02);
     fprintf(fp,"%s",s03);
     fprintf(fp,"%s",s04);
     fprintf(fp,"%s",s05);
 
     fprintf(fp,s06,0,0,epn0n(0) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,1,-1,epn0n(-1) );
     fprintf(fp,s06,1, 0,epn1n( 0) );
     fprintf(fp,s06,1, 1,epn0n( 1) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,2,-2,epn0n(-2) );
     fprintf(fp,s06,2,-1,epn1n(-1) );
     fprintf(fp,s06,2, 0,epn2n( 0) );
     fprintf(fp,s06,2, 1,epn1n( 1) );
     fprintf(fp,s06,2, 2,epn0n( 2) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,3,-3,epn0n(-3) );
     fprintf(fp,s06,3,-2,epn1n(-2) );
     fprintf(fp,s06,3,-1,epn2n(-1) );
     fprintf(fp,s06,3, 0,epn3n( 0) );
     fprintf(fp,s06,3, 1,epn2n( 1) );
     fprintf(fp,s06,3, 2,epn1n( 2) );
     fprintf(fp,s06,3, 3,epn0n( 3) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,4,-4,epn0n(-4) );
     fprintf(fp,s06,4,-3,epn1n(-3) );
     fprintf(fp,s06,4,-2,epn2n(-2) );
     fprintf(fp,s06,4,-1,epn3n(-1) );
     fprintf(fp,s06,4, 0,epn4n( 0) );
     fprintf(fp,s06,4, 1,epn3n( 1) );
     fprintf(fp,s06,4, 2,epn2n( 2) );
     fprintf(fp,s06,4, 3,epn1n( 3) );
     fprintf(fp,s06,4, 4,epn0n( 4) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,5,-5,epn0n(-5) );
     fprintf(fp,s06,5,-4,epn1n(-4) );
     fprintf(fp,s06,5,-3,epn2n(-3) );
     fprintf(fp,s06,5,-2,epn3n(-2) );
     fprintf(fp,s06,5,-1,epn4n(-1) );
     fprintf(fp,s06,5, 0,epn5n( 0) );
     fprintf(fp,s06,5, 1,epn4n( 1) );
     fprintf(fp,s06,5, 2,epn3n( 2) );
     fprintf(fp,s06,5, 3,epn2n( 3) );
     fprintf(fp,s06,5, 4,epn1n( 4) );
     fprintf(fp,s06,5, 5,epn0n( 5) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,6,-6,epn0n(-6) );
     fprintf(fp,s06,6,-5,epn1n(-5) );
     fprintf(fp,s06,6,-4,epn2n(-4) );
     fprintf(fp,s06,6,-3,epn3n(-3) );
     fprintf(fp,s06,6,-2,epn4n(-2) );
     fprintf(fp,s06,6,-1,epn5n(-1) );
     fprintf(fp,s06,6, 0,epn6n( 0) );
     fprintf(fp,s06,6, 1,epn5n( 1) );
     fprintf(fp,s06,6, 2,epn4n( 2) );
     fprintf(fp,s06,6, 3,epn3n( 3) );
     fprintf(fp,s06,6, 4,epn2n( 4) );
     fprintf(fp,s06,6, 5,epn1n( 5) );
     fprintf(fp,s06,6, 6,epn0n( 6) );
     fprintf(fp,"%s",s07);
 
     fclose(fp);
}
 
/*----------------------------------------------------------------------------
                               info_hamilton()
-----------------------------------------------------------------------------*/
void info_hamilton() /* Liste der Auswahlregeln fuer die Vkq im Hamiltonian*/
{
    CHAR *name = "results/hamilton.info";
    FILE *fp,*fopen(const char*,const char*);
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11,*t12,*t13,*t14,*t15,*t16,*t17,*t18,*t19,*t20;
    CHAR *t21,*t22,*t23,*t24,*t25,*t26,*t27,*t28,*t29,*t30;
    CHAR *t31,*t32;

   
 
    fp = fopen_errchk(name,"w");
    clearscreen;
    printf("Information given in the File %s .\n",name);
 
t01 = "=======================================================================";
t02 = "|Independant non-vanishing Vkq for the different         |";
t03 = "|Symmetry numbers   ( v := or )                                     |";
t04 = "=======================================================================";
t05 = "";
t06 = "=======================================================================";
t07 = "|   |V20|V21|V22|**|V40|V41|V42|V43|V44|**|V60|V61|V62|V63|V64|V65|V66|";
t08 = "|===|===|===|===|**|===|===|===|===|===|**|===|===|===|===|===|===|===|";
t09 = "| 0 | R |RvI| C |**| R | C | C | C | C |**| R | C | C | C | C | C | C |";
t10 = "|---|---|---|---|**|---|---|---|---|---|**|---|---|---|---|---|---|---|";
t11 = "| 1 | R |   |RvI|**| R |   | C |   | C |**| R |   | C |   | C |   | C |";
t12 = "|---|---|---|---|**|---|---|---|---|---|**|---|---|---|---|---|---|---|";
t13 = "| 2 | R |   |RvI|**| R |   |RvI|   |RvI|**| R |   |RvI|   |RvI|   |RvI|";
t14 = "|---|---|---|---|**|---|---|---|---|---|**|---|---|---|---|---|---|---|";
t15 = "| 3 | R |   |   |**| R |   |   |   |RvI|**| R |   |   |   | C |   |   |";
t16 = "|---|---|---|---|**|---|---|---|---|---|**|---|---|---|---|---|---|---|";
t17 = "| 4 | R |   |   |**| R |   |   |   |RvI|**| R |   |   |   |RvI|   |   |";
t18 = "|---|---|---|---|**|---|---|---|---|---|**|---|---|---|---|---|---|---|";
t19 = "| 5 | R |   |   |**| R |   |   |RvI|   |**| R |   |   | C |   |   | C |";
t20 = "|---|---|---|---|**|---|---|---|---|---|**|---|---|---|---|---|---|---|";
t21 = "| 6 | R |   |   |**| R |   |   |RvI|   |**| R |   |   |RvI|   |   |RvI|";
t22 = "|---|---|---|---|**|---|---|---|---|---|**|---|---|---|---|---|---|---|";
t23 = "| 7 | R |   |   |**| R |   |   |   |   |**| R |   |   |   |   |   |RvI|";
t24 = "|---|---|---|---|**|---|---|---|---|---|**|---|---|---|---|---|---|---|";
t25 = "|8*)|   |   |   |**| R |   |   |   | R |**| R |   |   |   | R |   |   |";
t26 = "|===|===|===|===|**|===|===|===|===|===|**|===|===|===|===|===|===|===|";
t27 = "|   |V20|V21|V22|**|V40|V41|V42|V43|V44|**|V60|V61|V62|V63|V64|V65|V66|";
t28 = "=======================================================================";
t29 = "|*)In the cubic case: V44 = 5/2*V40 und  V64 = -21/2*V60 |";
t30 = "|----------------------------------------------------------------------";
t31 = "| C := complex Parameter  R/I := Real/Imaginary part   |";
t32 = "-----------------------------------------------------------------------";
 
 
    fprintf(fp,"%s\n",t01);
    fprintf(fp,"%s\n",t02);
    fprintf(fp,"%s\n",t03);
    fprintf(fp,"%s\n",t04);
    fprintf(fp,"%s\n",t05);
    fprintf(fp,"%s\n",t06);
    fprintf(fp,"%s\n",t07);
    fprintf(fp,"%s\n",t08);
    fprintf(fp,"%s\n",t09);
    fprintf(fp,"%s\n",t10);
 
    fprintf(fp,"%s\n",t11);
    fprintf(fp,"%s\n",t12);
    fprintf(fp,"%s\n",t13);
    fprintf(fp,"%s\n",t14);
    fprintf(fp,"%s\n",t15);
    fprintf(fp,"%s\n",t16);
    fprintf(fp,"%s\n",t17);
    fprintf(fp,"%s\n",t18);
    fprintf(fp,"%s\n",t19);
    fprintf(fp,"%s\n",t20);
 
    fprintf(fp,"%s\n",t21);
    fprintf(fp,"%s\n",t22);
    fprintf(fp,"%s\n",t23);
    fprintf(fp,"%s\n",t24);
    fprintf(fp,"%s\n",t25);
    fprintf(fp,"%s\n",t26);
    fprintf(fp,"%s\n",t27);
    fprintf(fp,"%s\n",t28);
    fprintf(fp,"%s\n",t29);
    fprintf(fp,"%s\n",t30);
    fprintf(fp,"%s\n",t31);
    fprintf(fp,"%s\n",t32);
 
     fclose(fp);
}
 
/*----------------------------------------------------------------------------
                               info_omegakq()
-----------------------------------------------------------------------------*/
void info_omegakq(INT k,INT q)   /* Liste der Faktoren omegakq */
{
    CHAR *name = "results/omegakq.info";
    FILE *fp,*fopen(const char*,const char*);
    CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
    CHAR *t11,*t12,*t13,*t14,*t15,*t16,*t17,*t18,*t19,*t20;
    CHAR *t21,*t22,*t23,*t24,*t25,*t26,*t27,*t28,*t29,*t30;
    CHAR *t31,*t32,*t33,*t34,*t35,*t36,*t37,*t38,*t39,*t40;
    CHAR *t41,*t42,*t43,*t44,*t45,*t46,*t47,*t48,*t49,*t50;
    CHAR *t51,*t52,*t53,*t54,*t55;
    CHAR *s01,*s02,*s03,*s04,*s05,*s06,*s07;
    DOUBLE omegan0n(INT);
    DOUBLE omegan1n(INT);
    DOUBLE omegan2n(INT);
    DOUBLE omegan3n(INT);
    DOUBLE omegan4n(INT);
    DOUBLE omegan5n(INT);
    DOUBLE omegan6n(INT),z;
 
   switch( ABS( k-ABS(q) ) ){
       case 0 : z=omegan0n(ABS(q));break;
       case 1 : z=omegan1n(ABS(q));break;
       case 2 : z=omegan2n(ABS(q));break;
       case 3 : z=omegan3n(ABS(q));break;
       case 4 : z=omegan4n(ABS(q));break;
       case 5 : z=omegan5n(ABS(q));break;
       case 6 : z=omegan6n(ABS(q));break;
       default: fprintf(stderr,"info_omegakq: Invalid k=%i\n",k); z=0.0; break;
   }
 
   printf("\n");
   printf("===================\n");
   printf("|  k  q | omega   |\n");
   printf("|       |      kq |\n");
   printf("===================\n");
   printf("| %2d %2d | %3.0f     |\n",k,ABS(q),z);
   printf("===================\n\n");
 
 
 
 
    fp = fopen_errchk(name,"w");
    printf("Information given in the File %s .\n",name);
 
    t01 = "===========================================================";
    t02 = "|Table of the factors Thetakq             |";
    t03 = "===========================================================";
    t04 = "|                                                         |";
    t05 = "| Es gilt :                                               |";
    t06 = "|                                                         |";
    t07 = "|    P  ( J )  =  omega   *  O   ( J )                    |";
    t08 = "|     kq  -            kq      kq  -                      |";
    t09 = "|---------------------------------------------------------|";
    t10 = "|                                                         |";
    t11 = "|  wobei :                                                |";
    t12 = "|                                                         |";
    t13 = "|   STEV    =  (  O    +  O     ) / 2                     |";
    t14 = "|       kq         kq      k,-q                           |";
    t15 = "|---------------------------------------------------------|";
    t16 = "|                                                         |";
    t17 = "|                                                         |";
    t18 = "| omega       =  omega       =  omega       =   1         |";
    t19 = "|      n,n            n+1,n          n+2,n                |";
    t20 = "|                                                         |";
    t21 = "|                                                         |";
    t22 = "|---------------------------------------------------------|";
    t23 = "|                                                         |";
    t24 = "| omega       =  1 + 2*D                                  |";
    t25 = "|       n+3,n           n;{2,5,8,11,..}                   |";
    t26 = "|                                                         |";
    t27 = "|                                                         |";
    t28 = "|---------------------------------------------------------|";
    t29 = "|                                                         |";
    t30 = "| omega       =  1 + 2*D             + 2*D                |";
    t31 = "|       n+4,n           n;{1,4,7,..}      n;{2,5,8,..}    |";
    t32 = "|                                                         |";
    t33 = "|                                                         |";
    t34 = "|---------------------------------------------------------|";
    t35 = "|                                                         |";
    t36 = "| omega       =  (  1 + 2*D             )                 |";
    t37 = "|       n+5,n              n;{1,4,7,..}                   |";
    t38 = "|                                                         |";
    t39 = "|                                                         |";
    t40 = "|             *  (  1 + 4*D          + 4*D           )    |";
    t41 = "|                          n;{3,8,..}     n;{4,9,..}      |";
    t42 = "|                                                         |";
    t43 = "|                                                         |";
    t44 = "|---------------------------------------------------------|";
    t45 = "|                                                         |";
    t46 = "| omega       =  1+4*D         +4*D         +4*D          |";
    t47 = "|      n+6,n          n;{2,7,..}   n;{3,8,..}   n;{4,9,..}|";
    t48 = "|                                                         |";
    t49 = "-----------------------------------------------------------";
    t50 = "|                                                         |";
    t51 = "|             /  1  , if n is an element of m               |";
    t52 = "|  D     :=   |                                           |";
    t53 = "|   n;M       \\  0  ,if n is not an element of m          |";
    t54 = "|                                                         |";
    t55 = "-----------------------------------------------------------";
 
 
    fprintf(fp,"%s\n",t01);
    fprintf(fp,"%s\n",t02);
    fprintf(fp,"%s\n",t03);
    fprintf(fp,"%s\n",t04);
    fprintf(fp,"%s\n",t05);
    fprintf(fp,"%s\n",t06);
    fprintf(fp,"%s\n",t07);
    fprintf(fp,"%s\n",t08);
    fprintf(fp,"%s\n",t09);
    fprintf(fp,"%s\n",t10);
 
    fprintf(fp,"%s\n",t11);
    fprintf(fp,"%s\n",t12);
    fprintf(fp,"%s\n",t13);
    fprintf(fp,"%s\n",t14);
    fprintf(fp,"%s\n",t15);
    fprintf(fp,"%s\n",t16);
    fprintf(fp,"%s\n",t17);
    fprintf(fp,"%s\n",t18);
    fprintf(fp,"%s\n",t19);
    fprintf(fp,"%s\n",t20);
 
    fprintf(fp,"%s\n",t21);
    fprintf(fp,"%s\n",t22);
    fprintf(fp,"%s\n",t23);
    fprintf(fp,"%s\n",t24);
    fprintf(fp,"%s\n",t25);
    fprintf(fp,"%s\n",t26);
    fprintf(fp,"%s\n",t27);
    fprintf(fp,"%s\n",t28);
    fprintf(fp,"%s\n",t29);
    fprintf(fp,"%s\n",t30);
 
    fprintf(fp,"%s\n",t31);
    fprintf(fp,"%s\n",t32);
    fprintf(fp,"%s\n",t33);
    fprintf(fp,"%s\n",t34);
    fprintf(fp,"%s\n",t35);
    fprintf(fp,"%s\n",t36);
    fprintf(fp,"%s\n",t37);
    fprintf(fp,"%s\n",t38);
    fprintf(fp,"%s\n",t39);
    fprintf(fp,"%s\n",t40);
 
    fprintf(fp,"%s\n",t41);
    fprintf(fp,"%s\n",t42);
    fprintf(fp,"%s\n",t43);
    fprintf(fp,"%s\n",t44);
    fprintf(fp,"%s\n",t45);
    fprintf(fp,"%s\n",t46);
    fprintf(fp,"%s\n",t47);
    fprintf(fp,"%s\n",t48);
    fprintf(fp,"%s\n",t49);
    fprintf(fp,"%s\n",t50);
 
    fprintf(fp,"%s\n",t51);
    fprintf(fp,"%s\n",t52);
    fprintf(fp,"%s\n",t53);
    fprintf(fp,"%s\n",t54);
    fprintf(fp,"%s\n",t55);
 
    s01 = "======================= \n";
    s02 = "�       �             � \n";
    s03 = "�  k  q �   omega     � \n";
    s04 = "�       �        kq   � \n";
    s05 = "======================= \n";
    s06 = "� %2d %2d �    %3.0f      �  \n";
    s07 = "----------------------- \n";
 
 
     fprintf(fp,"\n\n");
 
     fprintf(fp,"%s",s01);
     fprintf(fp,"%s",s02);
     fprintf(fp,"%s",s03);
     fprintf(fp,"%s",s04);
     fprintf(fp,"%s",s05);
 
     fprintf(fp,s06,0,0,omegan0n(0) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,1,-1,omegan0n( 1) );
     fprintf(fp,s06,1, 0,omegan1n( 0) );
     fprintf(fp,s06,1, 1,omegan0n( 1) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,2,-2,omegan0n( 2) );
     fprintf(fp,s06,2,-1,omegan1n( 1) );
     fprintf(fp,s06,2, 0,omegan2n( 0) );
     fprintf(fp,s06,2, 1,omegan1n( 1) );
     fprintf(fp,s06,2, 2,omegan0n( 2) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,3,-3,omegan0n( 3) );
     fprintf(fp,s06,3,-2,omegan1n( 2) );
     fprintf(fp,s06,3,-1,omegan2n( 1) );
     fprintf(fp,s06,3, 0,omegan3n( 0) );
     fprintf(fp,s06,3, 1,omegan2n( 1) );
     fprintf(fp,s06,3, 2,omegan1n( 2) );
     fprintf(fp,s06,3, 3,omegan0n( 3) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,4,-4,omegan0n( 4) );
     fprintf(fp,s06,4,-3,omegan1n( 3) );
     fprintf(fp,s06,4,-2,omegan2n( 2) );
     fprintf(fp,s06,4,-1,omegan3n( 1) );
     fprintf(fp,s06,4, 0,omegan4n( 0) );
     fprintf(fp,s06,4, 1,omegan3n( 1) );
     fprintf(fp,s06,4, 2,omegan2n( 2) );
     fprintf(fp,s06,4, 3,omegan1n( 3) );
     fprintf(fp,s06,4, 4,omegan0n( 4) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,5,-5,omegan0n(-5) );
     fprintf(fp,s06,5,-4,omegan1n( 4) );
     fprintf(fp,s06,5,-3,omegan2n( 3) );
     fprintf(fp,s06,5,-2,omegan3n( 2) );
     fprintf(fp,s06,5,-1,omegan4n( 1) );
     fprintf(fp,s06,5, 0,omegan5n( 0) );
     fprintf(fp,s06,5, 1,omegan4n( 1) );
     fprintf(fp,s06,5, 2,omegan3n( 2) );
     fprintf(fp,s06,5, 3,omegan2n( 3) );
     fprintf(fp,s06,5, 4,omegan1n( 4) );
     fprintf(fp,s06,5, 5,omegan0n( 5) );
     fprintf(fp,"%s",s07);
 
     fprintf(fp,s06,6,-6,omegan0n( 6) );
     fprintf(fp,s06,6,-5,omegan1n( 5) );
     fprintf(fp,s06,6,-4,omegan2n( 4) );
     fprintf(fp,s06,6,-3,omegan3n( 3) );
     fprintf(fp,s06,6,-2,omegan4n( 2) );
     fprintf(fp,s06,6,-1,omegan5n( 1) );
     fprintf(fp,s06,6, 0,omegan6n( 0) );
     fprintf(fp,s06,6, 1,omegan5n( 1) );
     fprintf(fp,s06,6, 2,omegan4n( 2) );
     fprintf(fp,s06,6, 3,omegan3n( 3) );
     fprintf(fp,s06,6, 4,omegan2n( 4) );
     fprintf(fp,s06,6, 5,omegan1n( 5) );
     fprintf(fp,s06,6, 6,omegan0n( 6) );
     fprintf(fp,"%s",s07);
 
     fclose(fp);
}
 
/*------------------------------------------------------------------------------
                           c_single_error()
------------------------------------------------------------------------------*/
void c_single_error()   /* Eingabefehler beim Versuch von  -c -s */
{
    clearscreen;
    printf("\n");
    printf("The syntax is :\n");
    printf("\n");
    printf("%s -c[reate] -s[ingleion] ",PROGRAMMNAME);
    printf("filename.filetype rk|sk|pk nn ion temp\n");
    printf("\n");
    printf("rk : If the next neighbours are in the right angled co-ordination\n");
    printf("     system .\n");
    printf("sk :If the next neighbours are in the spherical co-ordination\n");
    printf("     system.\n");
    printf("pk : If the next neighbours are in the polar co-ordination\n");
    printf("     system.\n");
    printf("nn  : number of nearest neighbours > 0 \n");
    printf("ion : Name of the ion in the surrounding enviroment\n");
    printf("temp: Temperature of the sample\n");
    printf("\n");
 
    exit(1);
}
/*------------------------------------------------------------------------------
                           i_kq_error()
------------------------------------------------------------------------------*/
void i_kq_error(CHAR c)   /* Eingabefehler beim Versuch von  -i e 'c' */
{
   clearscreen;
   printf("\n");
   printf("Die Eingabesyntax ist :\n");
   printf("\n");
   if( c != 'x' &&  c != 'X' )
    printf("%s -i[nfo] ew[problem] %c[kq] [symmetrienummer] \n",PROGRAMMNAME,c);
   else
    printf("%s -i[nfo] ew[problem]  %c[W]                   \n",PROGRAMMNAME,c);
   printf("\n");
 
    exit(1);
}
/*------------------------------------------------------------------------------
                           r_kq_error()
------------------------------------------------------------------------------*/
void r_kq_error(CHAR c,CHAR cs)   /* Eingabefehler beim Versuch von  -r -'c', -s -'c' */
{
    CHAR *s;
 
    clearscreen;
    printf("\n");
    printf("Die Eingabesyntax ist :\n");
    printf("\n");
 
 if( c=='r'|| c=='R'){
    if( cs!='x' && cs!='X' )
      printf("%s -r[ead] -%c[kq] [symmetrienummer] \n",PROGRAMMNAME,cs);
    else
      printf("%s -r[ead] -%c[W ]                   \n",PROGRAMMNAME,cs);
    printf("\n");
 }
 
 if( c=='o'|| c=='O'){
    if( cs!='x' && cs!='X' )
      printf("%s -o[rtho] -%c[kq] \n",PROGRAMMNAME,cs);
    else
      printf("%s -o[rtho] -%c[W ] \n",PROGRAMMNAME,cs);
    printf("\n");
 }
 
 if( c=='b'|| c=='B'){
    if( cs!='x' && cs!='X' )
      printf("%s -b[eside] -%c[kq] \n",PROGRAMMNAME,cs);
    else
      printf("%s -b[eside] -%c[W ] \n",PROGRAMMNAME,cs);
    printf("\n");
 }
 
 if( c=='s'|| c=='S'){
   s = PROGRAMMNAME;
   if( cs!='x' && cs!='X' )
printf("%s -s[uszept] -%c[kq] anf_temp end_temp lambda [theta]||[filename.type]\n",s,cs);
   else
printf("%s -s[uszept] -%c[W ] anf_temp end_temp lambda [theta]||[filename.type]\n",s,cs);
   printf("\n");
 }
 if( c=='m'|| c=='M'){
   s = PROGRAMMNAME;
   if( cs!='x' && cs!='X' )
printf("%s -m[oment] -%c[kq] [symmetrienr] anf_feld end_feld temp\n",s,cs);
   else
printf("%s -m[oment] -%c[W ] anf_feld end_feld temp\n",s,cs);
   printf("\n");
 }
 if( c=='k'|| c=='K'){
   s = PROGRAMMNAME;
   printf("%s -k[poly ] -%c[W ] anf_feld end_feld temp\n",s,cs);
   printf("\n");
 }
    exit(1);
}
/*------------------------------------------------------------------------------
                           c_kq_error()
------------------------------------------------------------------------------*/
void c_kq_error(CHAR c)   /* Eingabefehler beim Versuch von  -c -'c' */
{
    INT  i;
    /*CHAR *s="-c[reate] -c[kq]";*/
 
    clearscreen;
    printf("\n");
    printf("The syntax is :\n");
    printf("\n");
    if( c!='x' && c!='X' )
       printf("%s -c[reate] -%c[kq] ",PROGRAMMNAME,c);
    else
       printf("%s -c[reate] -%c[W ] ",PROGRAMMNAME,c);
 
    printf("E_P E_E ion ");
    if( c!='x' && c!='X' )  printf("symmetry number ");
    printf("temp [rk|sk|pk] \n\n");
 
 
    printf("e_p := Unit in which the parameters ");
    printf("are given.\n");
    printf("e_e := Unit in which the energy eigenvalues ");
    printf("are calculated.\n\n");
    printf("ion symmetry number: ");
    for(i=0;i<(int)ANZ_INFO_EINHEIT-1;++i)
       printf("%s,",INFO_EINHEIT[ i ].info);
    printf("%s\n\n",INFO_EINHEIT[ ANZ_INFO_EINHEIT-1 ].info);
 
    if( c!='x' && c!='X' )
        printf("symmetry number : between 0 und 8\n");
    printf("ion         : Name of the ion\n");
    printf("temp        : Temperature of the sample in Kelvin\n\n");
    printf("rk : If the magnetic field is given in right-angled ");
    printf("co-ordinates.\n");
    printf("sk : If the magnetic field is given in spherical");
    printf("co-ordinates.\n");
    printf("pk : If the magnetic field is given in polar");
    printf("co-ordinates.\n");
    printf("\n");
 
    exit(1);
}
/*------------------------------------------------------------------------------
                           create_error()
------------------------------------------------------------------------------*/
void create_error()   /* Eingabefehler beim Versuch von  -c */
{
    clearscreen;
    printf("\n");
    printf("The commands :\n");
    printf("\n");
    printf("%s -c[reate] -s[ingleion] \n",PROGRAMMNAME);
    printf("%s -c[reate] -A[kq] (see Hutchings p.255)\n",PROGRAMMNAME);
    printf("%s -c[reate] -B[kq] (see Hutchings p.255)\n",PROGRAMMNAME);
    printf("%s -c[reate] -D[kq]       \n",PROGRAMMNAME);
    printf("%s -c[reate] -L[kq] (Wybourne-Parameter) \n",PROGRAMMNAME);
    printf("%s -c[reate] -V[kq] (NOT Stevens Parameter) \n",PROGRAMMNAME);
    printf("%s -c[reate] -W[kq]       \n",PROGRAMMNAME);
    printf("%s -c[reate] -x[W ]       \n",PROGRAMMNAME);
    printf("\n");
 
    exit(1);
}
/*------------------------------------------------------------------------------
                                 leftcopy()
------------------------------------------------------------------------------*/
CHAR *leftcopy(CHAR*string,INT bufferlen)
{
    CHAR *buffer;
    INT  i,len;
 
    buffer = STRING_ALLOC(bufferlen+2);
 
    len  = strlen_own(string);
    for( i=0 ; i<=len; ++i)
         VALUE(buffer,i) = VALUE(string,i);
    for( i=len ; i<=bufferlen; ++i)
         VALUE(buffer,i) = ' ';
    VALUE(buffer,bufferlen+1) ='\0';
 
    return( buffer );
}
/*------------------------------------------------------------------------------
                                  strlen_own()
------------------------------------------------------------------------------*/
INT strlen_own(CHAR*s)
{
    INT len=0;
    while(    *(s+len++) !='\0'    );
    return(--len);
}
/*------------------------------------------------------------------------------
                               read_error()
------------------------------------------------------------------------------*/
void read_error(INT nr,FILE*fp,CHAR*name) /* Eingabefehler beim Versuch von  -r -s*/
{
    printf("\n");
    clearscreen;
 
    switch(nr){
    case 1: printf("The syntax is :\n");
            printf("\n");
            printf("%s -r[ead] -s[ingleion] ",PROGRAMMNAME);
            printf("filename.filetype symmetrienr\n");
            printf("\n");
           printf("symmetry number                       : Zahl zw 1 und 8\n");
           printf("Table of the symmetry numbers: %s -i s\n",PROGRAMMNAME);
            printf("\n");
            exit(1);
 
    case 2: printf("Error !\n");
            printf("The File : %s doesnt exist !\n",name);
            printf("\n");
            exit(1);
 
    case 3: printf("Error !\n");
            printf("cannot in File %s  determine the Mode rk|sk|pk \n",name);
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 4: printf("Error !\n");
            printf("Mistake in %s !\n",name);
            printf("\n");
            fclose(fp);
            exit(1);
 
 
    case 6: printf("Error !\n");
            printf("Wrong Symmetry number !\n");
            printf("It must be: 0<= Symmetry number <= 8 !\n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 7: printf("The syntax is :\n");
            printf("\n");
            printf("%s -i[nfo] c[lm] l m theta phi\n",PROGRAMMNAME);
            printf("\n");
            exit(1);
 
    case 8: printf("Error !\n");
            printf("The angluar momentum l must be >= 0 !\n");
            printf("\n");
            exit(1);
 
    case 9: printf("Error !\n");
            printf("The absolute values ");
            printf("of the angular momentum must be <= %s \n",name);
            printf("\n");
            exit(1);
 
    case 10:printf("Error !\n");
            printf("It must: 0<= Polar angle theta <= 180 !\n");
            printf("\n");
            exit(1);
 
    case 11:printf("Error !\n");
            printf("It must: 0<= Azimutal angle phi <= 360 !\n");
            printf("\n");
            exit(1);
 
    case 12:printf("The syntax is :\n");
            printf("\n");
            printf("%s -i[nfo] m[agnetfeld] J  Bx By Bz\n",PROGRAMMNAME);
            printf("The total angluar momentum J >= 0 ist.\n");
            printf("\n");
            exit(1);
 
    case 13:printf("Error !\n");
            printf("The total angluar momentum J must be >= 0 !\n");
            printf("\n");
            exit(1);
 
    case 14:printf("The syntax is :\n");
            printf("\n");
printf("%s -i[nfo] v[lm] eV|cm-1 filename.filetype symmetrienr\n",PROGRAMMNAME);
printf("\n");
printf("symmetry number                       : Zahl zw 1 und 8\n");
printf("Table of the Symmetry numbers: %s -i s\n",PROGRAMMNAME);
printf("\n");
exit(1);
 
 
    case 15:printf("The syntax is :\n");
            printf("\n");
            printf("%s -i[nfo] P[kq(J)] G[GT]|n[orm] k q J\n",PROGRAMMNAME);
            printf("\n");
            printf("GGT : <JM'|Pkq(J)|JM> by the biggest ");
            printf("common denominator\n");
            printf("      The Matrix elements <JM'|Pkq(J)|JM> .\n");
            printf("norm: <JM'|Pkq(J)|JM> in standard unit ");
            printf("of Pkq.\n");
            printf("\n");
            exit(1);
 
    case 16:printf("Error !\n");
            printf("The angular momentum k must be >= 0 !\n");
            printf("\n");
            exit(1);
 
    case 17:printf("Error !\n");
            printf("It must be | k-|q| |<=6 !\n");
            printf("\n");
            exit(1);
 
    case 18:printf("The syntax is :\n");
            printf("\n");
            printf("%s -i[nfo] S[TEVkq(J)] G[GT]|n[orm] k q J\n",PROGRAMMNAME);
            printf("\n");
            printf("GGT : <JM'|STEVkq(J)|JM> by the biggest ");
            printf("common denominator \n");
            printf("      The Matrix element <JM'|STEVkq(J)|JM> .\n");
            printf("norm: <JM'|STEVkq(J)|JM> in standard unit ");
            printf("of STEVkq.\n");
            printf("\n");
            exit(1);
 
    case 20:printf("The syntax is :\n");
            printf("\n");
            printf("%s -i[nfo] o[megakq] k q \n",PROGRAMMNAME);
            printf("\n");
            exit(1);
 
    case 21:printf("Error in %s !\n",name);
            printf("Energy unit in %s not recognised .\n",name);
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 22:printf("Error in %s !\n",name);
            printf("It must be 2*V44 =  5*V40  .\n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 23:printf("Error in %s !\n",name);
            printf("it must be  2*V64 = - 21*V60  .\n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 24:printf("Error in %s !\n",name);
            printf("The Magnetfield mode ");
            printf("not recognised.\n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 25:printf("Error in %s !\n",name);
            printf("It must be  B44 =  5*B40  .\n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 26:printf("Error in %s !\n",name);
            printf("It must be  B64 = - 21*B60 .\n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 27:printf("Error !\n");
            printf("The Temperature  is <0 !\n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 28:printf("Error in %s !\n",name);
            printf("It must be  2*W44 =  5*W40  .\n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 29:printf("Error in %s !\n",name);
            printf("It must  2*W64 = - 21*W60  be .\n");
            printf("\n");
            fclose(fp);
            exit(1);
 
 
    case 30:printf("Error !\n");
            printf("It must  2*V44 =  5*V40  be .\n");
            printf("Point charge model does not fulfill \n");
            printf("cubic symmetry!\n");
            fclose(fp);
            exit(1);
 
    case 31:printf("Error in %s !\n",name);
            printf("It must be  2*V64 = - 21*V60  .\n");
            printf("Point charge model does not fulfill \n");
            printf("cubic symmetry!\n");
            fclose(fp);
            exit(1);
 
    case 32:printf("The syntax is :\n");
            printf("\n");
            printf("%s -i[nfo] ew[problem] -s[ingleion] ",PROGRAMMNAME);
            printf("filename.filetype symmetrienr\n");
            printf("\n");
           printf("symmetry number                       : Zahl zw 1 und 8\n");
           printf("Table of the symmetry numbers : %s -i s\n",PROGRAMMNAME);
            printf("\n");
            exit(1);
 
 
    case 35:printf("The syntax is :\n");
            printf("\n");
            printf("%s -i[nfo] F[kq(J)] k q J\n",PROGRAMMNAME);
            printf("\n");
            exit(1);
 
    case 36:printf("The syntax is :\n");
            printf("\n");
            printf("%s -i[nfo] G[kq(J)] k q J\n",PROGRAMMNAME);
            printf("\n");
            exit(1);
 
 
    case 42:printf("Error in %s !\n",name);
            printf("It must be 2*epsilon  * D   = 5*epsilon  * D  .\n");
            printf("                 44   44            40   40      \n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 43:printf("Error in %s !\n",name);
            printf("It must be 2*epsilon  * D   = -21*epsilon  *D  .\n");
            printf("                 64   64              60  60      \n");
            printf("\n");
            fclose(fp);
            exit(1);
 
 
    case 52:printf("Error in %s !\n",name);
            printf("It must be 2*epsilon  * L   = 5*epsilon  * L  .\n");
            printf("                 44   44            40   40      \n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 53:printf("Error in %s !\n",name);
            printf("It must be 2*epsilon  * L   = -21*epsilon  *L  .\n");
            printf("                 64   64              60  60      \n");
            printf("\n");
            fclose(fp);
            exit(1);
 
    case 60:printf("Error !\n");
            printf("The start temperature must be > 0 !\n");
            printf("\n");
            exit(1);
 
    case 61:printf("Error !\n");
            printf("The end temperature must be > 0 !\n");
            printf("\n");
            exit(1);
 
    case 62:printf("Error !\n");
       printf("The start and end temperature must be different\n");
       printf("\n");
       exit(1);
 
    case 63:printf("Error !\n");
       printf("The temperature must be bigger than zero!\n");
       printf("\n");
       exit(1);
 
    case 64:printf("Error !\n");
            printf("The starting Field  must be >= 0  !\n");
            printf("\n");
            exit(1);
 
    case 65:printf("Error !\n");
            printf("The end field must be >= 0 !\n");
            printf("\n");
            exit(1);
 
    }
 
}
/*------------------------------------------------------------------------------
                               Bkq_error()
------------------------------------------------------------------------------*/
void Bkq_error(CHAR* filename,CHAR*ion,INT symmetrienr )
{
  clearscreen;
  printf("Error in %s !\n",filename);
  printf("\n");
  printf("The Ion %s is in gerneral for the ",ion);
  printf("Symmetry number %d \n",symmetrienr);
  printf("complex crystal fieldparameters .\n");
  printf("\n");
  printf("Take the Vkq-, Wkq-, Dkq- or the Lkq- Parameter !\n");
  printf("\n");
  exit(1);
}
/*------------------------------------------------------------------------------
                               Bkq_tip()
------------------------------------------------------------------------------*/
void Bkq_tip(CHAR* ion,INT symmetrienr )
 {
  clearscreen;
  printf("The Ion %s is in gerneral for the ",ion);
  printf("Symmetry number %d \n",symmetrienr);
  printf("only real crystal field parameter .\n");
  printf("\n");
  if( symmetrienr == 8 ){
     printf("You can also use the x, W parameters !\n");
     printf("\n");
     return;
  }
  printf("You can also use the Akq- or Bkq- Parameter !\n");
  printf("\n");
}
/*------------------------------------------------------------------------------
                               i_error()
------------------------------------------------------------------------------*/
void i_error() /* Eingabefehler beim Versuch von  -i */
{
    clearscreen;
    printf("\n");
    printf("Info of the Eigenvalue problem din :\n")  ;
    printf("\n");
    printf("%s -i[nfo] ew[problem] s[ingleion] \n",PROGRAMMNAME);
    printf("%s -i[nfo] ew[problem] A[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -i[nfo] ew[problem] B[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -i[nfo] ew[problem] D[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -i[nfo] ew[problem] L[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -i[nfo] ew[problem] V[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -i[nfo] ew[problem] W[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -i[nfo] ew[problem] x[W]                    \n",PROGRAMMNAME);
    printf("%s -i[nfo] ew[problem] ",PROGRAMMNAME);
    printf("s[ingleion] filename.type symmetrienr\n");
    printf("\n");
 
    exit(1);
}
/*------------------------------------------------------------------------------
                               r_error()
------------------------------------------------------------------------------*/
void r_error(CHAR c) /* Eingabefehler beim Versuch von  -r,-s */
{
 CHAR *s;
    clearscreen;
 
if(c=='r'||c=='R'){
    printf("\n");
    printf("The commands to read in a file are :\n")  ;
    printf("\n");
    printf("%s -r[ead] -s[ingleion] \n",PROGRAMMNAME);
    printf("%s -r[ead] -A[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -r[ead] -B[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -r[ead] -D[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -r[ead] -L[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -r[ead] -V[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -r[ead] -W[kq] [symmetrienummer] \n",PROGRAMMNAME);
    printf("%s -r[ead] -x[W ]                   \n",PROGRAMMNAME);
    printf("\n");
}
 
if(c=='o'||c=='O'){
    printf("\n");
    printf("The commmands for a rhombic \n")  ;
    printf("crystal field fit  ( Symmetrienummer = 2 ) are :\n")  ;
    printf("\n");
    printf("%s -o[rtho] -s[ingleion] \n",PROGRAMMNAME);
    printf("%s -o[rtho] -A[kq]  \n",PROGRAMMNAME);
    printf("%s -o[rtho] -B[kq]  \n",PROGRAMMNAME);
    printf("%s -o[rtho] -D[kq]  \n",PROGRAMMNAME);
    printf("%s -o[rtho] -L[kq]  \n",PROGRAMMNAME);
    printf("%s -o[rtho] -V[kq]  \n",PROGRAMMNAME);
    printf("%s -o[rtho] -W[kq]  \n",PROGRAMMNAME);
    printf("%s -o[rtho] -x[W ]  \n",PROGRAMMNAME);
    printf("\n");
}
 
if(c=='b'||c=='B'){
    printf("\n");
    printf("The commands for\n")  ;
    printf("the crystal field fit are :\n")  ;
    printf("\n");
    printf("%s -b[eside] -s[ingleion] \n",PROGRAMMNAME);
    printf("%s -b[eside] -A[kq]  \n",PROGRAMMNAME);
    printf("%s -b[eside] -B[kq]  \n",PROGRAMMNAME);
    printf("%s -b[eside] -D[kq]  \n",PROGRAMMNAME);
    printf("%s -b[eside] -L[kq]  \n",PROGRAMMNAME);
    printf("%s -b[eside] -V[kq]  \n",PROGRAMMNAME);
    printf("%s -b[eside] -W[kq]  \n",PROGRAMMNAME);
    printf("%s -b[eside] -x[W ]  \n",PROGRAMMNAME);
    printf("\n");
}
if(c=='s'||c=='S'){
    s = PROGRAMMNAME;
    printf("\n");
  printf("commmands to calculate the inverse susceptability-\n");
  printf("as a function of temperature (in Kelvin) and in\n");
  printf("small external fields. The Temperature interval  goes from \n");
  printf("anf_temp bis end_temp.                                           \n");
  printf("\n");
  printf("%s -s[uszept] -s[ingleion] \n",PROGRAMMNAME);
  printf("%s -s[uszept] -A[kq] anf_temp end_temp lambda [theta]||[filename.type]\n",s);
  printf("%s -s[uszept] -B[kq] anf_temp end_temp lambda [theta]||[filename.type]\n",s);
  printf("%s -s[uszept] -D[kq] anf_temp end_temp lambda [theta]||[filename.type]\n",s);
  printf("%s -s[uszept] -L[kq] anf_temp end_temp lambda [theta]||[filename.type]\n",s);
  printf("%s -s[uszept] -V[kq] anf_temp end_temp lambda [theta]||[filename.type]\n",s);
  printf("%s -s[uszept] -W[kq] anf_temp end_temp lambda [theta]||[filename.type]\n",s);
  printf("%s -s[uszept] -x[W ] end_temp lambda [theta]||[filename.type]\n",s);
  printf("\n");
  printf("\n");
  printf("lambda is the molecular field constant (1/X=1/X_cf-lambda)\n");
  printf("theta  is the Curie temperature with T>theta 1/X proportional T+theta)\n");
  printf("In <filename.type> it is possible to give a temperature\n");
  printf("dependant theta. (In RTPLOT-convention).\n");
}
 
if(c=='m'||c=='M'){
    s = PROGRAMMNAME;
    printf("\n");
  printf("commands for the calculation of the RE magnetic Moment\n");
  printf("at a given Temperature (in Kelvin) in a field  \n");
  printf("range : anf_feld until end_feld (in Tesla).                     \n");
  printf("\n");
  printf("%s -m[oment] -s[ingleion] \n",PROGRAMMNAME);
  printf("%s -m[oment] -A[kq] [symmetrienr] anf_feld end_feld temp\n",s);
  printf("%s -m[oment] -B[kq] [symmetrienr] anf_feld end_feld temp\n",s);
  printf("%s -m[oment] -D[kq] [symmetrienr] anf_feld end_feld temp\n",s);
  printf("%s -m[oment] -L[kq] [symmetrienr] anf_feld end_feld temp\n",s);
  printf("%s -m[oment] -V[kq] [symmetrienr] anf_feld end_feld temp\n",s);
  printf("%s -m[oment] -W[kq] [symmetrienr] anf_feld end_feld temp\n",s);
  printf("%s -m[oment] -x[W ] anf_feld end_feld temp\n",s);
  printf("\n");
  printf("\n");
}
 
if(c=='k'||c=='K'){
    s = PROGRAMMNAME;
    printf("\n");
  printf("commands for the calculation of a cubic average  \n");
  printf("of the magnetic moment of the RE ion at a given \n");
  printf("temperature (in Kelvin) in a field range of    \n");
  printf("anf_feld until end_feld(in Tesla) :                               \n");
  printf("\n");
  printf("%s -k[poly ] -x[W ] anf_feld end_feld temp\n",s);
  printf("\n");
  printf("\n");
}
 
    exit(1);
}
/*------------------------------------------------------------------------------
                              isinlimits()
------------------------------------------------------------------------------*/
void isinlimits(  FILE *fp,  CHAR *filename,INT nrion,DOUBLE x1,DOUBLE x2,DOUBLE x3,CHAR modus) /* x1,x2,x3 in seinen Grenzen ? */
{
    INT i=0;
 
 
    switch(modus){
         case 's': if( x1< 0.0             ) i=coor_error(fp,filename,nrion ,1);
                   if( x2< 0.0 || x2 >360.0) i=coor_error(fp,filename,nrion ,2);
                   if( x3< 0.0 || x3 >180.0) i=coor_error(fp,filename,nrion ,3);
                   break;
 
         case 'p': if( x1< 0.0             ) i=coor_error(fp,filename,nrion ,1);
                   if( x2< 0.0 || x2 >360.0) i=coor_error(fp,filename,nrion ,3);
                   break;
    }
    if( i!=0 ) { fclose(fp);exit(1);}
}
/*------------------------------------------------------------------------------
                              coor_error()
------------------------------------------------------------------------------*/
INT coor_error(FILE *fp,CHAR *filename,INT nrion,INT errornr)
{
    UNUSED_PARAMETER(fp);
 
    printf("\n");
    printf("mistake in File %s !\n",filename);
 
    switch(errornr){
         case 1 :  switch(nrion){
                        case 0 :  printf("The Magnetic field ");
                                  printf("is negative !\n");
                                  break;
                        default:  printf("The  %d-ten ",nrion);
                                  printf("surrrounding ion is negative !\n");
                                  break;
                   }
                   break;
 
         case 2 :  switch(nrion){
                        case 0 :  printf("The Polar angle of the Magnetic");
                                  printf("field is outside the ");
                                  printf("allowable area !\n");
                                  break;
                        default:  printf("The Polar angle of the %d-ten ",nrion);
                                  printf("surrounding ion is outside the\n");
                                  printf("allowed area !\n");
                                  break;
                   }
                   break;
 
         case 3 :  switch(nrion){
                        case 0 :  printf("The Azimutal angle of the Magnetic");
                                  printf("field is outide the ");
                                  printf("allowed area !\n");
                                  break;
                        default:  printf("The Azimutal angle of the %d-ten ",nrion);
                                  printf("surrounding ion is outside the\n");
                                  printf("allowed area !\n");
                                  break;
                   }
                   break;
 
         case 4 : printf("The %d-te surrounding ion sits ",nrion);
                  printf("in the origin .\n please select a ");
                  printf("different place !\n");
                  break;
 
 
 
    }
    return(1);
}
/*------------------------------------------------------------------------------
                              isimplementiert()
------------------------------------------------------------------------------*/
INT isimplementiert(CHAR *ion)/* fragt,ob ion implementiert ,falls nein exit*/
             /* falls ja : nummer des implementierten ions */
{
    INT i;
 
    for( i=0 ; i< ANZ_IONEN ; ++i )
              if(  equal( IONENIMP[i].ionname ,ion)  ) return( i );
 
    clearscreen;
    printf("Error !\n");
    printf("Ion %s is not implemented !\n",ion);
    printf("\n");
    printf("Implemented ions are :\n");
    for( i=0 ; i< ANZ_IONEN ; ++i )
       printf("%s\n",IONENIMP[i].ionname);
    printf("\n");
 
    exit(1);
}
/*------------------------------------------------------------------------------
                                  a_toi()
------------------------------------------------------------------------------*/
INT a_toi(CHAR *string,INT anfang,INT ende) /* Integer aus String bestimmen */
{
    CHAR *buffer,*a_tos(CHAR*,INT,INT);
    INT  zahl;
 
    buffer = a_tos( string,anfang,ende);
    zahl   = atoi( buffer );
    free_( buffer );
 
    return( zahl );
}
/*------------------------------------------------------------------------------
                                  a_tof()
------------------------------------------------------------------------------*/
DOUBLE a_tof( CHAR *string,INT anfang,INT ende) /* Reelle Zahl im String bestimmen */
{
    CHAR   *buffer,*a_tos(CHAR*,INT,INT);
    DOUBLE zahl;
 
 
    buffer = a_tos( string,anfang,ende);
    zahl   = atof( buffer );
    free_( buffer );
 
    return( zahl );
}
/*------------------------------------------------------------------------------
                                  a_tos()
------------------------------------------------------------------------------*/
CHAR  *a_tos(CHAR *string,INT anfang,INT ende)
{    INT i;
    CHAR   *buffer;
 
    buffer = STRING_ALLOC(ende-anfang+2);
    for( i=anfang ; i<=ende ; ++i )
         VALUE(buffer,i-anfang) = VALUE(string,i);
    VALUE(buffer,ende-anfang+2) = '\0';
 
    return( buffer );
}
/*------------------------------------------------------------------------------
                                  equal()
------------------------------------------------------------------------------*/
INT equal( CHAR *s ,CHAR *t ) /* enthaelt String t den String s ?*/
{
    INT  i,len;
    CHAR up(CHAR);
 
 
    len = strlen_own(s);
    if( strlen_own(t)  < strlen_own(s) )  return(0);
 
 
    for( i=0; i<len ; ++i){
         if(!(*(s+i)==*(t+i)||*(s+i)==up(*(t+i))||*(t+i)==up(*(s+i)) ))
              return(0);
    }
    return(1);
}
/*------------------------------------------------------------------------------
                                    up()
------------------------------------------------------------------------------*/
CHAR up(CHAR c )
{
    if(isupper(c) || isdigit(c) || c=='+'|| c=='-' )  return(c);
 
    return( c + 'A' - 'a' );
}
/*------------------------------------------------------------------------------
ENDEMODUL    C F I E L D    C
------------------------------------------------------------------------------*/
