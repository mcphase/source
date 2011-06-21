/*-----------------------------------------------------------------------------
 
                                M   O  D  U  L
 
                                  DIAHERMX.C
 
-------------------------------------------------------------------------------
 
Aufgabe               :  Diagonalisieren einer hermiteschen Matrix
                         ---                   ---          -    -
-------------------------------------------------------------------------------
 
Benutzte Literatur    :  J.H Wilkinson   C.Reinsch   "Linear Algebra",
                         Springer Verlag  Berlin - Heidelberg - New York
                         Man findet hierin eine kurze Beschreibung der
                         theoretischen Grundlagen sowie ein kommentiertes
                         ALGOL-Programm ,welches nach C konvertiert wurde.
 
                         Dieses Buch liegt in der kleinen Bibliothek des
                         IFF's der KFA Juelich aus.
 
-------------------------------------------------------------------------------
 
Definierte Funktionen :
 
-----------------------
 
Die nachfolgenden Funktionen wurden aus ALGOL konvertiert :
 
comhes()   :  komplexe quatratische Matrix auf Hessenbergform transformieren
comlr2()   :  transformierte Matrix diagonalisieren
combak()   :  die so diagonalisierte Matrix zuruecktransformieren
 
hermhes()  :  hermitesche Matrix auf Hessenbergform transformieren
hermlr2()  :  transformierte Matrix diagonalisieren
 
diagonalisiere() : hermhes() + hermlr2()
info_ewproblem() : Informationen ueber die Diagonalisierung ausgeben
 
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
Extern definierte Funktionen
-----------------------------------------------------------------------------*/
extern MATRIX  *mx_copy(); /* Matrix kopieren */
extern MATRIX  *mx_mult(); /* Matrizen multiplizieren */
extern MATRIX  *mx_ewev();
extern MATRIX  *mx_sub() ; /* Matrizen subtrahieren   */
extern MATRIX  *mx_alloc();/* Speicher fuer Matrix holen */
extern KOMPLEX *cdiv();    /* komplexe Zahlen dividieren */
extern KOMPLEX *csqroot();   /* Quadratwurzel komplexer Zahlen */
extern KOMPLEX *vskalar(); /*         <a|b>                  */
extern DOUBLE  cnorm();    /* Euklidnorm eines komplexen Vektors */
 
extern VEKTOR  *vr_normalisieren(); /* Vektor  normalisieren */
extern VEKTOR  *cv_mult();          /* |b> = c * |a>        */
extern VEKTOR  *vr_copy();          /* |b> =     |a>        */
extern VEKTOR  *_vr_copy();         /* |c> =     |a>        */
extern VEKTOR  *vr_sub();           /* |a> = |a>- |b>       */
 
extern VEKTOR  *vr_alloc();/* Speicher fuer komplexen Vektor holen */
extern DOUBLE  sqrt();     /* Quadratwurzel  */
extern DOUBLE  pow_();
extern INT     a_toi();
extern DOUBLE  a_tof();
extern INT     free_vr();  /* Vektor-Speicher wieder freigeben */
extern DOUBLE  log();      /* Logarithmus               */

extern FILE *fopen_errchk();         /* definiert in EINGABE.C*/
 
INT test_nullmatrix();
void write_title();
INT gi_entartung();
/*----------------------------------------------------------------------------
                               diagonalisiere()
------------------------------------------------------------------------------*/
/* hermitesche Matrix diagonalisieren   */
EWPROBLEM  *diagonalisiere(ewproblem,matrix,overwrite,setup)
    EWPROBLEM *ewproblem;
    MATRIX    *matrix;
    INT       overwrite;
    SETUP     *setup;
{
    DOUBLE    macheps,accuracy();
    COMHES    *hessenbergform , *hermhes();
    EWPROBLEM *dia_hessenberg , *hermlr2();
    EWPROBLEM *number_ev();
    EWPROBLEM *orthonormalisieren();
    EWPROBLEM *entartung();
    EWPROBLEM *ordnen_ew();
    EWPROBLEM *set_ewproblem();
 
/*[1]*/  macheps = (MACHEPSFACT*accuracy());
 
/*[0]*/  if(test_nullmatrix(matrix)==NEIN){
/*[2]*/     hessenbergform = hermhes( ewproblem, matrix ,overwrite);
/*[3]*/     dia_hessenberg = hermlr2( ewproblem,hessenbergform,macheps );
         }
         else{
          ewproblem = set_ewproblem( setup,ewproblem,matrix,macheps,
                                     overwrite );
          dia_hessenberg = ewproblem;
         }
 
         if( setup != (SETUP*)0 )
              dia_hessenberg->eps_setup       = MACHEPS(setup);
         dia_hessenberg->eps_machine          = macheps;
/*[4]*/  ewproblem        = number_ev( dia_hessenberg );
/*[5]*/  ewproblem        = ordnen_ew( ewproblem,overwrite );
/*[6]*/  ewproblem        = entartung( ewproblem,overwrite );
/*[7]*/  ewproblem        = orthonormalisieren( ewproblem ); 
         return( ewproblem );
 
/* [0] Rechnergenauigkeit bestimmen                                 */
/* [1] Festellen ob Matrix Nullmatrix ist                           */
/* [2] Matrix auf Hessenbergform transformieren                     */
/* [3] die auf Hessenbergform transformierte Matrix diagonalisieren */
/* [4] die Zahl der Eigenvektoren bestimmen                         */
/* [5] die EW'e der groesse nach ordnen : E1<=E2<= ... <= En        */
/* [6] die Entartung der Eigenwerte bestimmen und EW ordnen         */
/* [7] die Eigenvektoren orthonormalisieren                         */
}
/*----------------------------------------------------------------------------
                               set_ewproblem()
------------------------------------------------------------------------------*/
EWPROBLEM *set_ewproblem(setup,ewproblem,matrix,macheps,overwrite)
  SETUP     *setup;
  EWPROBLEM *ewproblem;
  MATRIX    *matrix;
  DOUBLE    macheps;
  INT       overwrite;
{
  INT n,i,j;
  MATRIX *mx_alloc(),*ev;
  VEKTOR *vr_alloc(),*ew;
  COMHES *comhes;
 
  n = MXDIM(matrix);
 
  if( ewproblem == (EWPROBLEM*)0 ){
      comhes              = COMHES_ALLOC(1);
      comhes -> overwrite = overwrite;
      comhes -> matrix    = matrix;
      comhes -> hesse     = matrix;
      comhes -> swap      = INT_ALLOC(n);
      comhes -> swap_dim  = n;
      ewproblem = COMLR2_ALLOC(1);
      ewproblem -> eigenvektoren   = mx_alloc(n,n);
      ewproblem -> eigenwerte      = vr_alloc(n);
      ewproblem -> iterationsteps  = 0;
      ewproblem -> iteration       = SUCCESSFUL;
      ewproblem -> eps_machine     = macheps;
      ewproblem -> eps_setup       = MACHEPS(setup);
      ewproblem -> comhes          = comhes;
  }
 
    for( i=1 ; i<=n ; ++i){      /* matrix 'matrix' auf Null setzen */
       for( j=1 ; j<=n ; ++j){
           R(matrix,i,j) = 0.0;
           I(matrix,i,j) = 0.0;
       }
    }
 
 
    ev = ewproblem -> eigenvektoren;
    for( i=1 ; i<=n ; ++i){
       for( j=1 ; j<=n ; ++j){
           R(ev,i,j) = 0.0;
           I(ev,i,j) = 0.0;
       }
       R(ev,i,i) = 1.0;
    }
 
    ew = ewproblem -> eigenwerte;
    for( i=1 ; i<=n ; ++i){
           RV(ew,i) = 0.0;
           IV(ew,i) = 0.0;
    }
 
  return( ewproblem );
 
}
/*----------------------------------------------------------------------------
                               cfield_setup()
------------------------------------------------------------------------------*/
SETUP *cfield_setup() /* Setupfile lesen falls vorhanden sonst erzeugen*/
{
   FILE  *fopen(),*fp;
   SETUP *read_setup();
   SETUP *create_setup();
   SETUP *setup;
 
   setup = SETUP_ALLOC(1);
 
   if( (fp=fopen(CFIELDSETUP,"rb")) == (FILE*)0 )
        setup = create_setup(setup);
   else setup = read_setup(fp,setup);
 
   return( setup );
}
/*----------------------------------------------------------------------------
                               read_setup()
------------------------------------------------------------------------------*/
SETUP *read_setup(fp,setup)  /* setupfile lesen    */
    FILE *fp;
   SETUP *setup;
{
   INT    buffer_size=81,stelle;
   INT    a_toi();
   CHAR   *string,*line,*fgets(),c;
   DOUBLE a_tof(),pow_(),versionsnr,dummy;
   SETUP  *create_setup();
   DOUBLE accuracy();
 
   string = STRING_ALLOC( buffer_size );
 
   while(  *(line=fgets( string , buffer_size , fp )) != '|'  );/* 1.| */
   c = VALUE(line,2);
   versionsnr = a_tof(line,11,23);
   printf("#versionsnr =%6.2f\n",versionsnr);
   if( c != 'V' || versionsnr < VERSION  ){
        setup = create_setup(setup);
        return( setup );
   }
 
   while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 1.==== */
   while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 2.==== */
   while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 3.==== */
 
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   KOMMASTELLE(setup) = a_toi(line,26,41);
   if(KOMMASTELLE(setup) < 1) KOMMASTELLE(setup) = 1;  /* da EW auf x.xx */
 
   while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 4.==== */
   while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 5.==== */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   dummy               = a_tof(line,26,41);
   NUMMERIERUNG(setup) = NEIN;
   if( dummy != 0.0 ) NUMMERIERUNG(setup)=JA;
 
 
   while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 4.==== */
   while(  *(line=fgets( string , buffer_size , fp )) != '='  );/* 5.==== */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   line=fgets( string , buffer_size , fp ); /*                     */
   ANZ_DATENPUNKTE(setup) = a_toi(line,26,41);
   if( ANZ_DATENPUNKTE(setup) < 2) ANZ_DATENPUNKTE(setup)=2;
 
   stelle = (INT)(-log(accuracy())/log(10.0)+1.0 );
   if( KOMMASTELLE(setup) > stelle || KOMMASTELLE(setup)<0 )
          setup = create_setup(setup);
 
   MACHEPS(setup) = 1/pow_(10.0, KOMMASTELLE(setup) );
 
   return(setup);
}
/*----------------------------------------------------------------------------
                               create_setup()
------------------------------------------------------------------------------*/
SETUP *create_setup(setup)  /* setupfile erzeugen */
   SETUP *setup;
{
   FILE *fopen(),*fp;
   CHAR *t01,*t02,*t03,*t04,*t05,*t06,*t07,*t08,*t09,*t10;
   CHAR *t11,*t12,*t13;/*,*t14,*t15,*t16,*t17,*t18,*t19,*t20;*/
   DOUBLE accuracy(),log();
   INT stelle;
 
   fp=fopen_errchk(CFIELDSETUP,"w");
   write_title(fp);
 
 
t01=" ______________________________________________________________ \n";
t02="|                  C F I E L D / S O 1 I O N                   |\n";
t03="|                        SETUP - FILE                          |\n";
t04="|                                                              |\n";
t05=" ______________________________________________________________ \n";
t06="\n";
fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);
fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);
 
t07="================================================================\n";
t08="|                                                              |\n";
t09="| Computational accuracy  =  %5.0e                             |\n";
t10="|                                                              |\n";
t11="================================================================\n";
t12="\n";
fprintf(fp,"%s",t07);fprintf(fp,"%s",t08);fprintf(fp,t09,accuracy());
fprintf(fp,"%s",t10);fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);
 
t01="================================================================\n";
t02="|                                                              |\n";
t03="| the inacuracy of the calculation  is %3d. decimal places     |\n";
t04="| This agrees when two numbers                                 |\n";
t05="| are the same for the computer. Because of rounding erro      |\n";
t06="| This decimal place should be put too:                        |\n";
t07="|                                                              |\n";
t08="|  Recommended          ==================                     |\n";
t09="|  decimal place:       |       3        |                     |\n";
t10="|                       ==================                     |\n";
t11="|                                                              |\n";
t12="================================================================\n";
t13="\n";
 
KOMMASTELLE(setup) = 3;
stelle = (INT)(-log(accuracy())/log(10.0)+1.0 );
fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,t03,stelle);
fprintf(fp,"%s",t04);fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);
fprintf(fp,"%s",t08);fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);fprintf(fp,"%s",t13);
 
t01="================================================================\n";
t02="| The program outputs record into file which are NOT           |\n";
t03="| numbered in the files.                                       |\n";
t04="| If you would like the program to number the record in a file |\n";
t05="| The put a different number in the box other than zero        |\n";
t06="|                                                              |\n";
t07="|                                                              |\n";
t08="|                       ==================                     |\n";
t09="|                       |       0        |                     |\n";
t10="|                       ==================                     |\n";
t11="|                                                              |\n";
t12="================================================================\n";
t13="\n";
 
fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);
fprintf(fp,"%s",t04);fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);
fprintf(fp,"%s",t08);fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);fprintf(fp,"%s",t13);
ANZ_DATENPUNKTE(setup) = 20;
 
t01="================================================================\n";
t02="|                                                              |\n";
t03="| The program writes to output fiels                           |\n";
t04="| Please give the                                              |\n";
t05="| number of datapoints that are supposed to be calculated      |\n";
t06="| in the lower box.                                            |\n";
t07="|                                                              |\n";
t08="|                       ==================                     |\n";
t09="|                       |       20       |                     |\n";
t10="|                       ==================                     |\n";
t11="|                                                              |\n";
t12="================================================================\n";
t13="\n";
 
fprintf(fp,"%s",t01);fprintf(fp,"%s",t02);fprintf(fp,"%s",t03);
fprintf(fp,"%s",t04);fprintf(fp,"%s",t05);fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);
fprintf(fp,"%s",t08);fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);
fprintf(fp,"%s",t11);fprintf(fp,"%s",t12);fprintf(fp,"%s",t13);
ANZ_DATENPUNKTE(setup) = 20;
fclose(fp);
return(setup);
}
/*----------------------------------------------------------------------------
                               write_title()
------------------------------------------------------------------------------*/
void write_title(fp)  /* setupfile erzeugen */
   FILE *fp;
{
   CHAR *t01,*t02,*t03,*t04,*t05/*,*t06,*t07*/,*t08/*,*t09,*t10*/;
   CHAR/* *t11,*/*t12,*t13,*t14/*,*t15*/,*t16;/*,*t17,*t18,*t19,*t20;*/
 
t01=" -----------------------\n";
t02="| VERSION : %6.2f      |\n";
t03=" -----------------------\n";
t04="\n";
/*fprintf(fp,"%s",t01);*/
fprintf(fp,t02,VERSION);/*fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);*/
 
t01=" -------------------------------------------------------------- \n";
t02="|                                                              |\n";
t03="|                  C F I E L D / S O 1 I O N                   |\n";
t04="|                                                              |\n";
t05="|                    A Crystal field program                   |\n"; /*
t06="|                                                              |\n";
t07=" -------------------------------------------------------------- \n"; */
t08="                __________________________________              \n"; /*
t09="               |                                  |             \n";
t10="               |          Programmautor:          |             \n";
t11="               |                                  |             \n"; */
t12="               |         Peter  Hoffmann          |             \n";
t13="               |    Forschungszentrum Juelich     |             \n";
t14="               |Institut fuer Festkoerperforschung|             \n"; /*
t15="               |Tel.: 02461-616896                |             \n"; */
t16="                __________________________________              \n";
 
 
fprintf(fp,"%s",t01);/*fprintf(fp,"%s",t02);*/fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);fprintf(fp,"%s",t05);
/*fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);*/fprintf(fp,"%s",t08);/*fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);*/
/*fprintf(fp,"%s",t11);*/fprintf(fp,"%s",t12);fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);/*fprintf(fp,"%s",t15);*/
fprintf(fp,"%s",t16);
 
 
}
/*----------------------------------------------------------------------------
                               write_titlecom()
------------------------------------------------------------------------------*/
void write_titlecom(fp)  /* setupfile erzeugen */
   FILE *fp;
{
   CHAR *t01/*,*t02*/,*t03/*,*t04*/,*t05/*,*t06,*t07*/,*t08/*,*t09,*t10*/;
   CHAR/* *t11,*/*t12,*t13,*t14/*,*t15*/,*t16;/*,*t17,*t18,*t19,*t20;

t01="#{-----------------------\n";
t02="#{VERSION : %6.2f      |\n";
t03="#-----------------------\n";
t04="#\n"; */
/*fprintf(fp,"%s",t01);*/
/*fprintf(fp,t02,VERSION);fprintf(fp,"%s",t03);fprintf(fp,"%s",t04);*/
 
t01="#{------------------------------------------------------------- \n"; /*
t02="#                                                              |\n"; */
t03="#                  C F I E L D / S O 1 I O N  %6.2f           |\n";  /*
t04="#                                                              |\n"; */
t05="#                    A crystal field program                   |\n"; /*
t06="#                                                              |\n";
t07="#-------------------------------------------------------------- \n"; */
t08="#               __________________________________              \n"; /*
t09="#              |                                  |             \n";
t10="#              |          Programmautor:          |             \n";
t11="#              |                                  |             \n"; */
t12="#              |         Peter  Hoffmann          |             \n";
t13="#              |    Forschungszentrum Juelich     |             \n";
t14="#              |Institut fuer Festkoerperforschung|             \n"; /*
t15="#              |Tel.: 02461-616896                |             \n"; */
t16="#               __________________________________              \n";
 
 
fprintf(fp,"%s",t01);/*fprintf(fp,"%s",t02);*/fprintf(fp,t03,VERSION);/*fprintf(fp,"%s",t04);*/fprintf(fp,"%s",t05);
/*fprintf(fp,"%s",t06);fprintf(fp,"%s",t07);*/fprintf(fp,"%s",t08);/*fprintf(fp,"%s",t09);fprintf(fp,"%s",t10);*/
/*fprintf(fp,"%s",t11);*/fprintf(fp,"%s",t12);fprintf(fp,"%s",t13);fprintf(fp,"%s",t14);/*fprintf(fp,"%s",t15);*/
fprintf(fp,"%s",t16);
 
 
}

/*----------------------------------------------------------------------------
                               is_equal()
------------------------------------------------------------------------------*/
INT is_equal(a,b,macheps) /* ist |a-b|< macheps ? */
   DOUBLE a,b,macheps;
{
   if( ABSD( a-b ) < macheps ) return( JA );
 
   return( NEIN );
}
/*----------------------------------------------------------------------------
                               is_null()
------------------------------------------------------------------------------*/
DOUBLE is_null(a,macheps) /* ist |a|< macheps ? */
   DOUBLE a,macheps;
{
   if( is_equal(a,0.0,macheps ) == JA  )  return( 0.0 );
 
   return( a );
}
/*----------------------------------------------------------------------------
                                 null()
------------------------------------------------------------------------------*/
INT    null(a,macheps) /* ist |a|< macheps ? */
   DOUBLE a,macheps;
{
 
   return( is_equal(a,0.0,macheps )  );
 
}
/*----------------------------------------------------------------------------
                               null_mx_error()
------------------------------------------------------------------------------*/
void null_mx_error()
{
  clearscreen;
  printf("\nThe total Hamiltonian is zero .\n");
  printf("therefore no observed splitting of the energy levels.\n\n");
 
  exit(1);
}
/*----------------------------------------------------------------------------
                               test_nullmatrix();
------------------------------------------------------------------------------*/
INT test_nullmatrix(matrix)  /* Ist die Matrix matrix eine Nullmatrix? */
     MATRIX *matrix;
{
 
   INT zeile,spalte,anz_nullen;
   anz_nullen = 0;
   for( zeile=1;zeile<=MXDIM(matrix);++zeile)
       for( spalte=1;spalte<=MXDIM(matrix);++spalte){
 
            if( R(matrix,zeile,spalte) == 0.0 ) ++anz_nullen;
            if( I(matrix,zeile,spalte) == 0.0 ) ++anz_nullen;
   }
 
   if( anz_nullen == 2*MXDIM(matrix)*MXDIM(matrix) ) return( JA );
 
   return( NEIN );
}
/*----------------------------------------------------------------------------
                                   accuracy()
------------------------------------------------------------------------------*/
DOUBLE accuracy() /* Rechnergenauigkeit bestimmen */
{
    DOUBLE work();
    DOUBLE accuracy;
    DOUBLE digit = 9.0;
 
    accuracy = work(1.0);/* Zehnerpotenz der Rechnergenauigkeit bestimmen */
 
    do{
        if( 1.0+digit*accuracy == 1.0 ) return( digit*accuracy );
        --digit;
    }while( digit >= 0 );
 
    return( accuracy );
}
/*                                -------------
                                     work()
                                  -------------                               */
DOUBLE work(zahl)
    DOUBLE zahl;
{
    return( (1.0+1/zahl==1.0) ? 1/zahl : work(10*zahl) );
}
/*----------------------------------------------------------------------------
                                   hermhes()
------------------------------------------------------------------------------*/
COMHES *hermhes(ewproblem,a,overwrite)/*transformiert die hermitesche */
                                      /* Matrix a auf Hessenbergform  */
    EWPROBLEM *ewproblem;
    MATRIX *a;         /* Pointer auf die Matrix a */
    INT    overwrite;  /* soll a ueberschrieben werden ?    */
{
    INT    low,upp;
    COMHES *comhes();
 
    low = 1;       /* bei hermiteschen Matrizen ist kein "balance" noetig*/
    upp = MXDIM(a);/* daher sind low und upp so zu setzen */
    return( comhes(ewproblem,low,upp,a,overwrite) );
}
/*----------------------------------------------------------------------------
                                   hermlr2()
------------------------------------------------------------------------------*/
                                 /* diagonalisert HESSENBERGFORM einer       */
                                 /* hermiteschen Matrix innerhalb der        */
                                 /* Rechnergenauigkeit macheps               */
COMLR2 *hermlr2(ewproblem,hermhes,macheps)
    EWPROBLEM *ewproblem;
    COMHES *hermhes;
    DOUBLE macheps;
{
    INT    low,upp;
    COMLR2 *comlr2();
 
    low = 1 ;
    upp = MXDIM(hermhes->hesse);
 
    return( comlr2(ewproblem,low,upp,hermhes,macheps ) );
}
/*----------------------------------------------------------------------------
                                  number_ev()
------------------------------------------------------------------------------*/
EWPROBLEM *number_ev(ewproblem)   /* Anzahl der Eigenvektoren bestimmen */
    EWPROBLEM *ewproblem;
{
    INT    nullvektor=0;
    INT    anz_sp,sp,ze;
    VEKTOR *v;
    MATRIX *m;
    COMHES *comhes;
 
    m  = ewproblem->eigenvektoren;
    sp = anz_sp = ANZ_SP(m);
 
    comhes = ewproblem -> comhes;
    if(test_nullmatrix(comhes->matrix)==JA){
        ANZ_SP(ewproblem->eigenvektoren) = ANZ_SP(comhes->matrix);
        return( ewproblem );
    }
 
    do{
         v = MXSP(m,sp);
         for( ze=1 ; ze<= anz_sp ; ++ze )
              if( RV(v,ze)!=0.0  || IV(v,ze)!=0.0 )
                        nullvektor = 1;
         --sp;
 
    }while(!nullvektor);
 
    ANZ_SP(ewproblem->eigenvektoren) = sp+1;
 
    return(ewproblem);
}
/*----------------------------------------------------------------------------
                               ordnen_ew()
------------------------------------------------------------------------------*/
/* Energieeigenwerte aufsteigend sortieren */
/* nach Shell's Algorithmus */
EWPROBLEM *ordnen_ew( ewproblem,overwrite )
      EWPROBLEM *ewproblem;
      INT       overwrite;
{
   INT *nummer,level,i,*sort();
   VEKTOR *ew;
   DOUBLE *werte;
   COMHES *comhes;
 
 
   ew     = ewproblem->eigenwerte;
 
   if( overwrite == NOSPACE ) nummer = ewproblem->energie_nummer;
   else                       nummer = INT_ALLOC( VRDIM(ew) );
 
   for(level=1; level<=VRDIM(ew) ; ++level )
          VALUE(nummer,level) = level;
 
   comhes = ewproblem -> comhes;
   if(test_nullmatrix(comhes->matrix)==JA){
      ewproblem->energie_nummer = nummer;
      return( ewproblem );
   }
 
   werte = DOUBLE_ALLOC(VRDIM(ew));
   for(i=1;i<=VRDIM(ew);++i)
      VALUE(werte,i) = RV(ew,i);
   nummer = sort(werte,nummer,VRDIM(ew));
   free_(werte);
 
   ewproblem->energie_nummer = nummer;
   return( ewproblem );
}
/*----------------------------------------------------------------------------
                               sort()
------------------------------------------------------------------------------*/
INT *sort(werte,nummer,anz)
   DOUBLE *werte;
   INT    *nummer,anz;
{
   INT level,i,j,gap,temp;
 
   for(level=1; level<=anz ; ++level )
          VALUE(nummer,level) = level;
 
   for( gap=anz/2 ; gap>0 ; gap /= 2 )
    for( i=gap ; i<=anz ; ++i )
     for(j=i-gap;j>0 && VALUE(werte,VALUE(nummer,j))
                       >VALUE(werte,VALUE(nummer,j+gap)); j-=gap)
     {
       temp                = VALUE(nummer,j);
       VALUE(nummer,j)     = VALUE(nummer,j+gap);
       VALUE(nummer,j+gap) = temp;
     }
 
   return(nummer);
}
/*----------------------------------------------------------------------------
                               entartung()
------------------------------------------------------------------------------*/
/*
  entartung ist am Anfang eine VRDIM(ew) X VRDIM(ew) Nullmatrix.
  In einer Zeile der Matrix entartung stehen spaeter die Energienummern
  zu gleichen Energieniveaus.
 
*/
/* Entartung der Eigenwerte bestimmen */
/* Energienievaus shiften ,sodass kleinster*/
/* Eigenwert Null ist */
 
EWPROBLEM *entartung( ewproblem,overwrite )
    EWPROBLEM *ewproblem;
    INT       overwrite;
{
    MATRIX    *entartung,*eintragen();
    VEKTOR    *ew;
    INT       *nummer,*gi,level,i,zeile,spalte;
    DOUBLE    eps,shift;
 
    ew        =  ewproblem->eigenwerte;
    nummer    =  ewproblem->energie_nummer;
    eps       =  ewproblem->eps_setup;
 
   if( overwrite == NOSPACE ) entartung = ewproblem->entartung;
   else                       entartung =  mx_alloc( VRDIM(ew),VRDIM(ew));
 
   for( zeile=1;zeile<=VRDIM(ew);++zeile)
      for( spalte=1;spalte<=VRDIM(ew);++spalte){
          R(entartung,zeile,spalte) = 0.0;
          I(entartung,zeile,spalte) = 0.0;
      }
 
 
    zeile = 1; /* verschiedene Energieniveaus bestimmen*/
    for( level=1 ; level<=VRDIM(ew) ; ++level )
     if(level==1) entartung = eintragen(entartung ,zeile,VALUE(nummer,level));
     else{
       if(is_equal(RV(ew,VALUE(nummer,level)),RV(ew,VALUE(nummer,level-1)),eps))
            entartung = eintragen( entartung ,zeile,VALUE(nummer,level));
       else
            entartung = eintragen( entartung ,++zeile,VALUE(nummer,level));
           }
 
     /* die Anzahl der verschiedenen Energieniveaus merken */
     ANZ_ZE(entartung) = zeile;
 
    /* Entartungen der verschiedenen Energieniveaus merken */
 
   if( overwrite == NOSPACE )  gi = ewproblem->gi;
   else                        gi = INT_ALLOC( ANZ_ZE(entartung) );
 
    for( i=1 ; i<=ANZ_ZE(entartung) ; ++i )
          VALUE(gi,i) = gi_entartung( entartung,i, VRDIM(ew) );
 
     /* energieniveaus shiften und shift merken*/
    shift = RV(ew,(INT)R(entartung,1,1));
    for( level=1 ; level<=VRDIM(ew) ; ++level )
       RV(ew,level) -= shift;
 
    ewproblem->gi        = gi;
    ewproblem->shift     = shift;
    ewproblem->entartung = entartung;
 
 
    return( ewproblem );
 
}
/*----------------------------------------------------------------------------
                               gi_entartung()
------------------------------------------------------------------------------*/
INT gi_entartung( entartung ,zeile, anz_sp) /* Entartung eines Energieniveaus*/
   MATRIX *entartung;                       /* gleich Anzahl der von Null    */
   INT     zeile; /* von Matrix entartung *//* verschiedenen Eintraege in der*/
   INT    anz_sp; /* von Matrix entartung *//* Zeile von entartung */
{
  INT spalte;
 
  for( spalte=1 ; spalte <anz_sp ; ++spalte )
      if( R(entartung,zeile,spalte) == 0.0 )
                return( --spalte );
 
 if( R(entartung,zeile,anz_sp) == 0.0 ) return(anz_sp-1);
 else                                   return(anz_sp);
}
/*----------------------------------------------------------------------------
                               eintragen()
------------------------------------------------------------------------------*/
MATRIX *eintragen( matrix, zeile ,energieniveau_nr ) /* Energieniveau_nr in  */
     MATRIX *matrix;                                 /* die Zeile der Matrix */
     INT    zeile,energieniveau_nr;                  /* matrix eintragen und */
{                                                    /* zwar in die Spalte   */
     INT spalte;                                     /* wo noch keine Nummer */
                                                     /* steht.               */
     for( spalte=1;spalte<=MXDIM(matrix);++spalte)
         if( R(matrix,zeile,spalte)==0.0 ){
             R(matrix,zeile,spalte)=(DOUBLE)energieniveau_nr;
             return(matrix);
         }
 
    printf("\n");
    printf("Error in Function : eintragen() in DIAHERMX.C !\n");
    printf("Energy level cannot register degeneracy in matrix.\n");
    printf("\n");
    exit(0);
}
/*----------------------------------------------------------------------------
                          orthonormalisieren()
------------------------------------------------------------------------------*/
EWPROBLEM *orthonormalisieren( ewproblem )/* Eigenvektoren */
    EWPROBLEM *ewproblem;                 /* orthonormalisieren */
{                                         /* nach Kram-Schmidt */
    MATRIX *ev;
/*  MATRIX *entartung; */
    VEKTOR *v,*w,*c,*d,*ortho;
    VEKTOR *vr_normalisieren();
    KOMPLEX *vw;
    INT    spalte,k/*,*gi*/;
    COMHES *comhes;
 
    comhes = ewproblem -> comhes;

    if(test_nullmatrix(comhes->matrix)==JA) return( ewproblem );
 
    ev        = ewproblem -> eigenvektoren;
/*  entartung = ewproblem -> entartung; */
/*  gi        = ewproblem -> gi; */
 
    for( spalte=1; spalte<=ANZ_SP( ev ) ; ++spalte){
 
         v = MXSP( ev , spalte );
 
         if( spalte == 1 )
                v = vr_normalisieren(v);
         else{
                ortho = vr_copy(v);
                for( k=1 ; k<=spalte-1 ; ++k ){
 
                     w  = MXSP( ev , k );
 
                     vw = vskalar(w,v);
                     c  = cv_mult(vw,w);
                     free_(vw);
                     d  = vr_sub(ortho,c);
                     free_vr(c);
                     free_vr(ortho);
                     ortho = d;
                 }
                 ortho = vr_normalisieren(ortho);
                 v     =_vr_copy(v,ortho);
                 free_vr( ortho );
         }
 
    }
 
 
    return( ewproblem );
}
/*----------------------------------------------------------------------------
                                    comhes()
------------------------------------------------------------------------------*/
COMHES *comhes(ewproblem,k,l,mx,overwrite)  /* transformiert eine komplexe */
                                            /* Matrix auf Hessenbergform   */
    EWPROBLEM *ewproblem;
    INT k;
    INT l;
    MATRIX *mx;   /* Pointer auf die Matrix mx */
    INT overwrite;
{
    INT    i,j,la,m,n;
    DOUBLE xr,xi,yr,yi;
 
    INT     *swap,swap_dim;
    KOMPLEX *cdiv();
    KOMPLEX *z;
    COMHES  *comhes;
    MATRIX  *mx_copy(),*alt;
 
    swap_dim = l;
    if( overwrite == NOSPACE){
         comhes = ewproblem -> comhes;
         swap   = comhes -> swap;
    }
    else swap = INT_ALLOC(swap_dim);
 
    for( i=1 ; i<=swap_dim ; ++i )
         VALUE(swap,i) = i;
 
    n  = MXDIM(mx);
    alt = mx;
    if( overwrite==NEIN )   alt = mx_copy( mx );
 
    la = l - 1;
    for( m=k+1 ; m<=la ; ++m ){
         i  = m;
         xr = xi = 0.0;
         for( j=m ; j<=l ;++j )
              if( ABSDR(mx,j,m-1)+ABSDI(mx,j,m-1) > ABSD(xr)+ABSD(xi)  ){
                   xr = R(mx,j,m-1);
                   xi = I(mx,j,m-1);
                   i  = j;
              }
         VALUE(swap,m) = i;
         if( i != m ){/* interchange rows and columns of matrix mx */
              for( j=m-1 ; j<= n ; ++j){
                   yr        = R(mx,i,j);
                   R(mx,i,j) = R(mx,m,j);
                   R(mx,m,j) = yr;
                   yi        = I(mx,i,j);
                   I(mx,i,j) = I(mx,m,j);
                   I(mx,m,j) = yi;
              }
              for( j=1 ; j<=l ; ++j){
                   yr        = R(mx,j,i);
                   R(mx,j,i) = R(mx,j,m);
                   R(mx,j,m) = yr;
                   yi        = I(mx,j,i);
                   I(mx,j,i) = I(mx,j,m);
                   I(mx,j,m) = yi;
              }
         }
         if( xr != 0.0   ||  xi != 0.0 )
              for( i=m+1 ; i<=l ; ++i ){
                   yr = R(mx,i,m-1);
                   yi = I(mx,i,m-1);
                   if( yr != 0.0  ||  yi != 0.0 ){
 
                        z  = cdiv(yr,yi,xr,xi);
                        yr = RT(z);
                        yi = IT(z);
                        free_(z);
 
                        R(mx,i,m-1) = yr;
                        I(mx,i,m-1) = yi;
                        for( j=m ; j<=n ; ++j){
                             R(mx,i,j) += -yr*R(mx,m,j)+yi*I(mx,m,j);
                             I(mx,i,j) += -yr*I(mx,m,j)-yi*R(mx,m,j);
                        }
                        for( j=1 ; j<=l ; ++j){
                             R(mx,j,m) +=  yr*R(mx,j,i)-yi*I(mx,j,i);
                             I(mx,j,m) +=  yr*I(mx,j,i)+yi*R(mx,j,i);
                        }
                   }
              }
    }
    if( overwrite == NOSPACE ) comhes = ewproblem -> comhes;
    else                       comhes = COMHES_ALLOC(1);
    comhes -> overwrite= overwrite;
    comhes -> matrix   = alt;
    comhes -> hesse    = mx;
    comhes -> swap     = swap;
    comhes -> swap_dim = swap_dim;
    return(comhes);
}
/*----------------------------------------------------------------------------
                                    comlr2()
------------------------------------------------------------------------------*/
COMLR2 *comlr2(ewproblem,low,upp,comhes,macheps)/*diagonalisiertHESSENBERGFORM*/
    EWPROBLEM *ewproblem;
    INT     low,upp;                      /* einer komplexen Matrix innerhalb */
    COMHES  *comhes;
    DOUBLE  macheps;                    /* der Rechnergenauigkeit macheps   */
{
    INT overwrite;
    INT    *swap,n;
    INT    iterationsteps=0;  /* Anzahl der Iterationsschritte */
    INT    its=0;
    INT    i,j,k,m,en;
    DOUBLE sr,si,tr,ti,xr,xi,yr,yi,zr,zi,norm,e1,e2;
    KOMPLEX *dummy,*cdiv(),*csqroot();
    MATRIX *ev;                                /* Eigenvektoren */
    VEKTOR *ew;                                /* Eigenwerte    */
    MATRIX *mx;
    COMLR2 *comlr2;
    CHAR   *text;
 
    swap= comhes -> swap;
    mx  = comhes -> hesse;   /* Hessenbergform der zu diagonalisierenen Matrix*/
    n   = MXDIM(mx);
    overwrite = comhes -> overwrite;
 
 if( overwrite == NOSPACE ){
     ev = ewproblem -> eigenvektoren;
     ew = ewproblem -> eigenwerte;
 }
 else{
    ev  = mx_alloc(n,n); /* Speicher fuer die Eigenvektoren holen ,welche */
                         /* hier spaeter als Spalten auftauchen */
    ew  = vr_alloc(n);   /* Speicher fuer die Eigenwerte */
     }
 
    tr   = ti = 0.0;
    text = SUCCESSFUL;/* man geht erstmal davon aus,dass die Iteration k lappt*/
 
    for( i=1 ; i<=n ; ++i){
       for( j=1 ; j<=n ; ++j){
           R(ev,i,j) = 0.0;
           I(ev,i,j) = 0.0;
       }
       R(ev,i,i) = 1.0;
    }
 
    for( i=upp-1 ; i>=low+1 ; --i ){
 
         j = VALUE(swap,i);
 
         for( k=i+1 ; k<=upp ; ++k ){
              R(ev,k,i) = R(mx,k,i-1);
              I(ev,k,i) = I(mx,k,i-1);
         }
 
         if( i != j ){
 
              for( k=i ; k<=upp ; ++k ){
 
                   R(ev,i,k) = R(ev,j,k);
                   I(ev,i,k) = I(ev,j,k);
                   R(ev,j,k) = I(ev,j,k) = 0.0;
              }
              R(ev,j,i) = 1.0;
         }
    }
 
    for( i=1 ; i<=low-1 ; ++i ){
         RV(ew,i) = R(mx,i,i);
         IV(ew,i) = I(mx,i,i);
    }
 
    for( i=upp+1 ; i<=n ; ++i ){
         RV(ew,i) = R(mx,i,i);
         IV(ew,i) = I(mx,i,i);
    }
 
    en = upp;
 
 
nextw:
    if( en < low ) goto fin;
    its = 0;
 
 
nextit:
    /* look for single small sub-diagonal element */
    for( k=en ; k>=low+1 ; --k){
         e1 = ABSDR(mx,k,k-1) + ABSDI(mx,k,k-1);
         e2 = macheps;
 e2 *= ABSDR(mx,k-1,k-1)+ABSDI(mx,k-1,k-1)+ABSDR(mx,k,k)+ABSDI(mx,k,k);
         if( e1 <= e2 )  goto cont1;
    }
    k  = low;
 
 
cont1:
    if( k  == en ) goto root;
    if( its== 30 ){
         text = FAILED;
         goto fail;
    }
 
    /* form shift */
    if( its==10 || its==20 ){
         sr = ABSDR(mx,en,en-1) + ABSDR(mx,en-1,en-2);
         si = ABSDI(mx,en,en-1) + ABSDI(mx,en-1,en-2);
    }
    else{
         sr = R(mx,en,en);
         si = I(mx,en,en);
 
         xr = R(mx,en-1,en)*R(mx,en,en-1) - I(mx,en-1,en)*I(mx,en,en-1);
         xi = R(mx,en-1,en)*I(mx,en,en-1) + I(mx,en-1,en)*R(mx,en,en-1);
 
         if( xr!=0.0  ||  xi!=0.0 ){
              yr    = ( R(mx,en-1,en-1)-sr  )/2.0;
              yi    = ( I(mx,en-1,en-1)-si  )/2.0;
 
              dummy = csqroot( yr*yr-yi*yi+xr , 2.0*yr*yi+xi );
              zr = RT(dummy);
              zi = IT(dummy);
              free_(dummy);
 
              if( yr*zr+yi*zi < 0.0 ){
                   zr *= -1.0;
                   zi *= -1.0;
              }
 
              dummy = cdiv( xr , xi , yr+zr , yi+zi );
              xr    = RT(dummy);
              xi    = IT(dummy);
              free_(dummy);
 
              sr   -= xr;
              si   -= xi;
         }
        }
 
    for( i=low ; i<=en ; ++i ){
         R(mx,i,i) -= sr;
         I(mx,i,i) -= si;
    }
 
    tr  += sr;
    ti  += si;
    its +=  1;
    iterationsteps += 1;
    j    = k+1;
 
    /* look for two consecutive small sub-diagonal elements */
    xr   = ABSDR(mx,en-1,en-1) + ABSDI(mx,en-1,en-1);
    yr   = ABSDR(mx,en  ,en-1) + ABSDI(mx,en  ,en-1);
    zr   = ABSDR(mx,en  ,en  ) + ABSDI(mx,en  ,en  );
 
    for( m=en-1 ; m>=j ; --m){
         yi = yr;
         yr = ABSDR(mx,m,m-1) + ABSDI(mx,m,m-1);
 
         xi = zr;
         zr = xr;
         xr = ABSDR(mx,m-1,m-1) + ABSDI(mx,m-1,m-1);
 
         if( yr <= macheps*zr/yi*(zr+xr+xi) ) goto cont2;
    }
 
    m = k;
 
 
cont2:
    /* triangular decomposition   H = L x R */
    for( i=m+1 ; i<=en ; ++i ){
         xr = R(mx,i-1,i-1);
         xi = I(mx,i-1,i-1);
         yr = R(mx,i  ,i-1);
         yi = I(mx,i  ,i-1);
 
         if( ABSD(xr)+ABSD(xi) < ABSD(yr)+ABSD(yi) ){
              /* interchange rows of mx */
              for( j=i-1 ; j<=n ; ++j){
                   zr          = R(mx,i-1,j);
                   R(mx,i-1,j) = R(mx,i  ,j);
                   R(mx,i  ,j) = zr;
 
                   zi          = I(mx,i-1,j);
                   I(mx,i-1,j) = I(mx,i  ,j);
                   I(mx,i  ,j) = zi;
              }
              dummy = cdiv( xr , xi , yr , yi );
              zr           = RT(dummy);
              zi           = IT(dummy);
              free_(dummy);
 
              RV(ew,i) = 1.0;
         }
         else{
              dummy = cdiv( yr , yi , xr , xi );
              zr           = RT(dummy);
              zi           = IT(dummy);
              free_(dummy);
 
              RV(ew,i) = -1.0;
             }
 
         R(mx,i,i-1) = zr;
         I(mx,i,i-1) = zi;
 
         for( j=i ; j<=n ; ++j ){
              R(mx,i,j) += -zr*R(mx,i-1,j) + zi*I(mx,i-1,j);
              I(mx,i,j) += -zr*I(mx,i-1,j) - zi*R(mx,i-1,j);
         }
    }
 
/* composition   R x L = H */
    for( j=m+1 ; j<=en ; ++j ){
         xr          = R(mx,j,j-1);
         xi          = I(mx,j,j-1);
         R(mx,j,j-1) = I(mx,j,j-1)  =  0.0;
 
         /* interchange columns of mx and ev ,if necessary */
         if( RV(ew,j) > 0.0 ){
              for( i=1 ; i<=j ; ++i){
                   zr          = R(mx,i,j-1);
                   R(mx,i,j-1) = R(mx,i,j  );
                   R(mx,i,j  ) = zr;
 
                   zi          = I(mx,i,j-1);
                   I(mx,i,j-1) = I(mx,i,j  );
                   I(mx,i,j  ) = zi;
              }
 
              for( i=low ; i<=upp ; ++i){
                   zr          = R(ev,i,j-1);
                   R(ev,i,j-1) = R(ev,i,j  );
                   R(ev,i,j  ) = zr;
 
                   zi          = I(ev,i,j-1);
                   I(ev,i,j-1) = I(ev,i,j  );
                   I(ev,i,j  ) = zi;
              }
         }
 
         for( i=1 ; i<=j ; ++i){
              R(mx,i,j-1) +=  xr*R(mx,i,j) - xi*I(mx,i,j);
              I(mx,i,j-1) +=  xr*I(mx,i,j) + xi*R(mx,i,j);
         }
 
         for( i=low ; i<=upp ; ++i){
              R(ev,i,j-1) +=  xr*R(ev,i,j) - xi*I(ev,i,j);
              I(ev,i,j-1) +=  xr*I(ev,i,j) + xi*R(ev,i,j);
         }
    }
 
    goto nextit;
 
 
root:
    /* a root found */
    RV(ew,en) = R(mx,en,en) + tr;
    IV(ew,en) = I(mx,en,en) + ti;
    --en;
    goto nextw;
 
 
fin:
    /* all roots found */
    norm = 0.0;
 
    for( i=1 ; i<=n ; ++i){
         norm += ABSDRV(ew,i) + ABSDIV(ew,i) ;
         for( j=i+1 ; j<=n ; ++j)
              norm += ABSDR(mx,i,j) + ABSDI(mx,i,j);
    }
 
    /* Backsubstitute to find vectors of upper triangular form */
    for( en=n ; en>=2 ; --en){
         xr = RV(ew,en);
         xi = IV(ew,en);
 
         for( i=en-1 ; i>=1 ; --i ){
              zr = R(mx,i,en);
              zi = I(mx,i,en);
 
              for( j=i+1 ; j<=en-1 ; ++j ){
                   zr += R(mx,i,j)*R(mx,j,en) - I(mx,i,j)*I(mx,j,en);
                   zi += R(mx,i,j)*I(mx,j,en) + I(mx,i,j)*R(mx,j,en);
              }
 
              yr = xr - RV(ew,i);
              yi = xi - IV(ew,i);
 
              if( yr==0.0 && yi==0.0 )
                   yr = macheps*norm;
 
              dummy = cdiv( zr , zi , yr , yi );
              R(mx,i,en)   = RT(dummy);
              I(mx,i,en)   = IT(dummy);
              free_(dummy);
         }
    }
 
    /* multiply by transformation matrix to give vectors of orginal      full*/
    /* matrix */
    for( i=1 ; i<=low-1 ; ++i )
         for( j=i+1 ; j<=n ; ++j ){
              R(ev,i,j) = R(mx,i,j);
              I(ev,i,j) = I(mx,i,j);
         }
 
 
    for( i=upp+1 ; i<=n ; ++i )
         for( j=i+1 ; j<=n ; ++j ){
              R(ev,i,j) = R(mx,i,j);
              I(ev,i,j) = I(mx,i,j);
         }
 
    for( j=n ; j>=low ; --j )
         for( i=low ; i<=upp ; ++i ){
              zr = R(ev,i,j);
              zi = I(ev,i,j);
 
              m  = (  (upp<j) ? upp : j-1  );
 
              for( k=low ; k<=m ; ++k){
                   zr += R(ev,i,k)*R(mx,k,j) - I(ev,i,k)*I(mx,k,j);
                   zi += R(ev,i,k)*I(mx,k,j) + I(ev,i,k)*R(mx,k,j);
              }
 
              R(ev,i,j) = zr;
              I(ev,i,j) = zi;
         }
 
 
fail:
    if( overwrite == NOSPACE ) comlr2 = ewproblem;
    else                       comlr2 = COMLR2_ALLOC(1);
    comlr2 -> comhes          = comhes;
    comlr2 -> eigenvektoren   = ev;
    comlr2 -> eigenwerte      = ew;
    comlr2 -> iterationsteps  = iterationsteps;
    comlr2 -> iteration       = text;
    comlr2 -> eps_machine     = macheps;
 
    return( comlr2 );
}
/*----------------------------------------------------------------------------
                                 info_ewproblem()
------------------------------------------------------------------------------*/
void info_ewproblem( ewproblem ) /* Diagonalisierungsroutine testen */
    EWPROBLEM *ewproblem;
{
    KOMPLEX *dummy;
    KOMPLEX *vskalar();
    COMHES *comhes;
    MATRIX *mx,*ev,*c,*d,*e;
    MATRIX *mx_mult(),*mx_sub(),*mx_ewev();
    MATRIX *entartung;
    VEKTOR *ew,*v,*vir,*vis;
    DOUBLE re,im,macheps,shift;
    DOUBLE is_null();
    INT    z,s,dim,its,ze,sp,*nummer;
    INT    *gi,zeile,spalte1,spalte2;
    CHAR   *text;
    CHAR   *textmx  = "RM(%2d,%2d) = %20.12e     IM(%2d,%2d) = %20.12e\n";
    CHAR   *textev1 = "R(ev(%2d) %2d) = %20.12e  I(ev(%2d) %2d) = %20.12e\n";
    CHAR   *textevi = "R(       %2d) = %20.12e  I(       %2d) = %20.12e\n";
    CHAR   *textdf1 = "R(df(%2d) %2d) = %20.12e  I(df(%2d) %2d) = %20.12e\n";
    CHAR   *textew  = "R( ew(%2d) ) = %20.12e     I( ew(%2d) ) = %20.12e\n";
    CHAR   *tdelta  = "< %2d , %2d | %2d , %2d > =  %20.12e + i * %20.12e\n";
    FILE   *fopen(),*out;
 
    out    = fopen_errchk("results/egnwert.info","w");
 
    comhes = ewproblem->comhes;
    mx     = comhes   ->matrix;        /* urspruengliche Matrix holen       */
    ev     = ewproblem->eigenvektoren; /* Eigenvektoren als Spalten in ev   */
    ew     = ewproblem->eigenwerte;    /* Eigenwerte als Zeilen im Spalten- */
                                       /* vektor ew                         */
    text      = ewproblem->iteration;
    its       = ewproblem->iterationsteps;
    macheps   = ewproblem->eps_machine;
    shift     = ewproblem->shift;
    nummer    = ewproblem->energie_nummer;
    entartung = ewproblem->entartung;
    gi        = ewproblem->gi;
 
    dim    = MXDIM(mx);
    fprintf(out,"--------------------------------------------\n");
    fprintf(out,"|Diagonalising a %2d x %2d  Matrix M |\n",dim,dim);
    fprintf(out,"--------------------------------------------\n");
    fprintf(out,"\n");
    fprintf(out,"---------------------------------------------------------\n");
    fprintf(out,"| Problem : diagonalised matix M                         |\n");
    fprintf(out,"---------------------------------------------------------\n");
    fprintf(out,"\n");
    for( s=1 ; s<=dim ; ++s ){
         fprintf(out,"\n");
         for( z=1 ; z<=dim ; ++z ){
              re = is_null( R(mx,z,s),macheps );
              im = is_null( I(mx,z,s),macheps );
              fprintf(out,textmx,z,s,re,z,s,im);
         }
    }
 
    fprintf(out,"\n");
    fprintf(out,"---------------------------------------------------------\n");
    fprintf(out,"| Answer :                                             |\n");
    fprintf(out,"---------------------------------------------------------\n");
    fprintf(out,"\n");
    if(strcmp(text,FAILED)){ fprintf(out," Iteration has failed .\n");
                        exit(1);fclose(out);}
    fprintf(out,"The interation has finsihed .\n");
    fprintf(out,"The calculated accuracy is                 : %6.1e.\n",
                accuracy());
    fprintf(out,"Accuracy in comparison : %6.1e.\n",
                macheps);
    if( its!=1 )
         fprintf(out,"Required %3d Iteration step.\n",its );
    else fprintf(out,"1 interation step is required.\n");
 
    fprintf(out,"\n");
 
if( comhes->overwrite == NEIN ){
   fprintf(out,"-----------------------------------------------------------\n");
   fprintf(out,"|The orthonormal Eigenvectors follow from                 M|\n");
   fprintf(out,"-----------------------------------------------------------\n");
   fprintf(out,"\n");
    for( s=1 ; s<=dim ; ++s ){
         v = MXSP(ev,VALUE(nummer,s));
         for( z=1 ; z<=VRDIM(v) ; ++z ){
              re = is_null( RV(v,z),macheps );
              im = is_null( IV(v,z),macheps );
              if( z==1 )  fprintf(out,textev1,s,z,re,s,z,im);
              else        fprintf(out,textevi,z,re,z,im);
         }
         fprintf(out,"\n");
 
    }
 
 
    fprintf(out,"\n");
    fprintf(out,"---------------------------------------------------------\n");
    fprintf(out,"|                                                       |\n");
    fprintf(out,"| Where H | i r >  = ( E  + E    ) | i r > , r = 1 .. n |\n");
    fprintf(out,"|                       i    shift                     i|\n");
    fprintf(out,"|                                                       |\n");
    fprintf(out,"| TEST    :  < i r | i s > = D                          |\n");
    fprintf(out,"|                             rs                        |\n");
    fprintf(out,"---------------------------------------------------------\n");
    fprintf(out,"\n");
    for( zeile=1 ; zeile<=ANZ_ZE(entartung) ; ++zeile ){
       for( spalte1=1 ; spalte1<=VALUE(gi,zeile); ++spalte1 )
          for( spalte2=spalte1; spalte2<=VALUE(gi,zeile) ; ++spalte2 ){
               vir = MXSP(ev,(INT)R(entartung,zeile,spalte1) );
               vis = MXSP(ev,(INT)R(entartung,zeile,spalte2) );
 
               dummy = vskalar( vir , vis );
               re    = is_null( RT(dummy),macheps );
               im    = is_null( IT(dummy),macheps );
 
               fprintf(out,tdelta,zeile,spalte1,zeile,spalte2,re,im);
          }
       fprintf(out,"\n");
    }
}
    fprintf(out,"\n");
    fprintf(out,"---------------------------------------------------------\n");
    fprintf(out,"|Solution: the intrinsic values of the matrix follow M  |\n");
    fprintf(out,"---------------------------------------------------------\n");
    fprintf(out,"\n");
    fprintf(out,"Energy shift = %20.12e\n",shift);
    fprintf(out,"(The 'true' Energy eigenvalues one gets with the ");
    fprintf(out,"Addition of\n Energy shift with the represented");
    fprintf(out,"Eigenvalue)\n");
    fprintf(out,"\n");
    for( z=1 ; z<=VRDIM(ew) ; ++z ){
         re = is_null( RV(ew,VALUE(nummer,z)),macheps );
         im = is_null( IV(ew,VALUE(nummer,z)),macheps );
         fprintf(out,textew,z,re,z,im);
    }
 
 
if( comhes->overwrite == NEIN ){
    fprintf(out,"\n");
    fprintf(out,"---------------------------------------------------------\n");
    fprintf(out,"| TEST    :  df(i) = ( M - ew(i) ) * ev(i)              |\n");
    fprintf(out,"---------------------------------------------------------\n");
    fprintf(out,"\n");
 
    c = mx_mult(mx,ev);
    d = mx_ewev(ew,ev,shift);
    e = mx_sub((DOUBLE)(1.0),c,d);    /* 1.0*(c-d) */
 
    for( sp=1 ; sp<=MXDIM(e) ; ++sp ){
         for( ze=1 ; ze<=MXDIM(e) ; ++ze ){
              re = is_null( R(e,ze,VALUE(nummer,sp)),macheps );
              im = is_null( I(e,ze,VALUE(nummer,sp)),macheps );
              if( ze==1 )  fprintf(out,textdf1,sp,ze,re,sp,ze,im);
              else        fprintf(out,textevi,ze,re,ze,im);
         }
         fprintf(out,"\n");
    }
 
 
}
 
fclose(out);
}
/*------------------------------------------------------------------------------
ENDEMODUL    D I A H E R M X    C
------------------------------------------------------------------------------*/
