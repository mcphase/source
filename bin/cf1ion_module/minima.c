 
/*-----------------------------------------------------------------------------
 
                                M   O  D  U  L
 
                                   MINIMA C
 
 
Minimiering einer N-dimensionalen Funktion
 
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
 
extern MATRIX    *mx_alloc();        /* definiert in MATRIX.C   */
extern VEKTOR    *vr_alloc();        /* definiert in MATRIX.C   */
extern INT       free_vr();          /* definiert in MATRIX.C   */
extern INT       menue();            /* definiert in ORTHO.C    */
extern INT       if_stabil();        /* definiert in ORTHO.C    */
extern INT       strategie();        /* definiert in ORTHO.C    */
extern INT       is_equal();         /* definiert in DIAHERMX.C */
extern DOUBLE    accuracy();         /* definiert in DIAHERMX.C */
extern INT       VA05A();            /* definiert in VA05A.C    */
 
/*----------------------------------------------------------------------------
Globales Arrays
-----------------------------------------------------------------------------*/
 
extern _CALFUN CF;  /* globale struktur fuer routine CALFUN, */
 
FIT FITIMP[]={
/* 0 */   {  0 , "              " } ,
/* 1 */   {  1 , "Downhill      " } ,
/* 2 */   {  2 , "Powell        " } ,
/* 3 */   {  3 , "Gradient      " } ,
/* 4 */   {  4 , "VA05A         " } ,
/* 5 */   {  5 , "DUMMY         " }
};
#define ANZ_FITIMP  ( sizeof(FITIMP) / sizeof(FIT) )
 
INT flag=0;
 
void CALFUN(m,n,f,x)   /* wird spaeter von der C-routine VA05A()  */
 INT     m, n;         /* aufgerufen                              */
 DOUBLE *f, *x;
{
   DOUBLE   (*funk)();
   SETUP     *setup;
   EWPROBLEM *ewproblem;
   ITERATION *iteration;
   VEKTOR    *p0,*p;
   INT i;
++flag;
 
   funk      = CF.funk;
   setup     = CF.setup;
   ewproblem = CF.ewproblem;
   iteration = CF.iteration;
   p         = CF.vektor;
   p0        = CF.start;
 
   for(i=1; i<= n; ++i)
       RV(p ,i) = *(x+i) * RV(p0,i);
 
 
  *(f+1)=(*funk)(setup,ewproblem,iteration,p );
  *(f+1)= sqrt( *(f+1) );
}
/*----------------------------------------------------------------------------
                                va05a_()
-----------------------------------------------------------------------------*/
MINIMUM *va05a_(_fitnr,setup,ewproblem,iteration,p0,funk)
      MINIMUM   *_fitnr;
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      VEKTOR *p0;
      DOUBLE (*funk)();
{
      DOUBLE macheps;
      VEKTOR *xi,*pp,*vr_alloc();
 
      /* variablendefinition fuer VA05A() */
      INT     m,n,maxfun,iprint,i;
      DOUBLE *f,*x,h,dmax,acc,*w;
 
      /* globale struktur CF setzen */
      CF.funk      = funk;
      CF.setup     = setup;
      CF.ewproblem = ewproblem;
      CF.iteration = iteration;
      CF.vektor    = vr_alloc( VRDIM(p0) );
      CF.start     = vr_alloc( VRDIM(p0) );
 
      maxfun  = FITMAX(iteration);
      macheps = sqrt(ewproblem->eps_machine);
 
 
      n = VRDIM(p0);
      x = DOUBLE_ALLOC(n);
      for(i=1; i<=n; ++i){
          *(x+i) = (DOUBLE)1.;
          RV(CF.start,i) = RV(p0,i);
      }
      m = 1;
      f = DOUBLE_ALLOC(m);
 
      h      = 0.0001;
      acc    = 0.01;
      dmax   = 1.;
      w      = DOUBLE_ALLOC(2*m*n+2*n*n+2*m+5*n);
      iprint = 0;
      VA05A(&m,&n,f,x,&h,&dmax,&acc,&maxfun,&iprint,w);
      switch( iprint ){
      case 1 :
      case 2 :   TEXT1(_fitnr)  = FAILED;
                 TEXT2(_fitnr)  = STABLE;
                 break;
      default:   TEXT1(_fitnr)  = SUCCESSFUL;
                 TEXT2(_fitnr)  = STABLE;
                 break;
      }
 
      xi = vr_alloc( n );
      for(i=1; i<=n; ++i)
         RV(xi,i) = *(w+i);
 
      pp = vr_alloc( n );
      for(i=1; i<=n; ++i)
         RV(pp,i) = *(x+i) * RV(p0,i);
 
      ITER_STEPS( _fitnr  ) = flag;
      P_VEKTOR(   _fitnr  ) = pp;
      XI_VEKTOR(  _fitnr  ) = xi;
      ITERATION(  _fitnr  ) = iteration;
      FRET(       _fitnr  ) = (*(f+1))*(*(f+1)); /* chi_2 */
 
      free_(x);free_(w);free_(f);
      free_vr(CF.vektor);
      free_vr(CF.start );
      return(_fitnr);
}
/*----------------------------------------------------------------------------
                                fitnr5()
-----------------------------------------------------------------------------*/
MINIMUM *fitnr5(_fitnr ,setup,ewproblem,iteration,p0,funk)
      MINIMUM   *_fitnr;
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      VEKTOR *p0;
      DOUBLE (*funk)();
{
      DOUBLE chi2,macheps;
      INT    itmax,iter;
 
      iter    = 0;
      itmax   = FITMAX(iteration);
      macheps = sqrt(ewproblem->eps_machine);
 
      TEXT1(_fitnr)  = SUCCESSFUL;
      TEXT2(_fitnr)  = UNSTABLE;
 
      chi2    = (*funk)(setup,ewproblem,iteration,p0);
 
      ITER_STEPS( _fitnr  ) = iter;
      P_VEKTOR(   _fitnr  ) = p0;
      ITERATION(  _fitnr  ) = iteration;
      FRET(       _fitnr  ) = chi2;      /* func(p) */
      return(_fitnr);
}
/*----------------------------------------------------------------------------
                                amoeba()
 
Siehe "Numerical Recipes", Seite 292-293
-----------------------------------------------------------------------------*/
MINIMUM *amoeba(_amoeba,setup,ewproblem,iteration,p0,funk)
      MINIMUM   *_amoeba;
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      VEKTOR *p0;
      DOUBLE (*funk)();
{
      DOUBLE alpha = 1.0;  /* nach ref: 1.0 */
      DOUBLE beta  = 0.5;  /*           0.5 */
      DOUBLE gamma = 2.0;  /*           2.0 */
      INT    itmax;
 
      MATRIX *p, *mx_alloc();
      VEKTOR *y, *vr_alloc(), *pr, *prr, *pbar ;
      INT    ndim,iter,ze,sp,mpts;
      INT    ilo,ihi,inhi,i,j,mw;
      DOUBLE rtol,ypr,yprr,macheps,accuracy(),ftol;
 
      itmax          = FITMAX(iteration);
      iter           = 0;
 
 
      macheps = sqrt(ewproblem->eps_machine);
      ftol    = sqrt(macheps);
      ndim    = VRDIM(p0);
      mw      = MAX_WIEDER(_amoeba);
 
      p = mx_alloc(ndim+1, ndim);
 
      for(sp=1; sp<=ndim; ++sp)       /* erste Zeile von p ist p0 */
             R(p,1,sp) = RV(p0,sp);
 
      for(ze=2; ze<=ndim+1; ++ze)     /*  rest von p loeschen */
         for(sp=1; sp<=ndim; ++sp)
             R(p,ze,sp) = 0.0;
 
      for(ze=2; ze<=ndim+1; ++ze)     /* restliche Zeilen Einheitsvektoren*/
             R(p,ze,ze-1) = 1.0;
 
 
      y = vr_alloc(ndim+1);
 
      pr    = vr_alloc(ndim+1);
      for(ze=1; ze<=ndim+1; ++ze)    /* y=( f(p0),f(e1),...,f(e9) ) */
          if( ze==1 )
             RV(y,1) = (*funk)(setup,ewproblem,iteration,p0);
          else{
                for(i=1; i<=ndim; ++i)
                    RV(pr,i) = 0.0;
                RV(pr,ze-1)  = 1.0;
                RV(y, ze)    = (*funk)(setup,ewproblem,iteration,pr);
              }
 
 
 
 
      prr    = vr_alloc( ndim+1 );
      pbar   = vr_alloc( ndim+1 );
 
 
      mpts           = ndim+1;
      TEXT1(_amoeba)  = SUCCESSFUL;
      TEXT2(_amoeba)  = UNSTABLE;
 
anfang:
 
      ilo = 1;
      if( RV(y,1) > RV(y,2) ) { ihi=1; inhi=2; }
      else                    { ihi=2; inhi=1; }
 
      for(i=1; i<=mpts; ++i){
         if( RV(y,i) < RV(y,ilo) )   ilo =i;
         if( RV(y,i) > RV(y,ihi) ){  inhi=ihi; ihi=i; }
         else  if( RV(y,i) > RV(y,inhi) ){ if(i!=ihi) inhi=i; }
      }
      rtol  = 2.0*ABSD(RV(y,ihi) - RV(y,ilo));
 
      if( !is_equal( ABSD(RV(y,ihi))+ABSD(RV(y,ilo)), 0.0,macheps ) )
               rtol /=     ABSD(RV(y,ihi))+ABSD(RV(y,ilo));
      else{ TEXT1(_amoeba)=FAILED;
            TEXT2(_amoeba)=UNSTABLE; goto ende;}
      if( rtol <= ftol)                          goto ende;
      if( iter==itmax) {TEXT1(_amoeba)=FAILED;
                        TEXT2(_amoeba)=UNSTABLE;  goto ende;}
      if( *(TEXT2(_amoeba))== *(STABLE) && _amoeba->anz_wiederholung ==mw)
          goto ende;
 
      ++iter;
 
      for(j=1; j<=ndim; ++j)
         RV(pbar,j) = 0.0;
 
 
      for(i=1; i<=mpts; ++i)
         if( i!= ihi)
             for(j=1; j<=ndim; ++j)
                RV(pbar,j) += R(p,i,j);
 
      for(j=1; j<=ndim; ++j){
         RV(pbar,j) /= ndim;
         RV(pr,  j)  = (1.0+alpha)*RV(pbar,j) - alpha*R(p,ihi,j);
      }
 
 
      ypr = (*funk)(setup,ewproblem,iteration,pr);
 
      if_stabil(_amoeba,ewproblem,iteration,pr,iter,ypr);
 
      if( IS_MENUE(iteration) )
          if(iter%MENUE(iteration) ==0)
              menue(_amoeba,setup,ewproblem,iteration,pr,iter,ypr);
 
 /*   if(iter%10 ==0) strategie(ewproblem,iteration,pr);  */
 
 
      if( ypr <= RV(y,ilo) ){
          for(j=1; j<=ndim; ++j)
             RV(prr,j) = gamma*RV(pr,j) + (1.0-gamma)*RV(pbar,j);
 
          yprr = (*funk)(setup,ewproblem,iteration,prr);
 
          if( yprr < RV(y,ilo) ){
              for(j=1; j<=ndim; ++j)
                  R(p,ihi,j) = RV(prr,j);
              RV(y,ihi) = yprr;
          }
          else{
                for(j=1; j<=ndim; ++j)
                    R(p,ihi,j) = RV(pr,j);
                RV(y,ihi) = ypr;
              }
      }
      else if( ypr >= RV(y,inhi) ){
 
              if( ypr < RV(y,ihi) ){
                  for(j=1; j<=ndim; ++j)
                      R(p,ihi,j) = RV(pr,j);
                  RV(y,ihi) = ypr;
              }
 
              for(j=1; j<=ndim; ++j)
                  RV(prr,j) = beta*R(p,ihi,j) + (1.0-beta)*RV(pbar,j);
 
              yprr = (*funk)(setup,ewproblem,iteration,prr);
 
              if( yprr < RV(y,ihi) ){
                  for(j=1; j<=ndim; ++j)
                      R(p,ihi,j) = RV(prr,j);
                  RV(y,ihi) = yprr;
              }
              else{
                     for(i=1; i<=mpts; ++i)
                         if( i!=ilo){
                             for(j=1; j<=ndim; ++j){
                                 RV(pr,j) = 0.5 *(R(p,i,j)+R(p,ilo,j));
                                 R(p,i,j) = RV(pr,j);
                             }
                            RV(y,i) = (*funk)(setup,ewproblem,iteration,pr);
                         }
                  }
          }
          else{
                  for(j=1; j<=ndim; ++j)
                      R(p,ihi,j) = RV(pr,j);
                  RV(y,ihi) = ypr;
              }
 
          goto anfang;
ende:
 
      free_vr(pr);
      free_vr(prr);
      free_vr(pbar);
 
      ITER_STEPS(_amoeba ) = iter;
      MATRIX(    _amoeba ) = p;
      VEKTOR(    _amoeba ) = y;
      ITERATION( _amoeba ) = iteration;
 
      return( _amoeba );
}
/*----------------------------------------------------------------------------
                                powell()
 
Siehe "Numerical Recipes", Seite 299-300
-----------------------------------------------------------------------------*/
MINIMUM *powell(_powell,setup,ewproblem,iteration,pv,func)
      MINIMUM   *_powell;
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      VEKTOR    *pv;
      DOUBLE    (*func)();
{
 
      MATRIX *xi, *mx_alloc();
      VEKTOR *vr_alloc(),  *pt, *ptt, *xit, *p;
      INT    n,iter,i,j,ibig,ze,sp,mw;
      DOUBLE macheps,ftol,fret;
      DOUBLE fp,fptt,del,t;
      LINMIN *linmin(),*_linmin;
      MINIMUM *brent();
 
      INT    itmax;
 
      itmax   = FITMAX(iteration);
 
      mw      = MAX_WIEDER(_powell);
      macheps = ewproblem->eps_machine;
 
      ftol    = sqrt(macheps);
      n       = VRDIM(pv);
      p       = vr_alloc( n );
      for(sp=1; sp<=n; ++sp)
          RV(p,sp) = RV(pv,sp);
 
 
      xi = mx_alloc(n,n);
      for(ze=1; ze<=n; ++ze)      /* Matrix xi loeschen */
         for(sp=1; sp<=n; ++sp)
             R(xi,ze,sp) = 0.0;
 
      for(sp=1; sp<=n; ++sp)     /* in die Spalten von xi   */
          R(xi,sp,sp) = 1.0;     /* Einheitsvektoren setzen */
 
 
 
 
      pt             = vr_alloc(n);
      ptt            = vr_alloc(n);
      xit            = vr_alloc(n);
      TEXT1(_powell)  = SUCCESSFUL;
      TEXT2(_powell)  = UNSTABLE;
 
      fret = (*func)(setup,ewproblem,iteration,p);
 
      for(j=1; j<=n; ++j)
          RV(pt,j) = RV(p,j);
 
      iter = 0;
anfang:
       ++iter;
 
       fp   = fret;
       ibig = 0;
       del  = 0.0;
 
       for(i=1; i<=n; ++i){
 
          for(j=1; j<=n; ++j)
              RV(xit,j) = R(xi,j,i);
 
          _linmin= LINMIN_ALLOC(1);
          _linmin= linmin(_linmin,brent,setup,ewproblem,iteration,p,xit,func);
            p    = P_VEKTOR(  _linmin);
            xit  = XI_VEKTOR( _linmin);
            fret = FRET(      _linmin);
          free_(_linmin);
 
          if( ABSD(fp-fret) > del ) {del=ABSD(fp-fret); ibig=i;}
 
       }
 
       if(2*ABSD(fp-fret)<=ftol*(ABSD(fp)+ABSD(fret))) goto ende;
 
       if( iter==itmax) {TEXT1(_powell)=FAILED; TEXT2(_powell)=UNSTABLE;goto ende;}
       if( *(TEXT2(_powell)) == *(STABLE) && _powell->anz_wiederholung==mw)
           goto ende;
 
      if_stabil(_powell,ewproblem,iteration,p,iter,fret);
      if( IS_MENUE(iteration) )
          if(iter%MENUE(iteration) ==0)
              menue(_powell, setup,ewproblem,iteration,p,iter,fret);
 
 /*   if(iter%10 ==0) strategie(ewproblem,iteration,p);  */
 
       for(j=1; j<=n; ++j){
           RV(ptt,j) = 2.0*RV(p,j)-RV(pt,j);
           RV(xit,j) =     RV(p,j)-RV(pt,j);
           RV(pt ,j) =     RV(p,j);
       }
 
       fptt = (*func)(setup,ewproblem,iteration,ptt);
 
       if( fptt >= fp ) goto anfang;
 
       t = 2.0*(fp-2.0*fret+fptt)*(fp-fret-del)*(fp-fret-del)
              -del*(fp-fptt)*(fp-fptt);
 
       if( t    >= 0.0 ) goto anfang;
 
       _linmin = LINMIN_ALLOC(1);
       _linmin = linmin(_linmin,brent,setup,ewproblem,iteration,p,xit,func);
          p    = P_VEKTOR(  _linmin);
          xit  = XI_VEKTOR( _linmin);
          fret = FRET(      _linmin);
       free_(_linmin);
 
       for(j=1; j<=n; ++j)
           R(xi,j,ibig) = RV(xit,j);
 
       goto anfang;
ende  :
 
      free_vr(pt);
      free_vr(ptt);
      free_vr(xit);
 
      ITER_STEPS( _powell ) = iter;
      MATRIX(     _powell ) = xi;
      VEKTOR(     _powell ) = p;
      ITERATION(  _powell ) = iteration;
      FRET(       _powell ) = fret;      /* func(p) */
 
      return( _powell );
}
/*----------------------------------------------------------------------------
                                frprmn()
 
Siehe "Numerical Recipes", Seite 305-306
-----------------------------------------------------------------------------*/
MINIMUM *frprmn(_frprmn,setup,ewproblem,iteration,pv,func)
      MINIMUM   *_frprmn;
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      VEKTOR    *pv;
      DOUBLE    (*func)();
{
 
      VEKTOR *vr_alloc(),  *xi, *g, *h, *p, *gradient();
      INT    n,iter,i,j,ze,sp,mw;
      DOUBLE macheps,ftol,fret,eps,gg,dgg,gam,fp;
      LINMIN *linmin(),*_linmin;
      MINIMUM *dbrent();
 
      INT    itmax,its;
 
      itmax   = FITMAX(iteration);
      mw      = MAX_WIEDER(_frprmn);
      macheps = ewproblem->eps_machine;
 
      ftol    = sqrt(macheps);
      eps     = sqrt(macheps);
 
      n       = VRDIM(pv);
      p       = vr_alloc( n );
      for(sp=1; sp<=n; ++sp)
          RV(p,sp) = RV(pv,sp);
 
      g              = vr_alloc(n);
      h              = vr_alloc(n);
      TEXT1(_frprmn) = SUCCESSFUL;
      TEXT2(_frprmn) = UNSTABLE;
 
      fp = (*func)(setup,ewproblem,iteration,p);
      xi = gradient(func,setup,ewproblem,iteration,p);
 
      for(j=1; j<=n; ++j){
          RV(g,j)=-RV(xi,j);  RV(h,j)=RV(g,j);  RV(xi,j)=RV(h,j);
      }
 
      for(its=1; its<=itmax; ++its){
          iter = its;
          _linmin= LINMIN_ALLOC(1);
          _linmin=  linmin(_linmin,dbrent,setup,ewproblem,iteration,p,xi,func);
            p    = P_VEKTOR(  _linmin);
            xi   = XI_VEKTOR( _linmin);
            fret = FRET(      _linmin);
          free_(_linmin);
 
          if( 2.0*ABSD(fret-fp)<=ftol*(ABSD(fret)+ABSD(fp)+eps) ) goto ende;
          if( iter==itmax) {TEXT1(_frprmn)=FAILED; TEXT2(_frprmn)=UNSTABLE;goto ende;}
          if( *(TEXT2(_frprmn))==*(STABLE)&& _frprmn->anz_wiederholung==mw) goto ende;
          if_stabil(_frprmn,ewproblem,iteration,p,iter,fret);
 
          if( IS_MENUE(iteration) )
             if(iter%MENUE(iteration) ==0)
                menue(_frprmn,setup,ewproblem,iteration,p,iter,fret);
 
 /*   if(iter%10 ==0) strategie(ewproblem,iteration,p);  */
 
          fp = (*func)(setup,ewproblem,iteration,p);
          xi = gradient(func,setup,ewproblem,iteration,p);
 
          gg=0.0; dgg=0.0;
          for(j=1; j<=n; ++j){
              gg  += RV(g ,j)*RV(g ,j);
             dgg  += RV(xi,j)*RV(xi,j);
             dgg  += (RV(g,j)+RV(xi,j))*RV(xi,j);
          }
 
          if( is_equal( gg, 0.0, macheps ) ) goto ende;
 
          gam = dgg/gg;
 
          for(j=1; j<=n; ++j){
              RV(g ,j) = -RV(xi,j);
              RV(h ,j) =  RV(g ,j) + gam*RV(h,j);
              RV(xi,j) =  RV(h ,j);
          }
 
      }
 
ende:
      free_vr(g);
      free_vr(h);
 
      ITER_STEPS( _frprmn ) = iter;
      P_VEKTOR(   _frprmn ) = p;
      XI_VEKTOR(  _frprmn ) = xi;
      ITERATION(  _frprmn ) = iteration;
      FRET(       _frprmn ) = fret;      /* func(p) */
 
      return( _frprmn );
}
/*----------------------------------------------------------------------------
                                linmin()
 
Siehe "Numerical Recipes", Seite 300-301
-----------------------------------------------------------------------------*/
LINMIN *linmin(_linmin,brent,setup,ewproblem,iteration,p,xi,func)
      LINMIN    *_linmin;
      MINIMUM   *(*brent)();
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      VEKTOR    *p,*xi;
      DOUBLE    (*func)();
{
      INT j;
      DOUBLE  ax,xx,bx,fa,fx,fb,fret,xmin;
      MNBRAK  *_mnbrak,*mnbrak();
      MINIMUM *_brent;
 
      ax = 0.0;
      xx = 1.0;
      _mnbrak = MNBRAK_ALLOC(1);
      _mnbrak = mnbrak(_mnbrak,setup,ewproblem,iteration,ax,xx,func,p,xi );
           ax = AX(_mnbrak);
           xx = BX(_mnbrak);
           bx = CX(_mnbrak);
      free_(_mnbrak);
 
      _brent  = MINIMUM_ALLOC(1);
      _brent  = (*brent)(_brent,setup,ewproblem,iteration,ax,xx,bx,func,p,xi );
        xmin  = XMIN( _brent);
        fret  = BRENT(_brent);
      free_(_brent);
 
      for( j=1; j<= VRDIM(p); ++j){
           RV(xi,j) = xmin * RV(xi,j);
           RV(p ,j) = RV(p,j) + RV(xi,j);
      }
 
      P_VEKTOR( _linmin) = p;
      XI_VEKTOR(_linmin) = xi;
      FRET(     _linmin) = fret;
 
      return(_linmin);
}
/*----------------------------------------------------------------------------
                                gradient()
 
 
-----------------------------------------------------------------------------*/
VEKTOR *gradient(func,setup,ewproblem,iteration,xt)
      DOUBLE    (*func)();
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      VEKTOR    *xt;
{
      VEKTOR *df,*vr_alloc();
      INT i;
      DOUBLE dfx;
      DOUBLE eps;
 
      eps = ewproblem->eps_machine;
      eps = sqrt(eps);
 
      df  = vr_alloc(  VRDIM(xt)  );
      for( i=1; i<= VRDIM(df); ++i){
           RV(xt,i) += eps; dfx = (*func)(setup,ewproblem,iteration,xt);
           RV(xt,i) -= eps; dfx-= (*func)(setup,ewproblem,iteration,xt);
           RV(df,i) = dfx/eps;
      }
      return(df);
}
/*----------------------------------------------------------------------------
                                f1dim()
 
Siehe "Numerical Recipes", Seite 301
-----------------------------------------------------------------------------*/
DOUBLE f1dim(setup,ewproblem,iteration,func,p,xi,x)
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      DOUBLE    (*func)();
      VEKTOR    *p,*xi;
      DOUBLE    x;
{
      INT     j;
      VEKTOR  *xt,*vr_alloc();
      DOUBLE  fxt;
 
 
      xt = vr_alloc( VRDIM(p) );
      for( j=1; j<=VRDIM(p); ++j)
           RV(xt,j) = RV(p,j) + x * RV(xi,j);
 
      fxt = (*func)(setup,ewproblem,iteration,xt);
 
      free_vr(xt);
      return( fxt );
}
/*----------------------------------------------------------------------------
                               df1dim()
 
Siehe "Numerical Recipes", Seite 306
-----------------------------------------------------------------------------*/
DOUBLE df1dim(setup,ewproblem,iteration,func,p,xi,x)
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      DOUBLE    (*func)();
      VEKTOR    *p,*xi;
      DOUBLE    x;
{
      INT     j;
      VEKTOR  *xt,*vr_alloc(),*gradient(),*df;
      DOUBLE  df1;
 
 
      xt = vr_alloc( VRDIM(p) );
      for( j=1; j<=VRDIM(p); ++j)
           RV(xt,j) = RV(p,j) + x * RV(xi,j);
 
      df  = gradient(func,setup,ewproblem,iteration,xt);
 
      df1=0.0;
      for( j=1; j<=VRDIM(p); ++j)
           df1 += RV(df,j)*RV(xi,j);
 
      free_vr(xt);
      free_vr(df);
      return( df1 );
}
/*----------------------------------------------------------------------------
                                mnbrak()
 
Siehe "Numerical Recipes", Seite 281-282
-----------------------------------------------------------------------------*/
MNBRAK *mnbrak(_mnbrak,setup,ewproblem,iteration,ax,bx,func,p,xi)
      MNBRAK    *_mnbrak;
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      DOUBLE    ax,bx,(*func)();
      VEKTOR    *p,*xi;
{
      DOUBLE gold;
      DOUBLE glimit = 100.0;
      DOUBLE tiny;
 
      DOUBLE ulim,u,q,r,cx,fa,fb,fc,fu,dum;
      DOUBLE f1dim();
 
      gold = 0.5*(1.0 + sqrt(5.0));
      tiny = sqrt(ewproblem->eps_machine);
 
      fa = f1dim(setup,ewproblem,iteration,func,p,xi,ax);
      fb = f1dim(setup,ewproblem,iteration,func,p,xi,bx);
      if( fb > fa ){
         dum = ax; ax  = bx; bx  = dum;
         dum = fb; fb  = fa; fa  = dum;
      }
 
      cx = bx + gold*(bx-ax);
      fc = f1dim(setup,ewproblem,iteration,func,p,xi,cx);
 
anfang:
 
      if( fb >= fc ){
          r = (bx-ax)*(fb-fc);
          q = (bx-cx)*(fb-fa);
          u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*DDSIGN(MAX(ABSD(q-r),tiny),q-r));
          ulim = bx + glimit*(cx-bx);
          if( (bx-u)*(u-cx)> 0.0 ){
              fu = f1dim(setup,ewproblem,iteration,func,p,xi,u);
              if( fu < fc ){
                  ax = bx; fa = fb;
                  bx = u;  fb = fu;
                  goto anfang;
              }
              else if( fu > fb ){
                       cx = u; fc = fu;
                       goto anfang;
                   }
              u  = cx + gold*(cx-bx);
              fu = f1dim(setup,ewproblem,iteration,func,p,xi,u);
 
          }
          else if( (cx-u)*(u-ulim) > 0.0 ){
              fu = f1dim(setup,ewproblem,iteration,func,p,xi,u);
              if( fu < fc ){
                  bx = cx; cx = u;
                  u  = cx + gold*(cx-bx);
 
                  fb = fc; fc = fu;
                  fu = f1dim(setup,ewproblem,iteration,func,p,xi,u);
              }
          }
          else if( (u-ulim)*(ulim-cx) >= 0.0 ){
                  u  = ulim;
                  fu = f1dim(setup,ewproblem,iteration,func,p,xi,u);
          }
          else{   u  = cx + gold*(cx-bx);
                  fu = f1dim(setup,ewproblem,iteration,func,p,xi,u);
          }
 
          ax = bx; bx = cx; cx = u;
          fa = fb; fb = fc; fc = fu;
          goto anfang;
 
      }
 
      AX(_mnbrak) = ax;
      BX(_mnbrak) = bx;
      CX(_mnbrak) = cx;
 
      return(_mnbrak);
}
/*----------------------------------------------------------------------------
                                brent()
 
Siehe "Numerical Recipes", Seite 284-286
-----------------------------------------------------------------------------*/
MINIMUM *brent(_brent,setup,ewproblem,iteration,ax,bx,cx,func,pv,xi)
      MINIMUM   *_brent;
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      DOUBLE    ax,bx,cx,(*func)();
      VEKTOR    *pv,*xi;
{
      INT    iter,itmax  = 100;
      DOUBLE cgold;
      DOUBLE zeps;
      DOUBLE f1dim();
      DOUBLE macheps;
 
      DOUBLE x,p,q,r,v,w,e,a,b,d,xm,tol,tol1,tol2,fx,fv,fw;
      DOUBLE etemp,u,fu;
 
      macheps         = ewproblem -> eps_machine;
      cgold           = 0.5*(3.0 - sqrt(5.0));
      TEXT1(_brent )  = SUCCESSFUL;
      TEXT2(_brent )  = STABLE;
 
      zeps = sqrt(macheps);
      tol  = sqrt(macheps);
 
      a = MIN( ax,cx );
      b = MAX( ax,cx );
 
      v = bx;
      w = v;
      x = v;
      e = 0.0;
 
      fx = f1dim(setup,ewproblem,iteration,func,pv,xi,x);
      fv = fx;
      fw = fx;
 
      for( iter=1; iter<=itmax; ++iter){
           xm   = 0.5*(a+b);
           tol1 = tol*ABSD(x)+zeps;
           tol2 = 2.0*tol1;
 
           if( ABSD(x-xm) <= (tol2-0.5*(b-a)) ) goto drei;
 
           if( ABSD(e) > tol1 ){
               r = (x-w)*(fx-fv);
               q = (x-v)*(fx-fw);
               p = (x-v)*q - (x-w)*r;
               q = 2.0*(q-r);
 
               if( q > 0.0 ) p *= -1.0;
               q = ABSD(q);
 
               etemp = e; e = d;
 
               if( (ABSD(p)>=ABSD(0.5*q*etemp))||(p<=q*(a-x))||
                   (p>=q*(b-x))) goto eins;
 
               d = p/q;
               u = x+d;
 
               if( ((u-a)<tol2)||((b-u)<tol2) ) d = DDSIGN(tol1,xm-x);
 
               goto zwei;
           }
eins:
           if( x >= xm ) e = a - x;
           else          e = b - x;
           d = cgold*e;
zwei:
           if( ABSD(d) >= tol1 ) u = x + d;
           else                  u = x + DDSIGN(tol1,d);
 
           fu = f1dim(setup,ewproblem,iteration,func,pv,xi,u);
 
           if( fu <= fx ){
               if( u >= x ) a = x;
               else         b = x;
 
               v = w; fv = fw;
               w = x; fw = fx;
               x = u; fx = fu;
           }
           else{  if( u < x )  a = u;
                  else         b = u;
 
                  if( (fu <= fw)|| (w == x) ){
                      v = w; fv = fw;
                      w = u; fw = fu;
                  }
                  else if( (fu<=fv) || (v==x) ||(v==w) )
                       v = u; fv = fu;
 
               }
 
      }/* end for(...iter... */
 
      TEXT1(_brent )=FAILED;
      TEXT2(_brent )=UNSTABLE;
 
drei:
      XMIN( _brent)  = x;
      BRENT(_brent)  = fx;
 
      return(_brent);
}
/*----------------------------------------------------------------------------
                               dbrent()
 
Siehe "Numerical Recipes", Seite 287-289
-----------------------------------------------------------------------------*/
MINIMUM *dbrent(_brent,setup,ewproblem,iteration,ax,bx,cx,func,pv,xi)
      MINIMUM   *_brent;
      SETUP     *setup;
      EWPROBLEM *ewproblem;
      ITERATION *iteration;
      DOUBLE    ax,bx,cx,(*func)();
      VEKTOR    *pv,*xi;
{
      INT    iter,itmax  =  100;
      INT    ok1, ok2;
      DOUBLE zeps;
      DOUBLE f1dim();
      DOUBLE macheps;
      DOUBLE df1dim();
 
 
      DOUBLE v,w,x,e,xm,a,b,tol,tol1,tol2,fx,fv,fw,fu,dx,dv,dw,du;
      DOUBLE etemp,u,d1,d2,u1,u2,d,olde;
 
      macheps         = ewproblem -> eps_machine;
      TEXT1(_brent )  = SUCCESSFUL;
      TEXT2(_brent )  = STABLE;
 
      zeps = sqrt(macheps);
      tol  = sqrt(macheps);
 
      a = MIN( ax,cx );
      b = MAX( ax,cx );
 
      v = bx;
      w = v;
      x = v;
      e = 0.0;
 
      fx = f1dim(setup,ewproblem,iteration,func,pv,xi,x);
      fv = fx;
      fw = fx;
 
      dx = df1dim(setup,ewproblem,iteration,func,pv,xi,x);
      dv = dx;
      dw = dx;
 
      for( iter=1; iter<=itmax; ++iter){
           xm   = 0.5*(a+b);
           tol1 = tol*ABSD(x)+zeps;
           tol2 = 2.0*tol1;
 
           if( ABSD(x-xm) <= (tol2-0.5*(b-a)) ) goto drei;
 
           if( ABSD(e) > tol1 ){
               d1 = 2.0*(b-a);
               d2 = d1;
               if( dw != dx ) d1=(w-x)*dx/(dx-dw);
               if( dv != dx ) d2=(v-x)*dx/(dx-dv);
               u1=x+d1; u2=x+d2;
               ok1=((a-u1)*(u1-b)>0.0)&&(dx*d1<=0.0);
               ok2=((a-u2)*(u2-b)>0.0)&&(dx*d2<=0.0);
               olde=e; e=d;
               if( !(ok1||ok2) ) goto eins;
               else
                     if( ok1 && ok2 ){
                        if( ABSD(d1)<ABSD(d2) ) d=d1;
                        else                    d=d2;
                     }
 
 
               else{
                     if(ok1) d=d1;
                     else    d=d2;
                   }
 
               if( ABSD(d) > ABSD(0.5*olde) ) goto eins;
 
               u=x+d;
 
               if(u-a<tol2 || b-u<tol2) d=DDSIGN(tol1,xm-x);
 
               goto zwei;
 
           }/* end if( ABSD(e) > tol1 ){  */
 
eins:     if( dx >= 0.0 )  e=a-x;
          else             e=b-x;
 
          d=0.5*e;
 
zwei:     if( ABSD(d) >= tol1 ){
              u  = x+d;
              fu = f1dim(setup,ewproblem,iteration,func,pv,xi,u);
          }
          else{
                u  = x + DDSIGN(tol1,d);
                fu = f1dim(setup,ewproblem,iteration,func,pv,xi,u);
                if( fu > fx ) goto drei;
              }
 
          du = df1dim(setup,ewproblem,iteration,func,pv,xi,u);
          if( fu <= fx ){
              if( u>=x ) a=x;
              else       b=x;
              v = w; fv=fw; dv=dw;
              w = x; fw=fx; dw=dx;
              x = u; fx=fu; dx=du;
          }
          else{
                if( u<x ) a=u;
                else      b=u;
                if( fu<=fw || w==x ){
                    v = w; fv=fw; dv=dw;
                    w = u; fw=fu; dw=du;
                }
                else{
                      if( fu<=fv || v==x || v==w ){
                          v = u; fv=fu; dv=du;
                      }
                    }
              }
 
      }/* end for(...iter... */
 
      TEXT1(_brent )=FAILED;
      TEXT2(_brent )=UNSTABLE;
 
drei:
      XMIN( _brent)  = x;
      BRENT(_brent)  = fx;
 
      return(_brent);
}
/*------------------------------------------------------------------------------
ENDEMODUL    M I N I M A    C
------------------------------------------------------------------------------*/
