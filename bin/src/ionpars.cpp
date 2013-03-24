// ionpars: class to load and store matrices for internal module cfield and so1ion

#include "ionpars.hpp"
#include "martin.h"
#include "ionpars.h"
#include "myev.h"
#include "perlparse.c"

#define NOF_OLM_MATRICES 48
#define NOF_RIXS_MATRICES 9

#define SMALL 1e-6   //!!! must match SMALL in mcdisp.c and ionpars.cpp !!!
                     // because it is used to decide wether for small transition
		     // energy the matrix Mijkl contains wn-wn' or wn/kT

 ionpars::ionpars (const ionpars & p) //copy constructor
 {J=p.J;so1ion=p.so1ion;
  Ja=p.Ja; Jb=p.Jb; Jc=p.Jc;Hcf=p.Hcf;
  Jaa=p.Jaa; Jbb=p.Jbb; Jcc=p.Jcc;
  gJ=p.gJ;nof_electrons=p.nof_electrons;
  alpha=p.alpha;beta=p.beta;gamma=p.gamma;
  r2=p.r2;r4=p.r4;r6=p.r6;
  sigma0=p.sigma0;sigma1=p.sigma1;sigma2=p.sigma2;
  Blm=p.Blm; // vector of crystal field parameters
  Llm=p.Llm; // vector of crystal field parameters
  // cnst is the Zlm constants - put them into the matrix ... (same code is reused in jjjpar.cpp, pointc.c)
   cnst=Matrix(0,6,-6,6);set_zlm_constants(cnst);
 
   int i;
   iontype = new char [strlen(p.iontype)+1];
   strcpy(iontype,p.iontype);
  
    Ri= new ComplexMatrix * [1+NOF_RIXS_MATRICES];
    for(i=1;i<=NOF_RIXS_MATRICES;++i){Ri[i]= new ComplexMatrix(1,(*p.Ri[i]).Rhi(),1,(*p.Ri[i]).Chi());
                   (*Ri[i])=(*p.Ri[i]);}
    Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
    OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 
    for(i=1;i<=NOF_OLM_MATRICES;++i){ Olm [i]= new Matrix(1,(*p.Olm[i]).Rhi(),1,(*p.Olm[i]).Chi()); 
                                      OOlm [i] = new ComplexMatrix(1,(*p.OOlm[i]).Rhi(),1,(*p.OOlm[i]).Chi()); 
                   (*Olm[i])=(*p.Olm[i]);
                   (*OOlm[i])=(*p.OOlm[i]);
    In= new ComplexMatrix * [1+IONPARS_MAXNOFCOMPONENTS+3+NOF_OLM_MATRICES+NOF_RIXS_MATRICES]; 
    for(i=1;i<=IONPARS_MAXNOFCOMPONENTS;++i){In[i]= new ComplexMatrix(1,(*p.In[i]).Rhi(),1,(*p.In[i]).Chi());
                   (*In[i])=(*p.In[i]);}
   }  

 }
ionpars::ionpars (int dimj) // constructor from dimj
 { 
   alpha=0;beta=0;gamma=0;r2=0;r4=0;r6=0;nof_electrons=0;

   so1ion=0;
  sigma0=0;sigma1=0;sigma2=0;
  J=((double)dimj-1)/2;
  Ja=Matrix(1,dimj,1,dimj);
  Jb=Matrix(1,dimj,1,dimj);
  Jc=Matrix(1,dimj,1,dimj);
  Hcf=Matrix(1,dimj,1,dimj);
  Jaa=ComplexMatrix(1,dimj,1,dimj);
  Jbb=ComplexMatrix(1,dimj,1,dimj);
  Jcc=ComplexMatrix(1,dimj,1,dimj);
// cnst is the Zlm constants - put them into the matrix ... (same code is reused in jjjpar.cpp, pointc.c)
   cnst=Matrix(0,6,-6,6);set_zlm_constants(cnst);
   iontype = new char [MAXNOFCHARINLINE];
   Blm=Vector(0,48);Blm=0; // vector of crystal field parameters
   Llm=Vector(0,45);Llm=0; // vector of crystal field parameters
 int i;   
   Ri = new ComplexMatrix * [1+NOF_RIXS_MATRICES];
  for(i=1;i<=NOF_RIXS_MATRICES;++i)Ri[i]= new ComplexMatrix(1,dimj,1,dimj);
   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 
  for(i=1;i<=NOF_OLM_MATRICES;++i){ Olm [i]= new Matrix(1,dimj,1,dimj); 
                                    OOlm [i] = new ComplexMatrix(1,dimj,1,dimj); 
                                  }  
   In = new ComplexMatrix * [1+IONPARS_MAXNOFCOMPONENTS+3+NOF_OLM_MATRICES+NOF_RIXS_MATRICES];
  for(i=1;i<=IONPARS_MAXNOFCOMPONENTS;++i)In[i]= new ComplexMatrix(1,dimj,1,dimj);
  
}

ionpars::ionpars (char * ion) // constructor from iontype (mind:no matrices filled with values !)
 {int dimj;
  getpar(ion, &dimj, &alpha, &beta, &gamma, &gJ,&r2, &r4,&r6, &nof_electrons );
   iontype = new char [strlen(ion)+1];
   strcpy(iontype,ion);

    so1ion=0;
   sigma0=0;sigma1=0;sigma2=0;
  J=((double)dimj-1)/2;
  Ja=Matrix(1,dimj,1,dimj);
  Jb=Matrix(1,dimj,1,dimj);
  Jc=Matrix(1,dimj,1,dimj);
  Hcf=Matrix(1,dimj,1,dimj);
  Jaa=ComplexMatrix(1,dimj,1,dimj);
  Jbb=ComplexMatrix(1,dimj,1,dimj);
  Jcc=ComplexMatrix(1,dimj,1,dimj);
// cnst is the Zlm constants - put them into the matrix ... (same code is reused in jjjpar.cpp, pointc.c)
cnst=Matrix(0,6,-6,6);set_zlm_constants(cnst);

   Blm=Vector(0,48);Blm=0; // vector of crystal field parameters
   Llm=Vector(0,45);Llm=0; // vector of crystal field parameters

 int i;   
   Ri = new ComplexMatrix * [1+NOF_RIXS_MATRICES];
  for(i=1;i<=NOF_RIXS_MATRICES;++i)Ri[i]= new ComplexMatrix(1,dimj,1,dimj);
   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 
  for(i=1;i<=NOF_OLM_MATRICES;++i){ Olm [i]= new Matrix(1,dimj,1,dimj); 
                                    OOlm [i] = new ComplexMatrix(1,dimj,1,dimj); 
                                  }  
   In = new ComplexMatrix * [1+IONPARS_MAXNOFCOMPONENTS+3+NOF_OLM_MATRICES+NOF_RIXS_MATRICES];
  for(i=1;i<=IONPARS_MAXNOFCOMPONENTS;++i)In[i]= new ComplexMatrix(1,dimj,1,dimj);
  
}

ionpars::~ionpars(){
 int i;
 delete []iontype;
 for (i=1;i<=NOF_RIXS_MATRICES;++i)delete Ri[i];
 for (i=1;i<=NOF_OLM_MATRICES;++i)  {delete Olm[i];delete OOlm[i];}
 for (i=1;i<=IONPARS_MAXNOFCOMPONENTS;++i)delete In[i];
   delete[] Olm;
   delete[] OOlm;
   delete[] Ri;
   delete[] In; 
 } //destructor

ionpars::ionpars(FILE * cf_file) 
//constructor with commands from file handle (filename of cf parameters etc)
{ static int pr=1;
  int dimj;complex<double> im(0,1);
  int i,j,l,m; 
  double alphar,betar,gammar,r2r,r4r,r6r,gJr;
  char instr[MAXNOFCHARINLINE];
  iontype= new char[MAXNOFCHARINLINE];
  char  moduletype[MAXNOFCHARINLINE];
   Blm=Vector(0,48);Blm=0; // vector of crystal field parameters
   Llm=Vector(0,45);Llm=0; // vector of crystal field parameters
   // cnst is the Zlm constants - put them into the matrix ... (same code is reused in jjjpar.cpp, pointc.c)
   cnst=Matrix(0,6,-6,6);set_zlm_constants(cnst);
   so1ion=0;strcpy(moduletype,"cfield");
   alpha=0;beta=0;gamma=0;r2=0;r4=0;r6=0;gJ=0;
   double s0r=0,s0i=0,s1r=0,s1i=0,s2r=0,s2i=0;
   fgets_errchk (instr, MAXNOFCHARINLINE, cf_file);
   // strip /r (dos line feed) from line if necessary
    char *token;while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}  

   if(!(strncmp(instr,"#!MODULE=cfield ",16)==0||
        strncmp(instr,"#!MODULE=cfield\n",16)==0||
        strncmp(instr,"#!cfield ",9)==0||
        strncmp(instr,"#!cfield\n",9)==0))
        {so1ion=1;strcpy(moduletype,"so1ion");
         if(!(strncmp(instr,"#!MODULE=so1ion ",16)==0||
              strncmp(instr,"#!MODULE=so1ion\n",16)==0||
              strncmp(instr,"#!so1ion ",9)==0||
              strncmp(instr,"#!so1ion\n",9)==0)){
         fprintf(stderr,"ERROR class ionpars - file does not start with #!MODULE=cfield or #!MODULE=so1ion\n");exit(EXIT_FAILURE);
         }
        }
const char kq[]="00\022S21S20\021\022\033S32S31S30\031\032\033\044S43S42S41S40\041\042\043\044\055S54S53S52S51S50\051\052\053\054\055\066S65S64S63S62S61S60\061\062\063\064\065\066\0";
char kq3[5];kq3[4]='\0';
  
// read in lines and get IONTYPE=  and CF parameters Blm
   while(feof(cf_file)==false)
  {fgets(instr, MAXNOFCHARINLINE, cf_file);
   if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
        extract(instr,"IONTYPE",iontype,(size_t)MAXNOFCHARINLINE);
        extract(instr,"nof_electrons",nof_electrons); //MR 120127

        kq3[0]='B';for(int i=0;i<=45;++i){strncpy(kq3+1,kq+i*3,3);        
                                          extract(instr,kq3,Blm(i));}
        extract(instr,"Dx2",Blm(46));
        extract(instr,"Dy2",Blm(47));
        extract(instr,"Dz2",Blm(48));

        kq3[0]='L';for(int i=0;i<=45;++i){strncpy(kq3+1,kq+i*3,3);        
                                          extract(instr,kq3,Llm(i));}

        extract(instr,"ALPHA",alphar);
        extract(instr,"BETA",betar);
        extract(instr,"GAMMA",gammar);

        extract(instr,"GJ",gJr);
       
        extract(instr,"R2",r2r);
        extract(instr,"R4",r4r);
        extract(instr,"R6",r6r);

        extract(instr,"SIGMA0r",s0r);
        extract(instr,"SIGMA1r",s1r);
        extract(instr,"SIGMA2r",s2r);
        extract(instr,"SIGMA0i",s0i);
        extract(instr,"SIGMA1i",s1i);
        extract(instr,"SIGMA2i",s2i);
	}}


  double Jxr[31*31],Jxi[31*31],Jyr[31*31],Jyi[31*31],Jzr[31*31],Jzi[31*31];

  double  mo22sr[31*31],mo22si[31*31];
  double  mo21sr[31*31],mo21si[31*31];
  double  mo20cr[31*31],mo20ci[31*31];
  double  mo21cr[31*31],mo21ci[31*31];
  double  mo22cr[31*31],mo22ci[31*31];

  double  mo33sr[31*31],mo33si[31*31];
  double  mo32sr[31*31],mo32si[31*31];
  double  mo31sr[31*31],mo31si[31*31];
  double  mo30cr[31*31],mo30ci[31*31];
  double  mo31cr[31*31],mo31ci[31*31];
  double  mo32cr[31*31],mo32ci[31*31];
  double  mo33cr[31*31],mo33ci[31*31];

  double  mo44sr[31*31],mo44si[31*31];
  double  mo43sr[31*31],mo43si[31*31];
  double  mo42sr[31*31],mo42si[31*31];
  double  mo41sr[31*31],mo41si[31*31];
  double  mo40cr[31*31],mo40ci[31*31];
  double  mo41cr[31*31],mo41ci[31*31];
  double  mo42cr[31*31],mo42ci[31*31];
  double  mo43cr[31*31],mo43ci[31*31];
  double  mo44cr[31*31],mo44ci[31*31];

  double  mo55sr[31*31],mo55si[31*31];
  double  mo54sr[31*31],mo54si[31*31];
  double  mo53sr[31*31],mo53si[31*31];
  double  mo52sr[31*31],mo52si[31*31];
  double  mo51sr[31*31],mo51si[31*31];
  double  mo50cr[31*31],mo50ci[31*31];
  double  mo51cr[31*31],mo51ci[31*31];
  double  mo52cr[31*31],mo52ci[31*31];
  double  mo53cr[31*31],mo53ci[31*31];
  double  mo54cr[31*31],mo54ci[31*31];
  double  mo55cr[31*31],mo55ci[31*31];

  double  mo66sr[31*31],mo66si[31*31];
  double  mo65sr[31*31],mo65si[31*31];
  double  mo64sr[31*31],mo64si[31*31];
  double  mo63sr[31*31],mo63si[31*31];
  double  mo62sr[31*31],mo62si[31*31];
  double  mo61sr[31*31],mo61si[31*31];
  double  mo60cr[31*31],mo60ci[31*31];
  double  mo61cr[31*31],mo61ci[31*31];
  double  mo62cr[31*31],mo62ci[31*31];
  double  mo63cr[31*31],mo63ci[31*31];
  double  mo64cr[31*31],mo64ci[31*31];
  double  mo65cr[31*31],mo65ci[31*31];
  double  mo66cr[31*31],mo66ci[31*31];

  double  modxcr[31*31],modxci[31*31];
  double  modycr[31*31],modyci[31*31];
  double  modzcr[31*31],modzci[31*31];
    
if (pr==1) {printf("#using %s ...\n",moduletype);
           }
  
  fprintf(stderr,"# module %s ... for ion %s\n",moduletype,iontype);
  cfield_mcphasnew(iontype,Jxr,Jxi,  Jyr, Jyi, Jzr, Jzi,
  mo22sr,mo22si,
  mo21sr,mo21si,
  mo20cr,mo20ci,
  mo21cr,mo21ci,
  mo22cr,mo22ci,

  mo33sr,mo33si,
  mo32sr,mo32si,
  mo31sr,mo31si,
  mo30cr,mo30ci,
  mo31cr,mo31ci,
  mo32cr,mo32ci,
  mo33cr,mo33ci,

  mo44sr,mo44si,
  mo43sr,mo43si,
  mo42sr,mo42si,
  mo41sr,mo41si,
  mo40cr,mo40ci,
  mo41cr,mo41ci,
  mo42cr,mo42ci,
  mo43cr,mo43ci,
  mo44cr,mo44ci,

  mo55sr,mo55si,
  mo54sr,mo54si,
  mo53sr,mo53si,
  mo52sr,mo52si,
  mo51sr,mo51si,
  mo50cr,mo50ci,
  mo51cr,mo51ci,
  mo52cr,mo52ci,
  mo53cr,mo53ci,
  mo54cr,mo54ci,
  mo55cr,mo55ci,

  mo66sr,mo66si,
  mo65sr,mo65si,
  mo64sr,mo64si,
  mo63sr,mo63si,
  mo62sr,mo62si,
  mo61sr,mo61si,
  mo60cr,mo60ci,
  mo61cr,mo61ci,
  mo62cr,mo62ci,
  mo63cr,mo63ci,
  mo64cr,mo64ci,
  mo65cr,mo65ci,
  mo66cr,mo66ci,

  modxcr,modxci,
  modycr,modyci,
  modzcr,modzci,

  &dimj,&alpha,&beta,&gamma,&gJ,&r2,&r4,&r6, &nof_electrons);

if(fabs(alphar-alpha)/fabs(alphar+1)>SMALL) {fprintf(stderr,"#Warning module %s internal value for Stevens Parameter (alpha=%g) different from input file (alpha=%g), using internal value\n",moduletype,alpha,alphar);}
if(fabs(betar-beta)/fabs(betar+1)>SMALL) {fprintf(stderr,"#Warning module %s internal value for Stevens Parameter (beta=%g) different from input file (beta=%g), using internal value\n",moduletype,beta,betar);}
if(fabs(gammar-gamma)/fabs(gammar+1)>SMALL) {fprintf(stderr,"#Warning module %s internal value for Stevens Parameter (gamma=%g) different from input file (gamma=%g), using internal value\n",moduletype,gamma,gammar);}
if(fabs(gJr-gJ)/fabs(gJr+1)>SMALL) {fprintf(stderr,"#Warning module %s internal value for Lande Factor (gJ=%g) different from input file (gJ=%g), using internal value\n",moduletype,gJ,gJr);}
if(fabs(r2r-r2)/fabs(r2r+1)>SMALL) {fprintf(stderr,"#Warning module %s internal value for radial Matrix element (<r2>=%g) different from input file (<r2>=%g), using internal value\n",moduletype,r2,r2r);}
if(fabs(r4r-r4)/fabs(r4r+1)>SMALL) {fprintf(stderr,"#Warning module %s internal value for radial Matrix element (<r4>=%g) different from input file (<r4>=%g), using internal value\n",moduletype,r4,r4r);}
if(fabs(r6r-r6)/fabs(r6r+1)>SMALL) {fprintf(stderr,"#Warning module %s internal value for radial Matrix element (<r6>=%g) different from input file (<r6>=%g), using internal value\n",moduletype,r6,r6r);}

if (pr==1) printf("#end using %s\n",moduletype);

   J=((double)dimj-1)/2; //momentum quantum number
   sigma0=complex<double>(s0r,s0i);
   sigma1=complex<double>(s1r,s1i);
   sigma2=complex<double>(s2r,s2i);
if (pr==1) printf("#J=%g\n",J);

   Ja = Matrix(1,dimj,1,dimj); 
   Jaa = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jaa(i,j)=im*Jxi[30*(i-1)+j-1]+Jxr[30*(i-1)+j-1];
    if(i<j){Ja(i,j)=Jxi[30*(j-1)+i-1];}else{Ja(i,j)=Jxr[30*(i-1)+j-1];}
   }

   Jb = Matrix(1,dimj,1,dimj); 
   Jbb = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jbb(i,j)=im*Jyi[30*(i-1)+j-1]+Jyr[30*(i-1)+j-1];
    if(i<j){Jb(i,j)=Jyi[30*(j-1)+i-1];}else{Jb(i,j)=Jyr[30*(i-1)+j-1];}
   }

   Jc = Matrix(1,dimj,1,dimj); 
   Jcc = ComplexMatrix(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jcc(i,j)=im*Jzi[30*(i-1)+j-1]+Jzr[30*(i-1)+j-1];
    if(i<j){Jc(i,j)=Jzi[30*(j-1)+i-1];}else{Jc(i,j)=Jzr[30*(i-1)+j-1];}
   }

//---------------------------------------------------------------------------
	Ri= new ComplexMatrix * [1+NOF_RIXS_MATRICES];
         for(i=1;i<=NOF_RIXS_MATRICES;++i)Ri[i]=new ComplexMatrix(1,dimj,1,dimj);
      // here fill the matrices Ri[1...9] with the 11 12 13 21 22 23 31 32 33
      // matrices of the RIXS scattering operator R
       complex<double> f1=sigma1/J;
       complex<double> f2=sigma2/(J*(2*J-1));
      // SIGMA0 contribution (haverkort PRL 105 (2010) 167404 equation (8)
       (*Ri[1])=sigma0-f2*0.6666*(Jaa*Jaa+Jbb*Jbb+Jcc*Jcc);
                        (*Ri[5])=(*Ri[1]);
                                         (*Ri[9])=(*Ri[1]);
      // SIGMA1 contribution (haverkort PRL 105 (2010) 167404 equation (8)
                        (*Ri[2])=f1*Jcc;(*Ri[3])=-f1*Jbb;
       (*Ri[4])=-f1*Jcc;                 (*Ri[6])=f1*Jaa;
       (*Ri[7])=f1*Jbb; (*Ri[8])=-f1*Jaa;
      // SIGMA2 contribution (haverkort PRL 105 (2010) 167404 equation (8)
       (*Ri[1])+=2.0*f2*Jaa*Jaa; (*Ri[2])+=f2*(Jaa*Jbb+Jbb*Jaa); (*Ri[3])+=f2*(Jaa*Jcc+Jcc*Jaa);
       (*Ri[4])+=f2*(Jbb*Jcc+Jcc*Jbb); (*Ri[5])+=2.0*f2*Jbb*Jbb; (*Ri[6])+=f2*Jaa;
       (*Ri[7])+=f2*(Jcc*Jaa+Jaa*Jcc); (*Ri[8])+=f2*(Jbb*Jcc+Jcc*Jbb); (*Ri[9])+=2.0*f2*Jcc*Jcc;

   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 
   
   for(i=1;i<=NOF_OLM_MATRICES;++i)   
    {   Olm[i]= new Matrix(1,dimj,1,dimj); 
 // define memory for all matrices 
        OOlm [i] = new ComplexMatrix(1,dimj,1,dimj); 
    }   
    
    
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {(*OOlm[1])(i,j)=im*mo22si[30*(i-1)+j-1]+mo22sr[30*(i-1)+j-1];
if(i<j){(*Olm[1])(i,j)=mo22si[30*(j-1)+i-1];}else{(*Olm[1])(i,j)=mo22sr[30*(i-1)+j-1];}
    (*OOlm[2])(i,j)=im*mo21si[30*(i-1)+j-1]+mo21sr[30*(i-1)+j-1];
if(i<j){(*Olm[2])(i,j)=mo21si[30*(j-1)+i-1];}else{(*Olm[2])(i,j)=mo21sr[30*(i-1)+j-1];}
    (*OOlm[3])(i,j)=im*mo20ci[30*(i-1)+j-1]+mo20cr[30*(i-1)+j-1];
if(i<j){(*Olm[3])(i,j)=mo20ci[30*(j-1)+i-1];}else{(*Olm[3])(i,j)=mo20cr[30*(i-1)+j-1];}
    (*OOlm[4])(i,j)=im*mo21ci[30*(i-1)+j-1]+mo21cr[30*(i-1)+j-1];
if(i<j){(*Olm[4])(i,j)=mo21ci[30*(j-1)+i-1];}else{(*Olm[4])(i,j)=mo21cr[30*(i-1)+j-1];}
    (*OOlm[5])(i,j)=im*mo22ci[30*(i-1)+j-1]+mo22cr[30*(i-1)+j-1];
if(i<j){(*Olm[5])(i,j)=mo22ci[30*(j-1)+i-1];}else{(*Olm[5])(i,j)=mo22cr[30*(i-1)+j-1];}
    
    (*OOlm[6])(i,j)=im*mo33si[30*(i-1)+j-1]+mo33sr[30*(i-1)+j-1];
if(i<j){(*Olm[6])(i,j)=mo33si[30*(j-1)+i-1];}else{(*Olm[6])(i,j)=mo33sr[30*(i-1)+j-1];}
    (*OOlm[7])(i,j)=im*mo32si[30*(i-1)+j-1]+mo32sr[30*(i-1)+j-1];
if(i<j){(*Olm[7])(i,j)=mo32si[30*(j-1)+i-1];}else{(*Olm[7])(i,j)=mo32sr[30*(i-1)+j-1];}
    (*OOlm[8])(i,j)=im*mo31si[30*(i-1)+j-1]+mo31sr[30*(i-1)+j-1];
if(i<j){(*Olm[8])(i,j)=mo31si[30*(j-1)+i-1];}else{(*Olm[8])(i,j)=mo31sr[30*(i-1)+j-1];}
    (*OOlm[9])(i,j)=im*mo30ci[30*(i-1)+j-1]+mo30cr[30*(i-1)+j-1];
if(i<j){(*Olm[9])(i,j)=mo30ci[30*(j-1)+i-1];}else{(*Olm[9])(i,j)=mo30cr[30*(i-1)+j-1];}
    (*OOlm[10])(i,j)=im*mo31ci[30*(i-1)+j-1]+mo31cr[30*(i-1)+j-1];
if(i<j){(*Olm[10])(i,j)=mo31ci[30*(j-1)+i-1];}else{(*Olm[10])(i,j)=mo31cr[30*(i-1)+j-1];}
    (*OOlm[11])(i,j)=im*mo32ci[30*(i-1)+j-1]+mo32cr[30*(i-1)+j-1];
if(i<j){(*Olm[11])(i,j)=mo32ci[30*(j-1)+i-1];}else{(*Olm[11])(i,j)=mo32cr[30*(i-1)+j-1];}
    (*OOlm[12])(i,j)=im*mo33ci[30*(i-1)+j-1]+mo33cr[30*(i-1)+j-1];
if(i<j){(*Olm[12])(i,j)=mo33ci[30*(j-1)+i-1];}else{(*Olm[12])(i,j)=mo33cr[30*(i-1)+j-1];}
    
    (*OOlm[13])(i,j)=im*mo44si[30*(i-1)+j-1]+mo44sr[30*(i-1)+j-1];
if(i<j){(*Olm[13])(i,j)=mo44si[30*(j-1)+i-1];}else{(*Olm[13])(i,j)=mo44sr[30*(i-1)+j-1];}
    (*OOlm[14])(i,j)=im*mo43si[30*(i-1)+j-1]+mo43sr[30*(i-1)+j-1];
if(i<j){(*Olm[14])(i,j)=mo43si[30*(j-1)+i-1];}else{(*Olm[14])(i,j)=mo43sr[30*(i-1)+j-1];}
    (*OOlm[15])(i,j)=im*mo42si[30*(i-1)+j-1]+mo42sr[30*(i-1)+j-1];
if(i<j){(*Olm[15])(i,j)=mo42si[30*(j-1)+i-1];}else{(*Olm[15])(i,j)=mo42sr[30*(i-1)+j-1];}
    (*OOlm[16])(i,j)=im*mo41si[30*(i-1)+j-1]+mo41sr[30*(i-1)+j-1];
if(i<j){(*Olm[16])(i,j)=mo41si[30*(j-1)+i-1];}else{(*Olm[16])(i,j)=mo41sr[30*(i-1)+j-1];}
    (*OOlm[17])(i,j)=im*mo40ci[30*(i-1)+j-1]+mo40cr[30*(i-1)+j-1];
if(i<j){(*Olm[17])(i,j)=mo40ci[30*(j-1)+i-1];}else{(*Olm[17])(i,j)=mo40cr[30*(i-1)+j-1];}
    (*OOlm[18])(i,j)=im*mo41ci[30*(i-1)+j-1]+mo41cr[30*(i-1)+j-1];
if(i<j){(*Olm[18])(i,j)=mo41ci[30*(j-1)+i-1];}else{(*Olm[18])(i,j)=mo41cr[30*(i-1)+j-1];}
    (*OOlm[19])(i,j)=im*mo42ci[30*(i-1)+j-1]+mo42cr[30*(i-1)+j-1];
if(i<j){(*Olm[19])(i,j)=mo42ci[30*(j-1)+i-1];}else{(*Olm[19])(i,j)=mo42cr[30*(i-1)+j-1];}
    (*OOlm[20])(i,j)=im*mo43ci[30*(i-1)+j-1]+mo43cr[30*(i-1)+j-1];
if(i<j){(*Olm[20])(i,j)=mo43ci[30*(j-1)+i-1];}else{(*Olm[20])(i,j)=mo43cr[30*(i-1)+j-1];}
    (*OOlm[21])(i,j)=im*mo44ci[30*(i-1)+j-1]+mo44cr[30*(i-1)+j-1];
if(i<j){(*Olm[21])(i,j)=mo44ci[30*(j-1)+i-1];}else{(*Olm[21])(i,j)=mo44cr[30*(i-1)+j-1];}
    
    (*OOlm[22])(i,j)=im*mo55si[30*(i-1)+j-1]+mo55sr[30*(i-1)+j-1];
if(i<j){(*Olm[22])(i,j)=mo55si[30*(j-1)+i-1];}else{(*Olm[22])(i,j)=mo55sr[30*(i-1)+j-1];}
    (*OOlm[23])(i,j)=im*mo54si[30*(i-1)+j-1]+mo54sr[30*(i-1)+j-1];
if(i<j){(*Olm[23])(i,j)=mo54si[30*(j-1)+i-1];}else{(*Olm[23])(i,j)=mo54sr[30*(i-1)+j-1];}
    (*OOlm[24])(i,j)=im*mo53si[30*(i-1)+j-1]+mo53sr[30*(i-1)+j-1];
if(i<j){(*Olm[24])(i,j)=mo53si[30*(j-1)+i-1];}else{(*Olm[24])(i,j)=mo53sr[30*(i-1)+j-1];}
    (*OOlm[25])(i,j)=im*mo52si[30*(i-1)+j-1]+mo52sr[30*(i-1)+j-1];
if(i<j){(*Olm[25])(i,j)=mo52si[30*(j-1)+i-1];}else{(*Olm[25])(i,j)=mo52sr[30*(i-1)+j-1];}
    (*OOlm[26])(i,j)=im*mo51si[30*(i-1)+j-1]+mo51sr[30*(i-1)+j-1];
if(i<j){(*Olm[26])(i,j)=mo51si[30*(j-1)+i-1];}else{(*Olm[26])(i,j)=mo51sr[30*(i-1)+j-1];}
    (*OOlm[27])(i,j)=im*mo50ci[30*(i-1)+j-1]+mo50cr[30*(i-1)+j-1];
if(i<j){(*Olm[27])(i,j)=mo50ci[30*(j-1)+i-1];}else{(*Olm[27])(i,j)=mo50cr[30*(i-1)+j-1];}
    (*OOlm[28])(i,j)=im*mo51ci[30*(i-1)+j-1]+mo51cr[30*(i-1)+j-1];
if(i<j){(*Olm[28])(i,j)=mo51ci[30*(j-1)+i-1];}else{(*Olm[28])(i,j)=mo51cr[30*(i-1)+j-1];}
    (*OOlm[29])(i,j)=im*mo52ci[30*(i-1)+j-1]+mo52cr[30*(i-1)+j-1];
if(i<j){(*Olm[29])(i,j)=mo52ci[30*(j-1)+i-1];}else{(*Olm[29])(i,j)=mo52cr[30*(i-1)+j-1];}
    (*OOlm[30])(i,j)=im*mo53ci[30*(i-1)+j-1]+mo53cr[30*(i-1)+j-1];
if(i<j){(*Olm[30])(i,j)=mo53ci[30*(j-1)+i-1];}else{(*Olm[30])(i,j)=mo53cr[30*(i-1)+j-1];}
    (*OOlm[31])(i,j)=im*mo54ci[30*(i-1)+j-1]+mo54cr[30*(i-1)+j-1];
if(i<j){(*Olm[31])(i,j)=mo54ci[30*(j-1)+i-1];}else{(*Olm[31])(i,j)=mo54cr[30*(i-1)+j-1];}
    (*OOlm[32])(i,j)=im*mo55ci[30*(i-1)+j-1]+mo55cr[30*(i-1)+j-1];
if(i<j){(*Olm[32])(i,j)=mo55ci[30*(j-1)+i-1];}else{(*Olm[32])(i,j)=mo55cr[30*(i-1)+j-1];}

    (*OOlm[33])(i,j)=im*mo66si[30*(i-1)+j-1]+mo66sr[30*(i-1)+j-1];
if(i<j){(*Olm[33])(i,j)=mo66si[30*(j-1)+i-1];}else{(*Olm[33])(i,j)=mo66sr[30*(i-1)+j-1];}
    (*OOlm[34])(i,j)=im*mo65si[30*(i-1)+j-1]+mo65sr[30*(i-1)+j-1];
if(i<j){(*Olm[34])(i,j)=mo65si[30*(j-1)+i-1];}else{(*Olm[34])(i,j)=mo65sr[30*(i-1)+j-1];}
    (*OOlm[35])(i,j)=im*mo64si[30*(i-1)+j-1]+mo64sr[30*(i-1)+j-1];
if(i<j){(*Olm[35])(i,j)=mo64si[30*(j-1)+i-1];}else{(*Olm[35])(i,j)=mo64sr[30*(i-1)+j-1];}
    (*OOlm[36])(i,j)=im*mo63si[30*(i-1)+j-1]+mo63sr[30*(i-1)+j-1];
if(i<j){(*Olm[36])(i,j)=mo63si[30*(j-1)+i-1];}else{(*Olm[36])(i,j)=mo63sr[30*(i-1)+j-1];}
    (*OOlm[37])(i,j)=im*mo62si[30*(i-1)+j-1]+mo62sr[30*(i-1)+j-1];
if(i<j){(*Olm[37])(i,j)=mo62si[30*(j-1)+i-1];}else{(*Olm[37])(i,j)=mo62sr[30*(i-1)+j-1];}
    (*OOlm[38])(i,j)=im*mo61si[30*(i-1)+j-1]+mo61sr[30*(i-1)+j-1];
if(i<j){(*Olm[38])(i,j)=mo61si[30*(j-1)+i-1];}else{(*Olm[38])(i,j)=mo61sr[30*(i-1)+j-1];}
    (*OOlm[39])(i,j)=im*mo60ci[30*(i-1)+j-1]+mo60cr[30*(i-1)+j-1];
if(i<j){(*Olm[39])(i,j)=mo60ci[30*(j-1)+i-1];}else{(*Olm[39])(i,j)=mo60cr[30*(i-1)+j-1];}
    (*OOlm[40])(i,j)=im*mo61ci[30*(i-1)+j-1]+mo61cr[30*(i-1)+j-1];
if(i<j){(*Olm[40])(i,j)=mo61ci[30*(j-1)+i-1];}else{(*Olm[40])(i,j)=mo61cr[30*(i-1)+j-1];}
    (*OOlm[41])(i,j)=im*mo62ci[30*(i-1)+j-1]+mo62cr[30*(i-1)+j-1];
if(i<j){(*Olm[41])(i,j)=mo62ci[30*(j-1)+i-1];}else{(*Olm[41])(i,j)=mo62cr[30*(i-1)+j-1];}
    (*OOlm[42])(i,j)=im*mo63ci[30*(i-1)+j-1]+mo63cr[30*(i-1)+j-1];
if(i<j){(*Olm[42])(i,j)=mo63ci[30*(j-1)+i-1];}else{(*Olm[42])(i,j)=mo63cr[30*(i-1)+j-1];}
    (*OOlm[43])(i,j)=im*mo64ci[30*(i-1)+j-1]+mo64cr[30*(i-1)+j-1];
if(i<j){(*Olm[43])(i,j)=mo64ci[30*(j-1)+i-1];}else{(*Olm[43])(i,j)=mo64cr[30*(i-1)+j-1];}
    (*OOlm[44])(i,j)=im*mo65ci[30*(i-1)+j-1]+mo65cr[30*(i-1)+j-1];
if(i<j){(*Olm[44])(i,j)=mo65ci[30*(j-1)+i-1];}else{(*Olm[44])(i,j)=mo65cr[30*(i-1)+j-1];}
    (*OOlm[45])(i,j)=im*mo66ci[30*(i-1)+j-1]+mo66cr[30*(i-1)+j-1];
if(i<j){(*Olm[45])(i,j)=mo66ci[30*(j-1)+i-1];}else{(*Olm[45])(i,j)=mo66cr[30*(i-1)+j-1];}

    (*OOlm[46])(i,j)=im*modxci[30*(i-1)+j-1]+modxcr[30*(i-1)+j-1];
if(i<j){(*Olm[46])(i,j)=modxci[30*(j-1)+i-1];}else{(*Olm[46])(i,j)=modxcr[30*(i-1)+j-1];}
    (*OOlm[47])(i,j)=im*modyci[30*(i-1)+j-1]+modycr[30*(i-1)+j-1];
if(i<j){(*Olm[47])(i,j)=modyci[30*(j-1)+i-1];}else{(*Olm[47])(i,j)=modycr[30*(i-1)+j-1];}
    (*OOlm[48])(i,j)=im*modzci[30*(i-1)+j-1]+modzcr[30*(i-1)+j-1];
if(i<j){(*Olm[48])(i,j)=modzci[30*(j-1)+i-1];}else{(*Olm[48])(i,j)=modzcr[30*(i-1)+j-1];}
    
   }}
//printf("%g\n",mo54sr[1][1]);

// ------------------------------------------------------------
// here transform the Llm (if present) to Blm ...
Vector thetaJ(0,6);thetaJ(0)=nof_electrons;thetaJ(2)=alpha;thetaJ(4)=beta;thetaJ(6)=gamma;

   fprintf(stderr,"#crystal field parameters:\n");  
   const char lm[]="B00 B22SB21SB20 B21 B22 B33SB32SB31SB30 B31 B32 B33 B44SB43SB42SB41SB40 B41 B42 B43 B44 B55SB54SB53SB52SB51SB50 B51 B52 B53 B54 B55 B66SB65SB64SB63SB62SB61SB60 B61 B62 B63 B64 B65 B66 Dx2 Dy2 Dz2 ";
   char lm4[5];lm4[4]='\0';
   for(i=0;i<=48;++i){strncpy(lm4,lm+i*4,4);l=lm4[1]-48;m=lm4[2]-48;if(lm4[3]=='S'){m=-m;}
                     if(i<=45&&Llm(i)!=0){if(l==3||l==5){lm4[0]='L';fprintf(stderr,"#Error internal module %s: wybourne parameter %s is not implemented\n",moduletype,lm4);
                                                  exit(EXIT_FAILURE);}
                                  double Blmcalc=Llm(i)*cnst(l,m)*sqrt(4.0*PI/(2*l+1))*thetaJ(l);if(m!=0){Blmcalc*=sqrt(2.0);}
                                  if((Blm(i)!=0)&(fabs(Blm(i)-Blmcalc)/(fabs(Blmcalc)+1e-14)>0.001)){fprintf(stderr,"#Warning internal module %s - reading %s=%12.6g meV is ignored, because Wybourne Parameter Llm=%12.6g meV does not correspond !\n Will use Blm=%12.6g calculated from Llm.\npresse enter to continue\n",moduletype,lm4,Blm(i),Llm(i),Blmcalc);getchar();}
                                  Blm(i)=Blmcalc;// here set the Blm as calculated from the Llm
                                  }
                     if(Blm(i)!=0){fprintf(stderr,"#! %s=%12.6g meV ",lm4,Blm(i));
                                   if(i<=45){if((l!=3)&(l!=5)){Llm(i)=Blm(i)/thetaJ(l)/cnst(l,m)/sqrt(4.0*PI/(2*l+1));if(m!=0){Llm(i)/=sqrt(2.0);}
                                                 lm4[0]='L';fprintf(stderr,"<-> %s=%12.6g meV",lm4,Llm(i));}
                                                else
                                                {lm4[0]='L';fprintf(stderr,"<-> %s=Wybourne parameter not implemented, ",lm4);}}
                                   fprintf(stderr,"\n");  
                                  }
                     }
// here set the crystal field Hamiltonian Matrix
  Hcf= Matrix(1,dimj,1,dimj); 
  Hcf=0;

   if(Hcf==(double)0.0){for(l=1;l<=48;++l){Hcf+=Blm(l)*(*Olm[l]);}}

if(so1ion==0)
  {//ATTENTION FOR cfield the AXES xyz are parallel to cab
   Matrix dummy(1,dimj,1,dimj); dummy=Jb;Jb=Jc;Jc=Ja;Ja=dummy;
   ComplexMatrix dummyc(1,dimj,1,dimj);dummyc=Jbb;Jbb=Jcc;Jcc=Jaa;Jaa=dummyc;

  if (pr==1) {printf("#Axis Convention using cfield as a module:  a||y b||z  c||x\n");
  printf("#xyz .... Coordinate system of the crystal field parameters used in cfield\n");
  printf("#abc .... Crystal axes\n");
  printf("#The interactions are described by the  PKQ Operators defined in cfield\n");
  printf("#O11(s) .... Ia=Jy\n#O10(c) .... Ib=Jz\n#O11(c) .... Ic=Jx\n");
  printf("#O22(s) .... Id\n#O21(s) .... Ie\n#O20(c) .... If\n#O21(c) .... Ig\nO22(c) .... Ih\n");
  printf("#O33(s) .... Ii\n#O32(s) .... Ij\n#O31(s) .... Ik\n#O30(c) .... Il\n#O31(c) .... Im\n");
  printf("# etc ... 45 moments up to l<=6\n\n");
              }
  }
  else
  {//ATTENTION FOR so1ion the AXES xyz are parallel to abc
  if (pr==1) {printf("#Axis Convention using so1ion as a module:  a||x b||y  c||z\n");
  printf("#xyz .... Coordinate system of the crystal field parameters used in so1ion\n");
  printf("#abc .... Crystal axes\n");
  printf("#The interactions are described by the  PKQ Operators defined in so1ion\n");
  printf("#O11(s) .... Ia=Jx\n#O10(c) .... Ib=Jy\n#O11(c) .... Ic=Jz\n");
  printf("#O22(s) .... Id\n#O21(s) .... Ie\n#O20(c) .... If\n#O21(c) .... Ig\nO22(c) .... Ih\n");
  printf("#O33(s) .... Ii\n#O32(s) .... Ij\n#O31(s) .... Ik\n#O30(c) .... Il\n#O31(c) .... Im\n");
  printf("# etc ... 45 moments up to l<=6\n\n");
             }
  }
 pr=0;
 // now fill interaction operators In with values
  In = new ComplexMatrix * [1+IONPARS_MAXNOFCOMPONENTS+3+NOF_OLM_MATRICES+NOF_RIXS_MATRICES];
  for(i=1;i<=IONPARS_MAXNOFCOMPONENTS;++i)In[i]= new ComplexMatrix(1,dimj,1,dimj);
  // standard operator sequence I1,....,I51
  // module so1ion: Jx Jy Jz O22S O21S O20 O21 O22 O33S O32S .... O66 Jx^2 Jy^2 Jz^2
  (*In[1])=Jaa;(*In[1])=Jbb;(*In[1])=Jcc;
  for(i=4;i<=IONPARS_MAXNOFCOMPONENTS;++i){(*In[i])=(*OOlm[i-3]);}
}


void ionpars::save(FILE * file) // save ion parameters to file 
{
  fprintf(file,"#-----------\nIONTYPE=%s\n#-----------\n\n",iontype);

  if(abs(Blm)>1e-10) {fprintf(file,"#--------------------------------------------------------------------------\n");
                      if(so1ion==0){
                      fprintf(file,"# Crystal Field parameters in Stevens Notation (coordinate system yzx||abc)\n");
                      }else{
                      fprintf(file,"# Crystal Field parameters in Stevens Notation (coordinate system xyz||abc)\n");
                      }
                      fprintf(file,"#--------------------------------------------------------------------------\n");
                      savBlm(file);fprintf(file,"\n");
                     }
  if(abs(Llm)>1e-10) {fprintf(file,"#---------------------------------------------------------------------------\n");
                       if(so1ion==0){
                      fprintf(file,"# Crystal Field parameters in Wybourne Notation (coordinate system yzx||abc)\n");
                      }else{
                      fprintf(file,"# Crystal Field parameters in Wybourne Notation (coordinate system xyz||abc)\n");
                      }
                      fprintf(file,"#---------------------------------------------------------------------------\n");
                      savLlm(file);fprintf(file,"\n");
                     }
   if(alpha*alpha+beta*beta+gamma*gamma>1e-10){fprintf(file,"#----------------\n# Stevens Factors\n#----------------\nALPHA=%g\nBETA=%g\nGAMMA=%g\n\n",alpha,beta,gamma);}
   if(r2+r4+r6>1e-10){fprintf(file,"#---------------------------------------------------------\n");
                      fprintf(file,"# Radial Matrix Elements (e.g. Abragam Bleaney 1971 p 399)\n");
                      fprintf(file,"#---------------------------------------------------------\n");
   if(r2>1e-10){fprintf(file,"#<r^2> in units of a0^2 a0=0.5292 Angstroem\nR2=%g\n",r2);}
   if(r4>1e-10){fprintf(file,"#<r^4> in units of a0^4 a0=0.5292 Angstroem\nR4=%g\n",r4);}
   if(r6>1e-10){fprintf(file,"#<r^6> in units of a0^6 a0=0.5292 Angstroem\nR6=%g\n",r6);}
                   fprintf(file,"\n");
                      }
}

void ionpars::savBlm(FILE * outfile)
{fprintf(outfile,"units=meV\n");
    char lm3[4];lm3[3]='\0';
   const char lm[]="00 22S21S20 21 22 33S32S31S30 31 32 33 44S43S42S41S40 41 42 43 44 55S54S53S52S51S50 51 52 53 54 55 66S65S64S63S62S61S60 61 62 63 64 65 66 ";
   for(int i=0;i<=45;++i){strncpy(lm3,lm+i*3,3);
   if(Blm(i)!=0){fprintf(outfile,"B%s=%g\n",lm3,myround(Blm(i)));}
                  }
   if(Blm(46)!=0){fprintf(outfile,"Dx2=%g\n",myround(Blm(46)));}
   if(Blm(47)!=0){fprintf(outfile,"Dy2=%g\n",myround(Blm(47)));}
   if(Blm(48)!=0){fprintf(outfile,"Dz2=%g\n",myround(Blm(48)));}
}

void ionpars::savLlm(FILE * outfile)
{fprintf(outfile,"units=meV\n");
   char lm3[4];lm3[3]='\0';
   const char lm[]="00 22S21S20 21 22 33S32S31S30 31 32 33 44S43S42S41S40 41 42 43 44 55S54S53S52S51S50 51 52 53 54 55 66S65S64S63S62S61S60 61 62 63 64 65 66 ";
   for(int i=0;i<=45;++i){strncpy(lm3,lm+i*3,3);
   if(Llm(i)!=0){fprintf(outfile,"L%s=%g\n",lm3,myround(Llm(i)));}
                         }
}

//------------------------------------------------------------------------------------------------
// ROUTINE CFIELD Icalc for full crystal field + higher order interactions
//------------------------------------------------------------------------------------------------
Vector & ionpars::cfield(double & T, Vector &  gjmbHxc,Vector & Hext, double & lnZs, double & U, ComplexMatrix & ests)
{//ABC not used !!!
    /*on input
    T		temperature[K]
    gJmbH	vector of effective field [meV]
    gJ          Lande factor
    ABC         single ion parameter values (A, B, C corresponding to <+|Ja|->,<-|Jb|->,<+|Jc|->/i
  on output    
    J		single ion momentum vector <J> (if T>0 thermal exp value <J>T 
                                                if T<0 the program asks for w_n and calculates
						       exp value <J>=sum_n w_n <n|J|n>
						       
    Z		single ion partition function
    U		single ion magnetic energy
*/
// check dimensions of vector


if(gjmbHxc.Hi()>48)
   {fprintf(stderr,"Error internal module cfield: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}
static Vector JJ(1,gjmbHxc.Hi());
cfieldJJ(JJ, T,  gjmbHxc,Hext, lnZs, U, ests);
return JJ;
}


void ionpars::cfieldJJ(Vector & JJ,double & T, Vector &  gjmbHxc,Vector & Hext, double & lnZs, double & U, ComplexMatrix & /*ests*/)
{   /*on input
    T		temperature[K]
    gJmbHxc	vector of exchange field [meV]
    Hext        external magnetic field [T]
  on output    
    JJ		single ion momentum vector <J> (if T>0 thermal exp value <J>T 
                                                if T<0 the program asks for w_n and calculates
						       exp value <J>=sum_n w_n <n|J|n>
						       
    Z		single ion partition function
    U		single ion magnetic energy
*/
Vector gjmbH(1,gjmbHxc.Hi());
gjmbH=gjmbHxc;
gjmbH(1)+=gJ*MU_B*Hext(1);
gjmbH(2)+=gJ*MU_B*Hext(2);
gjmbH(3)+=gJ*MU_B*Hext(3);
// check dimensions of vector
if(gjmbH.Hi()>48)
   {fprintf(stderr,"Error internal module cfield/so1ion: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

//  Driver routine to compute the  eigenvalues and normalized eigenvectors 
//  of a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. The real  and imaginary
//  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
// 			   int sort, int maxiter)


   // setup hamiltonian
   int dj,i,j;
 
   dj=Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);
//   Matrix Tam(1,dj,1,dj); // transformed Hamiltonian
   ComplexMatrix z(1,dj,1,dj);
   ComplexMatrix za(1,dj,1,dj);
   ComplexMatrix zb(1,dj,1,dj);
   ComplexMatrix zc(1,dj,1,dj);
   ComplexMatrix zolm(1,dj,1,dj);    
  

/*  int i1,j1; //printout matrix
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) printf ("%4.6g ",Ja(i1,j1));
    printf("\n");
    }
    printf ("\nH=%g %g %g\n",gjmbH(1),gjmbH(2),gjmbH(3));*/

   Ham=Hcf-gjmbH(1)*Ja-gjmbH(2)*Jb-gjmbH(3)*Jc;
 

// here the zeeman term is extended for multipolar fields
   for(j=4;j<=JJ.Hi();++j){Ham-=gjmbH(j)*(*Olm[j-3]);}

   // use old eigenstates ests to transform matrix to nearly diagonal form ... however we deleted this because it needs more time to transform than to solve the eigenvalue problem
/*   Tam=0;
   for(i=1;i<=dj;++i){for(j=1;j<=dj;++j){
   if(i<j){for(k=1;k<=dj;++k){for(l=1;l<=dj;++l){
           if(k<l){hkl=Ham(l,k);mukl=-Ham(k,l);}else{hkl=Ham(k,l);if(k==l){mukl=0;}else{mukl=Ham(l,k);}}           
           Tam(i,j)-=-imag(ests(k,i))*hkl*real(ests(l,j))+imag(ests(k,i))*mukl*imag(ests(l,j))+real(ests(k,i))*mukl*real(ests(l,j))+real(ests(k,i))*hkl*imag(ests(l,j));    
           }}
          }
   else   {for(k=1;k<=dj;++k){for(l=1;l<=dj;++l){
           if(k<l){hkl=Ham(l,k);mukl=-Ham(k,l);}else{hkl=Ham(k,l);if(k==l){mukl=0;}else{mukl=Ham(l,k);}}            
           Tam(i,j)+=real(ests(k,i))*hkl*real(ests(l,j))-real(ests(k,i))*mukl*imag(ests(l,j))+imag(ests(k,i))*mukl*real(ests(l,j))+imag(ests(k,i))*hkl*imag(ests(l,j));
           }}
          }
   }}         
  */  
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=0;int maxiter=1000000;
   if (T<0) sort=1;

   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);

   // calculate Z and wn (occupation probability)
     Vector wn(1,dj);
     double x,y;
     x=Min(En);
     double Zs;

     if (T>0)
     { for (i=1;i<=dj;++i)
       {if ((y=(En(i)-x)/KB/T)<600) wn[i]=exp(-y); 
        else wn[i]=0.0;
       }
       Zs=Sum(wn);wn/=Zs;
 
       lnZs=log(Zs)-x/KB/T;
     } 
     else
     { printf ("Temperature T<0: please choose probability distribution of states by hand\n");
                         printf ("Number   Energy     Excitation Energy\n");
     for (i=1;i<=dj;++i) printf ("%i    %4.4g meV   %4.4g meV\n",i,En(i),En(i)-x);
     char instr[MAXNOFCHARINLINE];
     for (i=1;i<=dj;++i)
      {printf("eigenstate %i: %4.4g meV %4.4g meV  - please enter probability w(%i):",i,En(i),En(i)-x,i);
       fgets(instr, MAXNOFCHARINLINE, stdin);
 
       wn(i)=strtod(instr,NULL);
      }
       Zs=Sum(wn);wn/=Zs;
 
       lnZs=log(Zs);
                         printf ("\n\nNumber   Energy     Excitation Energy   Probability\n");
     for (i=1;i<=dj;++i) printf ("%i    %4.4g meV   %4.4g meV %4.4g  \n",i,En(i),En(i)-x,wn(i));
     }

   // calculate U
     U=En*wn;
   // calculate <Ja>,<Jb>,<Jc>
//     z=ComplexMatrix(zr,zi);
//     z=ests(1,dj,1,dj)*z; // transform to original eigenstates ... however we deleted this because it needs more time to transform than to solve the eigenvalue problem
//     ests(1,dj,1,dj)=z;
//     for (i=1;i<=dj;++i) {ests(0,i)=complex <double> (En(i),wn(i));}
//     myPrintComplexMat(stdout,ests);     
//     myPrintComplexMat(stdout,z);

//   za=Jaa*z;zb=Jbb*z;zc=Jcc*z;
//     opZcol (i,za,Ja,zr,zi);JJ[1]+=wn(i)*real(z.Column(i)*za.Column(i));
//     opZcol (i,zb,Jb,zr,zi);JJ[2]+=wn(i)*real(z.Column(i)*zb.Column(i));
//     opZcol (i,zc,Jc,zr,zi);JJ[3]+=wn(i)*real(z.Column(i)*zc.Column(i));
//     zolm=(*OOlm[j-3])*z;
//    for (i=1;i<=dj;++i) JJ[j]+=wn(i)*real(z.Column(i)*zolm.Column(i));

     JJ=0;
//    ComplexVector ddd;
    for (i=1;i<=dj;++i)
    { if(wn(i)>1e-5){
  //     JJ[1]+=wn(i)*real(z.Column(i)*za.Column(i));
  //     JJ[2]+=wn(i)*real(z.Column(i)*zb.Column(i));
  //     JJ[3]+=wn(i)*real(z.Column(i)*zc.Column(i));
     JJ[1]+=wn(i)*matelr(i,i,zr,zi,Ja);
     JJ[2]+=wn(i)*matelr(i,i,zr,zi,Jb);
     JJ[3]+=wn(i)*matelr(i,i,zr,zi,Jc);   
     // here the expectation values of the multipolar moments are calculated
       for(j=4;j<=JJ.Hi();++j)JJ[j]+=wn(i)*matelr(i,i,zr,zi,(*Olm[j-3]));
                    }
    }

}
/**************************************************************************/
void ionpars::cfeigenstates(ComplexMatrix *eigenstates,Vector &  gjmbHxc,Vector & Hext, double & T)
{   /*on input
    gJmbH	vector of effective field [meV]
      on output
    Matrix containing the eigenvalues and eigenfunctions of the crystalfield problem
    eigenvalues ares stored as real part of row zero
    boltzmann population numbers are stored as imaginary part of row zero
*/
Vector gjmbH(1,gjmbHxc.Hi());
gjmbH=gjmbHxc;
gjmbH(1)+=gJ*MU_B*Hext(1);
gjmbH(2)+=gJ*MU_B*Hext(2);
gjmbH(3)+=gJ*MU_B*Hext(3);
// check dimensions of vector
if(gjmbH.Hi()>48)
   {fprintf(stderr,"Error internal module cfield: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

//  Driver routine to compute the  eigenvalues and normalized eigenvectors 
//  of a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. The real  and imaginary
//  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
// 			   int sort, int maxiter)
   // setup hamiltonian
   int dj,i,j;
//   complex <double> imag(0,1);
   dj=Hcf.Rhi();
// static ComplexMatrix eigenstates(0,dj,1,dj);
   (*eigenstates) = ComplexMatrix(0,dj,1,dj);
   Matrix Ham(1,dj,1,dj);
   ComplexMatrix z(1,dj,1,dj);
   ComplexMatrix za(1,dj,1,dj);
   ComplexMatrix zb(1,dj,1,dj);
   ComplexMatrix zc(1,dj,1,dj);
   ComplexMatrix zolm(1,dj,1,dj);    

 
   Ham=Hcf-gjmbH(1)*Ja-gjmbH(2)*Jb-gjmbH(3)*Jc;

/* myPrintComplexMat(stdout,Jaa);
 myPrintComplexMat(stdout,Jbb);
 myPrintComplexMat(stdout,Jcc);*/

   for(j=4;j<=gjmbH.Hi();++j){Ham-=gjmbH(j)*(*Olm[j-3]);}

/*   int i1,j1; //printout matrix
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) printf ("%4.6g ",(*Olm[j])(i1,j1));
    printf ("\n");
    }*/
      
    
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=1;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);

   for(i=1;i<=dj;++i){//eigenstates(0,i)=complex <double> (En(i),0);
     for(j=1;j<=dj;++j){(*eigenstates)(i,j)=complex <double> (zr(i,j),zi(i,j));
   }}
    //calculate partition sum
     double zz=0;double KBT,E0;KBT=T*KB;E0=En(1);
      for(j=1;j<=dj;++j){zz+=exp(-((En(j)-E0)/KBT));}
        // put boltzmann population into row 0 of eigenstates...
        for(j=1;j<=dj;++j)
         {(*eigenstates)(0,j)=complex<double>(En(j),exp(-(En(j)-E0)/KBT)/zz);}
   
// return eigenstates;
}

/**************************************************************************/
// for mcdisp this routine is needed
int ionpars::cfielddm(int & tn,double & T,Vector &  gjmbHxc,Vector & Hext,ComplexVector & u1,float & delta,ComplexMatrix & ests)
{  /*on input
    tn      ... number of transition to be computed 
    sign(tn)... 1... without printout, -1 with extensive printout
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta-+	energy of transition [meV]
    u1(i)	<-|Ji|+> sqrt(n+-n-),  n+,n-
    .... occupation number of states (- to + transition chosen according to transitionnumber)
*/
Vector gjmbH(1,gjmbHxc.Hi());
gjmbH=gjmbHxc;
gjmbH(1)+=gJ*MU_B*Hext(1);
gjmbH(2)+=gJ*MU_B*Hext(2);
gjmbH(3)+=gJ*MU_B*Hext(3);

// check dimensions of vector
if(gjmbH.Hi()>48)
   {fprintf(stderr,"Error loadable module cfield.so: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

//  Driver routine to compute the  eigenvalues and normalized eigenvectors 
//  of a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. The real  and imaginary
//  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
// 			   int sort, int maxiter)

 Vector JJ(1,gjmbH.Hi());
double lnz,u;JJ=0;
if (T>0){cfieldJJ(JJ,T,gjmbHxc,Hext,lnz,u,ests);  //expectation values <J>
        }
        else
        {T=-T;}
   double ninit=u1[1].real();
   double pinit=u1[1].imag();
  int pr;
  pr=0;
  if (tn<0) {pr=1;tn*=-1;}
   // setup hamiltonian
   int dj,j;
   dj=Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);
    
   Ham=Hcf-gjmbH(1)*Ja-gjmbH(2)*Jb-gjmbH(3)*Jc;
 for(j=4;j<=gjmbH.Hi();++j){Ham-=gjmbH(j)*(*Olm[j-3]);
//double dd; dd=NormFro((*OOlm[j-3])-(*OOlm[j-3]).Conjugate().Transpose());
//   if (dd>1e-5) {printf("j=%i\n",j);myPrintComplexMatrix(stderr,(*OOlm[j-3]));}
}

/*   int i1,j1; //printout matrix
    printf ("\n");
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) {printf ("%4.6g ",
    real(((*OOlm[5])-Jcc*Jcc+Jaa*Jaa)(i1,j1)));}
//    real((Jcc*Jaa+Jaa*Jcc)(i1,j1)));}
//    real((*OOlm[1])(i1,j1)));}
    printf ("\n");
    }
    printf ("\n");
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) {printf ("%4.6g ",
    imag(((*OOlm[5])-Jcc*Jcc+Jaa*Jaa)(i1,j1)));}
//   imag((Jcc*Jaa+Jaa*Jcc)(i1,j1)));}
//   imag((*OOlm[1])(i1,j1)));}
    printf ("\n");
    }
exit(0);      
*/    
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=1;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);
    
   // calculate Z and wn (occupation probability)
     Vector wn(1,dj);double Zs;
     double x,y;int i,k,l;
     x=Min(En);
     for (i=1;i<=dj;++i)
     {if ((y=(En(i)-x)/KB/T)<700) wn[i]=exp(-y); 
      else wn[i]=0.0;
//      printf("%4.4g\n",En(i));
      }
     Zs=Sum(wn);wn/=Zs;  
     Zs*=exp(-x/KB/T);
   if (ninit>dj)ninit=dj;
   if (pinit<SMALL)pinit=SMALL;
   double zsum=0,zii;
   int noft=0;for(i=1;(i<=ninit)&((zii=exp(-(En(i)-x)/KB/T))>(pinit*zsum));++i){noft+=dj-i+1;zsum+=zii;}


   // calculate Ja,Jb,Jc
     ComplexMatrix z(1,dj,1,dj);
     ComplexMatrix ** zp;
     zp=new ComplexMatrix*[gjmbH.Hi()+1];
     for(l=1;l<=gjmbH.Hi();++l)
      {zp[l]= new ComplexMatrix(1,dj,1,dj);}
     z=ComplexMatrix(zr,zi);
     
     (*zp[1])=Jaa*z;
     (*zp[2])=Jbb*z;
     (*zp[3])=Jcc*z;

     
 for(j=4;j<=gjmbH.Hi();++j)
    {(*zp[j])=(*OOlm[j-3])*z;  
}
     
// calculate mat and delta for transition number tn
// 1. get i and j from tn
k=0;
for(i=1;i<=dj;++i){for(j=i;j<=dj;++j)
{++k;if(k==tn)break;
}if(k==tn)break;}

// 2. set delta
delta=En(j)-En(i);

if (delta<-0.000001){fprintf(stderr,"ERROR module so1ion or cfield.so - dchargedensity_coeff1calc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(j==i)delta=-SMALL; //if transition within the same level: take negative delta !!- this is needed in routine intcalc

// 3. set mat
for(l=1;l<=gjmbH.Hi();++l)
{if(i==j){//take into account thermal expectation values <Jl>
          u1(l)=(((*zp[l]).Column(j)*z.Column(i))-JJ(l));}
 else    {u1(l)=((*zp[l]).Column(j)*z.Column(i));}}
           // ... in complex vector scalar product a*b is defined as: a.conj(b) !!! (see cvector.cc)

if (delta>SMALL)
   { if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      printf(" |<%i|Ja|%i>|^2=%4.4g |<%i|Jb|%i>|^2=%4.4g |<%i|Jc|%i>|^2=%4.4g",i,j,abs(u1(1))*abs(u1(1)),i,j,abs(u1(2))*abs(u1(2)),i,j,abs(u1(3))*abs(u1(3)));
      printf(" n%i-n%i=%4.4g\n",i,j,wn(i)-wn(j));}
    u1*=sqrt(wn(i)-wn(j)); // occupation factor
     }else
   {// quasielastic scattering has not wi-wj but wj*epsilon/kT
     if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      printf(" |<%i|Ja-<Ja>|%i>|^2=%4.4g |<%i|Jb-<Jb>|%i>|^2=%4.4g |<%i|Jc-<Jc>|%i>|^2=%4.4g",i,j,abs(u1(1))*abs(u1(1)),i,j,abs(u1(2))*abs(u1(2)),i,j,abs(u1(3))*abs(u1(3)));
      printf(" n%i=%4.4g\n",i,wn(i));}
    u1*=sqrt(wn(i)/KB/T);
   }

//clean up memory
     for(l=1;l<=gjmbH.Hi();++l)
      {delete zp[l];}
     delete []zp;

// return number of all transitions     
// return (int)((J+1)*(2*J+1));
//printf("noft=%i dj=%i\n",noft,dj);
return noft;

}

int ionpars::cfielddn(int & tn,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,Vector & Zc,ComplexMatrix & est,double & T,ComplexVector & dMQ)
{/*on input
    tn      ... number of transition to be computed 
    sign(tn)... 1... without printout, -1 with extensive printout
    est		matrix with eigenstates, eigenvalues [meV], population numbers
    th ph  .... polar angles of the scattering vector with respect to xyz=cab coordinate system (cfield) or xyz=abc (so1ion)
on output    
    int   	total number of transitions
    dMQ(i)	-2<-|Q-<Q>|+> sqrt(n+-n-),  n+,n-
     // note that  <M(Q)>=-2<Q>_TH in units of mb
    .... occupation number of states (- to + transition chosen according to transitionnumber)
*/
  int pr;pr=0;if (tn<0) {pr=1;tn*=-1;}
  int i,j=1,k,l;
  int dj=(int)(2*J+1);
  double delta;
   double ninit=dMQ(1).real();
   double pinit=dMQ(1).imag();
   if (ninit>dj)ninit=dj;
   if (pinit<SMALL)pinit=SMALL;
   double zsum=0,zii;
   int noft=0;for(i=1;(i<=ninit)&((zii=exp(-(real(est(0,i))-real(est(0,1)))/KB/T))>(pinit*zsum));++i){noft+=dj-i;zsum+=zii;}
//printf("!!! ddMQ noft = %i ninit= %g pinit= %g zii=%g zsum=%g T=%g!!!!\n",noft,ninit,pinit,zii,zsum,T);
//noft=(int)((J+1)*(2*J+1));

// calculate nat for transition number tn
// 1. get i and j from tn (as in du1calc
k=0;
for(i=1;i<=dj;++i){for(j=i;j<=dj;++j)
{++k;if(k==tn)break;
}if(k==tn)break;}

// 2. set delta
delta=real(est(0,j))-real(est(0,i));

 
if (delta<-0.000001){fprintf(stderr,"ERROR module so1ion/cfield.so - ddMQcalc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(j==i)delta=-SMALL; //if transition within the same level: take negative delta !!- this is needed in routine intcalc

	 ComplexMatrix * MQMi[4];
         MQMi[1]=new ComplexMatrix(1,dj,1,dj);
         MQMi[2]=new ComplexMatrix(1,dj,1,dj);
         MQMi[3]=new ComplexMatrix(1,dj,1,dj);
        MQM((*MQMi[1]),(*MQMi[2]),(*MQMi[3]),th,ph,J0,J2,J4,J6,Zc);
        //      x           y         z   // this has been fixed for module so1ion now 3.4.10 MR
        //      a           b         c   // ... for module cfield a backtransformation in ddMQcalc has been introduced in jjjpar.cpp
      
// 3. set dMQ
         int K,M,Md;
         ComplexVector Malpha(1,3);Malpha=0;
          for(K=1;K<=3;++K){for(M=1;M<=dj;++M){for(Md=1;Md<=dj;++Md){
             Malpha(K)+=conj(est(M,i))*(*MQMi[K])(M,Md)*est(Md,j); 
            }}} 
if(i==j){//take into account thermal expectation values <Jl> //MR120120
         ComplexVector mm(1,3); mm=0;                        //MR120120
         for(K=1;K<=dj;++K){for(M=1;M<=dj;++M){for(Md=1;Md<=dj;++Md){  //MR120120
           mm(1)+=imag(est(0,K))*conj(est(M,K))*(*MQMi[1])(M,Md)*est(Md,K); //MR120120
           mm(2)+=imag(est(0,K))*conj(est(M,K))*(*MQMi[2])(M,Md)*est(Md,K); //MR120120
           mm(3)+=imag(est(0,K))*conj(est(M,K))*(*MQMi[3])(M,Md)*est(Md,K); //MR120120
         }}} // --> mm(1,..3)  thermal expextation values of M              //MR120120
          Malpha-=mm;// subtract thermal expectation values                 //MR120120
         }  //MR120120
         delete MQMi[1];delete MQMi[2]; delete MQMi[3];


       // set vector dMQ=2* <i|Ml|j>
       dMQ=0;
          for(l=1;l<=3;++l)
          {dMQ(l)=-2.0*Malpha(l);}

// multiply by occupation number difference ...

if (delta>SMALL)
   { if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      printf(" |<%i|MQa-<MQa>|%i>|^2=%4.4g |<%i|MQb-<MQb>|%i>|^2=%4.4g |<%i|MQc-<MQc>|%i>|^2=%4.4g",i,j,abs(dMQ(1))*abs(dMQ(1)),i,j,abs(dMQ(2))*abs(dMQ(2)),i,j,abs(dMQ(3))*abs(dMQ(3)));
      printf(" n%i-n%i=%4.4g\n",i,j,imag(est(0,i))-imag(est(0,j)));}
    dMQ*=sqrt(imag(est(0,i))-imag(est(0,j))); // occupation factor
     }else
   {// quasielastic scattering has not wi-wj but wj*epsilon/kT
     if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      printf(" |<%i|MQa-<MQa>|%i>|^2=%4.4g |<%i|MQb-<MQb>|%i>|^2=%4.4g |<%i|MQc-<MQc>|%i>|^2=%4.4g",i,j,abs(dMQ(1))*abs(dMQ(1)),i,j,abs(dMQ(2))*abs(dMQ(2)),i,j,abs(dMQ(3))*abs(dMQ(3)));
      printf(" n%i=%4.4g\n",i,imag(est(0,i)));}
    dMQ*=sqrt(imag(est(0,i))/KB/T);
   }


// return number of all transitions     

return noft;
}

int ionpars::cfielddrixs1(int & tn,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,Vector & Zc,ComplexMatrix & est,double & T,ComplexVector & drixs)
{/*on input
    tn      ... number of transition to be computed 
    sign(tn)... 1... without printout, -1 with extensive printout
    est		matrix with eigenstates, eigenvalues [meV], population numbers
    th ph  .... polar angles of the scattering vector with respect to xyz=cab coordinate system (cfield) or xyz=abc (so1ion)
on output    
    int   	total number of transitions
    drixs(1...9)	-2<-|Rij-<Rij>|+> sqrt(n+-n-),  n+,n-
     // Rij transition operator according to Haverkort PRL (9 components for different polarisation channels
    .... occupation number of states (- to + transition chosen according to transitionnumber)
*/
  int pr;pr=0;if (tn<0) {pr=1;tn*=-1;}
  int i,j=1,k,l;
  int dj=(int)(2*J+1);
  double delta;
   double ninit=drixs(1).real();
   double pinit=drixs(1).imag();
   if (ninit>dj)ninit=dj;
   if (pinit<SMALL)pinit=SMALL;
   double zsum=0,zii;
   int noft=0;for(i=1;(i<=ninit)&((zii=exp(-(real(est(0,i))-real(est(0,1)))/KB/T))>(pinit*zsum));++i){noft+=dj-i;zsum+=zii;}
//printf("!!! ddMQ noft = %i ninit= %g pinit= %g zii=%g zsum=%g T=%g!!!!\n",noft,ninit,pinit,zii,zsum,T);
//noft=(int)((J+1)*(2*J+1));

// calculate nat for transition number tn
// 1. get i and j from tn (as in du1calc
k=0;
for(i=1;i<=dj;++i){for(j=i;j<=dj;++j)
{++k;if(k==tn)break;
}if(k==tn)break;}

// 2. set delta
delta=real(est(0,j))-real(est(0,i));

 
if (delta<-0.000001){fprintf(stderr,"ERROR module so1ion/cfield.so - ddMQcalc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(j==i)delta=-SMALL; //if transition within the same level: take negative delta !!- this is needed in routine intcalc


      
// 3. set drixs
         int K,M,Md;
         drixs=0;
          for(K=1;K<=NOF_RIXS_MATRICES;++K){for(M=1;M<=dj;++M){for(Md=1;Md<=dj;++Md){
             drixs(K)+=conj(est(M,i))*(*Ri[K])(M,Md)*est(Md,j); 
            }}} 
if(i==j){//take into account thermal expectation values <Jl> 
         ComplexVector Rav(1,NOF_RIXS_MATRICES); Rav=0;                        
         for(K=1;K<=dj;++K){for(M=1;M<=dj;++M){for(Md=1;Md<=dj;++Md){  
           for(l=1;l<=NOF_RIXS_MATRICES;++l){
           Rav(l)+=imag(est(0,K))*conj(est(M,K))*(*Ri[l])(M,Md)*est(Md,K); }           
         }}} // --> Rav(1,..NOF_RIXS_MATRICES)  thermal expextation values of R              
          drixs-=Rav;// subtract thermal expectation values                
         }  //MR120120


// multiply by occupation number difference ...

if (delta>SMALL)
   { if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      for(int i1=1;i1<=NOF_RIXS_MATRICES;++i1)
       {printf(" |<%i|Rij%i-<Rij%i>|%i>|^2=%4.4g",i,i1,i1,j,abs(drixs(i1))*abs(drixs(i1)));
       if(i1%3)printf("\n");}
      printf(" n%i-n%i=%4.4g\n",i,j,imag(est(0,i))-imag(est(0,j)));}
    drixs*=sqrt(imag(est(0,i))-imag(est(0,j))); // occupation factor
     }else
   {// quasielastic scattering has not wi-wj but wj*epsilon/kT
     if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      for(int i1=1;i1<=NOF_RIXS_MATRICES;++i1)printf(" |<%i|Rij%i-<Rij%i>|%i>|^2=%4.4g",i,i1,i1,j,abs(drixs(i1))*abs(drixs(i1)));
      printf(" n%i=%4.4g\n",i,imag(est(0,i)));}
    drixs*=sqrt(imag(est(0,i))/KB/T);
   }


// return number of all transitions     

return noft;
}


//**********************************************************************/
// routine to calculate the charge density coefficients of Zlm() R(r)^2
// *********************************************************************
void ionpars::chargedensity_coeffcalc(Vector &mom, double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & parstorage)
{
    /*on input
    T		temperature[K]
    Hxc 	vector of exchange field [meV]
  on output    
    mom		chargedensity coefficients: (if T>0 thermal exp value <mom>T 
                                                if T<0 the program asks for w_n and calculates
						       exp value <mom>=sum_n w_n <n|mom|n>
						       
*/
Vector gjmbH(1,Hxc.Hi());
gjmbH=Hxc;
gjmbH(1)+=gJ*MU_B*Hext(1);
gjmbH(2)+=gJ*MU_B*Hext(2);
gjmbH(3)+=gJ*MU_B*Hext(3);
// check dimensions of vector
if(gjmbH.Hi()>48)
   {fprintf(stderr,"Error internal module so1ion/cfield: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

if(mom.Hi()!=28){fprintf(stderr,"Error internal module so1ion/cfield chargedensity_coeff: moment vector has not dimension 28 but %i\n",mom.Hi());
    exit(EXIT_FAILURE);}

//  Driver routine to compute the  eigenvalues and normalized eigenvectors 
//  of a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. The real  and imaginary
//  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
// 			   int sort, int maxiter)


   // setup hamiltonian
   int dj,i,j;
  
   dj=Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);
//   Matrix Tam(1,dj,1,dj); // transformed Hamiltonian
   ComplexMatrix z(1,dj,1,dj);
   ComplexMatrix za(1,dj,1,dj);
   ComplexMatrix zb(1,dj,1,dj);
   ComplexMatrix zc(1,dj,1,dj);
   ComplexMatrix zolm(1,dj,1,dj);    

   Ham=Hcf-gjmbH(1)*Ja-gjmbH(2)*Jb-gjmbH(3)*Jc;

// here the zeeman term is extended for multipolar fields
   for(j=4;j<=Hxc.Hi();++j){Ham-=gjmbH(j)*(*Olm[j-3]);}

// diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=0;int maxiter=1000000;
   if (T<0) sort=1;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);

   // calculate Z and wn (occupation probability)
     Vector wn(1,dj);
     double x,y;
     x=Min(En);
     double Zs;

     if (T>0)
     { for (i=1;i<=dj;++i)
       {if ((y=(En(i)-x)/KB/T)<600) wn[i]=exp(-y); 
        else wn[i]=0.0;
       }
       Zs=Sum(wn);wn/=Zs;
     } 
     else
     { printf ("Temperature T<0: please choose probability distribution of states by hand\n");
                         printf ("Number   Energy     Excitation Energy\n");
     for (i=1;i<=dj;++i) printf ("%i    %4.4g meV   %4.4g meV\n",i,En(i),En(i)-x);
     char instr[MAXNOFCHARINLINE];
     for (i=1;i<=dj;++i)
      {printf("eigenstate %i: %4.4g meV %4.4g meV  - please enter probability w(%i):",i,En(i),En(i)-x,i);
       fgets(instr, MAXNOFCHARINLINE, stdin);
 
       wn(i)=strtod(instr,NULL);
      }
       Zs=Sum(wn);wn/=Zs;
                         printf ("\n\nNumber   Energy     Excitation Energy   Probability\n");
     for (i=1;i<=dj;++i) printf ("%i    %4.4g meV   %4.4g meV %4.4g  \n",i,En(i),En(i)-x,wn(i));
     }
   z=ComplexMatrix(zr,zi);
   mom=0;
 if(nof_electrons==0){fprintf(stderr,"Error so1ion/cfield: nof_electrons=0 ... perhaps single ion property file does not contain the number of electrons in the shell: 'nof_electrons=...'\n");
     exit(EXIT_FAILURE);}
  mom(1) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
// a(0, 0) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
// Indices for chargedensity
//            0 not used
//          0 1  2  3 4 5 6  7  8  9 101112131415 16 17 18 19 20 2122232425262728 
int k[] = {-1,0, 2, 2,2,2,2, 4, 4, 4, 4,4,4,4,4,4, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-2,-1,0,1,2,-4,-3,-2,-1,0,1,2,3,4,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

// here the expectation values of the multipolar moments are calculated
 for(j=4;j<=8;++j){zolm=(*OOlm[j-3])*z;for (i=1;i<=dj;++i) mom[j-2]+=alpha*cnst(k[j-2],q[j-2])*wn(i)*real(z.Column(i)*zolm.Column(i));}
 for(j=16;j<=24;++j){zolm=(*OOlm[j-3])*z;for (i=1;i<=dj;++i) mom[j-9]+=beta*cnst(k[j-9],q[j-9])*wn(i)*real(z.Column(i)*zolm.Column(i));}
 for(j=36;j<=48;++j){zolm=(*OOlm[j-3])*z;for (i=1;i<=dj;++i) mom[j-20]+=gamma*cnst(k[j-20],q[j-20])*wn(i)*real(z.Column(i)*zolm.Column(i));}

     // theta_J*cnst(l,m)  are prefactors to get coefficients of Zlm*R(r)^2 
    //in case of module cfield and so1ion(stevens parameters tetan and zlm prefactors)

}


//**********************************************************************/
// routine to calculate the transition matrix elements of the charge density coefficients (of Zlm() R(r)^2)
// *********************************************************************
int ionpars::dchargedensity_coeff1calc(int & tn,double & T,Vector &  gjmbHxc,Vector & Hext, ComplexVector & cd1,float & delta,ComplexMatrix & ests)
{  /*on input
    tn      ... number of transition to be computed 
    sign(tn)... 1... without printout, -1 with extensive printout
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta-+	energy of transition [meV]
    cd1(i)	<-|theta_l plm Olm|+> sqrt(n+-n-),  n+,n-
    .... occupation number of states (- to + transition chosen according to transitionnumber)
*/
Vector gjmbH(1,gjmbHxc.Hi());
gjmbH=gjmbHxc;
gjmbH(1)+=gJ*MU_B*Hext(1);
gjmbH(2)+=gJ*MU_B*Hext(2);
gjmbH(3)+=gJ*MU_B*Hext(3);

// check dimensions of vector
if(gjmbH.Hi()>48)
   {fprintf(stderr,"Error loadable module cfield.so: wrong number of dimensions - check number of columns in file mcphas.j\n");
    exit(EXIT_FAILURE);}

//  Driver routine to compute the  eigenvalues and normalized eigenvectors 
//  of a complex Hermitian matrix z.The real parts of the elements must be
//  stored in the lower triangle of z,the imaginary parts (of the elements
//  corresponding to the lower triangle) in the positions
//  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
//  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
//  True, otherwise  not ordered for sort = False. The real  and imaginary
//  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
//  The storage requirement is 3*n*n + 4*n complex numbers. 
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
// void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
// 			   int sort, int maxiter)

 Vector JJ(1,28);JJ=0;
if (T>0){chargedensity_coeffcalc(JJ,T,gjmbHxc,Hext,ests); // expectation values for cd coeffs
        }
        else
        {T=-T;}
   double ninit=cd1[1].real();
   double pinit=cd1[1].imag();
  int pr;
  pr=0;
  if (tn<0) {pr=1;tn*=-1;}
   // setup hamiltonian
   int dj,j;
   dj=Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);
    
   Ham=Hcf-gjmbH(1)*Ja-gjmbH(2)*Jb-gjmbH(3)*Jc;
 for(j=4;j<=gjmbH.Hi();++j){Ham-=gjmbH(j)*(*Olm[j-3]);
//double dd; dd=NormFro((*OOlm[j-3])-(*OOlm[j-3]).Conjugate().Transpose());
//   if (dd>1e-5) {printf("j=%i\n",j);myPrintComplexMatrix(stderr,(*OOlm[j-3]));}
}

/*   int i1,j1; //printout matrix
    printf ("\n");
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) {printf ("%4.6g ",
    real(((*OOlm[5])-Jcc*Jcc+Jaa*Jaa)(i1,j1)));}
//    real((Jcc*Jaa+Jaa*Jcc)(i1,j1)));}
//    real((*OOlm[1])(i1,j1)));}
    printf ("\n");
    }
    printf ("\n");
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) {printf ("%4.6g ",
    imag(((*OOlm[5])-Jcc*Jcc+Jaa*Jaa)(i1,j1)));}
//   imag((Jcc*Jaa+Jaa*Jcc)(i1,j1)));}
//   imag((*OOlm[1])(i1,j1)));}
    printf ("\n");
    }
exit(0);      
*/    
   // diagonalize
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=1;int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);
    
   // calculate Z and wn (occupation probability)
     Vector wn(1,dj);double Zs;
     double x,y;int i,l;
     x=Min(En);
     for (i=1;i<=dj;++i)
     {if ((y=(En(i)-x)/KB/T)<700) wn[i]=exp(-y); 
      else wn[i]=0.0;
//      printf("%4.4g\n",En(i));
      }
     Zs=Sum(wn);wn/=Zs;  
     Zs*=exp(-x/KB/T);
   if (ninit>dj)ninit=dj;
   if (pinit<SMALL)pinit=SMALL;
   double zsum=0,zii;
   int noft=0;for(i=1;(i<=ninit)&((zii=exp(-(En(i)-x)/KB/T))>(pinit*zsum));++i){noft+=dj-i+1;zsum+=zii;}


   // calculate Ja,Jb,Jc
     ComplexMatrix z(1,dj,1,dj);
     ComplexMatrix ** zp;
     zp=new ComplexMatrix*[28+1];
     for(l=1;l<=28;++l)
      {zp[l]= new ComplexMatrix(1,dj,1,dj);}
     z=ComplexMatrix(zr,zi);
     (*zp[1])=0;
// a(0, 0) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
// Indices for chargedensity
//            0 not used
//          0 1  2  3 4 5 6  7  8  9 101112131415 16 17 18 19 20 2122232425262728 
int k[] = {-1,0, 2, 2,2,2,2, 4, 4, 4, 4,4,4,4,4,4, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-2,-1,0,1,2,-4,-3,-2,-1,0,1,2,3,4,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

     for(j=2;j<=6;++j){(*zp[j])=alpha*cnst(k[j],q[j])*(*OOlm[j-1])*z;}
     for(j=7;j<=15;++j){(*zp[j])=beta*cnst(k[j],q[j])*(*OOlm[j-1+7])*z;}
     for(j=16;j<=28;++j){(*zp[j])=gamma*cnst(k[j],q[j])*(*OOlm[j-1+7+11])*z;}

     // theta_J*cnst(l,m)  are prefactors to get coefficients of Zlm*R(r)^2 
    //in case of module cfield and so1ion(stevens parameters tetan and zlm prefactors)
  

   
// calculate mat and delta for transition number tn
// 1. get i and j from tn
int kk=0;
for(i=1;i<=dj;++i){for(j=i;j<=dj;++j)
{++kk;if(kk==tn)break;
}if(kk==tn)break;}

// 2. set delta
delta=En(j)-En(i);

if (delta<-0.000001){fprintf(stderr,"ERROR module so1ion or cfield.so - dchargedensity_coeff1calc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
if(j==i)delta=-SMALL; //if transition within the same level: take negative delta !!- this is needed in routine intcalc

// 3. set mat
for(l=1;l<=28;++l)
{if(i==j){//take into account thermal expectation values <Jl>
          cd1(l)=(((*zp[l]).Column(j)*z.Column(i))-JJ(l));}
 else    {cd1(l)=((*zp[l]).Column(j)*z.Column(i));}}
           // ... in complex vector scalar product a*b is defined as: a.conj(b) !!! (see cvector.cc)

if (delta>SMALL)
   { if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      for(l=1;l<=28;++l)printf(" |<%i|cd_coeff%i|%i>|^2=%4.4g ",i,l,j,abs(cd1(l))*abs(cd1(l)));
      printf(" n%i-n%i=%4.4g\n",i,j,wn(i)-wn(j));}
    cd1*=sqrt(wn(i)-wn(j)); // occupation factor
     }else
   {// quasielastic scattering has not wi-wj but wj*epsilon/kT
     if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      for(l=1;l<=28;++l)printf(" |<%i|cd_coeff%i-<cd_coeff%i>|%i>|^2=%4.4g ",i,l,l,j,abs(cd1(l))*abs(cd1(l)));
      printf(" n%i=%4.4g\n",i,wn(i));}
    cd1*=sqrt(wn(i)/KB/T);
   }

//clean up memory
     for(l=1;l<=28;++l)
      {delete zp[l];}
     delete []zp;

// return number of all transitions     
// return (int)((J+1)*(2*J+1));
//printf("noft=%i dj=%i\n",noft,dj);
return noft;
}   


//**********************************************************************/
// routine to calculate the scattering operator to go beyond dip approx
// *********************************************************************

// just another routine to calculakte Z(K)
double Z(int K, float J0, float J2, float J4, float J6, Vector Zc)
{// calculate Z(K)
 if (K==1) return Zc(1)*J0+Zc(2)*J2;
 if (K==3) return Zc(3)*J2+Zc(4)*J4;
 if (K==5) return Zc(5)*J4+Zc(6)*J6;
 if (K==7) return Zc(7)*J6;
 
 return 0;
}

void ionpars::MQM(ComplexMatrix & MQXM,ComplexMatrix & MQYM,ComplexMatrix & MQZM, double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc)
{complex <double> im(0,1);
         // .... calculate scattering operator ...(formula 11.141-143 in lovesey)
	    int K,Qd,M,Md;double MJ,MJd,PKdQd,thj,factor;
	    complex <double>bracketx,brackety,bracketz;
	    MQXM=0;MQYM=0;MQZM=0;
	    for(K=1;K<=7;K+=2){
			     factor=sqrt(4.0*PI)*Z(K,J0,J2,J4,J6,Zc)/K;
                               if (factor!=0){
					     thj=threej((float)K,J,J,0,J,-J);
					     for(Qd=-K;Qd<=K;Qd+=1){
                                                 bracketx=0;brackety=0;
					       if(K-1>=Qd+1&&K-1>=-Qd-1)
					       {bracketx+=SphericalHarmonicY (K-1,Qd+1,th,ph)*sqrt((double)(K-Qd)*(K-Qd-1));
					        brackety+=SphericalHarmonicY (K-1,Qd+1,th,ph)*sqrt((double)(K-Qd)*(K-Qd-1));
					       }
					       if(K-1>=Qd-1&&K-1>=-Qd+1)
                                                 {bracketx-=SphericalHarmonicY (K-1,Qd-1,th,ph)*sqrt((double)(K+Qd)*(K+Qd-1));
						        brackety+=SphericalHarmonicY (K-1,Qd-1,th,ph)*sqrt((double)(K+Qd)*(K+Qd-1));
					       }
					       if(K-1>=Qd&&K-1>=-Qd)
                                                 {bracketz=SphericalHarmonicY (K-1,Qd,th,ph)*sqrt((double)(K-Qd)*(K+Qd));
                                                 }else {bracketz=0;}

//.(1)..USE !		     ThreeJSymbolM	(J1,J2,J3,M1,&M2min,&M2max,*thrcof,ndim,errflag);
                                                 double thrj[30];int ndim=30; double MJdmin,MJdmax; int errflag;
                                                                      
                                                 ThreeJSymbolM ((float)K,J,J,-(float)Qd,MJdmin,MJdmax,thrj,ndim,errflag);
                                                 if (errflag!=0){fprintf(stderr,"ERROR mcdiff: threejsymbol error %i\n",errflag);exit(EXIT_FAILURE);}           
                                                 for (Md=int(MJdmin+1+J);Md<=int(MJdmax+1+J);++Md){
						                 MJd=(float)Md-1-J;
								 MJ=-Qd+MJd;M=int(MJ+1+J);
							         PKdQd=thrj[Md-int(MJdmin+1+J)]/thj; 
							         PKdQd*=odd(int(J-MJd)) ? -1 : 1;
							         MQXM(M,Md)+=0.5*factor*PKdQd*bracketx;
							         MQYM(M,Md)+=-im*brackety*0.5*factor*PKdQd;
							         MQZM(M,Md)+=factor*PKdQd*bracketz;
                                                                 }

/*
//.(2)..                     3jsymb=threej(J1,J2,J3,M1,M2,M3) 
//							       for(M=1;M<=dj;++M){MJ=(float)M-1-J;  
//							        for(Md=1;Md<=dj;++Md){MJd=(float)Md-1-J; 
//							         // according to 11.140 lovesey book        
//							         PKdQd=threej((float)K,J,J,-(float)Qd,MJd,-MJ)/thj; 
//							         PKdQd*=odd(int(J-MJd)) ? -1 : 1;
//							         MQXM(M,Md)+=im*0.5*factor*PKdQd*bracketx;
//							         MQYM(M,Md)+=-0.5*factor*PKdQd*brackety;
//							         MQZM(M,Md)+=factor*PKdQd*bracketz;
//							        }
//							       }
*/
                                             }
							      }
                                                             }

}
// calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb
// according to stored eigenstate matrix est
// calculates the scattering operator given the polar angles th, ph (with respect to the CEF coordinate 
// system xyz and the <jl(qr)> and the eigenstate matrix with eigenstates and thermal population numbers
ComplexVector & ionpars::MQ(double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc,ComplexMatrix & est)
    {  
       int dj=(int)(2*J+1);
       
	 ComplexMatrix MQXM(1,dj,1,dj),MQYM(1,dj,1,dj),MQZM(1,dj,1,dj);
         MQM(MQXM,MQYM,MQZM,th,ph,J0,J2,J4,J6,Zc);
							     // ... calculate thermal expectation values
							     // using the eigenstates and T 
							     // mom(1) = ....
          						     //
       static ComplexVector mm(1,3); mm=0;
       int K,M,Md;
       for(K=1;K<=dj;++K){for(M=1;M<=dj;++M){for(Md=1;Md<=dj;++Md){
         mm(1)+=imag(est(0,K))*conj(est(M,K))*MQXM(M,Md)*est(Md,K); 
         mm(2)+=imag(est(0,K))*conj(est(M,K))*MQYM(M,Md)*est(Md,K); 
         mm(3)+=imag(est(0,K))*conj(est(M,K))*MQZM(M,Md)*est(Md,K); 
       }}}
// myPrintComplexMatrix(stdout,MQXM);
// myPrintComplexMatrix(stdout,MQYM);
// myPrintComplexMatrix(stdout,MQZM);
// 					       myPrintComplexMatrix(stdout,est);}
					       
                     mm*=2; // this is now <M(Q)>=-2x<Q>_TH in units of mb
// myPrintComplexVector(stdout,mm);//equivalent to moment ...
    return mm;
    }

// for testing the code uncomment and make test and start ionpars.exe
/* int main(int argc, char **argv)
{FILE * cf_file;
cf_file = fopen_errchk (argv[1], "rb"); // reopen file

      ionpars iops(cf_file);
       myPrintComplexMatrix(stderr,(*iops.OOlm[26-3]));getchar();
      
      fclose(cf_file);cf_file = fopen_errchk (argv[1], "rb"); // reopen file
      ionpars iops1(cf_file);
       myPrintComplexMatrix(stderr,(*iops1.OOlm[26-3]));getchar();
      
      
return 1;
}*/
