// ionpars: class to load and store matrices for internal module cfield and so1ion

#include "ionpars.hpp"
#include "martin.h"
#include "ionpars.h"
#include "myev.h"
#include "perlparse.h"


 ionpars::ionpars (const ionpars & p) //copy constructor
 {J=p.J;so1ion=p.so1ion;
  Ja=p.Ja; Jb=p.Jb; Jc=p.Jc;Hcf=p.Hcf;
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
    for(i=1;i<=NOF_OLM_MATRICES;++i){ Olm [i]= new Matrix(1,(*p.Olm[i]).Rhi(),1,(*p.Olm[i]).Chi()); 
                   (*Olm[i])=(*p.Olm[i]);}
    In= new Matrix * [1+IONPARS_MAXNOFCOMPONENTS]; 
    for(i=0;i<=IONPARS_MAXNOFCOMPONENTS;++i){In[i]= new Matrix(1,(*p.In[i]).Rhi(),1,(*p.In[i]).Chi());
                   (*In[i])=(*p.In[i]);}  
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
// cnst is the Zlm constants - put them into the matrix ... (same code is reused in jjjpar.cpp, pointc.c)
   cnst=Matrix(0,6,-6,6);set_zlm_constants(cnst);
   iontype = new char [MAXNOFCHARINLINE];
   Blm=Vector(0,NOF_OLM_MATRICES);Blm=0; // vector of crystal field parameters
   Llm=Vector(0,45);Llm=0; // vector of crystal field parameters
 int i;   
   Ri = new ComplexMatrix * [1+NOF_RIXS_MATRICES];
  for(i=1;i<=NOF_RIXS_MATRICES;++i)Ri[i]= new ComplexMatrix(1,dimj,1,dimj);
   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
  for(i=1;i<=NOF_OLM_MATRICES;++i){ Olm [i]= new Matrix(1,dimj,1,dimj); }  
   In = new Matrix * [1+IONPARS_MAXNOFCOMPONENTS];
  for(i=0;i<=IONPARS_MAXNOFCOMPONENTS;++i)In[i]= new Matrix(1,dimj,1,dimj);  
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
// cnst is the Zlm constants - put them into the matrix ... (same code is reused in jjjpar.cpp, pointc.c)
cnst=Matrix(0,6,-6,6);set_zlm_constants(cnst);
   Blm=Vector(0,NOF_OLM_MATRICES);Blm=0; // vector of crystal field parameters
   Llm=Vector(0,45);Llm=0; // vector of crystal field parameters
 int i;   
   Ri = new ComplexMatrix * [1+NOF_RIXS_MATRICES];
  for(i=1;i<=NOF_RIXS_MATRICES;++i)Ri[i]= new ComplexMatrix(1,dimj,1,dimj);
   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
  for(i=1;i<=NOF_OLM_MATRICES;++i){ Olm [i]= new Matrix(1,dimj,1,dimj); }  
   In = new Matrix * [1+IONPARS_MAXNOFCOMPONENTS];
  for(i=0;i<=IONPARS_MAXNOFCOMPONENTS;++i)In[i]= new Matrix(1,dimj,1,dimj);
}

ionpars::~ionpars(){
 int i;
 delete []iontype;
 for (i=1;i<=NOF_RIXS_MATRICES;++i)delete Ri[i];
 for (i=1;i<=NOF_OLM_MATRICES;++i)  {delete Olm[i];}
 for (i=0;i<=IONPARS_MAXNOFCOMPONENTS;++i)delete In[i];
   delete[] Olm;
   delete[] Ri;
   delete[] In; 
 } //destructor

ionpars::ionpars(FILE * cf_file, char * cffilename) 
//constructor with commands from file handle (filename of cf parameters etc)
// MIND: this code has to be kept consistent with the code in cf1ion_module/eingabe.c (function read_new_format())
// ( reason: the stand alone c-progams so1ion, cfield should be consistent with the c++ programs
// mcphas, mcdisp etc (using this input routine)
{ static int pr=1;
  int dimj;complex<double> im(0,1);
  int i,j,l,m,nof_electronsr,inof=0,ir2=0,ir4=0,ir6=0,ialpha=0,ibeta=0,igamma=0,igj=0; 
  double alphar,betar,gammar,r2r,r4r,r6r,gJr;
  char instr[MAXNOFCHARINLINE];
  iontype= new char[MAXNOFCHARINLINE];
  char  moduletype[MAXNOFCHARINLINE];
   Blm=Vector(0,NOF_OLM_MATRICES);Blm=0; // vector of crystal field parameters
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
const char kq[]="00 22S21S20 21 22 33S32S31S30 31 32 33 44S43S42S41S40 41 42 43 44 55S54S53S52S51S50 51 52 53 54 55 66S65S64S63S62S61S60 61 62 63 64 65 66 x2 y2 z2 x4 y4 z4";
char kq3[5];kq3[4]='\0';kq3[3]='\0';
int perlp=1;//by default try perlparse  
// read in lines and get IONTYPE=  and CF parameters Blm
   while(feof(cf_file)==false)
  {fgets(instr, MAXNOFCHARINLINE, cf_file);
          // in IONTYPE line there is no #!perl ... do not perlparse
        if(!extract(instr,"IONTYPE",iontype,(size_t)MAXNOFCHARINLINE)&&strstr (instr,"#!perl")==NULL)perlp=0;
        inof+=1-extract(instr,"nof_electrons",nof_electrons); //MR 120127
        kq3[0]='B';for(int i=0;i<=45;++i){strncpy(kq3+1,kq+i*3,3);if(kq3[3]==' ')kq3[3]='\0';
                                          extract(instr,kq3,Blm(i));}
        extract(instr,"Dx2",Blm(46));
        extract(instr,"Dy2",Blm(47));
        extract(instr,"Dz2",Blm(48));
        extract(instr,"Dx4",Blm(49));
        extract(instr,"Dy4",Blm(50));
        extract(instr,"Dz4",Blm(51));
        kq3[0]='L';for(int i=0;i<=45;++i){strncpy(kq3+1,kq+i*3,3);if(kq3[3]==' ')kq3[3]='\0';        
                                          extract(instr,kq3,Llm(i));}
        ialpha+=1-extract(instr,"ALPHA",alpha);
        ibeta+=1-extract(instr,"BETA",beta);
        igamma+=1-extract(instr,"GAMMA",gamma);
        igj+=1-extract(instr,"GJ",gJ);       
        ir2+=1-extract(instr,"R2",r2);
        ir4+=1-extract(instr,"R4",r4);
        ir6+=1-extract(instr,"R6",r6);
        extract(instr,"SIGMA0r",s0r);
        extract(instr,"SIGMA1r",s1r);
        extract(instr,"SIGMA2r",s2r);
        extract(instr,"SIGMA0i",s0i);
        extract(instr,"SIGMA1i",s1i);
        extract(instr,"SIGMA2i",s2i);
 } //fclose(cf_file);
 fseek(cf_file,0,SEEK_SET);

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

  double  modxdxcr[31*31],modxdxci[31*31];
  double  modydycr[31*31],modydyci[31*31];
  double  modzdzcr[31*31],modzdzci[31*31];
    
if (pr==1) {printf("#using %s ...\n",moduletype);}
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

  modxdxcr,modxdxci,
  modydycr,modydyci,
  modzdzcr,modzdzci,

  &dimj,&alphar,&betar,&gammar,&gJr,&r2r,&r4r,&r6r, &nof_electronsr);


if (pr==1) printf("#end using %s\n",moduletype);

   J=((double)dimj-1)/2; //momentum quantum number
if (pr==1) printf("#J=%g\n",J);
   Ja = Matrix(1,dimj,1,dimj); 
   ComplexMatrix Jaa(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jaa(i,j)=im*Jxi[30*(i-1)+j-1]+Jxr[30*(i-1)+j-1];
    if(i<j){Ja(i,j)=Jxi[30*(j-1)+i-1];}else{Ja(i,j)=Jxr[30*(i-1)+j-1];}
   }
   Jb = Matrix(1,dimj,1,dimj); 
   ComplexMatrix Jbb(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jbb(i,j)=im*Jyi[30*(i-1)+j-1]+Jyr[30*(i-1)+j-1];
    if(i<j){Jb(i,j)=Jyi[30*(j-1)+i-1];}else{Jb(i,j)=Jyr[30*(i-1)+j-1];}
   }
   Jc = Matrix(1,dimj,1,dimj); 
   ComplexMatrix Jcc(1,dimj,1,dimj); 
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   {Jcc(i,j)=im*Jzi[30*(i-1)+j-1]+Jzr[30*(i-1)+j-1];
    if(i<j){Jc(i,j)=Jzi[30*(j-1)+i-1];}else{Jc(i,j)=Jzr[30*(i-1)+j-1];}
   }
if(so1ion==0)
  {//ATTENTION FOR cfield the AXES xyz are parallel to cab
   Matrix dummy(1,dimj,1,dimj); dummy=Jb;Jb=Jc;Jc=Ja;Ja=dummy;
   ComplexMatrix dummyc(1,dimj,1,dimj);dummyc=Jbb;Jbb=Jcc;Jcc=Jaa;Jaa=dummyc;

  if (pr==1) {printf("#Axis Convention for module cfield:  z'||b, x'||(a x b) and y' normal to x' and z'\n");
  printf("#x'y'z' .... Coordinate system of the crystal field parameters, magnetic moment\n");
  printf("#abc .... Crystal axes\n");
  printf("#Standard Interactions operator sequence  (Stevens Operators Olm imported from cfield code of Peter Hoffmann/Fabi)\n");
  printf("#I1,I2,I3,....= Jy(=O11(s)),Jz=(O10(c)),Jx(=O11(c)),O22(s),O21(s),O20(c),O21(c),O22(c),O33(s),...,O66(c)\n");  
              }
  }
  else
  {//ATTENTION FOR so1ion the AXES xyz are parallel to abc
  if (pr==1) {printf("#Axis Convention using so1ion as a module:  y||b, z||(a x b) and x normal to y and z\n");
  printf("#xyz .... Coordinate system of the crystal field parameters, magnetic moment\n");
  printf("#abc .... Crystal axes\n");
  printf("#Standard Interactions operator sequence  (Stevens Operators Olm imported from cfield code of Peter Hoffmann/Fabi)\n");
  printf("#I1,I2,I3,....= Jx(=O11(c)),Jy=(O11(s)),Jz(=O10(c)),O22(s),O21(s),O20(c),O21(c),O22(c),O33(s),...,O66(c)\n");  
             }
  }

//---------------------------------------------------------------------------
   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
      for(i=1;i<=NOF_OLM_MATRICES;++i)   
    {   Olm[i]= new Matrix(1,dimj,1,dimj); } 
   for(i=1;i<=dimj;++i){for(j=1;j<=dimj;++j)
   {
if(i<j){(*Olm[1])(i,j)=mo22si[30*(j-1)+i-1];}else{(*Olm[1])(i,j)=mo22sr[30*(i-1)+j-1];}
if(i<j){(*Olm[2])(i,j)=mo21si[30*(j-1)+i-1];}else{(*Olm[2])(i,j)=mo21sr[30*(i-1)+j-1];}
if(i<j){(*Olm[3])(i,j)=mo20ci[30*(j-1)+i-1];}else{(*Olm[3])(i,j)=mo20cr[30*(i-1)+j-1];}
if(i<j){(*Olm[4])(i,j)=mo21ci[30*(j-1)+i-1];}else{(*Olm[4])(i,j)=mo21cr[30*(i-1)+j-1];}
if(i<j){(*Olm[5])(i,j)=mo22ci[30*(j-1)+i-1];}else{(*Olm[5])(i,j)=mo22cr[30*(i-1)+j-1];}
    
if(i<j){(*Olm[6])(i,j)=mo33si[30*(j-1)+i-1];}else{(*Olm[6])(i,j)=mo33sr[30*(i-1)+j-1];}
if(i<j){(*Olm[7])(i,j)=mo32si[30*(j-1)+i-1];}else{(*Olm[7])(i,j)=mo32sr[30*(i-1)+j-1];}
if(i<j){(*Olm[8])(i,j)=mo31si[30*(j-1)+i-1];}else{(*Olm[8])(i,j)=mo31sr[30*(i-1)+j-1];}
if(i<j){(*Olm[9])(i,j)=mo30ci[30*(j-1)+i-1];}else{(*Olm[9])(i,j)=mo30cr[30*(i-1)+j-1];}
if(i<j){(*Olm[10])(i,j)=mo31ci[30*(j-1)+i-1];}else{(*Olm[10])(i,j)=mo31cr[30*(i-1)+j-1];}
if(i<j){(*Olm[11])(i,j)=mo32ci[30*(j-1)+i-1];}else{(*Olm[11])(i,j)=mo32cr[30*(i-1)+j-1];}
if(i<j){(*Olm[12])(i,j)=mo33ci[30*(j-1)+i-1];}else{(*Olm[12])(i,j)=mo33cr[30*(i-1)+j-1];}
    
if(i<j){(*Olm[13])(i,j)=mo44si[30*(j-1)+i-1];}else{(*Olm[13])(i,j)=mo44sr[30*(i-1)+j-1];}
if(i<j){(*Olm[14])(i,j)=mo43si[30*(j-1)+i-1];}else{(*Olm[14])(i,j)=mo43sr[30*(i-1)+j-1];}
if(i<j){(*Olm[15])(i,j)=mo42si[30*(j-1)+i-1];}else{(*Olm[15])(i,j)=mo42sr[30*(i-1)+j-1];}
if(i<j){(*Olm[16])(i,j)=mo41si[30*(j-1)+i-1];}else{(*Olm[16])(i,j)=mo41sr[30*(i-1)+j-1];}
if(i<j){(*Olm[17])(i,j)=mo40ci[30*(j-1)+i-1];}else{(*Olm[17])(i,j)=mo40cr[30*(i-1)+j-1];}
if(i<j){(*Olm[18])(i,j)=mo41ci[30*(j-1)+i-1];}else{(*Olm[18])(i,j)=mo41cr[30*(i-1)+j-1];}
if(i<j){(*Olm[19])(i,j)=mo42ci[30*(j-1)+i-1];}else{(*Olm[19])(i,j)=mo42cr[30*(i-1)+j-1];}
if(i<j){(*Olm[20])(i,j)=mo43ci[30*(j-1)+i-1];}else{(*Olm[20])(i,j)=mo43cr[30*(i-1)+j-1];}
if(i<j){(*Olm[21])(i,j)=mo44ci[30*(j-1)+i-1];}else{(*Olm[21])(i,j)=mo44cr[30*(i-1)+j-1];}
    
if(i<j){(*Olm[22])(i,j)=mo55si[30*(j-1)+i-1];}else{(*Olm[22])(i,j)=mo55sr[30*(i-1)+j-1];}
if(i<j){(*Olm[23])(i,j)=mo54si[30*(j-1)+i-1];}else{(*Olm[23])(i,j)=mo54sr[30*(i-1)+j-1];}
if(i<j){(*Olm[24])(i,j)=mo53si[30*(j-1)+i-1];}else{(*Olm[24])(i,j)=mo53sr[30*(i-1)+j-1];}
if(i<j){(*Olm[25])(i,j)=mo52si[30*(j-1)+i-1];}else{(*Olm[25])(i,j)=mo52sr[30*(i-1)+j-1];}
if(i<j){(*Olm[26])(i,j)=mo51si[30*(j-1)+i-1];}else{(*Olm[26])(i,j)=mo51sr[30*(i-1)+j-1];}
if(i<j){(*Olm[27])(i,j)=mo50ci[30*(j-1)+i-1];}else{(*Olm[27])(i,j)=mo50cr[30*(i-1)+j-1];}
if(i<j){(*Olm[28])(i,j)=mo51ci[30*(j-1)+i-1];}else{(*Olm[28])(i,j)=mo51cr[30*(i-1)+j-1];}
if(i<j){(*Olm[29])(i,j)=mo52ci[30*(j-1)+i-1];}else{(*Olm[29])(i,j)=mo52cr[30*(i-1)+j-1];}
if(i<j){(*Olm[30])(i,j)=mo53ci[30*(j-1)+i-1];}else{(*Olm[30])(i,j)=mo53cr[30*(i-1)+j-1];}
if(i<j){(*Olm[31])(i,j)=mo54ci[30*(j-1)+i-1];}else{(*Olm[31])(i,j)=mo54cr[30*(i-1)+j-1];}
if(i<j){(*Olm[32])(i,j)=mo55ci[30*(j-1)+i-1];}else{(*Olm[32])(i,j)=mo55cr[30*(i-1)+j-1];}

if(i<j){(*Olm[33])(i,j)=mo66si[30*(j-1)+i-1];}else{(*Olm[33])(i,j)=mo66sr[30*(i-1)+j-1];}
if(i<j){(*Olm[34])(i,j)=mo65si[30*(j-1)+i-1];}else{(*Olm[34])(i,j)=mo65sr[30*(i-1)+j-1];}
if(i<j){(*Olm[35])(i,j)=mo64si[30*(j-1)+i-1];}else{(*Olm[35])(i,j)=mo64sr[30*(i-1)+j-1];}
if(i<j){(*Olm[36])(i,j)=mo63si[30*(j-1)+i-1];}else{(*Olm[36])(i,j)=mo63sr[30*(i-1)+j-1];}
if(i<j){(*Olm[37])(i,j)=mo62si[30*(j-1)+i-1];}else{(*Olm[37])(i,j)=mo62sr[30*(i-1)+j-1];}
if(i<j){(*Olm[38])(i,j)=mo61si[30*(j-1)+i-1];}else{(*Olm[38])(i,j)=mo61sr[30*(i-1)+j-1];}
if(i<j){(*Olm[39])(i,j)=mo60ci[30*(j-1)+i-1];}else{(*Olm[39])(i,j)=mo60cr[30*(i-1)+j-1];}
if(i<j){(*Olm[40])(i,j)=mo61ci[30*(j-1)+i-1];}else{(*Olm[40])(i,j)=mo61cr[30*(i-1)+j-1];}
if(i<j){(*Olm[41])(i,j)=mo62ci[30*(j-1)+i-1];}else{(*Olm[41])(i,j)=mo62cr[30*(i-1)+j-1];}
if(i<j){(*Olm[42])(i,j)=mo63ci[30*(j-1)+i-1];}else{(*Olm[42])(i,j)=mo63cr[30*(i-1)+j-1];}
if(i<j){(*Olm[43])(i,j)=mo64ci[30*(j-1)+i-1];}else{(*Olm[43])(i,j)=mo64cr[30*(i-1)+j-1];}
if(i<j){(*Olm[44])(i,j)=mo65ci[30*(j-1)+i-1];}else{(*Olm[44])(i,j)=mo65cr[30*(i-1)+j-1];}
if(i<j){(*Olm[45])(i,j)=mo66ci[30*(j-1)+i-1];}else{(*Olm[45])(i,j)=mo66cr[30*(i-1)+j-1];}

if(i<j){(*Olm[46])(i,j)=modxci[30*(j-1)+i-1];}else{(*Olm[46])(i,j)=modxcr[30*(i-1)+j-1];}
if(i<j){(*Olm[47])(i,j)=modyci[30*(j-1)+i-1];}else{(*Olm[47])(i,j)=modycr[30*(i-1)+j-1];}
if(i<j){(*Olm[48])(i,j)=modzci[30*(j-1)+i-1];}else{(*Olm[48])(i,j)=modzcr[30*(i-1)+j-1];}

if(i<j){(*Olm[49])(i,j)=modxdxci[30*(j-1)+i-1];}else{(*Olm[49])(i,j)=modxdxcr[30*(i-1)+j-1];}
if(i<j){(*Olm[50])(i,j)=modydyci[30*(j-1)+i-1];}else{(*Olm[50])(i,j)=modydycr[30*(i-1)+j-1];}
if(i<j){(*Olm[51])(i,j)=modzdzci[30*(j-1)+i-1];}else{(*Olm[51])(i,j)=modzdzcr[30*(i-1)+j-1];}
   }}
 pr=0;
//---------------------------------------------------------------------------
	Ri= new ComplexMatrix * [1+NOF_RIXS_MATRICES];
         for(i=1;i<=NOF_RIXS_MATRICES;++i){Ri[i]=new ComplexMatrix(1,dimj,1,dimj);(*Ri[i])=0;}
  In = new Matrix * [1+IONPARS_MAXNOFCOMPONENTS];
  for(i=0;i<=IONPARS_MAXNOFCOMPONENTS;++i)In[i]= new Matrix(1,dimj,1,dimj);
//---------------------------------------------------------------------------
 if (perlp){
   ComplexMatrix * OOlm  [1+NOF_OLM_MATRICES]; 
   // define memory for all matrices 
   for(i=1;i<=NOF_OLM_MATRICES;++i)OOlm [i] = new ComplexMatrix(1,dimj,1,dimj); 
   // fill ComplexMatrix OOlm with values (for perl parsing)
   for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
    for(int k=1;k<=NOF_OLM_MATRICES;++k)if(i<j){(*OOlm[k])(i,j)=complex<double>((*Olm[k])(j,i),-(*Olm[k])(i,j));}
                                           else{(*OOlm[k])(i,j)=complex<double>((*Olm[k])(i,j),(*Olm[k])(j,i));}

   // now fill interaction operators IIn with values
   ComplexMatrix * IIn [1+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES+3];
   for(i=0;i<=IONPARS_MAXNOFCOMPONENTS;++i)IIn[i]= new ComplexMatrix(1,dimj,1,dimj);
   // standard operator sequence I1,....,I51
   // module so1ion: Jx Jy Jz O22S O21S O20 O21 O22 O33S O32S .... O66 Jx^2 Jy^2 Jz^2 Jx^4 Jy^4 Jz^4
   (*IIn[0])=1;(*IIn[1])=Jaa;(*IIn[2])=Jbb;(*IIn[3])=Jcc;
   for(i=4;i<=IONPARS_MAXNOFCOMPONENTS;++i)(*IIn[i])=(*OOlm[i-3]);
   for(i=1;i<=NOF_OLM_MATRICES;++i)IIn[i+IONPARS_MAXNOFCOMPONENTS]=OOlm[i]; // only give pointer 
   for(i=1;i<=NOF_RIXS_MATRICES;++i)IIn[i+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES]=Ri[i];
   IIn[1+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES]=&Jaa;
   IIn[2+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES]=&Jbb;
   IIn[3+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES]=&Jcc;
    char * operatornames[2+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES+3];
    char opnam[]="one\0";   
    operatornames[0]=opnam;
    for (i=1;i<=IONPARS_MAXNOFCOMPONENTS;++i){operatornames[i]=new char[4];sprintf(operatornames[i],"I%i",i);
                                             }
    for (i=1;i<=NOF_OLM_MATRICES;++i){operatornames[IONPARS_MAXNOFCOMPONENTS+i]=new char[5];
                                      operatornames[IONPARS_MAXNOFCOMPONENTS+i][0]='O';operatornames[IONPARS_MAXNOFCOMPONENTS+i][4]='\0';
                                      if(i>45)operatornames[IONPARS_MAXNOFCOMPONENTS+i][0]='J'; // Jx2 Jy2 Jz2 Jx4 Jy4 Jz4
                              strncpy(operatornames[IONPARS_MAXNOFCOMPONENTS+i]+1,kq+i*3,3);
                                   if(operatornames[IONPARS_MAXNOFCOMPONENTS+i][3]==' ')
                                      operatornames[IONPARS_MAXNOFCOMPONENTS+i][3]='\0';
                                     }
    char Rnam[]="R11\0R12\0R13\0R21\0R22\0R23\0R31\0R32\0R33\0";
    for (i=1;i<=NOF_RIXS_MATRICES;++i)operatornames[IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+i]=Rnam+4*(i-1);
    char Jnam[]="Jx\0Jy\0Jz\0";
    operatornames[IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+1]=Jnam;
    operatornames[IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+2]=Jnam+3;
    operatornames[IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+3]=Jnam+6;

    operatornames[IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES+3+1]=NULL;

                 char numnam[]="GJ     \0"
                               "ALPHA  \0BETA   \0GAMMA  \0"
                               "R2     \0R4     \0R6     \0"
                               "SIGMA0r\0SIGMA1r\0SIGMA2r\0"
                               "SIGMA0i\0SIGMA1i\0SIGMA2i\0"
                               "Dx2    \0Dy2    \0Dz2    \0"
                               "Dx4    \0Dy4    \0Dz4    \0"
                                 ;
                 char * numbernames[113];
                 for (i=0;i<19;++i)numbernames[i]=numnam+8*i;
                 for (i=0;i<=45;++i){ numbernames[19+i]=new char[5];
                                      numbernames[19+i][0]='B';numbernames[19+i][4]='\0';
                              strncpy(numbernames[19+i]+1,kq+i*3,3);
                                   if(numbernames[19+i][3]==' ')
                                      numbernames[19+i][3]='\0';
                                        }
                 for (i=0;i<=45;++i){ numbernames[65+i]=new char[5];
                                      numbernames[65+i][0]='L';numbernames[65+i][4]='\0';
                              strncpy(numbernames[65+i]+1,kq+i*3,3);
                                   if(numbernames[65+i][3]==' ')
                                      numbernames[65+i][3]='\0';
                                        }
                 numbernames[111]=NULL;
                 double *numbers[111];
                 numbers[0]=&gJ;
                 numbers[1]=&alpha;
                 numbers[2]=&beta;
                 numbers[3]=&gamma;
                 numbers[4]=&r2;
                 numbers[5]=&r4;
                 numbers[6]=&r6;
                 numbers[7]=&s0r;
                 numbers[8]=&s1r;
                 numbers[9]=&s2r;
                 numbers[10]=&s0i;
                 numbers[11]=&s1i;
                 numbers[12]=&s2i;
                 numbers[13]=&Blm[46];
                 numbers[14]=&Blm[47];
                 numbers[15]=&Blm[48];
                 numbers[16]=&Blm[49];
                 numbers[17]=&Blm[50];
                 numbers[18]=&Blm[51];
                 for(i=0;i<=45;++i){numbers[19+i]=&Blm[i];numbers[65+i]=&Llm[i];}

                 char * strings[3];char * stringnames[3];
                 char strnam[]="IONTYPE\0";   
                 stringnames[0]=strnam;stringnames[1]=NULL;
                 strings[0]=iontype;


                if(perlparse(cffilename,numbers,numbernames,strings,stringnames,IIn,operatornames)==false)
                     {printf("Error perl parsing sipf file %s\n",cffilename);exit(EXIT_FAILURE);}
  // now fill interaction operators In with values
  for(l=1;l<=IONPARS_MAXNOFCOMPONENTS;++l)for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j)
   if(i<j){(*In[l])(i,j)=imag((*IIn[l])(j,i));}else{(*In[l])(i,j)=real((*IIn[l])(i,j));}
  //perl parsing allows to redefine the following operators: In, Rij, Olm, Jx Jy Jz.
  // therefore reread these operators, too and put them into Matrix.
  for(i=1;i<=dimj;++i)for(j=1;j<=dimj;++j){
  for(l=1;l<=NOF_OLM_MATRICES;++l)
   if(i<j){(*Olm[l])(i,j)=imag((*IIn[IONPARS_MAXNOFCOMPONENTS+l])(j,i));}
      else{(*Olm[l])(i,j)=real((*IIn[IONPARS_MAXNOFCOMPONENTS+l])(i,j));}
  for(l=1;l<=NOF_RIXS_MATRICES;++l)
   if(i<j){(*Ri[l])(i,j)=imag((*IIn[IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+l])(j,i));}
      else{(*Ri[l])(i,j)=real((*IIn[IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+l])(i,j));}
  if(i<j){Ja (i,j)=imag((*IIn[1+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES])(j,i));}
     else{Ja(i,j)=real((*IIn[1+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES])(i,j));}
  if(i<j){Jb(i,j)=imag((*IIn[2+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES])(j,i));}
     else{Jb(i,j)=real((*IIn[2+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES])(i,j));}
  if(i<j){Jc(i,j)=imag((*IIn[3+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES])(j,i));}
     else{Jc(i,j)=real((*IIn[3+IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES+NOF_RIXS_MATRICES])(i,j));}
                                           }
  for (i=1;i<=NOF_OLM_MATRICES;++i)  {delete OOlm[i];}
  for(i=0;i<=IONPARS_MAXNOFCOMPONENTS;++i){delete IIn[i];}
  for (i=0;i<=45;++i){delete numbernames[16+i];
                      delete numbernames[62+i];
                       }
  for (i=1;i<=IONPARS_MAXNOFCOMPONENTS+NOF_OLM_MATRICES;++i)delete operatornames[i];
               } //if perlp
  else
  {
  // now fill interaction operators In with standard values
    (*In[1])=Ja;
    (*In[2])=Jb;
    (*In[3])=Jc;
    for(i=1;i<=NOF_OLM_MATRICES;++i){(*In[i+3])=(*Olm[i]);}
  }
//myPrintMatrix(stdout,(*In[51]));
if(ialpha&&fabs(alphar-alpha)/fabs(alphar+1)>SMALL_DEVIATION) {fprintf(stderr,"#Warning module %s internal value for Stevens Parameter (alpha=%g) different from input file (alpha=%g), using value from input file\n",moduletype,alphar,alpha);}
else{alpha=alphar;}
if(ibeta&&fabs(betar-beta)/fabs(betar+1)>SMALL_DEVIATION) {fprintf(stderr,"#Warning module %s internal value for Stevens Parameter (beta=%g) different from input file (beta=%g), using value from input file\n",moduletype,betar,beta);}
else{beta=betar;}
if(igamma&&fabs(gammar-gamma)/fabs(gammar+1)>SMALL_DEVIATION) {fprintf(stderr,"#Warning module %s internal value for Stevens Parameter (gamma=%g) different from input file (gamma=%g), using value from input file\n",moduletype,gammar,gamma);}
else{gamma=gammar;}
if(igj&&fabs(gJr-gJ)/fabs(gJr+1)>SMALL_DEVIATION) {fprintf(stderr,"#Warning module %s internal value for Lande Factor (gJ=%g) different from input file (gJ=%g), using value from input file\n",moduletype,gJr,gJ);}
else{gJ=gJr;}
if(ir2&&fabs(r2r-r2)/fabs(r2r+1)>SMALL_DEVIATION) {fprintf(stderr,"#Warning module %s internal value for radial Matrix element (<r2>=%g) different from input file (<r2>=%g), using value from input file\n",moduletype,r2r,r2);}
else{r2=r2r;}
if(ir4&&fabs(r4r-r4)/fabs(r4r+1)>SMALL_DEVIATION) {fprintf(stderr,"#Warning module %s internal value for radial Matrix element (<r4>=%g) different from input file (<r4>=%g), using value from input file\n",moduletype,r4r,r4);}
else{r4=r4r;}
if(ir6&&fabs(r6r-r6)/fabs(r6r+1)>SMALL_DEVIATION) {fprintf(stderr,"#Warning module %s internal value for radial Matrix element (<r6>=%g) different from input file (<r6>=%g), using value from input file\n",moduletype,r6r,r6);}
else{r6=r6r;}
if(inof&&nof_electronsr-nof_electrons!=0){fprintf(stderr,"#Warning module %s internal value for nr of electrons (=%i) different from input file (=%i), using value from input file\n",moduletype,nof_electronsr,nof_electrons);}
else{nof_electrons=nof_electronsr;}
//---------------------------------------------------------------------------
      // here fill the matrices Ri[1...9] with the 11 12 13 21 22 23 31 32 33
      // matrices of the RIXS scattering operator R using spherical symmetry (M Haverkort PRL)
      // (only if sigma is nonzero)
if(fabs(s0r)+fabs(s1r)+fabs(s2r)+fabs(s0i)+fabs(s1i)+fabs(s2i)>SMALL_DEVIATION)
{   sigma0=complex<double>(s0r,s0i);
   sigma1=complex<double>(s1r,s1i);
   sigma2=complex<double>(s2r,s2i);
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
}
//---------------------------------------------------------------------------
  // ------------------------------------------------------------
  // here transform the Llm (if present) to Blm ...
  Vector thetaJ(0,6);thetaJ(0)=nof_electrons;thetaJ(2)=alpha;thetaJ(4)=beta;thetaJ(6)=gamma;

   fprintf(stderr,"#crystal field parameters:\n");  
   const char lm[]="B00 B22SB21SB20 B21 B22 B33SB32SB31SB30 B31 B32 B33 B44SB43SB42SB41SB40 B41 B42 B43 B44 B55SB54SB53SB52SB51SB50 B51 B52 B53 B54 B55 B66SB65SB64SB63SB62SB61SB60 B61 B62 B63 B64 B65 B66 Dx2 Dy2 Dz2 Dx4 Dy4 Dz4";
   char lm4[5];lm4[4]='\0';
   for(i=0;i<=NOF_OLM_MATRICES;++i){strncpy(lm4,lm+i*4,4);l=lm4[1]-48;m=lm4[2]-48;if(lm4[3]=='S'){m=-m;}
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
  // ------------------------------------------------------------
  // here set the crystal field Hamiltonian Matrix
  Hcf= Matrix(1,dimj,1,dimj); 
  Hcf=0;if(Hcf==(double)0.0){for(l=1;l<=NOF_OLM_MATRICES;++l){Hcf+=Blm(l)*(*Olm[l]);}}

// for compatibility
 //cf_file=fopen_errchk(cffilename,"r");
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
  if(alpha*alpha+beta*beta+gamma*gamma>1e-10){   fprintf(file,"#----------------\n# Stevens Factors\n#----------------\nALPHA=%g\nBETA=%g\nGAMMA=%g\n\n",alpha,beta,gamma); }
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

   if(Blm(49)!=0){fprintf(outfile,"Dx4=%g\n",myround(Blm(49)));}
   if(Blm(50)!=0){fprintf(outfile,"Dy4=%g\n",myround(Blm(50)));}
   if(Blm(51)!=0){fprintf(outfile,"Dz4=%g\n",myround(Blm(51)));}
}

void ionpars::savLlm(FILE * outfile)
{fprintf(outfile,"units=meV\n");
   char lm3[4];lm3[3]='\0';
   const char lm[]="00 22S21S20 21 22 33S32S31S30 31 32 33 44S43S42S41S40 41 42 43 44 55S54S53S52S51S50 51 52 53 54 55 66S65S64S63S62S61S60 61 62 63 64 65 66 ";
   for(int i=0;i<=45;++i){strncpy(lm3,lm+i*3,3);
   if(Llm(i)!=0){fprintf(outfile,"L%s=%g\n",lm3,myround(Llm(i)));}
                         }
}
 // -------------------------------------------------------------------------------------
 // -------------------------------------------------------------------------------------
 //                        Module  Functions 
 // -------------------------------------------------------------------------------------
 // -------------------------------------------------------------------------------------

/**************************************************************************/
void ionpars::cfeigenstates(ComplexMatrix *eigenstates,Vector &  Hxc,Vector & Hext, double & T)
{   /*on input
    gJmbH	vector of effective field [meV]
      on output
    Matrix containing the eigenvalues and eigenfunctions of the crystalfield problem
    eigenvalues ares stored as real part of row zero
    boltzmann population numbers are stored as imaginary part of row zero
*/
   int i,j,sort=1,dj=Hcf.Rhi();
   (*eigenstates) = ComplexMatrix(0,dj,1,dj);
   Vector En(1,dj),wn(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   setup_and_solve_Hamiltonian(Hxc,Hext,En,zr,zi,sort);
   
   for(i=1;i<=dj;++i)for(j=1;j<=dj;++j)(*eigenstates)(i,j)=complex <double> (zr(i,j),zi(i,j));
 
    //calculate partition sum
    double zz;
    calculate_Z_wn(En,T,zz,wn);
     // put boltzmann population into row 0 of eigenstates...
     for(j=1;j<=dj;++j){(*eigenstates)(0,j)=complex<double>(En(j),wn(j));}
}

/**************************************************************************/
void ionpars::Icalc(Vector & I,double & T, Vector &  Hxc,Vector & Hext, double & lnZs, double & U, ComplexMatrix & /*ests*/)
{   /*on input
    T		temperature[K]
    Hxc	        vector of exchange field [meV]
    Hext        external magnetic field [T]
  on output    
    I		single ion momentum vector <I> (if T>0 thermal exp value <J>T 
                                                if T<0 the program asks for w_n and calculates
						       exp value <J>=sum_n w_n <n|I|n>
    Z		single ion partition function
    U		single ion magnetic energy
    */
   int dj=Hcf.Rhi(),sort=0;if (T<0) sort=1;
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   setup_and_solve_Hamiltonian(Hxc,Hext,En,zr,zi,sort);
   // calculate Z and wn (occupation probability)
   Vector wn(1,dj); double Zs;
   calculate_Z_wn(En,T,Zs,lnZs,wn);
   // calculate U
   U=En*wn;
   // calculate <I1>,<I2>,<I3>
   I=0;
    for (int i=1;i<=dj;++i)
    { if(wn(i)>SMALL_PROBABILITY){
     // here the expectation values of the multipolar moments are calculated
     for(int j=1;j<=I.Hi();++j)I[j]+=wn(i)*matelr(i,i,zr,zi,(*In[j]));
                                   }
    } //printf("%g %g\n",I[1],I[2]);
}


/**************************************************************************/
// for mcdisp this routine is needed
int ionpars::du1calc(int & tn,double & T,Vector &  Hxc,Vector & Hext,ComplexVector & u1,float & delta,int & n, int & nd,ComplexMatrix & ests)
{  /*on input
    tn      ... number of transition to be computed 
    sign(tn)... 1... without printout, -1 with extensive printout
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta-+	energy of transition [meV]
    u1(i)	<-|Ii|+> sqrt(n+-n-),  n+,n-
    .... occupation number of states (- to + transition chosen according to transitionnumber)
*/
   double ninit=u1[1].real();
   double pinit=u1[1].imag();
  int pr=0;if (tn<0) {pr=1;tn*=-1;}
  int i,j,dj=Hcf.Rhi();
  // set eigenvectors
   Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   for(i=1;i<=dj;++i)for(j=1;j<=dj;++j){zr(i,j)=real(ests(i,j));zi(i,j)=imag(ests(i,j));}
  // calculate mat and delta for transition number tn
  // 1. get i and j  delta=En(j)-En(i) from tn
  getijdelta_from_transitionnumber(i,j,delta,dj,tn,pr,ests);
  n=i;nd=j;
  char optype[5];
  for(int l=1;l<=Hxc.Hi();++l){sprintf(optype,"I%i",l);
  u1(l)=observable1(i,j,delta,zr,zi,T,ests,pr,optype,(*In[l]));}

// return number of all transitions     
     return noft(ests,T,pinit,ninit);
}

Matrix ionpars::opmat (int & n ,Vector &  Hxc,Vector & Hext)
   // on input
   // n		which operator 0=Hamiltonian, 1,2,3=J1,J2,J3
   // Hext,Hxc	vector of external and exchange field [meV]
   // on output    
   // operator matrix of Hamiltonian, I1, I2, I3 depending on n
{
if(n==0){// return Hamiltonian
Vector gjmbH(1,max(3,Hxc.Hi()));gjmbH=0;
for(int i=1;i<=Hxc.Hi();++i)gjmbH(i)=Hxc(i);
gjmbH(1)+=gJ*MU_B*Hext(1);
gjmbH(2)+=gJ*MU_B*Hext(2);
gjmbH(3)+=gJ*MU_B*Hext(3);
// check dimensions of vector
if(gjmbH.Hi()>NOF_OLM_MATRICES+3)
   {fprintf(stderr,"Error module so1ion/cfield: dimension of exchange field=%i > %i - check number of columns in file mcphas.j\n",gjmbH.Hi(),NOF_OLM_MATRICES+3);
    exit(EXIT_FAILURE);}
   (*In[n])=Hcf;for(int j=1;j<=gjmbH.Hi();++j){(*In[n])-=gjmbH(j)*(*In[j]);}
  }
else
 {if(n>NOF_OLM_MATRICES+3)
   {fprintf(stderr,"Error module so1ion/cfield: operatormatrix index=%i > %i - check number of columns in file mcphas.j\n",n,NOF_OLM_MATRICES+3);
    exit(EXIT_FAILURE);}
 }
return (*In[n]);
}


// *************************************************************************************************
// ************************* private helper functions *********************************************
// ************************************************************************************************
void ionpars::setup_and_solve_Hamiltonian(Vector &  Hxc,Vector & Hext,Vector & En,Matrix & zr,Matrix & zi,int sort)
{
Vector gjmbH(1,max(3,Hxc.Hi()));gjmbH=0;
for(int i=1;i<=Hxc.Hi();++i)gjmbH(i)=Hxc(i);
gjmbH(1)+=gJ*MU_B*Hext(1);
gjmbH(2)+=gJ*MU_B*Hext(2);
gjmbH(3)+=gJ*MU_B*Hext(3);

// check dimensions of vector
if(gjmbH.Hi()>NOF_OLM_MATRICES+3)
   {fprintf(stderr,"Error module so1ion/cfield: dimension of exchange field=%i > %i - check number of columns in file mcphas.j\n",gjmbH.Hi(),NOF_OLM_MATRICES+3);
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
   int j,dj=Hcf.Rhi();
   Matrix Ham(1,dj,1,dj);    
   Ham=Hcf;for(j=1;j<=gjmbH.Hi();++j){Ham-=gjmbH(j)*(*In[j]);}
   // diagonalize
   int maxiter=1000000;
   EigenSystemHermitean (Ham,En,zr,zi,sort,maxiter);

}

void ionpars::calculate_Z_wn(Vector & En,double & T,double & Zs,Vector & wn)
{double lnZs;
 calculate_Z_wn(En,T,Zs,lnZs,wn);
}
void ionpars::calculate_Z_wn(Vector & En,double & T,double & Zs,double & lnZs,Vector & wn)
{   // calculate Z and wn (occupation probability)
     int i,dj=wn.Hi();
     double y,x=Min(En);
     if (T>0)
     { for (i=1;i<=dj;++i)
       {if ((y=(En(i)-x)/KB/T)<600) wn[i]=exp(-y); 
        else wn[i]=0.0;
       }
       Zs=Sum(wn);wn/=Zs;
       lnZs=log(Zs)-x/KB/T;    
       Zs*=exp(-x/KB/T);
     } 
     if (T==0)
     { printf ("Temperature T==0: please choose probability distribution of states by hand\n");
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
     if (T<0){Zs=1;lnZs=0;wn=0;wn((int)(-T))=1;}     
}


int ionpars::noft(ComplexMatrix & est,double & T,double & pinit,double & ninit)
{// calculate number of transitions
   int dj=est.Rhi();
   double n=ninit;
   if (n>dj)n=dj;
   //if (pinit<SMALL_PROBABILITY)pinit=SMALL_PROBABILITY;
   double zsum=0,zii,x;
   int noft=0;
   if(T>0)for(int i=1;(i<=n)&((((x=(real(est(0,i))-real(est(0,1)))/KB/T)<200)? zii=exp(-x):zii=0)>=(pinit*zsum));++i)
   {noft+=dj-i+1;zsum+=zii;}
   
   if(T<0)for(int i=1;i<=n;++i)
   {noft+=dj-i+1;}
//printf("dj=%i  ninit=%g pinit=%g T=%g noft=%i\n",dj,ninit,pinit,T,noft);
   return noft;
}



void ionpars::popnr_diff(ComplexVector& dMQ,int &i,int &j,ComplexMatrix & est,float & delta,double & T,int &pr,const char * n)
{if (delta>SMALL_QUASIELASTIC_ENERGY)
   { if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      for(int l=1;l<=dMQ.Hi();++l)printf(" |<%i|%s%i-<%s%i>|%i>|^2=%4.4g",i,n,l,n,l,j,abs(dMQ(l))*abs(dMQ(l)));
      printf(" n%i-n%i=%4.4g\n",i,j,imag(est(0,i))-imag(est(0,j)));}
      dMQ*=sqrt(imag(est(0,i))-imag(est(0,j))); // occupation factor
   }else
   {// quasielastic scattering has not wi-wj but wj*epsilon/kT
     if(pr==1){
      printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
      for(int l=1;l<=dMQ.Hi();++l)printf(" |<%i|%s%i-<%s%i>|%i>|^2=%4.4g",i,n,l,n,l,j,abs(dMQ(l))*abs(dMQ(l)));
      printf(" n%i=%4.4g\n",i,imag(est(0,i)));
      }
    dMQ*=sqrt(imag(est(0,i))/KB/fabs(T));
   }
}
// ----------------------------------------------------------------------

void ionpars::getijdelta_from_transitionnumber(int & i,int & j,float & delta,int & dj,int & tn,int & pr,ComplexMatrix &ests)
{int  k=0;for(i=1;i<=dj;++i){for(j=i;j<=dj;++j){++k;if(k==tn)break;}if(k==tn)break;} 
  if(j==i){delta=-SMALL_QUASIELASTIC_ENERGY; //if transition within the same level: take negative delta !!- this is needed in routine intcalc
          }else{delta=real(ests(0,j))-real(ests(0,i));
                if (delta<-0.000001){fprintf(stderr,"ERROR module so1ion or cfield.so - "
                                     "d*calc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
           }
  if(pr==1){printf("delta(%i->%i)=%4.4gmeV",i,j,delta);
            if(delta>SMALL_QUASIELASTIC_ENERGY){printf(" n%i-n%i=%4.4g",i,j,imag(ests(0,i))-imag(ests(0,j)));
           }else{ printf(" n%i=%4.4g",i,imag(ests(0,i)));}
           }
}

// ----------------------------------------------------------------------
// calculates the transition matrix element <i|O-<O>|j>sqrt((ni-nj))
  complex<double> ionpars::observable1(int & i,int & j,float & delta,Matrix & zr,Matrix & zi,
                         double & T,ComplexMatrix&est,int & pr,const char *optype,Matrix & O)
  {static complex<double> ret;
   ret=complex<double>(matelr(i,j,zr,zi,O),mateli(i,j,zr,zi,O));
   double Oav=0;int dj=Hcf.Rhi();
   // calculate and subtract expectation value if T>0
   if(T>0&&i==j){for(int l=1;l<=dj;++l)if(imag(est(0,l))>SMALL_PROBABILITY)Oav+=imag(est(0,l))*matelr(l,l,zr,zi,O);
                 ret-=Oav; 
                }
   if(pr==1)printf(" |<%i|%s-<%s>|%i>|^2=%4.4g",i,optype,optype,j,abs(ret)*abs(ret));
   if (delta>SMALL_QUASIELASTIC_ENERGY){ret*=sqrt(imag(est(0,i))-imag(est(0,j))); // occupation factor
                                  }else{ret*=sqrt(imag(est(0,i))/KB/fabs(T));
                                       }
  return ret;
  }


/**************************************************************************/
//                          OBSERVABLES
/**************************************************************************/
void ionpars::Jcalc(Vector & JJ,double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & /*ests*/)
{   /*on input
    T		temperature[K]
    Hxc	vector of exchange field [meV]
    Hext        external magnetic field [T]
  on output    
    JJ		single ion momentum vector <J> (if T>0 thermal exp value <J>T 
                                                if T<0 the program asks for w_n and calculates
						       exp value <J>=sum_n w_n <n|J|n>
    */
   int dj=Hcf.Rhi(),sort=0;if (T<0) sort=1;
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   setup_and_solve_Hamiltonian(Hxc,Hext,En,zr,zi,sort);
   // calculate Z and wn (occupation probability)
   Vector wn(1,dj); double Zs,lnZs;
   calculate_Z_wn(En,T,Zs,lnZs,wn);
   // calculate <Ja>,<Jb>,<Jc>
   JJ=0;
    for (int i=1;i<=dj;++i)
    { if(wn(i)>SMALL_PROBABILITY){
     JJ(1)+=wn(i)*matelr(i,i,zr,zi,Ja);
     JJ(2)+=wn(i)*matelr(i,i,zr,zi,Jb);
     JJ(3)+=wn(i)*matelr(i,i,zr,zi,Jc);   
                                   }
    }
}
//**********************************************************************/
// routine to calculate the charge density coefficients of Zlm() R(r)^2
// *********************************************************************
void ionpars::chargedensity_coeffcalc(Vector &mom, double & T, Vector &  Hxc,Vector & Hext, ComplexMatrix & parstorage)
{/*on input
    T		temperature[K]
    Hxc 	vector of exchange field [meV]
  on output    
    mom		chargedensity coefficients: (if T>0 thermal exp value <mom>T 
                                                if T<0 the program asks for w_n and calculates
						       exp value <mom>=sum_n w_n <n|mom|n>						       
*/
   int i,j,dj=Hcf.Rhi();
   Vector En(1,dj);Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   int sort=0;if (T<0) sort=1;
   setup_and_solve_Hamiltonian(Hxc,Hext,En,zr,zi,sort); 
   Vector wn(1,dj);double Zs;
   calculate_Z_wn(En,T,Zs,wn);

   if(nof_electrons==0){fprintf(stderr,"Error so1ion/cfield: nof_electrons=0"
                              " ... perhaps single ion property file does not contain"
                              " the number of electrons in the shell: 'nof_electrons=...'\n");
                             exit(EXIT_FAILURE);}
   mom=0;
   mom(1) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
// a(0, 0) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
// Indices for chargedensity
//          0 not used
//          0 1  2  3 4 5 6  7  8  9 101112131415 16 17 18 19 20 2122232425262728 
int k[] = {-1,0, 2, 2,2,2,2, 4, 4, 4, 4,4,4,4,4,4, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-2,-1,0,1,2,-4,-3,-2,-1,0,1,2,3,4,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

// here the expectation values of the multipolar moments are calculated
for (i=1;i<=dj;++i)if(wn(i)>SMALL_PROBABILITY){
 for(j=4;j<=8;++j) mom[j-2]+=alpha*cnst(k[j-2],q[j-2])*wn(i)*matelr(i,i,zr,zi,(*Olm[j-3]));
 for(j=16;j<=24;++j)mom[j-9]+=beta*cnst(k[j-9],q[j-9])*wn(i)*matelr(i,i,zr,zi,(*Olm[j-3]));
 for(j=36;j<=48;++j)mom[j-20]+=gamma*cnst(k[j-20],q[j-20])*wn(i)*matelr(i,i,zr,zi,(*Olm[j-3]));
                                              }
     // theta_J*cnst(l,m)  are prefactors to get coefficients of Zlm*R(r)^2 
    //in case of module cfield and so1ion(stevens parameters tetan and zlm prefactors)
}

//**********************************************************************/
// routine to calculate the scattering operator to go beyond dip approx
// *********************************************************************
// calculate scattering operator <M(Q)>=-2x<Q>_TH in units of mb according to S Lovesey & Balcar
// according to stored eigenstate matrix est
// calculates the scattering operator given the polar angles th, ph (with respect to the CEF coordinate 
// system xyz and the <jl(qr)> and the eigenstate matrix with eigenstates and thermal population numbers
ComplexVector & ionpars::MQ(double th, double ph,double J0,double J2,double J4,double J6, Vector & Zc,ComplexMatrix & est)
    {   int dj=(int)(2*J+1);       
	 ComplexMatrix MQXM(1,dj,1,dj),MQYM(1,dj,1,dj),MQZM(1,dj,1,dj);
         MQM(MQXM,MQYM,MQZM,th,ph,J0,J2,J4,J6,Zc);
							     // ... calculate thermal expectation values
							     // using the eigenstates and T 
							     // mom(1) = ....
       static ComplexVector mm(1,3); mm=0;
       int K,M,Md;
       for(K=1;K<=dj;++K){for(M=1;M<=dj;++M){for(Md=1;Md<=dj;++Md){
         mm(1)+=imag(est(0,K))*conj(est(M,K))*MQXM(M,Md)*est(Md,K); 
         mm(2)+=imag(est(0,K))*conj(est(M,K))*MQYM(M,Md)*est(Md,K); 
         mm(3)+=imag(est(0,K))*conj(est(M,K))*MQZM(M,Md)*est(Md,K); 
       }}}					       
     mm*=2; // this is now <M(Q)>=-2x<Q>_TH in units of mb
// myPrintComplexMatrix(stdout,est);}
// myPrintComplexVector(stdout,mm);//equivalent to moment ...
    return mm;
    }
// --------------- some helper functions for that ---------------------------------------------
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
                                                                   } // for Qd
                                            } // factor != 0
                            } //for K
}


/**************************************************************************/
// for mcdisp this routine is needed
int ionpars::dJ1calc(int & tn,double & T,Vector &  Hxc,Vector & Hext,ComplexVector & J1,float & delta,ComplexMatrix & ests)
{  /*on input
    tn      ... number of transition to be computed 
    sign(tn)... 1... without printout, -1 with extensive printout
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta-+	energy of transition [meV]
    J1(i)	<-|Ji|+> sqrt(n+-n-),  n+,n-
    .... occupation number of states (- to + transition chosen according to transitionnumber)
*/
   double ninit=J1[1].real();
   double pinit=J1[1].imag();
  int pr=0;if (tn<0) {pr=1;tn*=-1;}
  int i,j,dj=Hcf.Rhi();
  // set eigenvectors
   Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   for(i=1;i<=dj;++i)for(j=1;j<=dj;++j){zr(i,j)=real(ests(i,j));zi(i,j)=imag(ests(i,j));}
  // calculate mat and delta for transition number tn
  // 1. get i and j  delta=En(j)-En(i) from tn
  getijdelta_from_transitionnumber(i,j,delta,dj,tn,pr,ests);
  J1(1)=observable1(i,j,delta,zr,zi,T,ests,pr,"Ja",Ja);
  J1(2)=observable1(i,j,delta,zr,zi,T,ests,pr,"Jb",Jb);
  J1(3)=observable1(i,j,delta,zr,zi,T,ests,pr,"Jc",Jc);



// return number of all transitions     
     return noft(ests,T,pinit,ninit);
}
//**********************************************************************/
// routine to calculate the transition matrix elements of the charge density coefficients (of Zlm() R(r)^2)
// *********************************************************************
int ionpars::dchargedensity_coeff1calc(int & tn,double & T,Vector &  Hxc,Vector & Hext, ComplexVector & cd1,float & delta,ComplexMatrix & ests)
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
   double ninit=cd1[1].real();
   double pinit=cd1[1].imag();
  int pr=0;if (tn<0) {pr=1;tn*=-1;}
  int i,j,l,dj=Hcf.Rhi();
  // set eigenvectors
   Matrix zr(1,dj,1,dj);Matrix zi(1,dj,1,dj);
   for(i=1;i<=dj;++i)for(j=1;j<=dj;++j){zr(i,j)=real(ests(i,j));zi(i,j)=imag(ests(i,j));}
  // calculate mat and delta for transition number tn
  // 1. get i and j  delta=En(j)-En(i) from tn
  getijdelta_from_transitionnumber(i,j,delta,dj,tn,pr,ests);

// a(0, 0) = nof_electrons / sqrt(4.0 * 3.1415); // nofelectrons 
// Indices for chargedensity
//            0 not used
//          0 1  2  3 4 5 6  7  8  9 101112131415 16 17 18 19 20 2122232425262728 
int k[] = {-1,0, 2, 2,2,2,2, 4, 4, 4, 4,4,4,4,4,4, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-2,-1,0,1,2,-4,-3,-2,-1,0,1,2,3,4,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

char optype[12];
cd1(1)=0;if(i==j){cd1(1)=-nof_electrons / sqrt(4.0 * 3.1415);}
for(l=2;l<=6;++l){sprintf(optype,"cd_coeff%i",l);
                  cd1(l)=alpha*cnst(k[l],q[l])*observable1(i,j,delta,zr,zi,T,ests,pr,optype,(*Olm[l-1]));}
for(l=7;l<=15;++l){sprintf(optype,"cd_coeff%i",l);
                   cd1(l)=beta*cnst(k[l],q[l])*observable1(i,j,delta,zr,zi,T,ests,pr,optype,(*Olm[l-1+7]));}
for(l=16;l<=28;++l){sprintf(optype,"cd_coeff%i",l);
                    cd1(l)=gamma*cnst(k[l],q[l])*observable1(i,j,delta,zr,zi,T,ests,pr,optype,(*Olm[l-1+7+11]));}
     // theta_J*cnst(l,m)  are prefactors to get coefficients of Zlm*R(r)^2 
    //in case of module cfield and so1ion(stevens parameters tetan and zlm prefactors)


// return number of all transitions     
     return noft(ests,T,pinit,ninit);
}   


int ionpars::dMQ1(int & tn,double & th,double & ph,double & J0,double & J2,double & J4,double & J6,Vector & Zc,ComplexMatrix & est,double & T,ComplexVector & dMQ)
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
  int i,j=1,l;
  int dj=(int)(2*J+1);
  float delta;
   double ninit=dMQ(1).real();
   double pinit=dMQ(1).imag();

 // calculate mat and delta for transition number tn
  // 1. get i and j  delta=En(j)-En(i) from tn
  getijdelta_from_transitionnumber(i,j,delta,dj,tn,pr,est);

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
       dMQ=0;for(l=1;l<=3;++l){dMQ(l)=-2.0*Malpha(l);}

       // multiply by occupation number difference ... and print matrix elements
       popnr_diff(dMQ,i,j,est,delta,T,pr,"MQ");


// return number of all transitions     
return noft(est,T,pinit,ninit);
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
  int i,j=1,l;
  int dj=(int)(2*J+1);
  float delta;
   double ninit=drixs(1).real();
   double pinit=drixs(1).imag();
   
 // calculate mat and delta for transition number tn
  // 1. get i and j  delta=En(j)-En(i) from tn
  getijdelta_from_transitionnumber(i,j,delta,dj,tn,pr,est);

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


       // multiply by occupation number difference ... and print matrix elements
       popnr_diff(drixs,i,j,est,delta,T,pr,"Rij");


// return number of all transitions     
return noft(est,T,pinit,ninit);
}



// ----------------------------------------------------------------------
// for testing the code uncomment and make test and start ionpars.exe
/* int main(int argc, char **argv)
{FILE * cf_file;
cf_file = fopen_errchk (argv[1], "rb"); // reopen file

      ionpars iops(cf_file);
      
      fclose(cf_file);cf_file = fopen_errchk (argv[1], "rb"); // reopen file
      ionpars iops1(cf_file);
            
return 1;
}*/
