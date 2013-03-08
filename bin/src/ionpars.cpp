// ionpars: class to load and store matrices for internal module cfield and so1ion

#include "ionpars.hpp"
#include "martin.h"
#include "ionpars.h"
#include "myev.h"

#define NOF_OLM_MATRICES 48

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
  Blm=p.Blm; // vector of crystal field parameters
  Llm=p.Llm; // vector of crystal field parameters
  // cnst is the Zlm constants - put them into the matrix ... (same code is reused in jjjpar.cpp, pointc.c)
cnst=Matrix(0,6,-6,6);int l,m;
cnst(0,0) = 0.28209479;
cnst(2,0) = 0.3153962;
cnst(2,1)=  1.092548;
cnst(2,2)=  0.5462823;
cnst(4,0)=  0.1057871;
cnst(4,1)=  0.6690465;
cnst(4,2)=  0.4730943;
cnst(4,3)=  1.77013;
cnst(4,4)=  0.625845;
cnst(6,0)=  0.06357014;
cnst(6,1)=  0.582621;
cnst(6,2)=  0.4606094;
cnst(6,3)=  0.921205;
cnst(6,4)=  0.5045723;
cnst(6,5)=  2.366619;
cnst(6,6)=  0.6831942;
for(l=2;l<=6;l+=2){for(m=0;m<=l;++m)cnst(l,-m)=cnst(l,m);}

  
   int i;
   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 
   iontype = new char [strlen(p.iontype)+1];
   strcpy(iontype,p.iontype);
  

 for(i=1;i<=NOF_OLM_MATRICES;++i)
 { Olm [i]= new Matrix(1,(*p.Olm[i]).Rhi(),1,(*p.Olm[i]).Chi()); 
 // define first matrix 
   OOlm [i] = new ComplexMatrix(1,(*p.OOlm[i]).Rhi(),1,(*p.OOlm[i]).Chi()); 
   (*Olm[i])=(*p.Olm[i]);
   (*OOlm[i])=(*p.OOlm[i]);
 }  

}
ionpars::ionpars (int dimj) // constructor from dimj
 {J=((double)dimj-1)/2;
  Ja=Matrix(1,dimj,1,dimj);
  Jb=Matrix(1,dimj,1,dimj);
  Jc=Matrix(1,dimj,1,dimj);
  Hcf=Matrix(1,dimj,1,dimj);
  Jaa=ComplexMatrix(1,dimj,1,dimj);
  Jbb=ComplexMatrix(1,dimj,1,dimj);
  Jcc=ComplexMatrix(1,dimj,1,dimj);
   so1ion=0;
   Blm=Vector(0,48);Blm=0; // vector of crystal field parameters
   Llm=Vector(0,45);Llm=0; // vector of crystal field parameters

   alpha=0;beta=0;gamma=0;r2=0;r4=0;r6=0;nof_electrons=0;
// cnst is the Zlm constants - put them into the matrix ... (same code is reused in jjjpar.cpp, pointc.c)
cnst=Matrix(0,6,-6,6);int l,m;
cnst(0,0) = 0.28209479;
cnst(2,0) = 0.3153962;
cnst(2,1)=  1.092548;
cnst(2,2)=  0.5462823;
cnst(4,0)=  0.1057871;
cnst(4,1)=  0.6690465;
cnst(4,2)=  0.4730943;
cnst(4,3)=  1.77013;
cnst(4,4)=  0.625845;
cnst(6,0)=  0.06357014;
cnst(6,1)=  0.582621;
cnst(6,2)=  0.4606094;
cnst(6,3)=  0.921205;
cnst(6,4)=  0.5045723;
cnst(6,5)=  2.366619;
cnst(6,6)=  0.6831942;
for(l=2;l<=6;l+=2){for(m=0;m<=l;++m)cnst(l,-m)=cnst(l,m);}
   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 
   iontype = new char [MAXNOFCHARINLINE];


 int i;   
 for(i=1;i<=NOF_OLM_MATRICES;++i)
 { Olm [i]= new Matrix(1,dimj,1,dimj); 
 // define first matrix 
   OOlm [i] = new ComplexMatrix(1,dimj,1,dimj); 
 }  

}

 
ionpars::ionpars (char * ion) // constructor from iontype (mind:no matrices filled with values !)
 {int dimj;
  getpar(ion, &dimj, &alpha, &beta, &gamma, &gJ,&r2, &r4,&r6, &nof_electrons );
   iontype = new char [strlen(ion)+1];
   strcpy(iontype,ion);
   so1ion=0;

  J=((double)dimj-1)/2;
  Ja=Matrix(1,dimj,1,dimj);
  Jb=Matrix(1,dimj,1,dimj);
  Jc=Matrix(1,dimj,1,dimj);
  Hcf=Matrix(1,dimj,1,dimj);
  Jaa=ComplexMatrix(1,dimj,1,dimj);
  Jbb=ComplexMatrix(1,dimj,1,dimj);
  Jcc=ComplexMatrix(1,dimj,1,dimj);
// cnst is the Zlm constants - put them into the matrix ... (same code is reused in jjjpar.cpp, pointc.c)
cnst=Matrix(0,6,-6,6);int l,m;
cnst(0,0) = 0.28209479;
cnst(2,0) = 0.3153962;
cnst(2,1)=  1.092548;
cnst(2,2)=  0.5462823;
cnst(4,0)=  0.1057871;
cnst(4,1)=  0.6690465;
cnst(4,2)=  0.4730943;
cnst(4,3)=  1.77013;
cnst(4,4)=  0.625845;
cnst(6,0)=  0.06357014;
cnst(6,1)=  0.582621;
cnst(6,2)=  0.4606094;
cnst(6,3)=  0.921205;
cnst(6,4)=  0.5045723;
cnst(6,5)=  2.366619;
cnst(6,6)=  0.6831942;
for(l=2;l<=6;l+=2){for(m=0;m<=l;++m)cnst(l,-m)=cnst(l,m);}

   Blm=Vector(0,48);Blm=0; // vector of crystal field parameters
   Llm=Vector(0,45);Llm=0; // vector of crystal field parameters

   Olm = new Matrix * [1+NOF_OLM_MATRICES];  // define array of pointers to our Olm matrices
   OOlm= new ComplexMatrix * [1+NOF_OLM_MATRICES]; 


 int i;   
 for(i=1;i<=NOF_OLM_MATRICES;++i)
 { Olm [i]= new Matrix(1,dimj,1,dimj); 
 // define first matrix 
   OOlm [i] = new ComplexMatrix(1,dimj,1,dimj); 
 }  

}

ionpars::~ionpars(){
 int i;
 delete []iontype;

 for (i=1;i<=NOF_OLM_MATRICES;++i)
  {delete Olm[i];delete OOlm[i];}
   delete[] Olm;
   delete[] OOlm;
  
 } //destructor


ionpars::ionpars(FILE * cf_file) 
//constructor with commands from file handle (filename of cf parameters etc)
{      
   static int pr=1;
//  FILE * tryfile;
  int dimj;complex<double> im(0,1);
  int i,j,l,m; //30 ... maximum number of 2j+1
  double alphar,betar,gammar,r2r,r4r,r6r,gJr;
  char instr[MAXNOFCHARINLINE];
  iontype= new char[MAXNOFCHARINLINE];
  char  moduletype[MAXNOFCHARINLINE];
   Blm=Vector(0,48);Blm=0; // vector of crystal field parameters
   Llm=Vector(0,45);Llm=0; // vector of crystal field parameters
   // cnst is the Zlm constants - put them into the matrix ... (same code is reused in jjjpar.cpp, pointc.c)
cnst=Matrix(0,6,-6,6);
cnst(0,0) = 0.28209479;
cnst(2,0) = 0.3153962;
cnst(2,1)=  1.092548;
cnst(2,2)=  0.5462823;
cnst(4,0)=  0.1057871;
cnst(4,1)=  0.6690465;
cnst(4,2)=  0.4730943;
cnst(4,3)=  1.77013;
cnst(4,4)=  0.625845;
cnst(6,0)=  0.06357014;
cnst(6,1)=  0.582621;
cnst(6,2)=  0.4606094;
cnst(6,3)=  0.921205;
cnst(6,4)=  0.5045723;
cnst(6,5)=  2.366619;
cnst(6,6)=  0.6831942;
for(l=2;l<=6;l+=2){for(m=0;m<=l;++m)cnst(l,-m)=cnst(l,m);}
so1ion=0;strcpy(moduletype,"cfield");
   alpha=0;beta=0;gamma=0;r2=0;r4=0;r6=0;gJ=0;
  fgets_errchk (instr, MAXNOFCHARINLINE, cf_file);
  // strip /r (dos line feed) from line if necessary
  char *token;  
  while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}  
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
  
// read in lines and get IONTYPE=  and CF parameters Blm
   while(feof(cf_file)==false)
  {fgets(instr, MAXNOFCHARINLINE, cf_file);
   if(instr[strspn(instr," \t")]!='#'){//unless the line is commented ...
        extract(instr,"IONTYPE",iontype,(size_t)MAXNOFCHARINLINE);
        extract(instr,"nof_electrons",nof_electrons); //MR 120127

        extract(instr,"ALPHA",alphar);
        extract(instr,"BETA",betar);
        extract(instr,"GAMMA",gammar);

        extract(instr,"GJ",gJr);
       
        extract(instr,"R2",r2r);
        extract(instr,"R4",r4r);
        extract(instr,"R6",r6r);

        extract(instr,"B00",Blm(0));

        extract(instr,"B22S",Blm(1));
        extract(instr,"B21S",Blm(2));
	extract(instr,"B20",Blm(3));
        extract(instr,"B21",Blm(4));
	extract(instr,"B22",Blm(5));
   
	extract(instr,"B33S",Blm(6));
	extract(instr,"B32S",Blm(7));
	extract(instr,"B31S",Blm(8));
	extract(instr,"B30",Blm(9));
   extract(instr,"B31",Blm(10));
   extract(instr,"B32",Blm(11));
   extract(instr,"B32",Blm(12));

   extract(instr,"B44S",Blm(13));
   extract(instr,"B43S",Blm(14));
   extract(instr,"B42S",Blm(15));
   extract(instr,"B41S",Blm(16));
   extract(instr,"B40",Blm(17));
   extract(instr,"B41",Blm(18));
   extract(instr,"B42",Blm(19));
   extract(instr,"B43",Blm(20));
   extract(instr,"B44",Blm(21));
  
   extract(instr,"B55S",Blm(22));
   extract(instr,"B54S",Blm(23));
   extract(instr,"B53S",Blm(24));
   extract(instr,"B52S",Blm(25));
   extract(instr,"B51S",Blm(26));
   extract(instr,"B50",Blm(27));
   extract(instr,"B51",Blm(28));
   extract(instr,"B52",Blm(29));
   extract(instr,"B53",Blm(30));
   extract(instr,"B54",Blm(31));
   extract(instr,"B55",Blm(32));
 
   extract(instr,"B66S",Blm(33));
   extract(instr,"B65S",Blm(34));
   extract(instr,"B64S",Blm(35));
   extract(instr,"B63S",Blm(36));
   extract(instr,"B62S",Blm(37));
   extract(instr,"B61S",Blm(38));
   extract(instr,"B60",Blm(39));
   extract(instr,"B61",Blm(40));
   extract(instr,"B62",Blm(41));
   extract(instr,"B63",Blm(42));
   extract(instr,"B64",Blm(43));
   extract(instr,"B65",Blm(44));
   extract(instr,"B66",Blm(45));
   extract(instr,"Dx2",Blm(46));
   extract(instr,"Dy2",Blm(47));
   extract(instr,"Dz2",Blm(48));

	extract(instr,"L00",Llm(0));
 
        extract(instr,"L22S",Llm(1));
        extract(instr,"L21S",Llm(2));
	extract(instr,"L20",Llm(3));
        extract(instr,"L21",Llm(4));
	extract(instr,"L22",Llm(5));
   
	extract(instr,"L33S",Llm(6));
	extract(instr,"L32S",Llm(7));
	extract(instr,"L31S",Llm(8));
	extract(instr,"L30",Llm(9));
   extract(instr,"L31",Llm(10));
   extract(instr,"L32",Llm(11));
   extract(instr,"L32",Llm(12));

   extract(instr,"L44S",Llm(13));
   extract(instr,"L43S",Llm(14));
   extract(instr,"L42S",Llm(15));
   extract(instr,"L41S",Llm(16));
   extract(instr,"L40",Llm(17));
   extract(instr,"L41",Llm(18));
   extract(instr,"L42",Llm(19));
   extract(instr,"L43",Llm(20));
   extract(instr,"L44",Llm(21));
  
   extract(instr,"L55S",Llm(22));
   extract(instr,"L54S",Llm(23));
   extract(instr,"L53S",Llm(24));
   extract(instr,"L52S",Llm(25));
   extract(instr,"L51S",Llm(26));
   extract(instr,"L50",Llm(27));
   extract(instr,"L51",Llm(28));
   extract(instr,"L52",Llm(29));
   extract(instr,"L53",Llm(30));
   extract(instr,"L54",Llm(31));
   extract(instr,"L55",Llm(32));
 
   extract(instr,"L66S",Llm(33));
   extract(instr,"L65S",Llm(34));
   extract(instr,"L64S",Llm(35));
   extract(instr,"L63S",Llm(36));
   extract(instr,"L62S",Llm(37));
   extract(instr,"L61S",Llm(38));
   extract(instr,"L60",Llm(39));
   extract(instr,"L61",Llm(40));
   extract(instr,"L62",Llm(41));
   extract(instr,"L63",Llm(42));
   extract(instr,"L64",Llm(43));
   extract(instr,"L65",Llm(44));
   extract(instr,"L66",Llm(45));


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
                                  double BlIcalc=Llm(i)*cnst(l,m)*sqrt(4.0*PI/(2*l+1))*thetaJ(l);if(m!=0){BlIcalc*=sqrt(2.0);}
                                  if((Blm(i)!=0)&(fabs(Blm(i)-BlIcalc)/(fabs(BlIcalc)+1e-14)>0.001)){fprintf(stderr,"#Warning internal module %s - reading %s=%12.6g meV is ignored, because Wybourne Parameter Llm=%12.6g meV does not correspond ! \npresse enter to continue\n",moduletype,lm4,Blm(i),Llm(i));getchar();}
                                  Blm(i)=BlIcalc;// here set the Blm as calculated from the Llm
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

   if(Hcf==(double)0.0){

   for(l=1;l<=48;++l){Hcf+=Blm(l)*(*Olm[l]);
//                   if(Blm(l)!=0){if(l<24){fprintf(stderr,"B%c=%g   ",l+99,Blm(l));}
//		                     else{fprintf(stderr,"B(z+%i)=%g   ",l-23,Blm(l));}
//		                }
                  }
   }

if(so1ion==0)
 {//ATTENTION FOR cfield the AXES xyz are parallel to cab
 Matrix dummy(1,dimj,1,dimj);
 dummy=Jb;Jb=Jc;Jc=Ja;Ja=dummy;
 ComplexMatrix dummyc(1,dimj,1,dimj);
 dummyc=Jbb;Jbb=Jcc;Jcc=Jaa;Jaa=dummyc;

 if (pr==1) {printf("#Axis Convention using cfield as a module:  a||y b||z  c||x\n");
 printf("#xyz .... Coordinate system of the crystal field parameters used in cfield\n");
 printf("#abc .... Crystal axes\n");
 printf("#The interactions are described by the  PKQ Operators defined in cfield\n");
 printf("#O11(s) .... Ia=Jy\n");
 printf("#O10(c) .... Ib=Jz\n");
 printf("#O11(c) .... Ic=Jx\n");
 printf("#O22(s) .... Id\n");
 printf("#O21(s) .... Ie\n");
 printf("#O20(c) .... If\n");
 printf("#O21(c) .... Ig\n");
 printf("#O22(c) .... Ih\n");
 printf("#O33(s) .... Ii\n");
 printf("#O32(s) .... Ij\n");
 printf("#O31(s) .... Ik\n");
 printf("#O30(c) .... Il\n");
 printf("#O31(c) .... Im\n");
 printf("# etc ... 45 moments up to l<=6\n");
 printf("#\n");
             }
 }
 else
 {//ATTENTION FOR so1ion the AXES xyz are parallel to abc
 if (pr==1) {printf("#Axis Convention using so1ion as a module:  a||x b||y  c||z\n");
 printf("#xyz .... Coordinate system of the crystal field parameters used in so1ion\n");
 printf("#abc .... Crystal axes\n");
 printf("#The interactions are described by the  PKQ Operators defined in so1ion\n");
 printf("#O11(s) .... Ia=Jx\n");
 printf("#O10(c) .... Ib=Jy\n");
 printf("#O11(c) .... Ic=Jz\n");
 printf("#O22(s) .... Id\n");
 printf("#O21(s) .... Ie\n");
 printf("#O20(c) .... If\n");
 printf("#O21(c) .... Ig\n");
 printf("#O22(c) .... Ih\n");
 printf("#O33(s) .... Ii\n");
 printf("#O32(s) .... Ij\n");
 printf("#O31(s) .... Ik\n");
 printf("#O30(c) .... Il\n");
 printf("#O31(c) .... Im\n");
 printf("# etc ... 45 moments up to l<=6\n");
 printf("#\n");
             }
 }
pr=0;
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
   if(Blm(0)!=0){fprintf(outfile,"B00=%g\n",myround(Blm(0)));}

   if(Blm(1)!=0){fprintf(outfile,"B22S=%g\n",myround(Blm(1)));}
   if(Blm(2)!=0){fprintf(outfile,"B21S=%g\n",myround(Blm(2)));}
   if(Blm(3)!=0){fprintf(outfile,"B20=%g\n",myround(Blm(3)));}
   if(Blm(4)!=0){fprintf(outfile,"B21=%g\n",myround(Blm(4)));}
   if(Blm(5)!=0){fprintf(outfile,"B22=%g\n",myround(Blm(5)));}
   
   if(Blm(6)!=0){fprintf(outfile,"B33S=%g\n",myround(Blm(6)));}
   if(Blm(7)!=0){fprintf(outfile,"B32S=%g\n",myround(Blm(7)));}
   if(Blm(8)!=0){fprintf(outfile,"B31S=%g\n",myround(Blm(8)));}
   if(Blm(9)!=0){fprintf(outfile,"B30=%g\n",myround(Blm(9)));}
   if(Blm(10)!=0){fprintf(outfile,"B31=%g\n",myround(Blm(10)));}
   if(Blm(11)!=0){fprintf(outfile,"B32=%g\n",myround(Blm(11)));}
   if(Blm(12)!=0){fprintf(outfile,"B32=%g\n",myround(Blm(12)));}

   if(Blm(13)!=0){fprintf(outfile,"B44S=%g\n",myround(Blm(13)));}
   if(Blm(14)!=0){fprintf(outfile,"B43S=%g\n",myround(Blm(14)));}
   if(Blm(15)!=0){fprintf(outfile,"B42S=%g\n",myround(Blm(15)));}
   if(Blm(16)!=0){fprintf(outfile,"B41S=%g\n",myround(Blm(16)));}
   if(Blm(17)!=0){fprintf(outfile,"B40=%g\n",myround(Blm(17)));}
   if(Blm(18)!=0){fprintf(outfile,"B41=%g\n",myround(Blm(18)));}
   if(Blm(19)!=0){fprintf(outfile,"B42=%g\n",myround(Blm(19)));}
   if(Blm(20)!=0){fprintf(outfile,"B43=%g\n",myround(Blm(20)));}
   if(Blm(21)!=0){fprintf(outfile,"B44=%g\n",myround(Blm(21)));}
  
   if(Blm(22)!=0){fprintf(outfile,"B55S=%g\n",myround(Blm(22)));}
   if(Blm(23)!=0){fprintf(outfile,"B54S=%g\n",myround(Blm(23)));}
   if(Blm(24)!=0){fprintf(outfile,"B53S=%g\n",myround(Blm(24)));}
   if(Blm(25)!=0){fprintf(outfile,"B52S=%g\n",myround(Blm(25)));}
   if(Blm(26)!=0){fprintf(outfile,"B51S=%g\n",myround(Blm(26)));}
   if(Blm(27)!=0){fprintf(outfile,"B50=%g\n",myround(Blm(27)));}
   if(Blm(28)!=0){fprintf(outfile,"B51=%g\n",myround(Blm(28)));}
   if(Blm(29)!=0){fprintf(outfile,"B52=%g\n",myround(Blm(29)));}
   if(Blm(30)!=0){fprintf(outfile,"B53=%g\n",myround(Blm(30)));}
   if(Blm(31)!=0){fprintf(outfile,"B54=%g\n",myround(Blm(31)));}
   if(Blm(32)!=0){fprintf(outfile,"B55=%g\n",myround(Blm(32)));}

   if(Blm(33)!=0){fprintf(outfile,"B66S=%g\n",myround(Blm(33)));}
   if(Blm(34)!=0){fprintf(outfile,"B65S=%g\n",myround(Blm(34)));}
   if(Blm(35)!=0){fprintf(outfile,"B64S=%g\n",myround(Blm(35)));}
   if(Blm(36)!=0){fprintf(outfile,"B63S=%g\n",myround(Blm(36)));}
   if(Blm(37)!=0){fprintf(outfile,"B62S=%g\n",myround(Blm(37)));}
   if(Blm(38)!=0){fprintf(outfile,"B61S=%g\n",myround(Blm(38)));}
   if(Blm(39)!=0){fprintf(outfile,"B60=%g\n",myround(Blm(39)));}
   if(Blm(40)!=0){fprintf(outfile,"B61=%g\n",myround(Blm(40)));}
   if(Blm(41)!=0){fprintf(outfile,"B62=%g\n",myround(Blm(41)));}
   if(Blm(42)!=0){fprintf(outfile,"B63=%g\n",myround(Blm(42)));}
   if(Blm(43)!=0){fprintf(outfile,"B64=%g\n",myround(Blm(43)));}
   if(Blm(44)!=0){fprintf(outfile,"B65=%g\n",myround(Blm(44)));}
   if(Blm(45)!=0){fprintf(outfile,"B66=%g\n",myround(Blm(45)));}
   if(Blm(46)!=0){fprintf(outfile,"Dx2=%g\n",myround(Blm(46)));}
   if(Blm(47)!=0){fprintf(outfile,"Dy2=%g\n",myround(Blm(47)));}
   if(Blm(48)!=0){fprintf(outfile,"Dz2=%g\n",myround(Blm(48)));}

}

void ionpars::savLlm(FILE * outfile)
{fprintf(outfile,"units=meV\n");
   if(Llm(0)!=0){fprintf(outfile,"L00=%g\n",myround(Llm(0)));}

   if(Llm(1)!=0){fprintf(outfile,"L22S=%g\n",myround(Llm(1)));}
   if(Llm(2)!=0){fprintf(outfile,"L21S=%g\n",myround(Llm(2)));}
   if(Llm(3)!=0){fprintf(outfile,"L20=%g\n",myround(Llm(3)));}
   if(Llm(4)!=0){fprintf(outfile,"L21=%g\n",myround(Llm(4)));}
   if(Llm(5)!=0){fprintf(outfile,"L22=%g\n",myround(Llm(5)));}
  
   if(Llm(6)!=0){fprintf(outfile,"L33S=%g\n",myround(Llm(6)));}
   if(Llm(7)!=0){fprintf(outfile,"L32S=%g\n",myround(Llm(7)));}
   if(Llm(8)!=0){fprintf(outfile,"L31S=%g\n",myround(Llm(8)));}
   if(Llm(9)!=0){fprintf(outfile,"L30=%g\n",myround(Llm(9)));}
   if(Llm(10)!=0){fprintf(outfile,"L31=%g\n",myround(Llm(10)));}
   if(Llm(11)!=0){fprintf(outfile,"L32=%g\n",myround(Llm(11)));}
   if(Llm(12)!=0){fprintf(outfile,"L32=%g\n",myround(Llm(12)));}

   if(Llm(13)!=0){fprintf(outfile,"L44S=%g\n",myround(Llm(13)));}
   if(Llm(14)!=0){fprintf(outfile,"L43S=%g\n",myround(Llm(14)));}
   if(Llm(15)!=0){fprintf(outfile,"L42S=%g\n",myround(Llm(15)));}
   if(Llm(16)!=0){fprintf(outfile,"L41S=%g\n",myround(Llm(16)));}
   if(Llm(17)!=0){fprintf(outfile,"L40=%g\n",myround(Llm(17)));}
   if(Llm(18)!=0){fprintf(outfile,"L41=%g\n",myround(Llm(18)));}
   if(Llm(19)!=0){fprintf(outfile,"L42=%g\n",myround(Llm(19)));}
   if(Llm(20)!=0){fprintf(outfile,"L43=%g\n",myround(Llm(20)));}
   if(Llm(21)!=0){fprintf(outfile,"L44=%g\n",myround(Llm(21)));}
 
   if(Llm(22)!=0){fprintf(outfile,"L55S=%g\n",myround(Llm(22)));}
   if(Llm(23)!=0){fprintf(outfile,"L54S=%g\n",myround(Llm(23)));}
   if(Llm(24)!=0){fprintf(outfile,"L53S=%g\n",myround(Llm(24)));}
   if(Llm(25)!=0){fprintf(outfile,"L52S=%g\n",myround(Llm(25)));}
   if(Llm(26)!=0){fprintf(outfile,"L51S=%g\n",myround(Llm(26)));}
   if(Llm(27)!=0){fprintf(outfile,"L50=%g\n",myround(Llm(27)));}
   if(Llm(28)!=0){fprintf(outfile,"L51=%g\n",myround(Llm(28)));}
   if(Llm(29)!=0){fprintf(outfile,"L52=%g\n",myround(Llm(29)));}
   if(Llm(30)!=0){fprintf(outfile,"L53=%g\n",myround(Llm(30)));}
   if(Llm(31)!=0){fprintf(outfile,"L54=%g\n",myround(Llm(31)));}
   if(Llm(32)!=0){fprintf(outfile,"L55=%g\n",myround(Llm(32)));}
 
   if(Llm(33)!=0){fprintf(outfile,"L66S=%g\n",myround(Llm(33)));}
   if(Llm(34)!=0){fprintf(outfile,"L65S=%g\n",myround(Llm(34)));}
   if(Llm(35)!=0){fprintf(outfile,"L64S=%g\n",myround(Llm(35)));}
   if(Llm(36)!=0){fprintf(outfile,"L63S=%g\n",myround(Llm(36)));}
   if(Llm(37)!=0){fprintf(outfile,"L62S=%g\n",myround(Llm(37)));}
   if(Llm(38)!=0){fprintf(outfile,"L61S=%g\n",myround(Llm(38)));}
   if(Llm(39)!=0){fprintf(outfile,"L60=%g\n",myround(Llm(39)));}
   if(Llm(40)!=0){fprintf(outfile,"L61=%g\n",myround(Llm(40)));}
   if(Llm(41)!=0){fprintf(outfile,"L62=%g\n",myround(Llm(41)));}
   if(Llm(42)!=0){fprintf(outfile,"L63=%g\n",myround(Llm(42)));}
   if(Llm(43)!=0){fprintf(outfile,"L64=%g\n",myround(Llm(43)));}
   if(Llm(44)!=0){fprintf(outfile,"L65=%g\n",myround(Llm(44)));}
   if(Llm(45)!=0){fprintf(outfile,"L66=%g\n",myround(Llm(45)));}

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
   for(j=4;j<=JJ.Hi();++j){Ham-=gjmbH(j)*(*Olm[j-3]);}

/*   int i1,j1; //printout matrix
   for (i1=1;i1<=dj;++i1){
    for (j1=1;j1<=dj;++j1) printf ("%4.6g ",(*Olm[j])(i1,j1));
    printf ("\n");
    }*/
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
     z=ComplexMatrix(zr,zi);
//     z=ests(1,dj,1,dj)*z; // transform to original eigenstates ... however we deleted this because it needs more time to transform than to solve the eigenvalue problem
//     ests(1,dj,1,dj)=z;
//     for (i=1;i<=dj;++i) {ests(0,i)=complex <double> (En(i),wn(i));}
//     myPrintComplexMat(stdout,ests);     
//     myPrintComplexMat(stdout,z);

   za=Jaa*z;zb=Jbb*z;zc=Jcc*z;
    
     JJ=0;
//    ComplexVector ddd;
    for (i=1;i<=dj;++i)
    {
     JJ[1]+=wn(i)*real(z.Column(i)*za.Column(i));
     JJ[2]+=wn(i)*real(z.Column(i)*zb.Column(i));
     JJ[3]+=wn(i)*real(z.Column(i)*zc.Column(i));
    }
     
// here the expectation values of the multipolar moments are calculated
   for(j=4;j<=JJ.Hi();++j)
   {
     zolm=(*OOlm[j-3])*z;
    for (i=1;i<=dj;++i) JJ[j]+=wn(i)*real(z.Column(i)*zolm.Column(i));
   };
  

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
double dd; dd=NormFro((*OOlm[j-3])-(*OOlm[j-3]).Conjugate().Transpose());
   if (dd>1e-5) {printf("j=%i\n",j);myPrintComplexMatrix(stderr,(*OOlm[j-3]));}
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
double dd; dd=NormFro((*OOlm[j-3])-(*OOlm[j-3]).Conjugate().Transpose());
   if (dd>1e-5) {printf("j=%i\n",j);myPrintComplexMatrix(stderr,(*OOlm[j-3]));}
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
