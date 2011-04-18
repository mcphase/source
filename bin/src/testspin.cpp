 // *************************************************************************
 // ************************ testspincf *************************************
 // *************************************************************************
// methods for class testspincf 
#include "testspin.hpp"
#include "../../version"



// add spinconfiguration sps to list
int testspincf::addspincf(spincf & sps)
{int i;
 if (n<maxn)
  {++n;configurations[n]= new spincf(1,1,1,nofatoms,nofcomponents);
   (*configurations[n])=sps; // store configuration n 
   return n;
   }
  else
  {for (i=1;i<maxn;++i)
   {if (configurations[i]->wasstable==0)
     {(*configurations[i])=sps;
      return i;
      }
    }
   } 
return 0; //configuration could not be added
}

// print out table of all configurations
void testspincf::save (const char * filemode)
{save(savfilename,filemode);}


// print out table of all configurations
void testspincf::save (const char*filename, const char * filemode)
{if(n>0)
 { int i;
 FILE * fout;
// we should print to a file all used configurations
 fout = fopen_errchk (filename,filemode);
 fprintf(fout,"# Test spin configurations for meanfield calculation - module %s\n#<!--mcphase.mcphas.tst-->\n",MCPHASVERSION);
 fprintf(fout,"#*********************************************************\n");
 fprintf(fout,"# mcphas - program to calculate static magnetic properties\n");
 fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
 fprintf(fout,"#**********************************************************\n"); 
 fprintf(fout,"#! nofatoms=%i\n",nofatoms);
 fprintf(fout,"#! nofcomponents=%i\n",nofcomponents);
 fprintf(fout,"# File Format: blocks of 3xnofatoms lines consisting of\n");
 fprintf(fout,"# ma[1] ma[2] ma[3] .... ma[n]\n");
 fprintf(fout,"# mb[1] mb[2] mb[3] .... mb[n]\n");
 fprintf(fout,"# mc[1] mc[2] mc[3] .... mc[n]\n");

 for (i=1;i<=n;++i)
   {// print  configuration i 
    fprintf(fout,"#{momentum configuration %i}\n",i);
    configurations[i]->print(fout);
    }
  fclose(fout);
  printf("file %s saved\n",filename);
 }
}

//constructor - read maximum of nofconf spinconfigurations from file file and 
// rememeber filename for lateron save of the list (savfile)
testspincf::testspincf (int nofconf, const char * file,const char * savfile,int na,int nm)
{ savfilename= new char [strlen(savfile)+1];
  strcpy(savfilename,savfile);
  FILE *fin_coq;
  FILE *fout;
  long int pos=0,pos_old=0,j;
  char instr[MAXNOFCHARINLINE];
  char text[MAXNOFCHARINLINE];
  int i;
  maxn=nofconf;
  configurations = new  spincf *  [nofconf+1];
  errno = 0;
//  nofatoms=na;
//  nofcomponents=nm;

  fin_coq = fopen(file, "rb");
  if (fin_coq == NULL)
    {
      fprintf (stderr, "No file %s: %s - no test spinconfigurations read - no problem - continuing ...\n",file, strerror (errno));
      n=0;nofatoms=na;nofcomponents=nm;
     }
  else
    {
  printf("reading file %s\n",file);
   // input file header ------------------------------------------------------------------
  instr[0]='#';
   while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
  {pos_old=pos;pos=ftell(fin_coq);
   if (pos==-1) 
       {fprintf(stderr,"Error reading file %s\n",savfile);exit (EXIT_FAILURE);}
   fgets(instr,MAXNOFCHARINLINE,fin_coq);
   extract(instr,"nofatoms",nofatoms);
   extract(instr,"nofcomponents",nofcomponents);
  }
if (nofatoms!=na||nofcomponents!=nm)
   {fprintf(stderr,"Error mcphas.tst: nofatoms=%i or nofcomponents=%i does not agree with main module, it should be nofatoms=%i and nofcomponents=%i\n",nofatoms,nofcomponents,na,nm);exit (EXIT_FAILURE);}

 // load spinconfigurations and check which one is nearest -------------------------------   
   j=fseek(fin_coq,pos,SEEK_SET); 
    if (j!=0){fprintf(stderr,"Error: wrong file format of mcphas.tst\n");exit (EXIT_FAILURE);}

 for (i=1;i<=nofconf&&feof(fin_coq)==0                      //end of file
//                    &&fgets(instr,MAXNOFCHARINLINE,fin_coq)!=NULL   //error in line reading
		    ;++i)
    {
     configurations[i]=new spincf (1,1,1,nofatoms,nofcomponents); // declare configuration i (with j spins)
     if ((*configurations[i]).load(fin_coq)==0)
        {delete configurations[i];
	break;
	}
      printf("n1=%i n2=%i n3=%i\n",(*configurations[i]).na(),(*configurations[i]).nb(),(*configurations[i]).nc());
      configurations[i]->reduce();
      fout = fopen_errchk ("./results/.spins.eps", "w");
      sprintf(text,"Testspinconfiguration No %i",i);
      configurations[i]->eps(fout,text);
    fclose (fout);
      }
 fclose (fin_coq);
 n=i-1;
   }
}

//kopier-konstruktor
testspincf::testspincf (const testspincf & p)
{ savfilename= new char [strlen(p.savfilename)+1];
  strcpy(savfilename,p.savfilename);
  int i;
  maxn=p.maxn;
  n=p.n;
  nofatoms=p.nofatoms;
  nofcomponents=p.nofcomponents;
  configurations = new spincf * [maxn+1];


 for (i=1;i<n;++i)
   {// declare configuration i (with j spins)
    configurations[i]=new spincf(p.configurations[i]->n(),nofatoms,nofcomponents);  
    (*configurations[i])=(*p.configurations[i]);
    }

}




//destruktor
testspincf::~testspincf ()
{int i;
 for (i=1;i<n;++i)
   {//  delete configuration i 
    delete configurations[i];
    }
delete []configurations;
delete []savfilename;
}

