#include "paramete.hpp"
#include <martin.h>

#define MAXNOFNUMBERSINLINE 20
#define MAXNOFBCPLANES 20
#define MAXNOFCHARINLINE 7024


 // *************************************************************************
 // ************************ parameters *************************************
 // *************************************************************************
 // methods for class parameters - the parameter file input will be 
 // changed into a more feasible format



// save parameters to file (jjj format)
int parameters::savjjj(FILE *file)
{ int i;
  errno = 0;
  fprintf(file,"%s",rems[1]);
  fprintf(file,"%s",rems[2]);
  fprintf(file,"a=%4.6g  b=%4.6g c=%4.6g [A] alpha=%4.6g beta=%4.6g gamma=%4.6g [deg]\n",a,b,c,alpha,beta,gamma);
  fprintf(file,"r1a=0.5 r2a=0 r3a=0\n");
  fprintf(file,"r1b=0.0 r2b=0.5 r3b=0   primitive basis vectors\n");
  fprintf(file,"r1c=0.5 r2c=0 r3c=1.0\n");
  fprintf(file,"nofatoms=1\n");

  fprintf(file,"****************************************************************************\n");
  fprintf(file,"atom=1 da=0 [a]  db=0 [b] dc=0 [c] nofneighbours=%i diagonalexchange=%i gJ=%4.6g sipffilename=mcphas.cf\n",paranz,diagonalexchange,gJ);
  fprintf(file,"da[a]     db[b]     dc[c]      Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");
  for (i = 1; i <= paranz; ++i)
    {
   fprintf(file,"%#+8.6f %#+8.6f %#+8.6f  %#+8.6f %#+8.6f %#+8.6f\n",
           dn[i](1),dn[i](2),dn[i](3),jij[i](1),jij[i](2),jij[i](3));
    }
  return 0;
}  

// save cf parameters to file (cf format)
int parameters::savcf(FILE *file)
{ fprintf(file,"%s",rems[1]);
  fprintf(file,"%s",rems[2]);
  fprintf(file,"%4.6g\n",J);
  fprintf(file,"%s",rems[4]);
  fprintf(file,"%4.6g %4.6g %4.6g %4.6g\n",tpara,dmag,tneel,hstar);
  fprintf(file,"%s",rems[5]);
  fprintf(file,"%4.6g %4.6g %4.6g %4.6g %4.6g %4.6g %4.6g %4.6g %4.6g %4.6g %4.6g\n",
                m0,alphaa,alphab,alphac,tferro,ex,m0xa,m0xb,m0xc,sta,staerr);
  fprintf(file,"%s",rems[6]);
  fprintf(file,"%4.6g %4.6g %4.6g\n",M,A,C);
  return 0;
}




// save parameters to file (coq format)
int parameters::save(FILE *file)
{ int i;
  errno = 0;
  fprintf(file,"%s",rems[1]);
  fprintf(file,"%s",rems[2]);
  fprintf(file,"%4.6g\n",J);
  fprintf(file,"%s",rems[4]);
  fprintf(file,"%4.6g %4.6g %4.6g %4.6g\n",tpara,dmag,tneel,hstar);
  fprintf(file,"%s",rems[5]);
  fprintf(file,"%4.6g %4.6g %4.6g %4.6g %4.6g %4.6g %4.6g %4.6g %4.6g %4.6g %4.6g\n",
                m0,alphaa,alphab,alphac,tferro,ex,m0xa,m0xb,m0xc,sta,staerr);
  fprintf(file,"%s",rems[6]);
  fprintf(file,"%4.6g %4.6g %4.6g\n",M,A,C);
  fprintf(file,"%s",rems[7]);
  fprintf(file,"%i\n",paranz);
  fprintf(file,"%s",rems[8]);
  for (i = 1; i <= paranz; ++i)
    {
     fprintf(file,"%4.6g %4.6g %4.6g\n",rij[i],jij[i](2),jijerr[i](2));
    }
  fprintf(file,"%s",rems[9]);
  for (i = 1; i <= paranz; ++i)
    {
     fprintf(file,"%4.6g %4.6g %4.6g\n",rij[i],jij[i](1),jijerr[i](1));
    }
  fprintf(file,"%s",rems[10]);
  for (i = 1; i <= paranz; ++i)
    {
     fprintf(file,"%4.6g %4.6g %4.6g\n",rij[i],jij[i](3),jijerr[i](3));
    }
  fprintf(file,"%s",rems[11]);
  for (i = 1; i <= paranz; ++i)
    {
     fprintf(file,"%4.6g %4.6g %4.6g %4.6g\n",rij[i],dn[i](1),dn[i](2),dn[i](3));
    }
   return 0;
}

//constructor
parameters::parameters (char *file)
{
  FILE *fin_coq;
  char instr[MAXNOFCHARINLINE];
  float nn[MAXNOFNUMBERSINLINE];
  nn[0]=MAXNOFNUMBERSINLINE;
  int i,ii,j,mult;
  errno = 0;
  
  fin_coq = fopen_errchk (file, "rb");
  fgets (instr, MAXNOFCHARINLINE, fin_coq);rems[1]=new char[strlen(instr)+1];strcpy(rems[1],instr);
  fgets (instr, MAXNOFCHARINLINE, fin_coq);rems[2]=new char[strlen(instr)+1];strcpy(rems[2],instr);
  inputline (fin_coq, nn);
  J = nn[1];
  fgets (instr, MAXNOFCHARINLINE, fin_coq);rems[4]=new char[strlen(instr)+1];strcpy(rems[4],instr);
  inputline (fin_coq, nn);
  tpara = nn[1];
  dmag  = nn[2];
  tneel = nn[3];
  hstar = nn[4];
  fgets (instr, MAXNOFCHARINLINE, fin_coq);rems[5]=new char[strlen(instr)+1];strcpy(rems[5],instr);
  inputline (fin_coq, nn);
  m0 = nn[1];
  alphaa = nn[2];
  alphab = nn[3];
  alphac = nn[4];
  tferro = nn[5];
  ex = nn[6];
  m0xa = nn[7];
  m0xb = nn[8];
  m0xc = nn[9];
  sta = nn[10];
  staerr = nn[11];
  fgets (instr, MAXNOFCHARINLINE, fin_coq);rems[6]=new char[strlen(instr)+1];strcpy(rems[6],instr);
  inputline (fin_coq, nn);
  M = nn[1];
  A = nn[2];
  C = nn[3];

// lande factor calculated from m0 and M
  gJ = m0 / M;

  fgets (instr, MAXNOFCHARINLINE, fin_coq);rems[7]=new char[strlen(instr)+1];strcpy(rems[7],instr);
  inputline (fin_coq, nn);
  paranz = (int)nn[1];

//dimension arrays
  rij = new float[paranz+1];
  if (rij == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    }
  jij = new Vector[paranz+1];for(i=0;i<=paranz;++i){jij[i]=Vector(1,3);}
  if (jij == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    }
  jijerr = new Vector[paranz+1];for(i=0;i<=paranz;++i){jijerr[i]=Vector(1,3);}
  if (jijerr == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    }
  dn = new Vector[paranz+1];for(i=0;i<=paranz;++i){dn[i]=Vector(1,3);}
  if (dn == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    }

  fgets (instr, MAXNOFCHARINLINE, fin_coq);rems[8]=new char[strlen(instr)+1];strcpy(rems[8],instr);
  for (i = 1; i <= paranz; ++i)
    {
      inputline (fin_coq, nn);
      rij[i] = nn[1];
      jij[i](2) = nn[2];
      jijerr[i](2) = nn[3];
    }
  fgets (instr, MAXNOFCHARINLINE, fin_coq);rems[9]=new char[strlen(instr)+1];strcpy(rems[9],instr);
  for (i = 1; i <= paranz; ++i)
    {
      inputline (fin_coq, nn);
      rij[i] = nn[1];
      jij[i](1) = nn[2];
      jijerr[i](1) = nn[3];
    }
  fgets (instr, MAXNOFCHARINLINE, fin_coq);rems[10]=new char[strlen(instr)+1];strcpy(rems[10],instr);
  for (i = 1; i <= paranz; ++i)
    {
      inputline (fin_coq, nn);
      rij[i] = nn[1];
      jij[i](3) = nn[2];
      jijerr[i](3) = nn[3];
    }
  fgets (instr, MAXNOFCHARINLINE, fin_coq);rems[11]=new char[strlen(instr)+1];strcpy(rems[11],instr);
  for (i = 1; i <= paranz; ++i)
    {
      inputline (fin_coq, nn);
      rij[i] = nn[1];
      dn[i](1) = nn[2];
      dn[i](2) = nn[3];
      dn[i](3) = nn[4];
      if (nn[2]==1.0&&nn[3]==0.0&&nn[4]==0.0) a=nn[1];
      if (nn[2]==0.0&&nn[3]==1.0&&nn[4]==0.0) b=nn[1];
      if (nn[2]==0.0&&nn[3]==0.0&&nn[4]==1.0) c=nn[1];
     }
      alpha=90;beta=90;gamma=90;diagonalexchange=1;
  fclose (fin_coq);
// here we have to calculate jj[] and nofbcplanes 
// by summing the couplinfs over all bc planes
	jj = new Vector[MAXNOFBCPLANES+1];for(i=0;i<=MAXNOFBCPLANES;++i){jj[i]=Vector(1,3);jj[i]=0;}
        if (jj == NULL)
	{ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}

  nofbcplanes=0;
  for (i = 1; i <= paranz; ++i)
      {if (nofbcplanes<(j=(int)rint(fabs(2*dn[i](1))))){nofbcplanes=j;}
       mult=1;// determination of multiplicity of parameter i
       if (fabs(dn[i](2))>0.00001){mult*=2;}
       if (fabs(dn[i](3))>0.00001&&
           fabs(fabs(2*dn[i](3))-rint(fabs(2*dn[i](3))))<0.00001)
	   {mult*=2;}
	// sum up parameter 
       for (ii=1;ii<=3;++ii)
       {jj[j](ii)+=mult*jij[i](ii);}
      } 
}

//kopier-konstruktor
parameters::parameters (const parameters & p)
{ int i;
  a=p.a;b=p.b;c=p.c;alpha=p.alpha;beta=p.beta;gamma=p.gamma;
  diagonalexchange=p.diagonalexchange;
  J = p.J;
  tpara = p.tpara;
  dmag = p.dmag;
  hstar = p.hstar;
  m0 = p.m0;
  alphaa = p.alphaa;
  alphab = p.alphab;
  alphac = p.alphac;
  tferro = p.tferro;
  ex = p.ex;
  m0xa = p.m0xa;
  m0xb = p.m0xb;
  m0xc = p.m0xc;
  sta = p.sta;
  staerr = p.staerr;
  M = p.M;
  A = p.A;
  C = p.C;
  gJ=p.gJ;
  paranz = p.paranz;

//dimension arrays
  for (i=1;i<=11;++i)
  {rems[i] = new char[strlen(p.rems[i])+1];
   strcpy(rems[i],p.rems[i]);}

  rij = new float[paranz+1];
  if (rij == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    }
  jij = new Vector[paranz+1];for(i=0;i<=paranz;++i){jij[i]=Vector(1,3);}
  if (jij == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    }
  jijerr = new Vector[paranz+1];for(i=0;i<=paranz;++i){jijerr[i]=Vector(1,3);}
  if (jijerr == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    }
  dn = new Vector[paranz+1];for(i=0;i<=paranz;++i){dn[i]=Vector(1,3);}
  if (rij == NULL)
    {
      fprintf (stderr, "Out of memory\n");
      exit (EXIT_FAILURE);
    }
  for (i=1;i<=paranz;++i)
  {rij[i]=p.rij[i];jij[i]=p.jij[i];jijerr[i]=p.jijerr[i];dn[i]=p.dn[i];}

  nofbcplanes=p.nofbcplanes;
	jj = new Vector[MAXNOFBCPLANES+1];for(i=0;i<=MAXNOFBCPLANES;++i){jj[i]=Vector(1,3);}
        if (jj == NULL)
	{ fprintf (stderr, "Out of memory\n");exit (EXIT_FAILURE);}
  for (i=1;i<=MAXNOFBCPLANES;++i)
  {jj[i]=p.jj[i];}

}


//destruktor
parameters::~parameters ()
{ int i;
  for(i=1;i<=11;++i)
  {delete []rems[i];}
  delete []rij;
  delete []jij;
  delete []jijerr;
  delete []dn;
  delete []jj;
}


