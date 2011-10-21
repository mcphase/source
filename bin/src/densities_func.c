/*****************************************************************/
// functions for programs charges, spindensity, orbmomdensity momdensity
//                        currdensities
/*****************************************************************/

void check_for_best(FILE *fin_coq,double Tin, double hain,double hbin, double hcin, spincf & savmf, double & T,double & ha, double& hb,double&hc,char*outstr)
{// load mfconfigurations and check which one is nearest -------------------------------

int n;
   double ddT,ddHa,ddHb,ddHc,dd,delta;
 float numbers[13];numbers[9]=1;numbers[10]=3;
 numbers[0]=13;

 for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

    { spincf spins(1,1,1,(int)numbers[9],(int)numbers[10]);
      spins.load(fin_coq);
      ddT=Tin-numbers[3];ddT*=ddT;
      ddHa=hain-numbers[5];ddHa*=ddHa;
      ddHb=hbin-numbers[6];ddHb*=ddHb;
      ddHc=hcin-numbers[7];ddHc*=ddHc;
      dd=sqrt(ddT+ddHa+ddHb+ddHc+0.000001);
      if (dd<delta)
       {delta=dd;
        sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g n=%g spins nofatoms=%i in primitive basis nofcomponents=%i",myround(numbers[3]),myround(numbers[5]),myround(numbers[6]),myround(numbers[7]),myround(numbers[8]),(int)numbers[9],(int)numbers[10]);
        savmf=spins;T=numbers[3];ha=numbers[5];hb=numbers[6];hc=numbers[7];
       }
    }
}
   
int headerinput(FILE * fin_coq,FILE* fout,graphic_parameters & gp,cryststruct & cs)
{ char instr[MAXNOFCHARINLINE];
 long int pos=0,j;int n=0;
cs.nofatoms=0;cs.nofcomponents=3;
char *token;cs.abc=0;
  instr[0]='#';
 while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring'
  { pos=ftell(fin_coq);
   if (pos==-1)
       {fprintf(stderr,"Error: wrong mf/sps file format\n");exit (EXIT_FAILURE);}
   fgets_errchk(instr,MAXNOFCHARINLINE,fin_coq);
   // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
   int i;
   for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
   // strip /r (dos line feed) from line if necessary
    while ((token=strchr(instr,'\r'))!=NULL){*token=' ';}

   if (instr[strspn(instr," \t")]=='#'){fprintf(fout,"%s",instr);}
   if(cs.abc[1]==0){extract(instr,"a",cs.abc[1]);extract(instr,"b",cs.abc[2]); extract(instr,"c",cs.abc[3]);
                 extract(instr,"alpha",cs.abc[4]);  extract(instr,"beta",cs.abc[5]);extract(instr,"gamma",cs.abc[6]);
   }
   extract(instr,"show_abc_unitcell",gp.show_abc_unitcell);
   extract(instr,"show_primitive_crystal_unitcell",gp.show_primitive_crystal_unitcell);
   extract(instr,"show_magnetic_unitcell",gp.show_magnetic_unitcell);
   extract(instr,"show_atoms",gp.show_atoms);
   extract(instr,"spins_scale_moment",gp.spins_scale_moment);
   extract(instr,"show_chargedensity",gp.show_density);

   extract(instr,"scale_view_1",gp.scale_view_1);
   extract(instr,"scale_view_2",gp.scale_view_2);
   extract(instr,"scale_view_3",gp.scale_view_3);

   extract(instr,"r1x",cs.r[1][1]);extract(instr,"r2x",cs.r[1][2]); extract(instr,"r3x",cs.r[1][3]);
   extract(instr,"r1y",cs.r[2][1]); extract(instr,"r2y",cs.r[2][2]); extract(instr,"r3y",cs.r[2][3]);
   extract(instr,"r1z",cs.r[3][1]); extract(instr,"r2z",cs.r[3][2]); extract(instr,"r3z",cs.r[3][3]);
   extract(instr,"r1a",cs.r[1][1]);extract(instr,"r2a",cs.r[1][2]); extract(instr,"r3a",cs.r[1][3]);
   extract(instr,"r1b",cs.r[2][1]); extract(instr,"r2b",cs.r[2][2]); extract(instr,"r3b",cs.r[2][3]);
   extract(instr,"r1c",cs.r[3][1]); extract(instr,"r2c",cs.r[3][2]); extract(instr,"r3c",cs.r[3][3]);
   extract(instr,"nofatoms",cs.nofatoms);    extract(instr,"nofcomponents",cs.nofcomponents);
   if ((cs.nofatoms>0)&&((extract(instr,"x",cs.x[n+1])+
                    extract(instr,"y",cs.y[n+1])+
  		       extract(instr,"z",cs.z[n+1])==0)||
		       (extract(instr,"da",cs.x[n+1])+
                   extract(instr,"db",cs.y[n+1])+
		       extract(instr,"dc",cs.z[n+1])==0)))
		  {++n;if(n>cs.nofatoms||cs.nofatoms>cs.maxnofatoms)
                    {fprintf(stderr,"ERROR reading file:maximum number of atoms in unit cell exceeded\n");exit(EXIT_FAILURE);}
                   cs.cffilenames[n]=new char[MAXNOFCHARINLINE];
                   extract(instr,"cffilename",cs.cffilenames[n],(size_t)MAXNOFCHARINLINE);
                   extract(instr,"gJ",cs.gJ[n]);
//		   printf("%s\n",cs.cffilenames[n]);
                  }
  }
    j=fseek(fin_coq,pos,SEEK_SET);
    if (j!=0){fprintf(stderr,"Error: wrong mf file format\n");exit (EXIT_FAILURE);}
return n;
}
