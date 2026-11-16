// functions for program spins option -q  - Eigenvector mode

void print_eigenvector_human_legible(par & p,spincf & ev_r,spincf & ev_i,char * outhklstr)
{printf("%s\n#sipf-file at (da[a] db[b] dc[c]) eigenvector\n",outhklstr);
 
  for(int ii=1;ii<=p.cs.nofatoms;++ii)
   for(int i=1;i<=ev_r.na();++i)
   for(int j=1;j<=ev_r.nb();++j)
   for(int k=1;k<=ev_r.nc();++k) //loop supercells
 { Vector dd3(1,3),dd4(1,3);
    dd3=ev_r.pos_dabc(i,j,k,ii, p.cs);
    dd4=ev_r.pos_dr123(i,j,k,ii, p.cs);
    printf("# %s at dabc=(%5.3f %5.3f %5.3f) dr1r2r3=(%5.3f %5.3f %5.3f) [ ",(*p.jjj[ii]).sipffilename,dd3(1),dd3(2),dd3(3),dd4(1),dd4(2),dd4(3));
    for(int a=1;a<=ev_r.nofcomponents;++a) {if( ev_r.m(i,j,k)((ii-1)*ev_r.nofcomponents+a)!=0.0||
                                                ev_i.m(i,j,k)((ii-1)*ev_r.nofcomponents+a)==0.0)
                                            printf("%5.4f",ev_r.m(i,j,k)((ii-1)*ev_r.nofcomponents+a)); 
                                           if( ev_i.m(i,j,k)((ii-1)*ev_r.nofcomponents+a)!=0.0)
                                            printf("%+5.4fi",ev_i.m(i,j,k)((ii-1)*ev_i.nofcomponents+a));
                                           printf(" ");
                                         }     
    printf("]\n");  
 }  
}


// apply symmetry operation sn from ciffile  and return character
// if verbose=true: print eigenvectors
complex<double> symop(par & p,spincf * scr,spincf * sci,int deg,char * outhklstr,int sn,std::string ciffile, bool polar,bool verbose)
{FILE * fin; fin = fopen_errchk (ciffile.c_str(),"rb");
 char instr[MAXNOFCHARINLINE];bool found=false,begin=false;
 while(feof(fin)==false&&found==false)
 {fgets(instr, MAXNOFCHARINLINE, fin);
  if(instr[strspn(instr," \t")]!='#')
  {if(begin==true)
    {if(sn==atoi(instr))found=true;
    }
   if(strncmp(instr,"_space_group_symop_operation_xyz",32)==0)begin=true;
  }
 }
 fclose(fin);
if(found==false){fprintf(stderr,"Error program spins: symmetry operation %i not found\n",sn);exit(EXIT_FAILURE); }

if(verbose)printf("#Applying symmetry operation %s\n",instr);
  
parser ob;
char * xd=mystrnstr(instr,"'",MAXNOFCHARINLINE);++xd;
char * yd=mystrnstr(xd,",",MAXNOFCHARINLINE);*yd='\0';++yd;
char * zd=mystrnstr(yd,",",MAXNOFCHARINLINE);*zd='\0';++zd;
char * d=mystrnstr(zd,"'",MAXNOFCHARINLINE);*d='\0';
//printf("%s\n%s\n%s\n",xd,yd,zd);
char f[MAXNOFCHARINLINE];
Vector t(1,3);
snprintf(f,MAXNOFCHARINLINE,"x=%g y=%g z=%g xd=%s yd=%s zd=%s",0.0,0.0,0.0,xd,yd,zd);
parseline(f,ob);extract(f,"xd",t(1),ob);extract(f,"yd",t(2),ob);extract(f,"zd",t(3),ob);
if(verbose)myPrintVector(t,"corresponding to x'=A*x+t and ev'=B*ev\n# translation t:");
Matrix A(1,3,1,3);
snprintf(f,MAXNOFCHARINLINE,"x=%g y=%g z=%g xd=%s yd=%s zd=%s",1.0,0.0,0.0,xd,yd,zd);
parseline(f,ob);extract(f,"xd",A(1,1),ob);extract(f,"yd",A(2,1),ob);extract(f,"zd",A(3,1),ob);
                A(1,1)-=t(1);A(2,1)-=t(2);A(3,1)-=t(3);
snprintf(f,MAXNOFCHARINLINE,"x=%g y=%g z=%g xd=%s yd=%s zd=%s",0.0,1.0,0.0,xd,yd,zd);
parseline(f,ob);extract(f,"xd",A(1,2),ob);extract(f,"yd",A(2,2),ob);extract(f,"zd",A(3,2),ob);
                A(1,2)-=t(1);A(2,2)-=t(2);A(3,2)-=t(3);
snprintf(f,MAXNOFCHARINLINE,"x=%g y=%g z=%g xd=%s yd=%s zd=%s",0.0,0.0,1.0,xd,yd,zd);
parseline(f,ob);extract(f,"xd",A(1,3),ob);extract(f,"yd",A(2,3),ob);extract(f,"zd",A(3,3),ob);
                A(1,3)-=t(1);A(2,3)-=t(2);A(3,3)-=t(3);
if(verbose)myPrintMatrix(A,"Roto(inversion) Matrix A:");
// determine transformation matrix for vectors (and pseudovectors) from
// spatial transformation matrix A: B should transform eigenvectors like a spatial vector
// in ijk coordinates.
// dabcd=A*dabc  
// dijk=abc_in_ijk*dabc 
// dijkd=abc_in_ijk*dabcd 
//  --> dijkd=abc_in_ijk*A*abc_in_ijk.Inverse()*dijk
// -->
Matrix B(1,3,1,3);Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,p.cs.abc);
B=abc_in_ijk*A*abc_in_ijk.Inverse(); 
if(polar==false&&Det(B)<0){B=-B;} // for axial vectors ignore the inversion symmetry operation
if(verbose)myPrintMatrix(B,"Observable transformation Matrix B:");

// go through all degenerate configurations, transform them and store original and transformed
// vectors in matrices evm and evm'
int dim=scr[1].nofcomponents*scr[1].nofatoms*scr[1].n();
ComplexMatrix evm(1,dim,1,deg),evmd(1,dim,1,deg);
for(int dg=1;dg<=deg;++dg)
{if(scr[dg].nofcomponents!=3){fprintf(stderr,"Error spins: nofcomponents=%i symmetry operations can currently only be performed for vector observables with nofcomponents=3\n",scr[dg].nofcomponents);exit(EXIT_FAILURE);}
spincf evd_r=scr[dg].symop(A,t,B,p.cs);
spincf evd_i=sci[dg].symop(A,t,B,p.cs);
if(verbose){printf("# Original eigenvector:\n");print_eigenvector_human_legible(p,scr[dg],sci[dg],outhklstr);}
if(verbose){printf("# Transformed eigenvector:\n");print_eigenvector_human_legible(p,evd_r,evd_i,outhklstr);}
// store as vector in order to determine character of representation ...
Vector vr=scr[dg].asVector();
Vector vi=sci[dg].asVector();
ComplexVector v(vr,vi);for(int i=v.Lo();i<=v.Hi();++i)evm(i,dg)=v(i);

Vector vdr=evd_r.asVector();
Vector vdi=evd_i.asVector();
ComplexVector vd(vdr,vdi);for(int i=vd.Lo();i<=vd.Hi();++i)evmd(i,dg)=vd(i);
}
// linear regression to determine deg x deg transformation matrix T
// evmd=evm*T
// evm.Transpose * evmd= evm.Transpose * evm * T = C * T
// T=C.inverse * evm.Transpose * evmd
ComplexMatrix C(1,deg,1,deg); C=evm.Hermitean()*evm;
ComplexMatrix T(1,deg,1,deg); T=C.Inverse()*(evm.Hermitean()*evmd);

complex<double> ch=Trace(T);
// detect how good is linear regression in order to see if something is odd
// evmd=evm*T should be tue. 
// evm*T-evmd=0 should be true - check it
double e=NormFro(evm*T-evmd)/dim/deg;
if(e>0.0001) {myPrintMatrix(evmd,"Matrix of transformed eigenvectors:");evmd=evm*T-evmd;
             myPrintMatrix(evmd,"Difference Matrix to Matrix of best fitting linear combination of original eigenvectors:");
             fprintf(stderr,"Error: spin configuration not invariant to symmetry operation %s |s-s'|=%g\n",instr,e);
             exit(EXIT_FAILURE); }

if(verbose)printf("ch%i=%5.3f%+5.3fi ",sn,ch.real(),ch.imag());
return ch;
}


bool identify_irrep(complex<double> * ch,int nofsymop,std::string ciffile)
// look if cif file contains point group character table and if yes find out irrep
//_point_group_name D4h
//_point_group_classes                   'E  EI C42 C42I 2C4 2C4I 2C2 2C2I 2C2d 2C2dI '
//_point_group_space_symmetry_operations '16 8  11  3    13  5    1   9    6    14'
//loop_
//_irrep_point_group
//_character_point_group
// A1g ' 1  1  1  1  1  1  1  1  1  1'
// A1u ' 1 -1  1 -1  1 -1  1 -1  1 -1'
{FILE * fin; fin = fopen_errchk (ciffile.c_str(),"rb");
 char instr[MAXNOFCHARINLINE];
 char pointgroup[MAXNOFCHARINLINE];
 char classes[MAXNOFCHARINLINE];
 float pso[nofsymop+1];pso[0]=nofsymop;
 bool found=false,begin=false,eot=false;int nofc=0,nd=0;char * irrep;
 while(feof(fin)==false&&found==false&&eot==false)
 {fgets(instr, MAXNOFCHARINLINE, fin);
  if(instr[strspn(instr," \t")]!='#')
  {char * t=instr+strspn(instr," \t");
   if(strncmp(t,"_point_group_name",17)==0)snprintf(pointgroup,MAXNOFCHARINLINE,"%s",t+18);
   if(strncmp(t,"_point_group_classes",20)==0)snprintf(classes,MAXNOFCHARINLINE,"%s",t+20);
   if(strncmp(t,"_point_group_space_symmetry_operations",38)==0)
       {t+=38;t+=strspn(t," \t'");
        if((nofc=splitstring(t,pso))>nofsymop)
        {std::cerr << "Error spins reading " << ciffile << ": _point_group_space_symmetry_operations > number of symmetry operations read " << nofsymop << "\n";
         exit(EXIT_FAILURE); }
       }
   if(begin==true) // we are in the table - compare irrep characters with irreps
    {if(instr[strspn(instr," \t")]=='_'||strncmp(instr+strspn(instr," \t"),"loop_",5)==0)eot=true;
     else{
     irrep=1+mystrnstr(instr,"'",MAXNOFCHARINLINE);
     float irrepch[nofsymop+1];irrepch[0]=nofsymop;
     int nofd=splitstring(irrep,irrepch);
     if(nofc!=nofd){fprintf(stderr,"Error - irrep contains %i characters not equal to %i in tag _point_group_space_symmetry_operations \n",nofc,nofd);exit(EXIT_FAILURE);}
     int i=0;
     while(i<nofc&&fabs(irrepch[i+1]-ch[(int)pso[i+1]].real())<1e-3){++i;//printf("%8.6f=%8.6f ",irrepch[i],ch[(int)pso[i]].real());
      }
     //printf("\n");
     if(i==nofc)found=true;
     //else printf("i=%i table: %8.6f  calc: %8.6f\n",(int)pso[i+1],irrepch[i+1],ch[(int)pso[i+1]].real());
     }}
   if(strncmp(instr,"_character_point_group",22)==0)begin=true;
   }
 ++nd;
 }
 fclose(fin);
 if(found==true){printf("irrep=%s",instr);return true;}
 else {std::cerr << "Error: irrep not found in " << ciffile << '\n';
       return false;}
}
                   

/**********************************************************************/
// programm fuer option -q (koennte auch eigenstaendig compiliert werden)
/**********************************************************************/
int spins_q(int argc, char **argv, LoadedInput & s)
{FILE * fin;
 
 graphic_parameters gp;bool deg=true;
 char prefix[MAXNOFCHARINLINE];prefix[0]='\0';
 par inputpars("mcphas.j",prefix); 
 spincf spinconf(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 spincf spinconfev_real(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 spincf spinconfev_imag(1,1,1,inputpars.cs.nofatoms,inputpars.cs.nofcomponents);
 Vector Hext(1,HEXT_DIMENSION);Hext=0;
 char outstr[MAXNOFCHARINLINE];
 char outhklstr[MAXNOFCHARINLINE];
 Vector hkl(1,3);double E;double T=0;char * sn=NULL;
 bool polar=true;
if (s.eigenvector.symmetryMode == SymmetryMode::Symoa)polar=false;

spincf  scr[MAXDEG];
spincf  sci[MAXDEG];
double Etable[MAXDEG];
Etable[0]=1;
if(!s.eigenvector.hasEnergy) // no energy is given - probe all energies
{fin = fopen_errchk (s.eigenvector.qFile.c_str(), "rb");
 char instr[MAXNOFCHARINLINE];instr[0]='#';
 while (feof(fin)==0)
 {fgets(instr,MAXNOFCHARINLINE,fin);
   if (instr[strspn(instr," \t")]!='#') // pointer to 'ltrimstring' 
    {if(0==extract(instr,"E",Etable[(int)Etable[0]]))Etable[0]++;
    }
 }Etable[0]--;
 fclose(fin);
}
printf("#* program spins ... reading from file %s\n",s.eigenvector.qFile.c_str()); 
for(int ie=1;ie<=Etable[0];++ie)
{E=Etable[ie];T=0;deg=true;
while(deg){fin = fopen_errchk (s.eigenvector.qFile.c_str(), "rb"); 
 // s.os has been correctly set for this in parseEigenvectorOptions()
 deg=check_for_best_excitation_and_close(fin,gp,"ddd",spinconf,spinconfev_real,spinconfev_imag,
                          0,0,T,Hext,outstr,outhklstr,"best",hkl,E,argc,argv,s.os,
                          inputpars.cs.abc,inputpars.cs.nofatoms,0,0.0,0.0,0.0,0); 
 T--;
 if(argc>s.os+8)print_eigenvector_human_legible(inputpars,spinconfev_real,spinconfev_imag,outhklstr);
 scr[-(int)T]=spincf(spinconfev_real);
 sci[-(int)T]=spincf(spinconfev_imag);
 }

 if(!s.eigenvector.cif.empty()) // if a cif file was given in the command line
                  {if(s.eigenvector.symmetryIndex<0) // label of symoperation in table)
                   {int nofsymop=-s.eigenvector.symmetryIndex;// label of symoperation in table;
                    complex<double> ch[nofsymop+1];
                    for(int i=1;i<=nofsymop;++i)
                    {
                     ch[i]=symop(inputpars,scr,sci,(int)-T,outhklstr,i,s.eigenvector.cif,polar,false);
                    }
                    printf("E=%g ",E);
                   // look if cif file contains point group character table and if yes find out irrep
                   if(identify_irrep(ch,nofsymop,s.eigenvector.cif)==false)
                    {for(int i=1;i<=nofsymop;++i)printf("ch%i=%5.3f%+5.3fi ",i,ch[i].real(),ch[i].imag());
                     printf("\n");
                    }
                   }
                   else symop(inputpars,scr,sci,(int)-T,outhklstr,s.eigenvector.symmetryIndex,s.eigenvector.cif,polar,true);
                  }
} 
 return 0; 
}
 
