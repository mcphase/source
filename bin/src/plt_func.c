/********************************************************************/
// functions for programs chrgplt, spindensplt, orbmomdensplt
/********************************************************************/


int read_pointcharge_parameters(graphic_parameters & gp,char ** sipffilenames,char * cffile,float * x, float *y, float*z,jjjpar & jjjps,Vector & abc)
{FILE * cf_file;int pchere;int nofpc=0;
 
cf_file = fopen_errchk (cffile, "rb");
 float par[100];par[0]=99;
while((pchere=inputparline("pointcharge",cf_file,par))==0&&feof(cf_file)==false){;}
while(pchere>3)
{if(jjjps.module_type==0||jjjps.module_type==4){printf("pointcharge %g |e| at xyz=%g %g %g mind: xyz=abc\n",par[1],par[2],par[3],par[4]);}
 if(jjjps.module_type==2){printf("pointcharge %g |e| at xyz=%g %g %g mind: xyz=cab\n",par[1],par[2],par[3],par[4]);}
 ++nofpc;if(nofpc>MAXNOFATOMS){fprintf(stderr,"Error chrgplt - too many pointcharges");exit(1);}
  sipffilenames[1+nofpc]=new char[MAXNOFCHARINLINE];
  snprintf(sipffilenames[1+nofpc],MAXNOFCHARINLINE,"pointcharge radius=%g",gp.scale_pointcharges*0.529177*signum(par[1])*pow((double)fabs(par[1]),0.3333));
  if(jjjps.module_type==0||jjjps.module_type==4){
  x[nofpc+1]=par[2]/abc(1)+x[1];// these are the positions in Angstroem
  y[nofpc+1]=par[3]/abc(2)+y[1];// however in order to be in line with the cfield xyz=cab
  z[nofpc+1]=par[4]/abc(3)+z[1];// and ic1ion xyz=abc we have to set these parameters
                          }
  if(jjjps.module_type==2){
  x[nofpc+1]=par[3]/abc(1)+x[1];// these are the positions in Angstroem (we set a=b=c=1A below)
  y[nofpc+1]=par[4]/abc(2)+y[1];// however in order to be in line with the cfield xyz=cab
  z[nofpc+1]=par[2]/abc(3)+z[1];// and ic1ion xyz=abc we have to set these parameters
                          }
 while((pchere=inputparline("pointcharge",cf_file,par))==0&&feof(cf_file)==false){}
 }
 fclose(cf_file);
return nofpc;
}
