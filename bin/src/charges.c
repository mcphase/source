/***************************************************
 * charges - display charges at given htpoint
 * Author: M. Rotter
 **************************************************/

#include "../../version"
#include "spincf.hpp"
#include "martin.h"
#include "myev.h"
#include <par.hpp>
#include <graphic_parameters.hpp>
#include <cryststruct.hpp>
#include "densities_func.c"


/**********************************************************************/
// main program
int main (int argc, char **argv)
{
printf("#***************************************************\n");
printf("# * charges - display charges at given htpoint\n");
printf("# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
printf("# * %s\n",MCPHASVERSION);
printf("# **************************************************\n");
// check command line
  if (argc < 6)
    { printf (" program charges - display charges at HT point\n\
                use as: charges threshhold T Ha Hb Hc  [file.mf]\n \
                            (default input file is results/mcphas.mf)\n \
                    or: charges threshhold T Ha Hb Hc h k l E \n \
                        (reads from results/mcphas.mf and results/mcdisp.qev\n\n \
                This program outputs a magnetic/charge structure (and magnetic/orbital excitation)\n \
                graphic/movie in the output files results/charges*.jvx (javaview)\n\n \
                the graphics output format can be fine tuned in .mf and .qev input files\n\n \
                jvx files can be viewed by: java javaview results/charges.jvx \n \
                                            java javaview \"model=results/charges.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n\n");
      exit (1);
    }
FILE * fin_coq, * fout;
 graphic_parameters gp;
  gp.threshhold=strtod(argv[1],NULL);

  cryststruct cs;
   spincf savmf;
   int i,n;
   char outstr[MAXNOFCHARINLINE];

// read input file with mf configuration
 if (argc==7)
 { fin_coq = fopen_errchk (argv[6], "rb");}
 else
 { fin_coq = fopen_errchk ("./results/mcphas.mf", "rb");}
  fout = fopen_errchk ("./results/charges.out", "w");
   // input file header and mfconf------------------------------------------------------------------
   n=headerinput(fin_coq,fout,gp,cs);
 
   int ext_nof_components[MAXNOFATOMS];
   for (i=1;i<=n;++i){ext_nof_components[i]=48;if (cs.gJ[i]==0){ext_nof_components[i]=51;}
                                       // here set for 3+48 components, module ic1ion
                                      }
   // check for spinfconfiguration which is nearest to the T/H values chosen by user in command line
   double T,ha,hb,hc;
   check_for_best(fin_coq,strtod(argv[2],NULL),strtod(argv[3],NULL),strtod(argv[4],NULL),strtod(argv[5],NULL),savmf,T,ha,hb,hc,outstr);
  fclose (fin_coq);
  
// create plot of spin+chargeconfiguration -----------------------------------------------------------
  int ii,nt,k,j;
  double lnz,u;
  

  par inputpars("./mcphas.j");
  
  Vector hh(1,savmf.nofcomponents*savmf.nofatoms);
  int max_ext_nof_components=51;

  spincf extendedspincf(savmf.na(),savmf.nb(),savmf.nc(),savmf.nofatoms,max_ext_nof_components);

          // the following is for the printout of charges.out ...........................
           fprintf(fout,"#!T=%g K Ha=%g T Hb= %g T Hc= %g T: nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",T,ha,hb,hc,savmf.na(),savmf.nb(),savmf.nc(),inputpars.nofatoms*savmf.na()*savmf.nb()*savmf.nc());
            fprintf(fout,"#J=value {atom-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ma> <Mb> <Mc> [mb] <Ja> <Jb> <Jc> ...\n");
          fprintf(fout,"#{corresponding effective fields heff [meV]- if passed to mcdiff only these are used for caculation (not the magnetic moments)}\n");

  // determine primitive magnetic unit cell
     Matrix p(1,3,1,3);Vector xyz(1,3),dd0(1,3);
     savmf.calc_prim_mag_unitcell(p,cs.abc,cs.r);
  // .............................................................................                                
	       
//  1. from the meanfieldconfiguration (savmf) the <Olm> have to be calculated for all l=2,4,6
// 1.a: the mcphas.j has to be used to determine the structure + single ione properties (copy something from singleion.c)
// 1.b: mcalc has to be used to calculate all the <Olm>.
hh=0;for(ii=1;ii<=inputpars.nofatoms;++ii)
{//(*inputpars.jjj[ii]).eigenstates(hh,T);} // initialize eigenstate matrices
 (*inputpars.jjj[ii]).mcalc_parameter_storage_init(hh,T);} // initialize mcalc module parameter storage

 for (i=1;i<=savmf.na();++i){for(j=1;j<=savmf.nb();++j){for(k=1;k<=savmf.nc();++k)
 {
    hh=savmf.m(i,j,k);
  extendedspincf.m(i,j,k)=0;
  for(ii=1;ii<=inputpars.nofatoms;++ii)
 {  Vector h(1,ext_nof_components[ii]);
    Vector moments(1,ext_nof_components[ii]);
    h=0;
   for(nt=1;nt<=savmf.nofcomponents;++nt){h(nt)=hh(nt+savmf.nofcomponents*(ii-1));}
            if((*inputpars.jjj[ii]).module_type!=1&&(*inputpars.jjj[ii]).module_type!=3)
            {(*inputpars.jjj[ii]).mcalc(moments,T,h,lnz,u,(*inputpars.jjj[ii]).mcalc_parstorage); // here we trigger single ion
                                                           // module to calculate all 48 (ext_nof_components)
                                                           // higher order moments 
            }
            else
            {moments=0;}

          // output atoms and moments in primitive unit cell to stdout
              Vector dd3(1,3);
              dd3=savmf.pos(i,j,k,ii, cs.abc, cs.r,cs.x,cs.y,cs.z);
              dd0=p.Inverse()*dd3;dd0(1)*=savmf.na();dd0(2)*=savmf.nb();dd0(3)*=savmf.nc();
              fprintf(fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
	              cs.cffilenames[ii],dd3(1)/cs.abc(1),dd3(2)/cs.abc(2),dd3(3)/cs.abc(3),dd0(1),dd0(2),dd0(3));
              printf("{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f \n",
	              cs.cffilenames[ii],dd3(1)/cs.abc(1),dd3(2)/cs.abc(2),dd3(3)/cs.abc(3),dd0(1),dd0(2),dd0(3));
                     for(nt=1;nt<=3;++nt){if(cs.gJ[ii]!=0){fprintf(fout," %4.4f",cs.gJ[ii]*moments(nt));}
                                          else         {fprintf(fout," %4.4f",2*moments(2*nt-1)+moments(2*nt));}
                                         }
                     for(nt=1;nt<=ext_nof_components[ii];++nt)                                                               // this else is when gJ=0: means intermediate coupling
		        {extendedspincf.m(i,j,k)(nt+max_ext_nof_components*(ii-1))=moments(nt);
                         fprintf(fout," %4.4f",extendedspincf.m(i,j,k)(nt+max_ext_nof_components*(ii-1)));
                        }
                         fprintf(fout,"\n");
                      fprintf(fout,"                  corresponding effective fields heff [meV]-->          ");
                      for(nt=1;nt<=savmf.nofcomponents;++nt)  // printout meanfields
                        {fprintf(fout," %4.4f",h(nt));}
                         fprintf(fout,"\n");

  }}}}

fclose(fout);
  gp.read();// read graphic parameters which are set by user in file results/graphic_parameters.set

//print out the long vector of moments 1-48
  printf("%s - spin configuration <Olm>(i)\n",outstr);
  extendedspincf.print(stdout);

             Vector hkl(1,3);hkl=0;
             spincf savev_real(extendedspincf*0.0);
             spincf savev_imag(extendedspincf*0.0);
  fout = fopen_errchk ("./results/charges.grid", "w");
     extendedspincf.cd(fout,cs,gp,savev_real,savev_imag,0.0,hkl,T,hh);
    fclose (fout);

  fout = fopen_errchk ("./results/charges.jvx", "w");
    gp.showprim=0;gp.spins_wave_amplitude=0;
     extendedspincf.jvx_cd(fout,outstr,cs,gp,
                  0.0,savev_real,savev_imag,hkl,T,hh);
    fclose (fout);

  fout = fopen_errchk ("./results/charges_prim.jvx", "w");
     gp.showprim=1;
    extendedspincf.jvx_cd(fout,outstr,cs,gp,
                  0.0,savev_real,savev_imag,hkl,T,hh);
    fclose (fout);


if (argc>=10){// try a spinwave picture
             double E;
             long int pos=0;
             int extended_eigenvector_dimension;
              char instr[MAXNOFCHARINLINE];
              float numbers[13];numbers[9]=1;numbers[10]=3;
              numbers[0]=13;
             gp.spins_wave_amplitude=1.0;gp.spins_show_ellipses=1.0;gp.spins_show_oscillation=1.0;
             fin_coq = fopen_errchk ("./results/mcdisp.qee", "rb");
             // input file header ------------------------------------------------------------------
             instr[0]='#';
              while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
              { pos=ftell(fin_coq); 
                if (pos==-1) 
                {fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
                fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
                for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
               // load evs and check which one is nearest -------------------------------   
               extract(instr,"spins_wave_amplitude",gp.spins_wave_amplitude);
               extract(instr,"spins_show_ellipses",gp.spins_show_ellipses);
               extract(instr,"spins_show_oscillation",gp.spins_show_oscillation);
               extract(instr,"extended_eigenvector_dimension",extended_eigenvector_dimension);
              }
//               if(extended_eigenvector_dimension!=48){fprintf(stderr,"Error program charges - extended_eigenvector_dimension in results/mcdisp.qee not equal to 48\n");exit(1);}
               j=fseek(fin_coq,pos,SEEK_SET); 
               if (j!=0){fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
   
               double delta,dd,ddT,ddHa,ddHb,ddHc,ddh,ddk,ddl,ddE;
               for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

               { fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                 spincf ev_real(extendedspincf.na(),extendedspincf.nb(),extendedspincf.nc(),extendedspincf.nofatoms,extended_eigenvector_dimension);
                 spincf ev_imag(extendedspincf.na(),extendedspincf.nb(),extendedspincf.nc(),extendedspincf.nofatoms,extended_eigenvector_dimension);
                 ev_real.load(fin_coq);ev_imag.load(fin_coq);
                 ddT=strtod(argv[2],NULL)-numbers[4];ddT*=ddT;
                 ddHa=strtod(argv[3],NULL)-numbers[1];ddHa*=ddHa;
                 ddHb=strtod(argv[4],NULL)-numbers[2];ddHb*=ddHb;
                 ddHc=strtod(argv[5],NULL)-numbers[3];ddHc*=ddHc;
                 ddh=strtod(argv[6],NULL)-numbers[5];ddh*=ddh;
                 ddk=strtod(argv[7],NULL)-numbers[6];ddk*=ddk;
                 ddl=strtod(argv[8],NULL)-numbers[7];ddl*=ddl;
                 ddE=strtod(argv[9],NULL)-numbers[9];ddE*=ddE;
                 
                 dd=sqrt(ddT+ddHa+ddHb+ddHc+ddh+ddk+ddl+ddE+0.000001);
                 if (dd<delta)
                 {delta=dd;
                  sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g h=%g k=%g l=%g E=%g",numbers[4],numbers[1],numbers[2],numbers[3],numbers[5],numbers[6],numbers[7],numbers[9]);
                  hkl(1)=numbers[5];hkl(2)=numbers[6];hkl(3)=numbers[7];E=numbers[9];                                        
                  savev_real=ev_real;
                  savev_imag=ev_imag;                  
                 }
               }
              fclose (fin_coq);
              fprintf(stdout,"#%s - eigenvector\n",outstr);
              fprintf(stdout,"#real\n");
              savev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              savev_imag.print(stdout);

              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              double phase;
              complex <double> im(0,1);
              for(i=0;i<16;++i)
              {phase=2*3.1415*i/15;
               printf("\n********************************************\n");
               printf(" calculating movie sequence %i(16)\n",i+1);
               printf("********************************************\n");
               char filename[MAXNOFCHARINLINE];
               sprintf(filename,"./results/charges.%i.jvx",i+1);
               fin_coq = fopen_errchk (filename, "w");gp.showprim=0;
                     extendedspincf.jvx_cd(fin_coq,outstr,cs,gp,
                                  phase,savev_real,savev_imag,hkl,T,hh);
               fclose (fin_coq);
               sprintf(filename,"./results/charges_prim.%i.jvx",i+1);
               fin_coq = fopen_errchk (filename, "w");gp.showprim=1;
                     extendedspincf.jvx_cd(fin_coq,outstr,cs,gp,
                                  phase,savev_real,savev_imag,hkl,T,hh);
               fclose (fin_coq);
              }
          printf("# %s\n",outstr);
          }

fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program charges\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/charges.jvx\n");
fprintf(stderr,"# * java javaview \"model=results/charges.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n");
fprintf(stderr,"# * saved density mesh in results/charges.grid\n");
fprintf(stderr,"# ************************************************************************\n");

  for(i=1;i<=cs.nofatoms;++i){  delete cs.cffilenames[i];}

  return 0;



}


