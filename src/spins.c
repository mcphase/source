/****************************************************
 * spins - display spinconfiguration at given htpoint
 * Author: Martin Rotter
 ****************************************************/


#define MAXNOFATOMS 100
#include "../../version"
#include "spincf.hpp"
#include "martin.h"
#include "graphic_parameters.hpp"
#include "cryststruct.hpp"
#include "densities_func.c"

/**********************************************************************/
// hauptprogramm
int main (int argc, char **argv)
{ spincf savspins;
 FILE * fin_coq, * fout;
 float delta,dd,ddT,ddHa,ddHb,ddHc;
double T,Ha,Hb,Hc;
 int i,n=0;
 cryststruct cs;
 long int pos=0,j;
 float numbers[13];numbers[9]=1;numbers[10]=3;
 numbers[0]=13;
 char instr[MAXNOFCHARINLINE];
 char outstr[MAXNOFCHARINLINE];
 char filename[MAXNOFCHARINLINE];
 graphic_parameters gp;
gp.show_abc_unitcell=1.0;
gp.show_primitive_crystal_unitcell=1.0;
gp.show_magnetic_unitcell=1.0;
gp.show_atoms=1.0;
gp.scale_view_1=1.0;
gp.scale_view_2=1.0;
gp.scale_view_3=1.0;
gp.spins_scale_moment=1.0;

 // check command line
printf("#****************************************************\n");
printf("# spins - display spinconfiguration at given htpoint\n");
printf("# Author: Martin Rotter %s\n",MCPHASVERSION);
printf("#****************************************************\n");
  if (argc < 5)
    { printf (" program spins - display spincf at HT point\n \
                use as: spins T Ha Hb Hc [file.sps||file.mf]\n \
                        (default input file is results/mcphas.sps)\n \
                    or: spins T Ha Hb Hc h k l E\n \
                        (reads from results/mcphas.sps and results/mcdisp.qev\n\n \
                This program outputs a magnetic structure (and magnetic excitation)\n \
                graphic/movie in the output files of different format:\n \
                                     results/spins*.eps (postscript), results/spins*.fst (fp_studio), \n \
                                     results/spins.out (ascii) and results/spins*.jvx (javaview)\n\n \
                the graphics output format can be fine tuned in .sps and .qev input files\n \
                by show_abc_unitcell, show_primitive_crystal_unitcell, spins_scale_moment\n \
                show_magnetic_unitcell, show_atoms, scale_view_1,scale_view_2, scale_view_3 ...\n\n \
                jvx files can be viewed by: java javaview results/spins.jvx \n \
                                            java javaview \"model=results/spins.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n \
               \n");
      exit (1);
    }

 if (argc==6) 
 { fin_coq = fopen_errchk (argv[5], "rb");}
 else
 { fin_coq = fopen_errchk ("./results/mcphas.sps", "rb");}
    
 fout = fopen_errchk ("./results/spins.out", "w");

  // input file header and mfconf------------------------------------------------------------------
   n=headerinput(fin_coq,fout,gp,cs);
gp.read();
gp.show_density=0;
// load spinsconfigurations and check which one is nearest -------------------------------   

     check_for_best(fin_coq,strtod(argv[1],NULL),strtod(argv[2],NULL),strtod(argv[3],NULL),strtod(argv[4],NULL),savspins,T,Ha,Hb,Hc,outstr);

  fclose (fin_coq);

  printf("%s - momentum configuration <J(i)>\n",outstr);
  fprintf(fout,"#! %s - momentum configuration <J(i)>\n",outstr);
  savspins.printall(fout,cs);
  savspins.print(stdout);
  fclose (fout);
  
// create plot of spinconfiguration -----------------------------------------------------------

    fin_coq = fopen_errchk ("./results/spins.eps", "w");
     savspins.eps(fin_coq,outstr);
    fclose (fin_coq);

// here the 3d file should be created
    fin_coq = fopen_errchk ("./results/spinsab.eps", "w");
Vector gJJ(1,n); for (i=1;i<=n;++i){gJJ(i)=cs.gJ[i];}
   savspins.eps3d(fin_coq,outstr,cs.abc,cs.r,cs.x,cs.y,cs.z,1,gJJ);
    fclose (fin_coq);
    fin_coq = fopen_errchk ("./results/spinsac.eps", "w");
     savspins.eps3d(fin_coq,outstr,cs.abc,cs.r,cs.x,cs.y,cs.z,2,gJJ);
    fclose (fin_coq);
    fin_coq = fopen_errchk ("./results/spinsbc.eps", "w");
     savspins.eps3d(fin_coq,outstr,cs.abc,cs.r,cs.x,cs.y,cs.z,3,gJJ);
    fclose (fin_coq);
    fin_coq = fopen_errchk ("./results/spins3dab.eps", "w");
     savspins.eps3d(fin_coq,outstr,cs.abc,cs.r,cs.x,cs.y,cs.z,4,gJJ);
    fclose (fin_coq);
    fin_coq = fopen_errchk ("./results/spins3dac.eps", "w");
     savspins.eps3d(fin_coq,outstr,cs.abc,cs.r,cs.x,cs.y,cs.z,5,gJJ);
    fclose (fin_coq);
    fin_coq = fopen_errchk ("./results/spins3dbc.eps", "w");
     savspins.eps3d(fin_coq,outstr,cs.abc,cs.r,cs.x,cs.y,cs.z,6,gJJ);
    fclose (fin_coq);

    fin_coq = fopen_errchk ("./results/spins.fst", "w");
     savspins.fst(fin_coq,outstr,cs.abc,cs.r,cs.x,cs.y,cs.z,gJJ);
    fclose (fin_coq);

    
   fin_coq = fopen_errchk ("./results/spins_prim.fst", "w");
     savspins.fstprim(fin_coq,outstr,cs.abc,cs.r,cs.x,cs.y,cs.z,gJJ);
    fclose (fin_coq);

             Vector hkl(1,3);hkl=0;
             Vector gjmbH(1,3);gjmbH=0;
             spincf savev_real(savspins*0.0);
             spincf savev_imag(savspins*0.0);
            // to do jvx output of static structure put zeros into these spinconfigurations

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   fin_coq = fopen_errchk ("./results/spins.jvx", "w");
    gp.showprim=0;gp.spins_wave_amplitude=0;
     savspins.jvx_cd(fin_coq,outstr,cs,gp,0.0,savev_real,savev_imag,hkl,T,gjmbH);
    fclose (fin_coq);

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   fin_coq = fopen_errchk ("./results/spins_prim.jvx", "w");
     gp.showprim=1;
     savspins.jvx_cd(fin_coq,outstr,cs,gp,0.0,savev_real,savev_imag,hkl,T,gjmbH);
    fclose (fin_coq);

if (argc>=9){// try a spinwave picture
             double E,ddh,ddk,ddl,ddE;
             gp.spins_wave_amplitude=1;gp.spins_show_ellipses=1;gp.spins_show_static_moment_direction=1;
             fin_coq = fopen_errchk ("./results/mcdisp.qev", "rb");
             // input file header ------------------------------------------------------------------
             instr[0]='#';
              while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
              { pos=ftell(fin_coq); 
                if (pos==-1) 
                {fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
                fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
                for(i=0;i<=(int)strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
               // load evs and check which one is nearest -------------------------------   
               extract(instr,"spins_wave_amplitude",gp.spins_wave_amplitude);
               extract(instr,"spins_show_ellipses",gp.spins_show_ellipses);
               extract(instr,"spins_show_static_moment_direction",gp.spins_show_static_moment_direction);
              }
              gp.read();gp.show_density=0;
               j=fseek(fin_coq,pos,SEEK_SET); 
               if (j!=0){fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
   
               for (delta=1000.0;feof(fin_coq)==0                      //end of file
                    &&(n=inputline(fin_coq,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

               { fgets(instr,MAXNOFCHARINLINE,fin_coq); 
                 spincf ev_real(savspins),ev_imag(savspins);
                 ev_real.load(fin_coq);ev_imag.load(fin_coq);
                 ddT=strtod(argv[1],NULL)-numbers[4];ddT*=ddT;
                 ddHa=strtod(argv[2],NULL)-numbers[1];ddHa*=ddHa;
                 ddHb=strtod(argv[3],NULL)-numbers[2];ddHb*=ddHb;
                 ddHc=strtod(argv[4],NULL)-numbers[3];ddHc*=ddHc;
                 ddh=strtod(argv[5],NULL)-numbers[5];ddh*=ddh;
                 ddk=strtod(argv[6],NULL)-numbers[6];ddk*=ddk;
                 ddl=strtod(argv[7],NULL)-numbers[7];ddl*=ddl;
                 ddE=strtod(argv[8],NULL)-numbers[9];ddE*=ddE;
                 
                 dd=sqrt(ddT+ddHa+ddHb+ddHc+ddh+ddk+ddl+ddE+0.000001);
                 if (dd<delta)
                 {delta=dd;
                  sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g h=%g k=%g l=%g E=%g",numbers[4],numbers[1],numbers[2],numbers[3],numbers[5],numbers[6],numbers[7],numbers[9]);
                  savev_real=ev_real;savev_imag=ev_imag;hkl(1)=numbers[5];hkl(2)=numbers[6];hkl(3)=numbers[7];E=numbers[9];        
                 }
               }
              fclose (fin_coq);
              fprintf(stdout,"#%s - eigenvector\n",outstr);
              fprintf(stdout,"#real\n");
              savev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              savev_imag.print(stdout);
              spincf spins(savspins);

              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              double phase;
              complex <double> im(0,1);
              for(i=0;i<16;++i)
              {phase=2*3.1415*i/15;
               printf("\n calculating movie sequence %i(16)\n",i+1);
               sprintf(filename,"./results/spins.%i.jvx",i+1);
               fin_coq = fopen_errchk (filename, "w");
               //printf("swa %g\n",gp.spins_wave_amplitude);
               gp.showprim=0;
                     savspins.jvx_cd(fin_coq,outstr,cs,gp,phase,savev_real,savev_imag,hkl,T,gjmbH);
               fclose (fin_coq);
               sprintf(filename,"./results/spins_prim.%i.jvx",i+1);
               fin_coq = fopen_errchk (filename, "w");
               gp.showprim=1;
                     savspins.jvx_cd(fin_coq,outstr,cs,gp,phase,savev_real,savev_imag,hkl,T,gjmbH);
               fclose (fin_coq);
              }
          printf("# %s\n",outstr);  
          }


  for(i=1;i<=cs.nofatoms;++i){  delete cs.cffilenames[i];}
  return 0;
}


