#include "graphic_parameters.hpp"

graphic_parameters::~graphic_parameters(){}


graphic_parameters::graphic_parameters()
{
snprintf(title,MAXNOFCHARINLINE,"chargedensity");
show_pointcharges=0;scale_pointcharges=1;
showprim=0;
show_abc_unitcell=1.0;
show_primitive_crystal_unitcell=1.0;
show_magnetic_unitcell=1.0;
show_atoms=1.0;
scale_view_1=1.0;
scale_view_2=1.0;
scale_view_3=1.0;
show_density=1.0;
scale_density_vectors=0;
spins_scale_moment=1.0;
spins_show_ellipses=1;
spins_colour=1;

spins_show_static_moment_direction=1;
spins_wave_amplitude=0;
spins_show_oscillation=1;
phonon_scale_static_displacements=1;
phonon_wave_amplitude=0;
threshhold=0.05;
density_dtheta=0.2;
density_dfi=0.2;
read();
gridi=30;
gridj=30;
gridk=30;
}

int graphic_parameters::read()
{FILE * fin;
fin=fopen("results/graphic_parameters.set","rb");
if (fin==NULL) return 0;

//std::fstream FILEIN; FILEIN.open("results/graphic_parameters.set", std::fstream::in);
//std::string strline;
//if(FILEIN.fail()==true) return 0;

char instr[MAXNOFCHARINLINE];
while(feof(fin)==false)
{
//while(!FILEIN.eof())
//{getline(FILEIN,strline);
if(fgets(instr,sizeof(instr),fin)!=NULL){
extract(instr,"show_pointcharges",show_pointcharges);
extract(instr,"scale_pointcharges",scale_pointcharges);
extract(instr,"showprim",showprim);
extract(instr,"show_abc_unitcell",show_abc_unitcell);
extract(instr,"show_primitive_crystal_unitcell",show_primitive_crystal_unitcell);
extract(instr,"show_magnetic_unitcell",show_magnetic_unitcell);
extract(instr,"show_atoms",show_atoms);
extract(instr,"scale_view_1",scale_view_1);
extract(instr,"scale_view_2",scale_view_2);
extract(instr,"scale_view_3",scale_view_3);
extract(instr,"gridi",gridi);
extract(instr,"gridj",gridj);
extract(instr,"gridk",gridk);
extract(instr,"show_density",show_density);
extract(instr,"scale_density_vectors",scale_density_vectors);
extract(instr,"density_dtheta",density_dtheta);
extract(instr,"density_dfi",density_dfi);
extract(instr,"spins_scale_moment",spins_scale_moment);
extract(instr,"spins_colour",spins_colour);
extract(instr,"spins_show_ellipses",spins_show_ellipses);

extract(instr,"spins_show_static_moment_direction",spins_show_static_moment_direction);
extract(instr,"spins_wave_amplitude",spins_wave_amplitude);
extract(instr,"spins_show_oscillation",spins_show_oscillation);
extract(instr,"phonon_wave_amplitude",phonon_wave_amplitude);
extract(instr,"phonon_scale_static_displacements",phonon_scale_static_displacements);
extract(instr,"density_threshhold",threshhold);
}
else
{*instr='\0';}
}

fclose(fin);
//FILEIN.close();
return 1;
}
