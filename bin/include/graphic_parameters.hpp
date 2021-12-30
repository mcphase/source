#ifndef GRAPHIC_PARAMETERS_H
#define GRAPHIC_PARAMETERS_H

#include<martin.h>
#include <fstream>

class graphic_parameters
{ public:
//POINTCHARGES
int show_pointcharges;
double scale_pointcharges;

//UNITCELL VIEW
double show_abc_unitcell,show_primitive_crystal_unitcell,show_magnetic_unitcell;
double show_atoms,scale_view_1,scale_view_2,scale_view_3;
int showprim;

//MOMENTS
double spins_scale_moment;
double spins_wave_amplitude;
double spins_show_ellipses;
double spins_show_static_moment_direction;
double spins_show_oscillation;
double spins_colour;

//PHONONS
double phonon_show_oscillation;
double phonon_wave_amplitude;


//DENSITY
double show_density,scale_density_vectors;
double density_dtheta,density_dfi;
double threshhold;
int gridi,gridj,gridk;

char title[100];

~graphic_parameters();


graphic_parameters();

int read();

};

#endif