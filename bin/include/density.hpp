//  class density for the storage of a charge/spin/orbmom...density surface
// used in program mcphase

#ifndef CHARGEDENSITY
#define CHARGEDENSITY

#include<martin.h>
#include<jjjpar.hpp>

class density
{ private:
  int nofpts;  
  Vector **rtetafi;
  int type;
   Vector pr;

Vector momentsx;
Vector momentsy;
Vector momentsz;
Vector momentlx;
Vector momently;
Vector momentlz;

  public:
  double dtheta;
  double dfi;
  Vector & rtf(int i); // pointer to rft (i=1 ... nofpoints)
  int nofpoints(); // returns nofpoints in densitysurface  
  
  // sub for calculation of charge density given a radiu R and polar angles teta, 
  // fi and expansion coeff. alm
  double denscalc(double & theta, double & fi,double & R,Vector & moments,jjjpar & ionpar,double & T, Vector &  gjmbH);
  void calc_cd_surface(Vector & moments,jjjpar & ionpar,double ccc,double & T, Vector &  gjmbH);
  //ccc ... surface value of density

  void moments_init(Vector & moments);
  //constructor
  density(char*title,double dtheta,double dfi);
  density(const density & p);//kopier-konstruktor
 ~density ();//destruktor
};

#endif

