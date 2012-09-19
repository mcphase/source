//  class testspincf for the storage of a set of test-spinconfigurations
// used in program mcphase

#ifndef TESTSPINS
#define TESTSPINS

#include<cstdio>
#include<cstring>
#include<cerrno>
#include<martin.h>
#include<spincf.hpp>


class testspincf
{ private:
  char * savfilename;
  
  public:
  int n,ninitial,maxn,nofatoms,nofcomponents; // number of configurations, maximum of number of configurations. nofatoms, nof moments
  spincf **configurations; // array of pointers to the different
                            //configurations

// add spinconfiguration sps to list
int  addspincf(spincf & sps);

// print out table of all configurations
void save(const char * filemode);
void save(const char*filename,const char * filemode);

//constructor - read maximum of nofconf spinconfigurations from file file and 
// remember filename for later on save of the list (savfile)
  testspincf (int nofconf, const char * file,const char * savfile, int nofatoms, int nofcomponents); //constructor
  testspincf (const testspincf & p);//kopier-konstruktor
 ~testspincf ();//destruktor
};

#endif

