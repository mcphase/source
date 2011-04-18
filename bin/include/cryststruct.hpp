// class to store a crystal structure


#ifndef CRYSTSTRUC
#define CRYSTSTRUC

#define MAXNOFATOMS 100

class cryststruct
{ public:
int nofatoms,nofcomponents,maxnofatoms;
   float x[MAXNOFATOMS],y[MAXNOFATOMS],z[MAXNOFATOMS];
   double gJ[MAXNOFATOMS];
   char * cffilenames[MAXNOFATOMS];
   Matrix r;
   Vector abc;

cryststruct()
{nofatoms=0;nofcomponents=3,maxnofatoms=MAXNOFATOMS;
 r=Matrix(1,3,1,3);
 abc=Vector(1,6);
 abc(4)=90;
 abc(5)=90;
 abc(6)=90;

}

};

#endif
