#ifndef PERLPARSE
#define PERLPARSE

void PrintComplexMatrix(FILE * fout, char * name, ComplexMatrix & mat);

int perlparse(char*sipffilename
              ,double ** numbers,char ** numbernames    // number 
              ,char **  strings,char ** stringnames    // string variables
              ,ComplexMatrix **  operators,char ** operatornames    // operators
              );

#endif
