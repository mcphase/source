#ifndef PERLPARSE
#define PERLPARSE

void PrintComplexMatrix(FILE * fout, char * name, ComplexMatrix & mat);

int perlparse(char*sipffilename
              ,double ** numbers,char ** numbernames    // number 
              ,char **  strings,char ** stringnames    // string variables
              ,ComplexMatrix **  operators,char ** operatornames    // operators
              );
int myparse(char*sipffilename
              ,double ** numbers,char ** numbernames    // number 
              ,char **  strings,char ** stringnames    // string variables
              ,ComplexMatrix **  operators,char ** operatornames    // operators
              ,int **sq2=NULL, double **csn=NULL, char **statements=NULL, int *iocs=NULL, int* oes=NULL
              );

int myparse_execute(ComplexMatrix **operators, 
              char **operatornames, int *seq2, double *constval, 
              char *statement, int ioc, int oe);
#endif
