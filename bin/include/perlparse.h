#ifndef PERLPARSE
#define PERLPARSE

#include "sparsecomplex.hpp"

void PrintComplexMatrix(FILE * fout, char * name, zsMat<double> & mat);
int perlparse(char*sipffilename
              ,double ** numbers,char ** numbernames    // number 
              ,char **  strings,char ** stringnames    // string variables
              ,zsMat<double> **  operators,char ** operatornames    // operators
              );
int myparse(char*sipffilename
              ,double ** numbers,char ** numbernames    // number 
              ,char **  strings,char ** stringnames    // string variables
              ,zsMat<double> **  operators,char ** operatornames    // operators
              ,int **sq2=NULL, double **csn=NULL, char **statements=NULL, int *iocs=NULL, int* oes=NULL
              );
int myparse_execute(zsMat<double> **operators, 
              char **operatornames, int *seq2, double *constval, 
              char *statement, int ioc, int oe);

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
