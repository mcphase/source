/* mms.cpp
 *
 * Loads and saves a matrix from a file in MatrixMarket format to a matrix of the sMat<> class
 *
 * Functions:
 *    void mm_gout(sMat<double> M, const char*filename, const char*comments="");// Outputs a matrix to MatrixMarket format
 *    sMat<double> mm_gin(const char *filename);                                // Reads in a MatrixMarket sparse matrix
 *    void mm_sout(sMat<double> M, const char*filename, const char*comments="");// Output to file a sparse symmetrix matrix
 *    sMat<double> mm_sin(const char *filename);                                // Reads in a sparse symmetrix matrix
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "maths.hpp"
#include <fstream>
#include <iomanip>
#include <cstring>

// --------------------------------------------------------------------------------------------------------------- //
// Function to output a matrix in Matrix Market coordinate real general format
// --------------------------------------------------------------------------------------------------------------- //
void mm_gout(sMat<double> M, const char *filename, const char *comments="")
{
   int i,r,c;
   std::vector< std::vector<int> > nz = M.find();
   int sz = (int)nz.size();

   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out | std::fstream::app);
   FILEOUT << "%%MatrixMarket matrix coordinate real general\n%" << comments << "\n";
   FILEOUT << M.nr() << "\t" << M.nc() << "\t" << sz << "\n";

   FILEOUT.precision(24);
   if(sz != 0)
   {
      for (i=0; i<sz; i++)
      {
         r = nz[i][0]; c = nz[i][1];
         FILEOUT << (r+1) << "\t" << (c+1) << "\t" << M(r,c) << "\n";
      }
   }
   
   FILEOUT.close();
}

// --------------------------------------------------------------------------------------------------------------- //
// Function to read in a matrix from a file in the Matrix Market coordinate real general format
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> mm_gin(const char *filename)
{
   char comments[1024];
   int i,r,c,sz,ppos;
   double elem;
   bool commentsflag = true;

   std::fstream FILEIN; FILEIN.open(filename, std::fstream::in);
   if(FILEIN.fail()==true) { sMat<double> retval; return retval; }

   FILEIN.getline(comments,1024);            // Gets the first line %%MatrixMarket matrix coordinate real general
   while(commentsflag)
   {
      ppos = FILEIN.tellg();
      FILEIN.getline(comments,1024);         // Gets the comments line
      if(comments[0]!=35 && comments[0]!=37) // 35==# 37==%
         commentsflag=false;
   }
   FILEIN.seekg(ppos);
   FILEIN >> r >> c >> sz;                   // Gets the number of rows, columns and non-zero elements

   FILEIN.precision(24);
   sMat<double> retval(r,c);

   for (i=0; i<sz; i++)
   {
      FILEIN >> r >> c >> elem; retval(r-1,c-1) = elem;
   }

   FILEIN.close();
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Function to output a matrix in Matrix Market coordinate real symmetric format
// --------------------------------------------------------------------------------------------------------------- //
void mm_sout(sMat<double> M, const char *filename, const char *comments="")
{
   if(!M.issymm()) return;

   int i,r,c;
   std::vector< std::vector<int> > nz = M.findupper();
   int sz = (int)nz.size();

   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out | std::fstream::app);
   FILEOUT << "%%MatrixMarket matrix coordinate real symmetric\n%" << comments << "\n";
   FILEOUT << M.nr() << "\t" << M.nc() << "\t" << sz << "\n";

   if(sz != 0)
   {
      for (i=0; i<sz; i++)
      {
         r = nz[i][0]; c = nz[i][1];
         FILEOUT << (r+1) << "\t" << (c+1) << "\t" << std::setprecision(16) << M(r,c) << "\n";
      }
   }
   
   FILEOUT.close();
}

// --------------------------------------------------------------------------------------------------------------- //
// Function to read in a matrix from a file in the Matrix Market coordinate real symmetric format
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> mm_sin(const char *filename)
{
   char comments[1024];
   int i,r,c,sz,ppos;
   double elem;
   bool commentsflag = true;

   std::fstream FILEIN; FILEIN.open(filename, std::fstream::in);
   if(FILEIN.fail()==true) { sMat<double> retval; return retval; }

   FILEIN.getline(comments,1024);            // Gets the first line %%MatrixMarket matrix coordinate real symmetric
   while(commentsflag)
   {
      ppos = FILEIN.tellg();
      FILEIN.getline(comments,1024);         // Gets the comments line
      if(comments[0]!=35 && comments[0]!=37) // 35==# 37==%
         commentsflag=false;
   }
   FILEIN.seekg(ppos);
   FILEIN >> r >> c >> sz;                   // Gets the number of rows, columns and non-zero elements

   sMat<double> retval(r,c);

   for (i=0; i<sz; i++)
   {
      FILEIN >> r >> c >> elem; retval(r-1,c-1) = elem;
      if(r!=c) retval(c-1,r-1) = elem;
   }

   FILEIN.close();
   return retval;
}
