/*-----------------------------------------------------------------------------*\
| function "IntMatrix MpMagicSquare (int order)"                 magicsquare.cc |
|                                                                               |
| Last change: Jul 14, 2000							|
|                                                                               |
| Matpack Library Release 1.5.1                                                 |
| Copyright (C) 1991-2000 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.						|
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.							|
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without prior consent of the author.			|
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration	|
| and installation of Matpack.							|
|                                                                               |
\*-----------------------------------------------------------------------------*/

#include "vector.h"

//-----------------------------------------------------------------------------//
// IntMatrix MpMagicSquare (int n)
//
// This function returns a magic square of order n in an integer valued
// matrix with row dimension 0 to n and column dimension 0 to n.
//
// The input argument must be within 3 and 1290, otherwise the 
// error handler Matpack.Error() is called, which terminates the program.
// The upper limit is determined by overflow condition on 32 bit machines.
// 
//-----------------------------------------------------------------------------//

static inline void NextPos (int &row, int &col, int current, int n)
{
  if ( (current % n) == 0 ) 
    row++;
  else {
    row--;
    col++;
  }
  if (row >= n) 
    row = 0;
  else if (row < 0) 
    row = n - 1;

  if (col >= n) 
    col = 0;
  else if (row < 0) 
    row = n - 1;
}

//-----------------------------------------------------------------------------//

static void MakeOddSquare (IntMatrix &M, int offset, int n, int srow, int scol)
{
  int count = offset, 
      row   = 0, 
      col   = (n - 1) / 2;

  for (int i = 0; i < n*n; ++i) {
    M(row+srow, col+scol) = count;
    NextPos(row, col, count, n);
    count++;
  }
}

//-----------------------------------------------------------------------------//

static void MakeSingleEvenSquare (IntMatrix &M)
{
  // Apply the De La Loubere algorithm in four evenly divided sectors of the matrix

  int n         = M.Rows(), 
      quadSize  = n / 2, 
      m         = (n - 2) / 4, 
      middleRow = (quadSize - 1) / 2, 
      targetRow = middleRow + quadSize;

  MakeOddSquare(M, 1, quadSize, 0, 0);
  MakeOddSquare(M, quadSize*quadSize + 1, quadSize, 0, quadSize);
  MakeOddSquare(M, quadSize*quadSize * 2 + 1, quadSize, quadSize, quadSize);
  MakeOddSquare(M, quadSize*quadSize * 3 + 1, quadSize, quadSize, 0);

  // swap middle row of quads A and D
  for (int i = 1; i <= m; i++) MpSwap( M(middleRow,i), M(targetRow,i) );

  // swap all other rows of quads A and D
  for (int j = 0; j < quadSize; j++)
    if (j != middleRow)
      for (int i = 0; i < m; i++) MpSwap( M(j,i), M(j+quadSize,i) );

  // swap all other rows of quads B and C
  for (int j = 0; j < quadSize; j++)
    for (int i = quadSize; i < quadSize + m + 2; i++) MpSwap( M(j,i), M(j+quadSize,i) );
}

//-----------------------------------------------------------------------------//

static inline void SkewSwap (IntMatrix &M, int m, int rowoff, int coloff)
{
  int n = M.Rows();
  for (int row = rowoff; row < rowoff + m; row++) 
    for (int col = coloff; col < coloff + m; col++) 
      MpSwap( M(row,col), M(n-row-1,n-col-1) );
}

//-----------------------------------------------------------------------------//

static void MakeDoubleEvenSquare (IntMatrix &M)
{
  int n     = M.Rows(), 
      m     = n / 4, 
      count = 0;

  // fill in consecutive numbers
  for (int row = 0; row < n; row++) 
    for (int col = 0; col < n; col++) {
      count++;
      M(row, col) = count;
    }
  
  // skew the matrix 
  SkewSwap(M, m, 0, m);
  SkewSwap(M, m, 0, 2*m);
  SkewSwap(M, m, m, 0);
  SkewSwap(M, m, m ,3*m);
}

//-----------------------------------------------------------------------------//

IntMatrix MpMagicSquare (int n)
{
  if (n < 3 || n > 1290) 
    Matpack.Error("MpMagicSquare: dimension must be within 3 to 1290");

  IntMatrix M(0,n-1,0,n-1, 0);

  // the magic sum = ( n * (n*n + 1) ) / 2;

  if ( ((n+1) % 2) == 0 ) 	    // magic square of odd order 2k + 1    
    MakeOddSquare(M, 1,n,0,0);
  else if  ( ((n-2) % 4) == 0 )     // magic square of the order 4k + 2
    MakeSingleEvenSquare(M);
  else                              // magic square of doubly even order 2k
    MakeDoubleEvenSquare(M);
  return M;
}

//-----------------------------------------------------------------------------//
