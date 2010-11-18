/*-----------------------------------------------------------------------------*\
| I/O functions for vectors and matrices                           loadstore.cc |
|                                                                               |
| Last change: Mar 7, 1998							|
|                                                                               |
| Matpack Library Release 1.3                                                   |
| Copyright (C) 1991-1998 by Berndt M. Gammel. All rights reserved.             |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.                                                |
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.                                                    |
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without prior consent of the author.                       |
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration |
| and installation of Matpack.                                                  |
|                                                                               |
\*-----------------------------------------------------------------------------*/

#include "matpack.h"

//-----------------------------------------------------------------------------//
// the following funcions are defined in this file:
//-----------------------------------------------------------------------------//

ostream& operator << (ostream& f, const Vector& x);
istream& operator >> (istream& f, Vector& x);

ostream& operator << (ostream& f, const ComplexVector& x);
istream& operator >> (istream& f, ComplexVector& x);

ostream& operator << (ostream& f, const IntVector& x);
istream& operator >> (istream& f, IntVector& x);

ostream& operator << (ostream& f, const Matrix& x);
istream& operator >> (istream& f, Matrix& x);

ostream& operator << (ostream& f, const ComplexMatrix& x);
istream& operator >> (istream& f, ComplexMatrix& x);

ostream& operator << (ostream& f, const IntMatrix& x);
istream& operator >> (istream& f, IntMatrix& x);

//-----------------------------------------------------------------------------//
// write real vector to ostream
//-----------------------------------------------------------------------------//

ostream& operator << (ostream& f, const Vector& x)
{
  int arraydim[2] = { x.Lo(), x.Hi() };
  int fmt = x.Format();

  MpWriteDataFormat (f,fmt,
		     MpRowMajor,
		     1,MpRealValued,MpScalar,
		     fmt & MpSingle,
		     fmt & MpSwapByteOrder,
		     arraydim);
  
  MpSave1dArrayOfRealScalar(f,x,fmt,
			    fmt & MpSingle,
			    fmt & MpSwapByteOrder);

  return f;
}

//-----------------------------------------------------------------------------//
// read real vector from istream
//-----------------------------------------------------------------------------//

istream& operator >> (istream& f, Vector& x)
{
  MpDataFormat d;

  // assert appropriate format
  if (! (MpReadDataFormat(f,d) && 
	 d.data == MpArray &&
	 d.array == 1 &&
	 d.element == MpScalar && 
	 d.number == MpRealValued &&
	 (x.Empty () ||   // either empty or conformant dimensions !
	  (d.arraydim[0] == x.Lo() && d.arraydim[1] == x.Hi())) &&
	 MpLoad1dArrayOfRealScalar(f,x,d.arraydim,d.format,d.precision,d.endian))
      )
    f.clear(ios::failbit); // set error flag for stream
  
  return f;
}

//-----------------------------------------------------------------------------//
// write complex vector to ostream
//-----------------------------------------------------------------------------//

ostream& operator << (ostream& f, const ComplexVector& x)
{
  int arraydim[2] = { x.Lo(), x.Hi() };
  int fmt = x.Format();

  MpWriteDataFormat (f,fmt,
		     MpRowMajor,
		     1,MpComplexValued,MpScalar,
		     fmt & MpSingle,
		     fmt & MpSwapByteOrder,
		     arraydim);
  
  MpSave1dArrayOfComplexScalar(f,x,fmt,
			       fmt & MpSingle,
			       fmt & MpSwapByteOrder);

  return f;
}

//-----------------------------------------------------------------------------//
// read complex vector from istream
//-----------------------------------------------------------------------------//

istream& operator >> (istream& f, ComplexVector& x)
{
  MpDataFormat d;

  // assert appropriate format
  if (! (MpReadDataFormat(f,d) && 
	 d.data == MpArray &&
	 d.array == 1 &&
	 d.element == MpScalar && 
	 d.number == MpComplexValued &&
	 (x.Empty () ||   // either empty or conformant dimensions !
	  (d.arraydim[0] == x.Lo() && d.arraydim[1] == x.Hi())) &&
	 MpLoad1dArrayOfComplexScalar(f,x,d.arraydim,d.format,d.precision,d.endian))
      )
    f.clear(ios::failbit); // set error flag for stream
  
  return f;
}

//-----------------------------------------------------------------------------//
// write real matrix to ostream
//-----------------------------------------------------------------------------//

ostream& operator << (ostream& f, const Matrix& x)
{
  int arraydim[4] = { x.Rlo(), x.Rhi(), x.Clo(), x.Chi() };
  int fmt = x.Format();

  MpWriteDataFormat (f,fmt,
		     fmt & MpColumnMajor,
		     2,MpRealValued,MpScalar,
		     fmt & MpSingle,
		     fmt & MpSwapByteOrder,
		     arraydim);
  
  MpSave2dArrayOfRealScalar(f,x,
			    fmt & MpColumnMajor,
			    fmt,
			    fmt & MpSingle,
			    fmt & MpSwapByteOrder);

  return f;
}

//-----------------------------------------------------------------------------//
// read real matrix from istream
//-----------------------------------------------------------------------------//

istream& operator >> (istream& f, Matrix& x)
{
  MpDataFormat d;

  // assert appropriate format
  if (! (MpReadDataFormat(f,d) && 
	 d.data == MpArray &&
	 d.array == 2 &&
	 d.element == MpScalar && 
	 d.number == MpRealValued &&
	 (x.Empty () ||   // either empty or conformant dimensions !
	  (d.arraydim[0] == x.Rlo() && d.arraydim[1] == x.Rhi() &&
	   d.arraydim[2] == x.Clo() && d.arraydim[3] == x.Chi())) &&
	 MpLoad2dArrayOfRealScalar(f,x,d.arraydim,d.order,
				   d.format,d.precision,d.endian))
      )
    f.clear(ios::failbit); // set error flag for stream
  
  return f;
}

//-----------------------------------------------------------------------------//
// write complex matrix to ostream
//-----------------------------------------------------------------------------//

ostream& operator << (ostream& f, const ComplexMatrix& x)
{
  int arraydim[4] = { x.Rlo(), x.Rhi(), x.Clo(), x.Chi() };
  int fmt = x.Format();

  MpWriteDataFormat (f,fmt,
		     fmt & MpColumnMajor,
		     2,MpComplexValued,MpScalar,
		     fmt & MpSingle,
		     fmt & MpSwapByteOrder,
		     arraydim);
  
  MpSave2dArrayOfComplexScalar(f,x,
			       fmt & MpColumnMajor,
			       fmt,
			       fmt & MpSingle,
			       fmt & MpSwapByteOrder);

  return f;
}

//-----------------------------------------------------------------------------//
// read complex matrix from istream
//-----------------------------------------------------------------------------//

istream& operator >> (istream& f, ComplexMatrix& x)
{
  MpDataFormat d;

  // assert appropriate format
  if (! (MpReadDataFormat(f,d) && 
	 d.data == MpArray &&
	 d.array == 2 &&
	 d.element == MpScalar && 
	 d.number == MpComplexValued &&
	 (x.Empty () ||   // either empty or conformant dimensions !
	  (d.arraydim[0] == x.Rlo() && d.arraydim[1] == x.Rhi() &&
	   d.arraydim[2] == x.Clo() && d.arraydim[3] == x.Chi())) &&
	 MpLoad2dArrayOfComplexScalar(f,x,d.arraydim,d.order,
				      d.format,d.precision,d.endian))
      )
    f.clear(ios::failbit); // set error flag for stream
  
  return f;
}

//-----------------------------------------------------------------------------//
// write integer vector to ostream
//-----------------------------------------------------------------------------//

ostream& operator << (ostream& f, const IntVector& x)
{
  int arraydim[2] = { x.Lo(), x.Hi() };
  int fmt = x.Format();

  MpWriteDataFormat (f,fmt,
		     MpRowMajor,
		     1,MpIntegerValued,MpScalar,
		     fmt & MpSingle,
		     fmt & MpSwapByteOrder,
		     arraydim);
  
  MpSave1dArrayOfIntegerScalar(f,x,fmt,
			       fmt & MpSingle,
			       fmt & MpSwapByteOrder);

  return f;
}

//-----------------------------------------------------------------------------//
// read integer vector from istream
//-----------------------------------------------------------------------------//

istream& operator >> (istream& f, IntVector& x)
{
  MpDataFormat d;

  // assert appropriate format
  if (! (MpReadDataFormat(f,d) && 
	 d.data == MpArray &&
	 d.array == 1 &&
	 d.element == MpScalar && 
	 d.number == MpIntegerValued &&
	 (x.Empty () ||   // either empty or conformant dimensions !
	  (d.arraydim[0] == x.Lo() && d.arraydim[1] == x.Hi())) &&
	 MpLoad1dArrayOfIntegerScalar(f,x,d.arraydim,d.format,d.precision,d.endian))
      )
    f.clear(ios::failbit); // set error flag for stream
  
  return f;
}

//-----------------------------------------------------------------------------//
// write integer matrix to ostream
//-----------------------------------------------------------------------------//

ostream& operator << (ostream& f, const IntMatrix& x)
{
  int arraydim[4] = { x.Rlo(), x.Rhi(), x.Clo(), x.Chi() };
  int fmt = x.Format();

  MpWriteDataFormat (f,fmt,
		     fmt & MpColumnMajor,
		     2,MpIntegerValued,MpScalar,
		     fmt & MpSingle,
		     fmt & MpSwapByteOrder,
		     arraydim);
  
  MpSave2dArrayOfIntegerScalar(f,x,
			       fmt & MpColumnMajor,
			       fmt,
			       fmt & MpSingle,
			       fmt & MpSwapByteOrder);

  return f;
}

//-----------------------------------------------------------------------------//
// read integer matrix from  istream
//-----------------------------------------------------------------------------//

istream& operator >> (istream& f, IntMatrix& x)
{
  MpDataFormat d;

  // assert appropriate format
  if (! (MpReadDataFormat(f,d) && 
	 d.data == MpArray &&
	 d.array == 2 &&
	 d.element == MpScalar && 
	 d.number == MpIntegerValued &&
	 (x.Empty () ||   // either empty or conformant dimensions !
	  (d.arraydim[0] == x.Rlo() && d.arraydim[1] == x.Rhi() &&
	   d.arraydim[2] == x.Clo() && d.arraydim[3] == x.Chi())) &&
	 MpLoad2dArrayOfIntegerScalar(f,x,d.arraydim,d.order,
				      d.format,d.precision,d.endian))
      )
    f.clear(ios::failbit); // set error flag for stream
  
  return f;
}

//-----------------------------------------------------------------------------//
