/*! \file Basso_Matrix.h

	\brief Defines the class Basso_Matrix

    A class that inherits from Vector but allows for 
    some numeric operations.

\author Jack Chessa, jfchessa@utep.edu
\date Sunday December 2, 2007

*/

#ifndef _BASSO_NARRAY2D_FEMO_H_
#define _BASSO_NARRAY2D_FEMO_H_

// std includes
#include <iostream>
#include <iomanip>

#include "Basso_defs.h"
#include "Basso_Array2D.h"
#include "Basso_BLAS_wrappers.h"

namespace Basso
{
/**
 \brief 

*/

template <class T>
class Basso_Matrix : public Basso_Array2D<T>
{

public:

	// PUBLIC DEFS
	typedef T value_type;
    typedef BASSO_ARRAY_INDEX indx_type;
	typedef Basso_Array2DIterator<T> iterator;

	// CONSTRUCTORS
	/** Void construtor */
	Basso_Matrix( void ) : Basso_Array2D<T>( ) { }
	Basso_Matrix( const string &s ) : Basso_Array2D<T>( s ) { }
	Basso_Matrix( indx_type NumRows, indx_type NumCols ) : Basso_Array2D<T>(NumRows,NumCols) { }
	Basso_Matrix( indx_type r, indx_type c, indx_type rmax, indx_type cmax ) : Basso_Array2D<T>(r,c,rmax,cmax) { }
	Basso_Matrix( const Basso_Matrix<T> &b ) : Basso_Array2D<T>(b) { }

	~Basso_Matrix() { }
    
	// basic linear algebra operations (mostly via BLAS)
	
	/** Sets all the values in the matrix to zero. */
	void Zero() { this->SetValue( 0.0 ); }
	
    // operator overloads    
    inline Basso_Matrix<T> &operator = ( const Basso_Matrix<T> &v ) { this->Copy(v); return *this; }
    inline Basso_Matrix<T> &operator = ( const T &s ) { this->Copy(s); return *this; }
    inline Basso_Matrix<T> &operator = ( const string &s ) { this->Copy(s); return *this; }
	
	virtual void Print( std::ostream &out=BASSO_STDOUT ) const 
	{
        out.precision(5);
		for  ( indx_type i=0; i<this->rowDim; ++i ) 
		{
			out << "\n";
			for  ( indx_type j=0; j<this->colDim; ++j  )
				out << scientific << setw(14) << (*this)(i,j) << " ";
		} 
		out << resetiosflags (ios_base::showbase);
	}
	  
private:
    	
};

} // end namespace

#endif




