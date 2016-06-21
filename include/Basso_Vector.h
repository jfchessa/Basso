/*! \file Basso_Vector.h

	\brief Defines the class Basso_Vector

    A class that inherits from Vector but allows for 
    some numeric operations.

\author Jack Chessa, jfchessa@utep.edu
\date Sunday December 2, 2007

*/

#ifndef _BASSO_NARRAY_FEMO_H_
#define _BASSO_NARRAY_FEMO_H_

// std includes
#include <iostream>
#include <iomanip>

#include "Basso_defs.h"
#include "Basso_Array.h"
#include "Basso_BLAS_wrappers.h"

namespace Basso
{
/**
 \brief 

*/

template <class T>
class Basso_Vector : public Basso_Array<T>
{

public:

	// PUBLIC DEFS
	typedef T value_type;
    typedef BASSO_ARRAY_INDEX indx_type;
	typedef Basso_ArrayIterator<T> iterator;

	// CONSTRUCTORS
	/** Void construtor for Vector */
	Basso_Vector( void ) : Basso_Array<T>( ) { }
	Basso_Vector( const string &s ) : Basso_Array<T>( s ) { }
	Basso_Vector( indx_type n ) : Basso_Array<T>(n) { }
	Basso_Vector( indx_type n, T v ) : Basso_Array<T>(n,v) { }
	Basso_Vector( const Basso_Vector<T> &b ) : Basso_Array<T>(b) { }
	
	 /** conversion operator; this allows for the point to convert to a Basso_Array  **/
    //operator Basso_Array<T> &() { return  *this; }

    // basic linear algebra operations (mostly via BLAS)
    void SCAL( const T a ) { scal( this->N(), a, this->Data(), 1 ); }
    
    T Norm( ) const { return nrm2( this->N(), this->Data(), 1 ); }
    
    T Dot( const T *x, int incx=1 ) const 
        { return dot( this->N(), this->Data(), 1, x, incx ); }
    
    void AXPY( const T a, const T *x, int incx=1 ) { axpy( this->N(), a, x, incx, this->Data(), 1 ); }
    
	/** sets all the values in the vector to zero */
	void Zero() { this->SetValue( 0.0 ); }
	
    // operator overloads    
    inline Basso_Vector<T> &operator = ( const Basso_Vector<T> &v ) { this->Copy(v); return *this; }
    inline Basso_Vector<T> &operator = ( const T &s ) { this->Copy(s); return *this; }
    inline Basso_Vector<T> &operator = ( const string &s ) { this->Copy(s); return *this; }
    
    inline Basso_Vector<T> &operator += ( const Basso_Vector<T> &v ) { AXPY( 1.0, v.Data() ); return *this; }
    inline Basso_Vector<T> &operator -= ( const Basso_Vector<T> &v ) { AXPY( -1.0, v.Data() ); return *this; }

	virtual void Print( std::ostream &out=BASSO_STDOUT ) const 
	{
		const T *ptr = this->fVptr;
		out.precision(5);
		for  ( indx_type i=0; i<this->dim; ++i, ++ptr ) 
			out << scientific << setw(14) << *ptr << " ";
        out << resetiosflags (ios_base::showbase);
	}
	
private:
    	
};

template <class T>
Basso_Vector<T> operator+(const Basso_Vector<T> &A, const Basso_Vector<T> &B)
{

	typedef typename Basso_Vector<T>::indx_type indx_type;
	indx_type n = A.N();

	if ( B.N() != n )
		return Basso_Vector<T>();

	else
	{
		Basso_Vector<T> C(n);

		for (BASSO_ARRAY_INDEX i=0; i<n; i++)
			C[i] = A[i] + B[i];

        return C; 
	}
}

template <class T>
Basso_Vector<T> operator-(const Basso_Vector<T> &A, const Basso_Vector<T> &B)
{
	typedef typename Basso_Vector<T>::indx_type indx_type;
	indx_type n = A.N();
	if (B.N() != n )
		return Basso_Vector<T>();

	else
	{
		Basso_Vector<T> C(n);

		for (indx_type i=0; i<n; i++)
		{
			C[i] = A[i] - B[i];
		}
		return C;
	}
}

template <class T>
Basso_Vector<T> operator * ( const Basso_Vector<T> &A, const T &s )
{
	typedef typename Basso_Vector<T>::indx_type indx_type;

	indx_type n = A.N();
	Basso_Vector<T> C(n);
	for (indx_type i=0; i<n; i++)
	{
		C[i] = s*A[i];
	}
	return C;
}
template <class T>
Basso_Vector<T> operator * ( const T &s, const Basso_Vector<T> &A )
    { return A*s; }

template <class T>
T dot( const Basso_Vector<T> &a, const Basso_Vector<T> &b ) // { return a.Dot(b.Data()); }
{ 
	typedef typename Basso_Vector<T>::indx_type indx_type;
	T dp=0.0;
	const T *aptr=a.Data(), *bptr=b.Data();
	for ( indx_type i=0; i<a.N(); ++i, ++aptr, ++bptr )
		dp += (*aptr)*(*bptr);
	return dp;  //ddot( a.N(), a.Data(), a.Inc(), b.Data(), b.Inc() );
}

template <class T>
T norm( const Basso_Vector<T> &a ) { return a.Norm(); }

} // end namespace

#endif




