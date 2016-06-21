/*! \file Basso_Array.h

	\brief Defines the class Basso_Array

\author Jack Chessa, jfchessa@utep.edu
\date Sunday December 2, 2007

*/

#ifndef _BASSO_ARRAY_H_
#define _BASSO_ARRAY_H_

// std includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <vector>
#include <math.h>

#include "Basso_defs.h"

namespace Basso
{

template < class T >
class Basso_ArrayIterator
{
public:
	typedef T value_type;
    
	Basso_ArrayIterator( ) { vptr_ = BASSO_NULL_PTR; }
	Basso_ArrayIterator( T *vptr ) { vptr_ = vptr; }
	
	// ++itr
	Basso_ArrayIterator &operator++ ( ) { ++vptr_; return *this; }
	// itr++
	Basso_ArrayIterator operator++( int i ) 
		{ Basso_ArrayIterator<T> itr(vptr_); ++vptr_; return itr; }
	// itr += i	
	Basso_ArrayIterator operator+=( int i ) 
		{ vptr_ += i; return *this; }
		
	// itr - i	
	Basso_ArrayIterator operator-( int i ) 
		{ Basso_ArrayIterator<T> itr(vptr_-i); return itr; }		
	// itr + i	
	Basso_ArrayIterator operator+( int i ) 
		{ Basso_ArrayIterator<T> itr(vptr_+i); return itr; }
		
	T &operator* ( ) { return *vptr_; }
	
	bool operator != ( const Basso_ArrayIterator &itr ) 
		{ 
			if ( vptr_ == itr.vptr_ )
				return false;
			return true;
		}
	bool operator == ( const Basso_ArrayIterator &itr ) 
		{ 
			if ( vptr_ == itr.vptr_ )
				return true;
			return false;
		}
	
private:
	T *vptr_;	
}; 

template < class T >
class Basso_ArrayConstIterator
{
public:
	typedef T value_type;
    
	Basso_ArrayConstIterator( ) { vptr_ = BASSO_NULL_PTR; }
	Basso_ArrayConstIterator( T *vptr ) { vptr_ = vptr; }
	
	// ++itr
	Basso_ArrayConstIterator &operator++ ( ) { ++vptr_; return *this; }
	// itr++
	Basso_ArrayConstIterator operator++( int i ) 
		{ Basso_ArrayIterator<T> itr(vptr_); ++vptr_; return itr; }
	// itr += i	
	Basso_ArrayConstIterator operator+=( int i ) 
		{ vptr_ += i; return *this; }
		
	// itr - i	
	Basso_ArrayConstIterator operator-( int i ) 
		{ Basso_ArrayConstIterator<T> itr(vptr_-i); return itr; }		
	// itr + i	
	Basso_ArrayConstIterator operator+( int i ) 
		{ Basso_ArrayConstIterator<T> itr(vptr_+i); return itr; }
		
	const T &operator* ( ) { return *vptr_; }
	
	bool operator != ( const Basso_ArrayConstIterator &itr ) 
		{ 
			if ( vptr_ == itr.vptr_ )
				return false;
			return true;
		}
	bool operator == ( const Basso_ArrayConstIterator &itr ) 
		{ 
			if ( vptr_ == itr.vptr_ )
				return true;
			return false;
		}
		
	
private:
	const T *vptr_;	
}; 



/**
 \brief Class that defines a vector class that allows for its shape to be
 redefined without really thrashing about any memory.  This is good for finite element
 matricies where the size might change from element to element with the topology

*/

template <class T>
class Basso_Array  
{

public:

	// PUBLIC DEFS
	typedef T value_type;
    typedef BASSO_ARRAY_INDEX indx_type;
    
	typedef Basso_ArrayIterator<T> iterator;
	typedef Basso_ArrayConstIterator<T> const_iterator;

	// CONSTRUCTORS
	/** Void constructor for Basso_Array */
	Basso_Array(  ) : dimMax(25), dim(0) { fVptr = new (std::nothrow) T[dimMax]; }
	
	/** Constructor for Basso_Array
	\param n the size and the maximum row size before memory reallocation
	**/
	Basso_Array( indx_type n ) : dimMax(n), dim(n) { fVptr = new (std::nothrow) T[dimMax]; }
	
	/**
	Constructor from a string
	example:  v = "1.0 2.0 3.14 -2.6";
	*/
    Basso_Array( const string &s ) : dimMax(25), dim(0) { fVptr = new (std::nothrow) T[dimMax]; Copy(s); }
	
	/** Constructor for Basso_Array
	\param r the size 
	\param rmax the maximum  size before memory reallocation
	\param a the value that the vector will be set to
	**/
	Basso_Array( indx_type n, const T &a ) : dimMax(n), dim(n) { fVptr = new (std::nothrow) T[dimMax]; this->SetValues(a); }
		
	/* * Constructor for Basso_Array
	\param r the size 
	\param rmax the maximum  size before memory reallocation
	* */
	//Basso_Array( indx_type r, indx_type rmax ) : dimMax(rmax), dim(r) { fVptr = new T[dimMax]; }

	
	
    /** Copy constructor
    */
    Basso_Array( const Basso_Array<T> &v ) : fVptr(BASSO_NULL_PTR), dimMax(0), dim(0) { Copy(v); }

	/** destructor **/
	virtual ~Basso_Array() 
	{ 	if ( fVptr!=BASSO_NULL_PTR ) 
		//std::cout << "deallocating Basso_Array " << this << " fVptr at " << fVptr << "\n";
			delete [] fVptr; 
	}
	
	/** Sets all values in the array to v */
	void SetValue( T v );	
	
	// Copy constructs
    void Copy( const T *b, indx_type incb=1 );
    void Copy( const T &b );
    void Copy( const Basso_Array<T> &b );
    
    /** copies the vector defined by a string into the existing vector.  This
    will resize if needed.  The string format is whitespace delineated (no commas)
      " 1.0 3.4  2.6 -3 "
    any leading or trailing white spaces are ignored.
    */
    void Copy( const string &b );

    // = overload  ( nonshallow copy )
    inline Basso_Array<T> &operator = ( const Basso_Array<T> &v ) { Copy(v); return *this; }
    inline Basso_Array<T> &operator = ( const T &s ) { Copy(s); return *this; }
    inline Basso_Array<T> &operator = ( const string &s ) { Copy(s); return *this; }
    //inline Basso_Array<T> &operator = ( T s ) { std::cout << "Hi there"; Copy(s); return *this; }

	// ACCESSORS
	/** Returns the jth entry in the vector **/
	inline T &operator [ ] ( indx_type j ) { return *(fVptr+j); }
	inline const T &operator [ ] ( indx_type j ) const { return *(fVptr+j); }
	
	inline T &operator ( ) ( indx_type j ) { return *(fVptr+j); }
	inline const T &operator ( ) ( indx_type j ) const { return *(fVptr+j); }
	
	/**
	Scatters an array of values into the Array.
	*/
	void Scatter( const T *aptr, const BASSO_IDTYPE *iis, indx_type n );
	
	/**
	Scatters an array of values into the Array.
	*/
	template< class NuMvEcToR, class InDxVeCt >
	void Scatter( const NuMvEcToR &v, const InDxVeCt &iis )
	{
		for ( indx_type i=0; i<iis.Length(); ++i )
			fVptr[ iis[i] ] += v[i];
	}
	
	/**
	Gathers the values in the array at the locations iis and
	puts them into aptr.
	\param aptr - on return has the gather data
	\param iis - the indcies to gather from
	\param n - the size of iis (and aptr)
	*/
    void Gather( T *aptr, const BASSO_IDTYPE *iis, indx_type n ) const;
	
    void PushBack( T val );
	
	// iterator accessors
	iterator Begin() { return Basso_ArrayIterator<T>(fVptr); }
	iterator End() { return Basso_ArrayIterator<T>( fVptr + dim ); }
	
	// constant iterator accessors
	const_iterator Begin() const { return Basso_ArrayConstIterator<T>(fVptr); }
	const_iterator End() const { return Basso_ArrayConstIterator<T>( fVptr + dim ); }
	
	/** returns a pointer to the matrix data **/
	T *Data() { return fVptr; }
	
	/** returns a pointer to the matrix data **/
	const T *Data() const { return fVptr; }
	
	// MEMBER FUNCTIONS
	/** returns the number of rows **/
	indx_type Length() const { return dim; }
	indx_type N() const { return dim; }
	indx_type M() const { return dim; }
	indx_type Size() const { return dim; }
	
    indx_type Inc() const { return 1; }
	
	/** returns the maximum number of rows that can be expanded to without reallocation of memory**/
    indx_type Capacity() const { return dimMax; }

	
protected:	
	/** 
	Resizes the vector and keeps the values 
	if memory needs to be reallocated then the capacity is set to \param nmax.  If the
	memory is not reallocated then n\param nmax is not used and the capacity stays the same.
	The old values in the vector are retained.
	**/
	void Resize( indx_type n, indx_type nmax );
	
	/** Acts same as Resize without copying the old values **/
    void Newsize( indx_type n, indx_type nmax );

public:
    /** 
	Resizes the vector and keeps the values 
	The old values in the vector are retained.
	**/
    void Resize( indx_type n ) { Resize(n,0); }  // check this
    
	/** Acts same as Resize without copying the old values **/
    void Newsize( indx_type n ) { Newsize(n,0); }
	
	/** Zeros out the matrix only in the active range **/
	void SetValues( const T &a )
	{
		T *aptr=fVptr;
		for  ( indx_type i=0; i<dim; ++i, ++aptr )
			*aptr = a;
	}

	virtual void Print( std::ostream &out=BASSO_STDOUT ) const 
	{
		const T *ptr = fVptr;
		for  ( indx_type i=0; i<dim; ++i, ++ptr ) 
			out << *ptr << " ";
	}

	/** Write to an ASCII matlab file */
	int WriteMatlab( const string &filename ) const;
	
protected:

	T *fVptr;
	indx_type dimMax, dim;
	
};
template <class T>
void Basso_Array<T>::SetValue( T v )
{
	T *aptr = fVptr;
	for ( indx_type i=0; i<dim; ++i, ++aptr )
    	*aptr = v;
}

template <class T>
int Basso_Array<T>::WriteMatlab( const string &filename ) const
{
	std::ofstream outfile;
	
    outfile.open( filename.c_str(), ios::out );
    if ( !outfile )
        return 1;
  
	const T *aptr = fVptr;
	for ( indx_type i=0; i<dim; ++i, ++aptr )
    	outfile << *aptr << "\n";

	
	outfile.close();
	
	return 0;
   
}

template< class T >
void Basso_Array<T>::Scatter( const T *aptr, const BASSO_IDTYPE *iis, indx_type n )
{
	const T *vptr=aptr;
	const BASSO_IDTYPE *iptr=iis;
	for ( indx_type i=0; i<n; ++i, ++iptr, ++vptr )
		fVptr[ *iptr ] += *vptr;
}

template< class T >
void Basso_Array<T>::Gather( T *aptr, const BASSO_IDTYPE *iis, indx_type n ) const
{
	T *vptr=aptr;
	const BASSO_IDTYPE *iptr=iis;
	for ( indx_type i=0; i<n; ++i, ++iptr, ++vptr )
		 *vptr = fVptr[ *iptr ];
}


template< class T >
void Basso_Array<T>::PushBack( T val )
{
    if ( N() < Capacity() )
    {
        fVptr[dim]=val;
        ++dim;
    }
    else
        Resize( dim, ceil(1.25)*dim );
}
    
template< class T >
void Basso_Array<T>::Copy( const Basso_Array<T> &v )
{
    if ( &v == this ) return;
    
    if ( v.N() != this->N() )
        this->Newsize( v.N() );
        
    this->Copy( v.Data() );     
}

template <class T>
void Basso_Array<T>::Copy( const T *b, indx_type incb )
{
    T *aptr = fVptr;
    const T *bptr = b;
    for ( indx_type i=0; i<dim; ++i, ++aptr, bptr+=incb )
        *aptr = *bptr;
}

template <class T>
void Basso_Array<T>::Copy( const string &b )
{
    indx_type n=0;

    T val;
    string s=string_trim(b);
    stringstream ss (s);
    ss << s;
    while ( !ss.eof() )
    {
        ss >> val;
        ++n;
    }    
    
    this->Newsize(n);
    
    ss.clear();
    ss << s;
    indx_type i=0;
    while ( !ss.eof() && i<n )
        ss >> fVptr[i++];
  
}
    
template <class T>
void Basso_Array<T>::Copy( const T &b )
{
    T *aptr = fVptr;
    for ( indx_type i=0; i<dim; ++i, ++aptr )
        *aptr = b;    
}

template <class T>
void Basso_Array<T>::Resize( indx_type n, indx_type nmax )
{
    if ( nmax < n )
        nmax=n;
    
	if ( n>dimMax )
	{
		T *swapPtr = fVptr;
		
		dimMax = nmax;
		fVptr = new (std::nothrow) T[dimMax];
		
		T *sptr=swapPtr, *vptr=fVptr;
		for ( indx_type i=0; i<dim; ++i, ++sptr, ++vptr )
			*vptr = *sptr;
		
		if ( swapPtr!=BASSO_NULL_PTR ) 
		{
			delete [] swapPtr; 
			swapPtr=BASSO_NULL_PTR;
		} 
                             
	}
	dim=n;
}

template <class T>
void Basso_Array<T>::Newsize( indx_type n, indx_type nmax )
{
    if ( nmax < n )
        nmax=n;
        
	if ( n>dimMax )
	{
		dimMax = nmax;
		if ( fVptr!=BASSO_NULL_PTR ) {
			
			delete[] fVptr; 
			fVptr=BASSO_NULL_PTR;
		}
		fVptr = new (std::nothrow) T[dimMax];
	}
	dim=n;
}

template <class T>
std::ostream &operator << ( std::ostream &out, const Basso_Array<T> &A )
{
	A.Print( out );
	return out;
}


} // end namespace

#endif




