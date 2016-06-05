/*! \file Basso_Array2D.h

	\brief Defines the class Basso_Array2D

\author Jack Chessa, jfchessa@utep.edu
\date Sunday December 2, 2007

*/

#ifndef _BASSO_ARRAY_2D_H_
#define _BASSO_ARRAY_2D_H_

// Basso defs
#include "Basso_defs.h"

// std includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace Basso
{

template < class T >
class Basso_Array2D;
	
template < class T >
class Basso_Array2DIterator
{
	
public:
	typedef T value_type; 
	Basso_Array2DIterator( ) { vptr_ = BASSO_NULL_PTR; }
	Basso_Array2DIterator( T *vptr,  BASSO_ARRAY_INDEX rowi, Basso_Array2D<T> *mptr ) 
		{ vptr_ = vptr; rowi_=rowi; mptr_=mptr; }
	
	// ++itr
	Basso_Array2DIterator &operator++ (  ) 
	{ 
		if ( rowi_ == mptr_->M()-1 )
		{
			vptr_ += mptr_->LDA() - mptr_->M() + 1;
			rowi_=0;
		}
		else
		{
			++vptr_;
			++rowi_; 
		}
		return *this; 
	}
		
	// itr++
	Basso_Array2DIterator operator++( int i ) 
	{ 
		Basso_Array2DIterator<T> itr(vptr_,rowi_,mptr_); 
		if ( rowi_ == mptr_->M()-1 )
		{
			vptr_ += mptr_->LDA() - mptr_->M() + 1;
			rowi_=0;
		}
		else
		{
			++vptr_;
			++rowi_; 
		}
		return itr; 
	}
			
	// itr+=
	Basso_Array2DIterator operator+=( int i ) 
	{ 
		Basso_Array2DIterator<T> itr(vptr_,rowi_,mptr_); 
		if ( rowi_ == mptr_->M()-i )
		{
			vptr_ += mptr_->LDA() - mptr_->M() + i;
			rowi_=0;
		}
		else
		{
			vptr_+=i;
			rowi_+=i; 
		}
		return itr; 
	}
	T &operator* ( ) { return *vptr_; }
	
	bool operator != ( const Basso_Array2DIterator &itr ) 
	{ 
		if ( vptr_ == itr.vptr_ )
			return false;
		return true;
	}
	bool operator == ( const Basso_Array2DIterator &itr ) 
	{ 
		if ( vptr_ == itr.vptr_ )
			return true;
		return false;
	}
		
	void Print( std::ostream &out=BASSO_STDOUT ) const { out << vptr_; }
	
protected:
	T *vptr_;	
	Basso_Array2D<T> *mptr_;
	int rowi_;
}; 

template<class T>
std::ostream &operator << ( std::ostream &out, const Basso_Array2DIterator<T> &itr )
{
	itr.Print( out );
	return out;
}


/*
template<class T>
class Basso_Array2DRowIterator : private Basso_Array2DIterator<T>
{
public:
	Basso_Array2DRowIterator( )  { }

	Basso_ArrayRow2DIterator &operator++ (  ) 
	{ 
	    (*this) += 
		return *this; 
	}
	
	    
};
*/
/**
 \brief Class that defines a matrix class that allows for its shape to be
 redefined without really thrashing about any memory.  This is good for finite element
 matricies where the size might change from element to element with the topology.

 This matrix uses Fortran column format.  So the data in the array is stored in an 
 array pointer as follows:
 
   A.Data() = [ A(0,0), A(1,0), A(2,0), ... A(LDA,0), A(0,1), A(1,1), A(2,1), ... A(LDA,1),
     A(2,0), A(2,1) ......  ]  
 
*/
template<class T>
class Basso_Array2D  {

public:

	// PUBLIC DEFS
	typedef T value_type;
    typedef BASSO_ARRAY_INDEX indx_type;
	typedef Basso_Array2DIterator<T> iterator;

	// CONSTRUCTORS
	/** Void construtor for Basso_Array2D */
	Basso_Array2D(  ) 
	{ 
		rowDimMax=0; colDimMax=0; rowDim=0; colDim=0; 
		fMptr = new (std::nothrow) T[colDimMax*rowDimMax];
	}
	
	/** Constructor for Basso_Array2D
	\param NumRows the row size and the maximum row size before memory reallocation
	\param NumCols the col size and the maximum col size before memory reallocation
	**/
	Basso_Array2D( indx_type NumRows, indx_type NumCols ) 
	{ 
		rowDimMax=NumRows; colDimMax=NumCols; rowDim=NumRows; colDim=NumCols; 
		fMptr = new (std::nothrow) T[colDimMax*rowDimMax];
	}
	
	/** Constructor for Basso_Array2D
	\param r the row size 
	\param rmax the maximum row size before memory reallocation
	\param c the col size 
	\param cmax the maximum col size before memory reallocation
	**/
	Basso_Array2D( indx_type r, indx_type c, indx_type rmax, indx_type cmax ) 
	{ 
		rowDimMax=rmax; colDimMax=cmax; rowDim=r; colDim=c; 
		fMptr = new (std::nothrow) T[colDimMax*rowDimMax];
	}

    /**
    Constructor from a string.  The format of the string is whitespace delineated (no commas) and 
    semicolon to delineate the rows.  If the rows are not of equal length then the largest
    is used as the column dimension and the others are padded with zeros.
      " 1 2 3; 4 5 6; 7 8 9"
    any leading or trailing white spaces are ignored.
    */
    Basso_Array2D( const string &b ) 
	{ 
		rowDimMax=0; colDimMax=0; rowDim=0; colDim=0; 
		fMptr = new (std::nothrow) T[colDimMax*rowDimMax];
        Copy(b);
	}
	
    /**
    Copy constructor
    */
    Basso_Array2D ( const Basso_Array2D<T> &m ) : fMptr(BASSO_NULL_PTR), rowDimMax(0), colDimMax(0), rowDim(0), colDim(0)  
        { Copy(m); }

	/** destructor **/
	virtual ~Basso_Array2D() 
	{ 
		if ( fMptr!=BASSO_NULL_PTR ) 
			delete[] fMptr; 
	}

    /**  Deep copy **/
    void Copy( const Basso_Array2D<T> &m );
    void Copy( const T &A );
    /** copies the 2d array defined by a string into the existing vector.  This
    will resize if needed.  The string format is whitespace delineated (no commas) and 
    semicolon to delineate the rows.  If the rows are not of equal length then the largest
    is used as the column dimension and the others are padded with zeros.
      " 1 2 3; 4 5 6; 7 8 9"
    any leading or trailing white spaces are ignored.
    */
    void Copy( const string &b );
    
    Basso_Array2D<T> &operator = ( const Basso_Array2D<T> &A ) { Copy(A); }
    Basso_Array2D<T> &operator = ( const T &A );
    Basso_Array2D<T> &operator = ( const string &s ) { Copy(s); return *this; }

	// ACCESSORS
	/** Returns the jth column vector.  So A[j][i] returns the ith row jth column 
	entry. In other words A[j][i] = A(i,j) 
	**/
	T *operator [ ] (indx_type j) { return fMptr+j*rowDimMax; }
	
	/** Returns the jth column vector.  So A[j][i] returns the ith row jth column 
	entry. In other words A[j][i] = A(i,j) 
	**/
	const T *operator [ ] (indx_type j)  const { return fMptr+j*rowDimMax; }
	
	/** Returns the  (i,j) value. Row i column j **/
	T &operator ( ) ( indx_type i, indx_type j )  { return fMptr[ j*rowDimMax + i ]; }
	
	/** Returns the  (i,j) value. Row i column j **/
	T operator ( ) ( indx_type i, indx_type j ) const { return fMptr[ j*rowDimMax + i ]; }
	
	/** returns a pointer to the matrix data.  Note the data is stored in column format **/
	T *Data() { return fMptr; }
	
	/** returns a pointer to the matrix data.  Note the data is stored in column format **/
	const T *Data() const { return fMptr; }
	
	/** Sets all values in the array to v */
	void SetValue( T v );
	
	// iterator accessors  
	iterator Begin() { return Basso_Array2DIterator<T>(fMptr,0,this); }
	iterator End() { return Basso_Array2DIterator<T>(fMptr+rowDimMax*colDim,rowDim,this); }
	
	// MEMBER FUNCTIONS
	/** returns the number of rows **/
	indx_type M() const { return rowDim; }
	indx_type NumRows() const { return rowDim; }
	
	/** returns the number of cols **/
	indx_type N() const { return colDim; }
	indx_type NumCols() const { return colDim; }
	
	/** returns the data size used **/
	indx_type Length() const { return colDim*rowDim; }
	
	/** returns the maximum number of rows that can be expanded to without reallocation of memory**/
	indx_type MMax() const { return rowDimMax; }
	
	/** returns the maximum number of cols that can be expanded to without reallocation of memory **/
	indx_type NMax() const { return colDimMax; }
	
	/** returns the allocated data size **/
	indx_type Capacity() const { return colDimMax*rowDimMax; }
	
	/** Returns the leading dimension of the matrix **/
	indx_type LDA() const { return rowDimMax; }

protected:	
	/** reshapes the matrix and sets all values to zero **/
	void Newsize( indx_type NumRows, indx_type NumCols, indx_type RowsMax, indx_type ColsMax );
	void Resize(  indx_type NumRows, indx_type RowsMax, indx_type NumCols, indx_type ColsMax );

public:	
    void Newsize( indx_type NumRows, indx_type NumCols ) {  Newsize( NumRows, NumCols, 0, 0 ); }
	void Resize(  indx_type NumRows, indx_type NumCols ) {  Resize( NumRows, NumCols, 0, 0 ); }
	
	/** Write to an ASCII matlab file */
	int WriteMatlab( const string &filename ) const;
	
	virtual void Print( std::ostream &out=BASSO_STDOUT ) const 
	{
		for  ( indx_type i=0; i<rowDim; ++i ) 
		{
			out << "\n";
			for  ( indx_type j=0; j<colDim; ++j  )
				out << setw(10) << (*this)(i,j) << " ";
		} 
	}

protected:

	indx_type rowDimMax, colDimMax, rowDim, colDim;
	T *fMptr;
	
};
template <class T>
void Basso_Array2D<T>::SetValue( T v )
{
	for ( indx_type i=0; i<rowDim; ++i )
		for ( indx_type j=0; j<colDim; ++j )
    		(*this)[j][i] = v;
}

template <class T>
int Basso_Array2D<T>::WriteMatlab( const string &filename ) const
{
	std::ofstream outfile;
	
    outfile.open( filename.c_str(), ios::out );
    if ( !outfile )
        return 1;

	for ( indx_type i=0; i<rowDim; ++i )
	{
		for ( indx_type j=0; j<colDim; ++j )
    		outfile << (*this)[j][i] << " ";
		outfile << "\n";
	}
   
	return 0;
}

template<class T>
void Basso_Array2D<T>::Copy( const string &b )  
{
	string row;  // determine the number of rows and columns and nnz
	istringstream ss(b);
	indx_type nr=0, nc=0, n=0;
	while ( getline(ss,row,';') ) // get the dimensions
	{
		++nr;
		istringstream rowss(row);
		T val;
		indx_type ncc=0;
		while ( !rowss.eof() )
		{
			rowss >> val;
			++n;
			++ncc;
		}
		if ( ncc > nc ) nc=ncc;
	}
	
	this->Newsize( nr, nc );  // allocate memory
	
	ss.clear();  // go back and fill
	ss.str(b);
	indx_type r=0;
	while ( getline(ss,row,';') ) 
	{
		istringstream rowss(row);
		indx_type c=0;
		while ( !rowss.eof() )
		{
			rowss >> fMptr[c*rowDimMax+r];
			++c;
		}
		++r;
	}
}

/* Old version that works but gives some initializations errors with teh .push_back()
template<class T>
void Basso_Array2D<T>::Copy( const string &b )  // JFC seems to be some issue here
{
    indx_type nr=0, nc=-1, c;

    vector< Basso_Array<T> > rmat;
    string temp, s(b);
    Basso_Array<T> row;
    
    while (s.find(";", 0) != std::string::npos)
    {  
        size_t  pos = s.find(";", 0);  
        temp = s.substr(0, pos);       
        s.erase(0, pos + 1);   
        ++nr;      
        
        row = temp;
        c = row.Length();
        if ( c>nc ) nc=c;
        rmat.push_back( row );     
    }
    row = s;
    c = row.Length();
    if ( c>nc ) nc=c;
    rmat.push_back( row );
    ++nr;
  
    this->Newsize(nr,nc);
    for ( indx_type i=0; i<nr; ++i )
    {
        for ( indx_type j=0; j<nc; ++j )
        {
            if ( j<rmat[i].Length() )
                (*this)(i,j) =  rmat[i][j];
            else
                (*this)(i,j) = static_cast<T>(0);
        }
    }
       
}
*/
template<class T>
void Basso_Array2D<T>::Copy( const T &A )
{
    T *iptr = fMptr;
    for ( indx_type j=0; j<N(); ++j, iptr += rowDimMax - rowDim )
        for ( indx_type i=0; i<M(); ++i, ++iptr )
            *iptr = A;
}

template<class T>
Basso_Array2D<T> &Basso_Array2D<T>::operator = ( const T &A )
{
    Copy(A);
    return *this;
}

template<class T>
void Basso_Array2D<T>::Copy( const Basso_Array2D<T> &m )
{
    if ( &m == this ) return;
    
    // if differing dimsniosn pass to newsize
    if ( m.M()!=M() || m.N()!=N() )
        Newsize( m.M(), m.N() );
      
    // then copy values
    for ( indx_type j=0; j<N(); ++j )
        for ( indx_type i=0; i<M(); ++i )
            fMptr[ i+j*rowDimMax ] = m.fMptr[ i+j*m.rowDimMax ];   
}

template<class T>
void Basso_Array2D<T>::Newsize( indx_type NumRows, indx_type NumCols, indx_type RowsMax, indx_type ColsMax  )
{
    if ( RowsMax<NumRows )  // set for default RowsMax and ColsMax
		RowsMax=max(NumRows,rowDimMax); // NumRows;
        
    if ( ColsMax<NumCols )
		ColsMax=max(NumCols,colDimMax); //NumCols;  

	if ( RowsMax > rowDimMax || ColsMax > colDimMax )
	//if ( NumRows*NumCols > Capacity() )
	{
		if ( fMptr!=NULL ) {
			delete [] fMptr; 
			fMptr=BASSO_NULL_PTR;
		}
		rowDimMax=RowsMax;
		colDimMax=ColsMax;
		fMptr = new (std::nothrow) T[colDimMax*rowDimMax];
	}
		
	rowDim=NumRows;
	colDim=NumCols;
}

template<class T>
void Basso_Array2D<T>::Resize( indx_type NumRows, indx_type NumCols, indx_type RowsMax, indx_type ColsMax )
{

    if ( RowsMax<NumRows )  // set for default RowsMax and ColsMax
		RowsMax=max(NumRows,rowDimMax); // NumRows;
        
    if ( ColsMax<NumCols )
		ColsMax=max(NumCols,colDimMax); //NumCols;    
	
	if ( RowsMax > rowDimMax || ColsMax > colDimMax )
//	if ( NumRows*NumCols > Capacity() )
	{
	    T *swapPtr = fMptr;
		
        indx_type lds = rowDimMax;  
    
		rowDimMax=RowsMax;
		colDimMax=ColsMax;
		fMptr = new (std::nothrow) T[colDimMax*rowDimMax];
		
		for ( indx_type j=0; j<colDim; ++j )
		    for ( indx_type i=0; i<rowDim; ++i )
            fMptr[ i + j*rowDimMax ] = swapPtr[ i + j*lds ];
		
		if ( swapPtr!=BASSO_NULL_PTR ) 
		{
			delete [] swapPtr;  
			swapPtr = BASSO_NULL_PTR;
		}

	}
	
	rowDim=NumRows;
	colDim=NumCols;
}

template<class T>
std::ostream &operator << ( std::ostream &out, const Basso_Array2D<T> &A )
{
	A.Print( out );
	return out;
}

}// end of namespace

#endif




