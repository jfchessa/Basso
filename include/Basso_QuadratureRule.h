/*! \file Basso_QuadratureRule.h

	\brief Defines the class Basso_QuadratureRule

\author Jack Chessa, jfchessa@utep.edu
\date Sunday April 21, 2012

*/

#ifndef _BASSO_QUADRATURE_RULE_H_
#define _BASSO_QUADRATURE_RULE_H_

// Basso defs
#include "Basso_defs.h"
#include "Basso_Numeric.h"
#include "Basso_Point.h"
#include "Basso_nVector.h"
#include "Basso_nMatrix.h"
#include "Basso_feaops.h"

namespace Basso
{
/*************************************/
/**
\brief a wrapper type quadrature point
This class acts really as a point for the quadrature rule iterators.  It is not
meant to work on its own.
**/
class Basso_QuadraturePoint
{
    friend class Basso_QuadratureIterator;
	
	/** void constructor, used only by friend classes **/
    Basso_QuadraturePoint( void ) { Set(NULL,NULL,0,0); /*pPtr_=NULL; wPtr_=NULL; dim_=0;*/ }
	
	/** constructor 
	\param p - a pointer (non-volitle) to the point
	\param w - a pointer (non-volitle) to the weight
	\param d - the dimension of the point (1,2,3)  The default is 3.
	**/
    Basso_QuadraturePoint( const Basso_Numeric *p, const Basso_Numeric *w, int d=3, int l=3 ) { Set(p,w,d,l); }

public:
    /** Returns a const pointer to the quadratuter point **/
	const Basso_Numeric *PtPtr( ) const { return pPtr_; }
	const Basso_Point  Pt( ) const { return Basso_Point(pPtr_,dim_); }
	
	/** Returns the quadrature weight **/
	Basso_Numeric Wt( ) const  { return *wPtr_; }
	
	/** returns the spacial dimension of the quadruature point **/
    int SDIM() const { return dim_; }
    
    /** Prints the point **/
	void Print( std::ostream &out=BASSO_STDOUT ) const
	{	
        out << "Quadrature point = { ";
        for  ( int i=0; i<dim_; ++i )
            cout << pPtr_[i] << " ";
		cout << "}, weight = " << *wPtr_ << "\n";	
	}
		
	/** returns true only it points to the same memory **/
	bool operator == ( const Basso_QuadraturePoint &pt )  const
	{
		if ( pt.pPtr_ == this->pPtr_  && pt.wPtr_ == this->wPtr_  && pt.dim_ == this->dim_ )
			return true;
		return false;
	}	
	
protected:		
	/** Moves the quadrature point to the next in memory.  Recalll this is 
	really only to be used by the iterators **/	
	Basso_QuadraturePoint &operator ++ (  ) 
	{ 
		pPtr_+=lda_; ++wPtr_;
		return  *this;
	}
	
	/** Moves the quadrature point up n in memory.  Recalll this is 
	really only to be used by the iterators 
	\param n - the number of points to move up
	**/	
	Basso_QuadraturePoint &operator += ( int n ) 
	{ 
		pPtr_+=lda_*n; wPtr_+=n;
		return  *this;
	}
    	
	void Set( const Basso_Numeric *p, const Basso_Numeric *w, int d, int l ) 
        { pPtr_=p; wPtr_=w; dim_=d; lda_=l; }
		
protected:
    int dim_, lda_;
	const Basso_Numeric *pPtr_;
	const Basso_Numeric *wPtr_;
		
};

/** overload for Basso_QuadraturePoint **/
std::ostream &operator << ( std::ostream &out, const Basso_QuadraturePoint &A )
{
	A.Print( out );
	return out;
}

/*************************************/
/**
\brief constant iterator for a quadrature rule
This class acts as a constant iterator for a quadrature rule class.
**/
class Basso_QuadratureIterator
{
    
public:
    /** void constructor **/
    Basso_QuadratureIterator( ) { }
    	
	/** constructor 
	\param p - a pointer (non-volitle) to the point
	\param w - a pointer (non-volitle) to the weight
	\param d - the dimension of the point (1,2,3)  The default is 3.
	**/
	Basso_QuadratureIterator(const Basso_Numeric *p, const Basso_Numeric *w, int d=3, int l=3 ) 
		: qpt_(p,w,d,l) { }
	
	/** reference operator **/
	const Basso_QuadraturePoint *operator -> () const { return &qpt_; }
	
	/** dereference operator **/
	const Basso_QuadraturePoint &operator * () const { return qpt_; }
	
	/** Unitary increment to the next quadrature point in the rule **/
	Basso_QuadratureIterator &operator ++ (  ) 
	{ 
		++qpt_;
		return  *this;
	}
	
	/** N increment to the  quadrature point n away in the rule 
	\param n- the number of increments
	**/
	Basso_QuadratureIterator &operator += ( int n ) 
	{ 
		qpt_ += n;
		return  *this;
	}
	
	/** == operator overload **/
	bool operator == ( const Basso_QuadratureIterator &itr )  const
	{ 
		if ( this->qpt_ == itr.qpt_  )
			return true;
		return false;
	}
	
	/** != operator overload **/
	bool operator != ( const Basso_QuadratureIterator &itr ) const 
	{ 
		if ( this->qpt_ == itr.qpt_  )
			return false;
		return true;
	}
	
protected:
	Basso_QuadraturePoint qpt_;
	
};

/*************************************/

/**
\brief A basic general purpose quadrature rule class
**/

enum Basso_QuadratureType { Basso_GAUSS1D_QUAD, Basso_GAUSS2D_QUAD, Basso_GAUSS3D_QUAD, 
                 Basso_TRIA_QUAD, Basso_TETRA_QUAD, Basso_GEN_QUAD };

class Basso_QuadratureRule
{
    
public:
	typedef Basso_QuadratureIterator Iterator;
	typedef Basso_QuadratureIterator ConstIterator;
	
    /** void constructor **/
	Basso_QuadratureRule( ) : qpts_(3,0), qwts_(0) { dim_=3; qtype_=Basso_GEN_QUAD; }
	
	/** constructor
	\param n -  the number of points in the rule
	\param d - the spacial dimension of the quadrature points 
	**/
	Basso_QuadratureRule( int n, int d=3 ) : qpts_(d,n), qwts_(n) { dim_=d; qpts_=0; qwts_=1.0/n; }
	
	/** constructor
	\param type -  The integration type
	\param order - The order polynomial that can be exactly integrated
	**/
	Basso_QuadratureRule( Basso_QuadratureType qtype, int order ) { SetQuadratureRule(qtype,order); }
	
	void SetQuadratureRule( Basso_QuadratureType qtype, int order );
	
	/** returns the polynomial order that can be integrated exctly **/
	//virtual int Order() const = 0;
	
	/** Returns the number of points in the quadrature rule **/
	int NumPoints( ) const { return qwts_.Length(); }
	
	/** Returns the number of points in the quadrature rule **/
	int Length( ) const { return qwts_.Length(); }
	
	/** Returens the spacial dimension of the qudarature rule **/
    int SDIM() const { return dim_; }
	
	/** Returns the quadrature type for the quadrature rule. */
	Basso_QuadratureType QuadratureType() const { return qtype_; }

	/** returns the leading dimension of the quadrature point matrix */
	int LDA() const { return qpts_.LDA(); }
	
	/**Returns a constant iterator to the first point in the quadrature rule**/
	Basso_QuadratureIterator Begin() const { return Basso_QuadratureIterator(qpts_.Data(),qwts_.Data(),dim_,qpts_.LDA()); }
	
	/**Returns a constant interator to the just past the end point in the quadrature rule. **/
	Basso_QuadratureIterator End() const 
	{ 
		Basso_QuadratureIterator itr(qpts_.Data(),qwts_.Data(),dim_,qpts_.LDA());
		itr += qwts_.Length();
		return itr; 
	}
	
	/** prints the quadrature rule **/
	void Print( std::ostream &out=BASSO_STDOUT ) const;

//protected:	
	/** Resizes the quadrature rule **/
	void Resize( int n, int d=3 ) { this->dim_=d; this->qpts_.Resize(d,n);  this->qwts_.Resize(n); }
	
	/** returns the pointer array for the quadrature points **/
	Basso_Numeric *PointData() { return this->qpts_.Data(); }
	
	/** Returns the pointer to the array of the quadrature weights **/
	Basso_Numeric *WeightData() { return this->qwts_.Data(); }
	
protected:
	int dim_;
	Basso_nMatrix qpts_;
	Basso_nVector qwts_;
	Basso_QuadratureType qtype_;
	
};


void Basso_QuadratureRule::SetQuadratureRule( Basso_QuadratureType qtype, int order ) 
{
	this->qtype_ = qtype;
	switch (qtype)
	{
		int npts, npt1d;
		case Basso_GAUSS1D_QUAD: 
			dim_ = 1;
			npt1d = ceil( 0.5*(order+1.0) );
			npts = npt1d;
			qpts_.Resize(dim_,npts);
			qwts_.Resize(npts);
			quadrature_gauss1d( qpts_.Data(), qwts_.Data(), npt1d, qpts_.LDA() );
			break;
			
		case Basso_GAUSS2D_QUAD:
			dim_ = 2;
			npt1d = ceil( 0.5*(order+1.0) );
			npts = npt1d*npt1d;
			qpts_.Resize(dim_,npts);
			qwts_.Resize(npts);
			quadrature_gauss2d( qpts_.Data(), qwts_.Data(), npt1d, qpts_.LDA() ); 
			break;
			
		case Basso_GAUSS3D_QUAD:
			dim_ = 3;
			npt1d = ceil( 0.5*(order+1.0) );
			npts = npt1d*npt1d*npt1d;
			qpts_.Resize(dim_,npts);
			qwts_.Resize(npts);
			quadrature_gauss2d( qpts_.Data(), qwts_.Data(), npt1d, qpts_.LDA() );
			break;
			
		case Basso_TRIA_QUAD:
			dim_ = 2;
			switch (order)
			{
				case 0:
				npts=1;
				break;

				case 1:
				npts=1;
				break;

				case 2:
				npts=3;
				break;	

				case 3:
				npts=4;
				break;

				case 4:
				npts=6;
				break;

				case 5:
				npts=7;
				break;

				case 6:
				npts=12;
				break;

				case 7:
				npts=13;
				break;	

				default:
				Basso_Warning("Basso_QuadratureRule( Basso_QuadratureType qtype, int order )","quadrature order >7");
				npts=13;
				break;

			}
			qpts_.Resize(dim_,npts);
			qwts_.Resize(npts);
			quadrature_tria( qpts_.Data(), qwts_.Data(), npts, qpts_.LDA() );
			break;
			
		case Basso_TETRA_QUAD:
			dim_ = 3;
			if ( order<=1 )
				npts=1;
			else if ( order<=2 )
				npts=4;
			else if ( order<=3 )
				npts=5;
			else {
				npts=5;
				Basso_Warning("Basso_QuadratureRule( Basso_QuadratureType qtype, int order )","quadrature order >3");
			}
			qpts_.Resize(dim_,npts);
			qwts_.Resize(npts);
			quadrature_tetra( qpts_.Data(), qwts_.Data(), npts, qpts_.LDA() );
			break;
			
/*		case Basso_GEN_QUAD:
		default:  */
	}
	
}

void Basso_QuadratureRule::Print( std::ostream &out ) const 
{
	out << "Quadrature Rule: maxdim=" << qpts_.MMax() << ", lda=" << qpts_.LDA() << "\n";
	for ( int i=0; i<qpts_.N(); ++i )
	{
		out << i << " : pt = { "; 
		int j;
		for ( j=0; j<dim_; ++j )
			out << qpts_[i][j] << " ";
		out << "},  wt = " << qwts_[i] << "\n";
	}
}

/** print stream operator overload **/
std::ostream &operator << ( std::ostream &out, const Basso_QuadratureRule &A )
{
	A.Print( out );
	return out;
}



 
} // end namespace
  
#endif

