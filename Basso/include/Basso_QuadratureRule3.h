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
	Basso_QuadratureRule( ) : qpts_(3,0), qwts_(0) { dim_=3; }
	
	/** constructor
	\param n -  the number of points in the rule
	\param d - the spacial dimension of the quadrature points 
	**/
	Basso_QuadratureRule( int n, int d=3 ) : qpts_(d,n), qwts_(n) { dim_=d; qpts_=0; qwts_=1.0/n; }
	
	/** returns the polynomial order that can be integrated exctly **/
	//virtual int Order() const = 0;
	
	/** Returns the number of points in the quadrature rule **/
	int NumPoints( ) const { return qwts_.Length(); }
	
	/** Returns the number of points in the quadrature rule **/
	int Length( ) const { return qwts_.Length(); }
	
	/** Returens the spacial dimension of the qudarature rule **/
    int SDIM() const { return dim_; }

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
	
	/** Resizes the quadrature rule **/
	void Resize( int n, int d=3 ) { this->dim_=d; this->qpts_.Resize(d,n);  this->qwts_.Resize(n); }
	
	/** returns the pointer array for the quadrature points **/
	Basso_Numeric *PointData() { return this->qpts_.Data(); }
	
	/** Returns the pointer to the array of the quadrature weights **/
	Basso_Numeric *WeightData() { return this->qwts_.Data(); }
	
	/** prints the quadrature rule **/
	void Print( std::ostream &out=BASSO_STDOUT ) const;

protected:
	int dim_;
	Basso_nMatrix qpts_;
	Basso_nVector qwts_;

};

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




//***********************************************************************************
/**
\brief class for 1D Gaussian quadrature rules 
*/

class Basso_QuadratureGauss1D : public Basso_QuadratureRule
{

public:
	typedef Basso_QuadratureIterator Iterator;
	typedef Basso_QuadratureIterator ConstIterator;
    
    /** Constructor
    \param n - the number of points in the n quadrature rule **/
    Basso_QuadratureGauss1D( int n ) : Basso_QuadratureRule(n,1)
    {
        if ( n>8 )
            Basso_Warning("Basso_QuadratureGauss1D( int n )","n must be 8 or less.");
            
        quadrature_gauss1d( this->qpts_.Data(), this->qwts_.Data(), n ); 
    }

    virtual int Order() const { return 2*NumPoints()-1; }
};


/**
\brief class for 2D Gaussian quadrature rules 
*/
class Basso_QuadratureGauss2D : public Basso_QuadratureRule
{

public:
	typedef Basso_QuadratureIterator Iterator;
	typedef Basso_QuadratureIterator ConstIterator;
    
    /** Constructor
    \param n - the number of points in the nxn quadrature rule 
    **/
    Basso_QuadratureGauss2D( int n ) : Basso_QuadratureRule(n*n,2)
    {
        n1_=n;
        if ( n>4 )
            Basso_Warning("Basso_QuadratureGauss2D( int n )","n must be 4 or less.");
            
        quadrature_gauss2d( this->qpts_.Data(), this->qwts_.Data(), n ); 
    }

    virtual int Order() const { return 2*n1_ - 1; }

    protected:
        int n1_;
};

/**
\brief class for 3D Gaussian quadrature rules 
*/
class Basso_QuadratureGauss3D : public Basso_QuadratureRule
{

public:
	typedef Basso_QuadratureIterator Iterator;
	typedef Basso_QuadratureIterator ConstIterator;
    
    /** Constructor
    \param n - the number of points in the nxnxn quadrature rule 
    **/
    Basso_QuadratureGauss3D( int n ) : Basso_QuadratureRule(n*n*n,3)
    {
        n1_=n;
        if ( n>3 )
            Basso_Warning("Basso_QuadratureGauss2D( int n )","n must be 3 or less.");
            
        quadrature_gauss3d( this->qpts_.Data(), this->qwts_.Data(), n ); 
    }

    virtual int Order() const { return 2*n1_ - 1; }

    protected:
        int n1_;
};

/**
\brief class for triangular quadrature rules 
*/
class Basso_QuadratureTria : public Basso_QuadratureRule
{

public:
	typedef Basso_QuadratureIterator Iterator;
	typedef Basso_QuadratureIterator ConstIterator;
    
    /** Constructor
    \param n - the number of points in the n-point triangular quadrature rule 
    **/
    Basso_QuadratureTria( int n ) : Basso_QuadratureRule(n,2)
    {
//        if ( n>3 )
//            Basso_Warning("Basso_QuadratureGauss2D( int n )","n must be 3 or less.");
            
        ord_ = quadrature_tria( this->qpts_.Data(), this->qwts_.Data(), n ); 
    }

    virtual int Order() const { return ord_; }

protected:
    int ord_;
    
};

/**
\brief class for tetrahedra quadrature rules 
*/
class Basso_QuadratureTetra : public Basso_QuadratureRule
{

public:
	typedef Basso_QuadratureIterator Iterator;
	typedef Basso_QuadratureIterator ConstIterator;
    
    /** Constructor
    \param n - the number of points in the n-point tetrahedral quadrature rule 
    **/
    Basso_QuadratureTetra( int n ) : Basso_QuadratureRule(n,3)
    {
//        if ( n>3 )
//            Basso_Warning("Basso_QuadratureGauss2D( int n )","n must be 3 or less.");
            
        ord_ = quadrature_tetra( this->qpts_.Data(), this->qwts_.Data(), n ); 
    }

    virtual int Order() const { return ord_; }

protected:
    int ord_;
    
};    
 
} // end namespace
   
    


#endif

