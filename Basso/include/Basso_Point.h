/*! \file  Basso_Point.h

Basso-ver_1.0

	\author Jack Chessa, jfchessa@utep.edu
	\date Wed Apr 6 2007

*/
#ifndef _BASSO_POINT_H_
#define _BASSO_POINT_H_

#include <iostream>

// Basso includes 
#include "Basso_defs.h"
#include "Basso_Numeric.h"
#include "Basso_Array.h"

namespace Basso {

/**
	\brief Geometric Basso_Point

		This class defines a geometric Basso_Point.  It can also be used to define a node

*/
class Basso_Point
{

public:
    
    /** shallow memory constructor
    \param p - pointer to the coord data
    \param dim - the dimenison of the point
    */
	Basso_Point( Basso_Numeric *p, int dim=3, bool shallow=true )
	{
		dim_=dim;
		if ( shallow )
		{
			xyz_ = p;
			localmem_=false;
		}
		else
		{
			xyz_ = new Basso_Numeric[dim];
			localmem_=true;
			const Basso_Numeric *ptr = p;
			Basso_Numeric *lptr=xyz_;
			for ( int i=0; i<dim_; ++i, ++ptr, ++lptr )
				*lptr = *ptr;
		}
	}
	
	/** deep copy constructor from a pointer */
	Basso_Point( const Basso_Numeric *p, int dim=3 )
	{
		dim_=dim;
		xyz_ = new Basso_Numeric[dim];
		localmem_=true;
		const Basso_Numeric *ptr = p;
		Basso_Numeric *lptr=xyz_;
		for ( int i=0; i<dim_; ++i, ++ptr, ++lptr )
			*lptr = *ptr;

	}
	
	/** Default constructor */
	Basso_Point (Basso_Numeric x=0.0, Basso_Numeric y=0.0, Basso_Numeric z=0.0)
	{
        dim_=3;
        xyz_ = new Basso_Numeric[dim_];
        localmem_=true;
        xyz_[0]=x; xyz_[1]=y, xyz_[2]=z;
	}
	
	
	Basso_Point ( int dim )
	{
        dim_= dim;
        xyz_ = new Basso_Numeric[dim_];
        localmem_=true;
		Basso_Numeric *ptr=xyz_;
		for ( int i=0; i<dim_; ++i, ++ptr )
			*ptr = 0.0;
	}
	
    ~Basso_Point() { if (localmem_) delete xyz_; }
    
    /** conversion operator **/
    operator const Basso_Numeric *() { return xyz_; }
	
	/** returns true if the memory allocation is local else false is returned */
    bool LocalMemory( ) const { return localmem_; }
	
	/** returns the dimension of the point in memory */
	int Dim() const { return dim_; }
	
	/** returns a constant pointer to the coordinates vector */
	const Basso_Numeric *Coords() const { return xyz_; }
	
	/** Returns the Ith coordinate of the Basso_Point with zero offset.  No dimension or out of range
		checking.	*/
    Basso_Numeric &x (int i) {  return xyz_[i];  }
    const Basso_Numeric &x (int i) const { return xyz_[i]; }
    
	/** Returns the [] operator overload. No range chancking. */
    Basso_Numeric &operator [] (int i) { return xyz_[i]; }
    const Basso_Numeric &operator [] (int i) const { return xyz_[i]; }
    
    /** Returns the x coordinate */
    //Basso_Numeric &x() { return x(0); }
    /** Returns the x coordinate */
    Basso_Numeric x() const { return x(0); }
    
    /** Returns the y coordinate */
    //Basso_Numeric &y() { return x(1); }
    /** Returns the y coordinate */
    Basso_Numeric y() const { if (dim_>1) return x(1); return 0.0; }
    
    /** Returns the z coordinate */
    //Basso_Numeric &z() { return x(2); }
    /** Returns the z coordinate */
    Basso_Numeric z() const { if (dim_>2) return x(2); return 0.0; }
    
	/** Computes the distance to the origin */
    Basso_Numeric Distance() const;
    
	/** Computes the distance to the Basso_Point b */
    Basso_Numeric Distance (const Basso_Point &p) const;
    
    /** assignment operator overlaod 
	If the point has shallow memory to start the copy will be shallow.
	If the point hsa local memory then the copy will be local.
	*/
	Basso_Point &operator = ( const Basso_Point &v ) 
	    { 
	        if ( localmem_ )
				for ( int i=0; i<dim_; ++i )
					xyz_[i] = v[i];
	        else
                xyz_ = v.xyz_;
	        
	        return *this; 
	    }

	
protected:
	Basso_Numeric *xyz_;
	int dim_;
    bool localmem_;
	
};

	Basso_Numeric Basso_Point::Distance( const Basso_Point &p ) const
	{
		return sqrt( pow(p.x()-this->x(),2) + pow(p.y()-this->y(),2) + pow(p.z()-this->z(),2) );
	}

	Basso_Numeric Basso_Point::Distance( ) const
	{
		return sqrt( pow(this->x(),2) + pow(this->y(),2) + pow(this->z(),2) );
	}

	/** Basso_Point from I/O stream */
	std::ostream &operator << ( std::ostream &out, const Basso_Point &p ) 
	{
        out.precision(3);
        out.setf(ios::fixed,ios::floatfield);
		
		out << "(" << p.x();
		for ( int i=1; i<p.Dim(); ++i )
			out << "," << p.x(i);
		out << ")";
	
		out.unsetf ( ios_base::showbase );
		return out;
	}
	
	
}  // end of namespace Basso

#endif

