#ifndef _BASSO_POINT_BC_SET_H_
#define  _BASSO_POINT_BC_SET_H_

#include <set>
#include <iostream>
#include <iomanip>

#include "Basso_defs.h"
#include "Basso_Numeric.h"

namespace Basso
{
/** \brief data structure like class for a point boundary conditions. 

This can hold both essential (Dirichelt) and natrual (Neumann) type boundary conditions.  
Typically, this is intended to use a locl node id number, but really can point to a global
node id as well.  
*/
class Basso_PointBC
{
public:
	/**
	Constructor
	\param lnid the node id number (typically local)
	\param ldof the local dof number
	\param value the spc value (default=0)
	*/
	Basso_PointBC( BASSO_IDTYPE lnid, int ldof, Basso_Numeric val=0.0 ) { lnid_=lnid; ldof_=ldof; value_=val; }
	
	/**Returns the node id.*/
	BASSO_IDTYPE LNID()   const { return lnid_; }
	/**Returns the lodal dof for the spc.*/
	int          LDOF()   const { return ldof_; }
	/**Returns the value of the PointBC.*/
	Basso_Numeric Value() const { return value_; }
	
	/**Overload of the rational operator <.*/
	bool operator < ( const Basso_PointBC &other ) const
	{
		if ( lnid_ < other.lnid_ ) 
			return true;
		if ( lnid_ > other.lnid_ )
			return false;
		// lnid are the same
		if ( ldof_ < other.ldof_ )
			return true;
		return false;
	}
	
	void Print( std::ostream &out=BASSO_STDOUT ) const { out << "PointBC (nid/" << lnid_ << "dof/" << ldof_ << ")->" << value_; }
	
protected:
	BASSO_IDTYPE lnid_;
	int ldof_;
	Basso_Numeric value_;
	
};	

std::ostream &operator << ( std::ostream &out, const Basso_PointBC &A )
{
	A.Print( out );
	return out;
}

//---------------------------------------------------------------------------------------------------------------------
/**\brief A set of single point constraints
*/
class Basso_PointBCSet
{
public:
	/**Constructor*/
	Basso_PointBCSet() {}
	
	/**A constant iterator over the set of spcs.  Points to a Basso_PointBC.*/
	typedef set<Basso_PointBC>::const_iterator  ConstIterator;
	
	/**Adds spcs for an array of nodes
	\param the node ids that the spcs are to be added.
	\param ldof local dof of the spc.
	\param val the value of the spc. The default value is 0.0.
	*/
	void AddPointBCs( const Basso_Array<BASSO_IDTYPE> &nids, int ldof, Basso_Numeric val=0.0 );
	
	/**Adds spcs for an array of nodes
	\param the node ids that the spcs are to be added.
	\param ldof local dof of the spc.
	\param val the value of the spc. The default value is 0.0.
	*/
	void AddPointBCs( const set<BASSO_IDTYPE> &nids, int ldof, Basso_Numeric val=0.0 );
	
	
	void AddPointBCs( const set<BASSO_IDTYPE> &nids, const Basso_Array<int> &ldofs, const Basso_Array<Basso_Numeric> &vals );
	
	void AddPointBCs( const set<BASSO_IDTYPE> &nids, const Basso_Array<int> &ldofs );
	
	void AddPointBC( BASSO_IDTYPE nid, int ldof, Basso_Numeric val=0.0 ) { spcs_.insert(Basso_PointBC(nid,ldof,val)); }
	
	/**Clears all the spcs in the set.*/
	void Clear() { spcs_.clear(); }
	
	/**Returns the total number of spcs in the set.*/
	BASSO_IDTYPE NumPointBCs() const { return spcs_.size(); }
	
	/**Returns a constant iterator to the begining of the spc set.*/
	ConstIterator Begin()  const { return spcs_.begin(); }
	
	/**Returns a constant iterator to the just past end of the spc set.*/
	ConstIterator End()    const { return spcs_.end(); }
	
	void Print( std::ostream &out=BASSO_STDOUT ) const;
	
protected:
	set<Basso_PointBC> spcs_;
	
};

void Basso_PointBCSet::AddPointBCs( const Basso_Array<BASSO_IDTYPE> &nids, int ldof, Basso_Numeric val )
{
	const BASSO_IDTYPE *nptr = nids.Data();
	for ( int i=0; i<nids.Length(); ++i, ++nptr )
		spcs_.insert( Basso_PointBC(*nptr,ldof,val) );
}

void Basso_PointBCSet::AddPointBCs( const set<BASSO_IDTYPE> &nids, int ldof, Basso_Numeric val )
{
	set<BASSO_IDTYPE>::const_iterator nitr;
	for ( nitr=nids.begin(); nitr!=nids.end(); ++nitr )
		spcs_.insert( Basso_PointBC(*nitr,ldof,val) );
}

void Basso_PointBCSet::AddPointBCs( const set<BASSO_IDTYPE> &nids, const Basso_Array<int> &ldofs, 
			const Basso_Array<Basso_Numeric> &vals )
{
	set<BASSO_IDTYPE>::const_iterator nitr;
	for ( nitr=nids.begin(); nitr!=nids.end(); ++nitr )
		for ( int i=0; i<ldofs.Length(); ++i )
			spcs_.insert( Basso_PointBC(*nitr,ldofs[i],vals[i]) );
}

void Basso_PointBCSet::AddPointBCs( const set<BASSO_IDTYPE> &nids, const Basso_Array<int> &ldofs )
{
	set<BASSO_IDTYPE>::const_iterator nitr;
	for ( nitr=nids.begin(); nitr!=nids.end(); ++nitr )
		for ( int i=0; i<ldofs.Length(); ++i )
			spcs_.insert( Basso_PointBC(*nitr,ldofs[i],0.0) );
}

void Basso_PointBCSet::Print( std::ostream &out ) const 
{ 
	set<Basso_PointBC>::const_iterator itr;
	out << "SPC set\n           lnid      ldof          value\n";
	for ( itr=spcs_.begin(); itr!=spcs_.end(); ++itr )
		out << setw(15) << itr->LNID() << setw(10) << itr->LDOF() << setw(15) << setprecision(4) << itr->Value() << "\n";
}

std::ostream &operator << ( std::ostream &out, const Basso_PointBCSet &A )
{
	A.Print( out );
	return out;
}


} // of namespace
#endif

