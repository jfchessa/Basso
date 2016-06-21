/*! \file Basso_FixedDOFMap.h

	\brief Defines the class Basso_FixedDOFMap

This class requires Trilinos Epetra	
	
\author Jack Chessa, jfchessa@utep.edu
\date Thursday, June 8, 2016

*/

#ifndef _BASSO_DOF_MAP_H_
#define _BASSO_DOF_MAP_H_

// std includes
#include <map>
#include <iomanip>

// trilinos includes

// basso includes
#include "Basso_defs.h"
#include "Basso_Array.h"
#include "Basso_DOFMap.h"

namespace Basso
{

	
/**
 \brief DOF map with a fixed number of active dofs per node.
	
*/
class Basso_FixedDOFMap : public Basso_DOFMap
{
public:
	/**
	Constructor
	\param gnids - an array of the global node ids that are supported by the processor.  This array needs to 
	map onto the rows of the local node coordinate matrix
	\param nldof the number of dofs active on each node.  Default value is 6.
	*/
	Basso_FixedDOFMap( const Basso_Array<BASSO_IDTYPE> &gnids, int nldof=6 ) 
				: gdof_map_(nldof,gnids.Length())
		{ fixed_=false; SetGNIDs(gnids); FixDOFMap(  ); }

	
	Basso_FixedDOFMap( int numnode, int nldof=6 ) : gdof_map_(nldof,numnode)
		{ fixed_=false; SetGNIDs(numnode); FixDOFMap(  ); }
	
		
	/** Deconstructor */
	~Basso_FixedDOFMap() {  }
		
	/** Sets the global node ids that are supported by the procesor.
	\param gnids - an array of the global node ids that are supported by the processor.  This array needs to 
	map onto the rows of the local node coordinate matrix
	*/	
	void SetGNIDs( const Basso_Array<BASSO_IDTYPE> &gnids ) { gnids_ = gnids; FixDOFMap(  ); }
	void SetGNIDs( int nn ) 
	{ 
		gdof_map_.Resize(NumDofPerNode(),nn);
		gnids_.Resize(nn); 
		for ( int I=0; I<nn; ++I )
			gnids_[I]=I;
		FixDOFMap(  );
	}
	
	/**
	Initialized the global dof map.  This takes the global node ids supported on each processor
	and determines a unique global id
	*/
	void FixDOFMap(  );
	
	/**Returns the number of active dofs for each node.*/
	int NumDofPerNode() const { return gdof_map_.NumRows(); }
	
	/**Returns the number of nodes supported on the local processor.*/
	int NumMyNodes()    const { return gnids_.Length(); }
	
	int NumDOF() const { return gdof_map_.Length(); }
	
	/**Initializes a scatter vector for the gdof map.  There is no bounds checking to expidite 
	the member function.  The form of the scatter vector is as follows:
		sctr = [ node0_dof0, node0_dof1, node0_dof2 ..., node1_dof0, node1_dof2, ... ] 
	\param econn a vector of the local node ids for which to compute the scatter vector. Do not use global node ids here.
	\param sctr on return has the scatter vector to the global dofs.  It should be of length econn.Length()*NumDofPerNode(). 
	*/
	void SetScatter( const Basso_Array<BASSO_IDTYPE> &econn, Basso_Array<BASSO_IDTYPE> &sctr ) const;
	
	/**Initializes a scatter vector for the gdof map.  There is no bounds checking to expidite 
	the member function.  The form of the scatter vector is as follows:
		sctr = [ node0_dof0, node0_dof1, node0_dof2 ..., node1_dof0, node1_dof2, ... ] 
	\param econn a pointer to an array of the local node ids for which to compute the scatter vector. Do not use global node ids here.
	\param nne the length of the array in econn
	\param sctr on return has the scatter vector to the global dofs.  It should be of length econn.Length()*NumDofPerNode(). 
	*/
	void SetScatter( const BASSO_IDTYPE *econn, int nne, Basso_Array<BASSO_IDTYPE> &sctr ) const;
	
	
	/**Initializes a scatter vector for the gdof map using only particular dofs.  There is no bounds checking to expidite 
	the member function.  The form of the scatter vector is as follows:
		sctr = [ node0_dof0, node0_dof1, node0_dof2 ..., node1_dof0, node1_dof2, ... ] 
	\param econn a pointer to an array of the local node ids for which to compute the scatter vector. Do not use global node ids here.
	\param nne the length of the array in econn
	\param sctr on return has the scatter vector to the global dofs.  It should be of length econn.Length()*NumDofPerNode(). 
	\param ldofs an array that has the local dofs to be used in forming the scatter vector.  Again there is no error checking here.
	*/
	void SetScatter( const BASSO_IDTYPE *econn, int nne, Basso_Array<BASSO_IDTYPE> &sctr, const Basso_Array<int>  &ldofs ) const;
	
	/**Returns the global dof given a lodal node id and a local dof
	 \param lnid the local node id
	\param ldof the local dof  for the lnid node to return the global dof.
	\return the globaldof corresponding to the ldof dof of node lnid.	 
	 */
	BASSO_IDTYPE GDOF( BASSO_IDTYPE lnid, int ldof ) const 
                { if ( !fixed_) return 0; return gdof_map_[lnid][ldof]; }

	
	/**Sets the dofs in the map 
	\param ldofs on return has the local dofs active in the map
	*/
	void LDOFs( Basso_Array<int> & ldofs ) const;
	
	virtual void Print( std::ostream &out=BASSO_STDOUT ) const;
	
protected:
	
	Basso_Array2D<BASSO_IDTYPE> gdof_map_;
	Basso_Array<BASSO_IDTYPE> gnids_;
	bool fixed_;
};

void Basso_FixedDOFMap::LDOFs( Basso_Array<int> & ldofs ) const
{
	ldofs.Resize( NumDofPerNode() );
	for ( int s=0; s<ldofs.Length(); ++ s )
		ldofs[s] = s;
}

void Basso_FixedDOFMap::Print( std::ostream &out ) const 
{
 	if ( !fixed_ ) 
	{
		out << "Basso_FixedDOFMap, not yet set\n";
		return;
	}
	
	int nn=NumMyNodes();
	if ( nn>1000 ) nn=1000;
	
	out << "Basso_FixedDOFMap " <<
	       "\n      lnid      gnid      ldof      gdof\n"; 
	for ( int i=0; i<nn; ++i )
		for ( int s=0; s<NumDofPerNode(); ++s )
			out << setw(10) << i << setw(10) << gnids_[i] << setw(10) << s << setw(10) << gdof_map_[i][s] << "\n"; 
}

void Basso_FixedDOFMap::SetScatter( const Basso_Array<BASSO_IDTYPE> &econn, Basso_Array<BASSO_IDTYPE> &sctr ) const
{
	SetScatter(econn.Data(), econn.Length(), sctr);		
}


void Basso_FixedDOFMap::SetScatter( const BASSO_IDTYPE *econn, int nne, Basso_Array<BASSO_IDTYPE> &sctr ) const
{
	const BASSO_IDTYPE *nptr = econn;
	BASSO_IDTYPE *cptr = sctr.Data(), I;
	int s;
	for ( I=0; I<nne; ++I, ++nptr )
		for ( s=0; s<NumDofPerNode(); ++s, ++cptr ) 
			*cptr = gdof_map_[ *nptr ][s];		
}

void Basso_FixedDOFMap::SetScatter( const BASSO_IDTYPE *econn, int nne, Basso_Array<BASSO_IDTYPE> &sctr, const Basso_Array<int>  &ldofs ) const
{
	const BASSO_IDTYPE *nptr = econn;
	BASSO_IDTYPE *cptr = sctr.Data(), I;
	int s;
	for ( I=0; I<nne; ++I, ++nptr )
	{	
		const int *sptr=ldofs.Data();
		for ( s=0; s<ldofs.Length(); ++s, ++cptr, ++sptr ) 
			*cptr = gdof_map_[ *nptr ][*sptr];
	}		
}

void Basso_FixedDOFMap::FixDOFMap( )
{
	int n=0;
	for ( int I=0; I<gnids_.Length(); ++I )
		for ( int i=0; i<gdof_map_.M(); ++i, ++n )
			gdof_map_[I][i]=n;
		
	// done! so set as fixed
	fixed_ = true;

}

std::ostream &operator << ( std::ostream &out, const Basso_FixedDOFMap &A )
{
	A.Print( out );
	return out;
}

} // end of namspace
#endif