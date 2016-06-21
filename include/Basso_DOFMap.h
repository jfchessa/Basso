/*! \file Basso_DOFMap.h

	\brief Defines the class Basso_DOFMap

	
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

namespace Basso
{

	
/**
 \brief DOF map with a fixed number of active dofs per node. Also assumes that the node numbering is 
 zero offset and continous.
	
*/
class Basso_DOFMap
{
public:
	/**
	Constructor
	\param gnids - an array of the global node ids that are supported by the processor.  This array needs to 
	map onto the rows of the local node coordinate matrix
	\param nldof the number of dofs active on each node.  Default value is 6.
	*/
	Basso_DOFMap( int nn=0, int nldof=6 ) : ndofpn_(nldof), numnode_(nn)  { }
	
	/** Deconstructor */
	~Basso_DOFMap() {  }
	
	
	/**Returns the number of active dofs for each node.*/
	int NumDofPerNode() const { return ndofpn_; }

	
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
	BASSO_IDTYPE GDOF( BASSO_IDTYPE lnid, int ldof ) const { return ndofpn_*lnid + ldof; }

	
	virtual void Print( std::ostream &out=BASSO_STDOUT ) const;
	
protected:
	
	int ndofpn_, numnode_;
};

void Basso_DOFMap::Print( std::ostream &out ) const 
{
	out << "Basso_DOFMap "; 
}

void Basso_DOFMap::SetScatter( const Basso_Array<BASSO_IDTYPE> &econn, Basso_Array<BASSO_IDTYPE> &sctr ) const
{
	SetScatter(econn.Data(), econn.Length(), sctr);		
}


void Basso_DOFMap::SetScatter( const BASSO_IDTYPE *econn, int nne, Basso_Array<BASSO_IDTYPE> &sctr ) const
{
	const BASSO_IDTYPE *nptr = econn;
	BASSO_IDTYPE *cptr = sctr.Data(), I;
	int s;
	for ( I=0; I<nne; ++I, ++nptr )
		for ( s=0; s<NumDofPerNode(); ++s, ++cptr ) 
			*cptr = ndofpn_*(*nptr) + s;		
}

void Basso_DOFMap::SetScatter( const BASSO_IDTYPE *econn, int nne, Basso_Array<BASSO_IDTYPE> &sctr, const Basso_Array<int>  &ldofs ) const
{
	const BASSO_IDTYPE *nptr = econn;
	BASSO_IDTYPE *cptr = sctr.Data(), I;
	int s;
	for ( I=0; I<nne; ++I, ++nptr )
	{	
		const int *sptr=ldofs.Data();
		for ( s=0; s<ldofs.Length(); ++s, ++cptr, ++sptr ) 
			*cptr = ndofpn_*(*nptr) + *sptr;	
	}		
}


std::ostream &operator << ( std::ostream &out, const Basso_DOFMap &A )
{
	A.Print( out );
	return out;
}

} // end of namspace
#endif