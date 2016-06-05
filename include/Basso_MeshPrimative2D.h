#ifndef _BASSO_MESH_PRIMATIVE_2D_H_
#define _BASSO_MESH_PRIMATIVE_2D_H_

#include "Basso_MeshPrimative.h"

namespace Basso
{

// ****************************************************************************************** //
/**
\brief VBC for esh primatives in 2d.
*/
template<class T>
class Basso_MeshPrimative2D : public Basso_MeshPrimative<T>
{
public:
	typedef T value_type;

    /**void constructor*/
	Basso_MeshPrimative2D() : Basso_MeshPrimative<T>() {}
	
	/**
	Sets the connectivity of line segments on a edge.  These are typically LINE2
	elements.  This member function will dynamically resize the data arrays passed
	to the member function as needed.
	\param conn - on return has the connectivity of the edge segments.
	\param id - an id the idenifies what edge that we wish to get the edge segments of.  
	This is typically an enumerated type in the derived class.
	*/
	virtual void EdgeSegs( Basso_Array2D<BASSO_IDTYPE> &conn, int id ) const = 0;
	
	/**
	Returns a set of node ids on an edge.
	\param nids- on return has the set on nodes ids on that edge
	\param id - an id the idenifies what edge that we wish to get the edge segments of.  
	This is typically an enumerated type in the derived class.
	*/
	virtual void EdgeNodes( set<BASSO_IDTYPE> &nids, int id ) const = 0;
	
	/**
	Returns a set of node ids on an edge.
	\param nids- on return has the set on nodes ids on that edge
	\param id - an id the idenifies what edge that we wish to get the edge segments of.  
	This is typically an enumerated type in the derived class.
	*/
	virtual void EdgeNodes( Basso_Array<BASSO_IDTYPE> &nids, int id ) const
	{
        set<BASSO_IDTYPE> idset;
        EdgeNodes(idset,id);
        nids.Newsize( idset.size() );
        set2array( idset, nids.Data() );
	}
	
	/** Returns the spacial dimension (2) */
    virtual int SpacialDimension() const { return 2; }
    
protected:
	

};


}// end of namespace
#endif

