#ifndef _BASSO_MESH_PRIMATIVE_3D_H_
#define _BASSO_MESH_PRIMATIVE_3D_H_

#include "Basso_MeshPrimative.h"

namespace Basso
{
	

// ****************************************************************************************** //
template<class T>
class Basso_MeshPrimative3D : public Basso_MeshPrimative<T>
{
public:
	typedef T value_type;

	Basso_MeshPrimative3D() : Basso_MeshPrimative<T>() {}
	
	virtual void SideSegs( Basso_Array2D<BASSO_IDTYPE> &conn, int id ) const = 0;
	virtual void SideNodes( set<BASSO_IDTYPE> &nids, int id ) const = 0;
	virtual void SideNodes( Basso_Array<BASSO_IDTYPE> &nids, int id ) const
	{
        set<BASSO_IDTYPE> idset;
        SideNodes(idset,id);
        nids.Newsize( idset.size() );
        set2array( idset, nids.Data() );
	}
	
	virtual void EdgeSegs( Basso_Array2D<BASSO_IDTYPE> &conn, int id ) const = 0;
	virtual void EdgeNodes( set<BASSO_IDTYPE> &nids, int id ) const = 0;
	virtual void EdgeNodes( Basso_Array<BASSO_IDTYPE> &nids, int id ) const
	{
        set<BASSO_IDTYPE> idset;
        EdgeNodes(idset,id);
        nids.Newsize( idset.size() );
        set2array( idset, nids.Data() );
	}
	
    virtual int SpacialDimension() const { return 3; }
    
protected:
	

};
	
	
	
}// end of namespace
#endif

