#ifndef _BASSO_MESH_PRIMATIVE_H_
#define _BASSO_MESH_PRIMATIVE_H_

#include "Basso_Array2D.h"
#include "Basso_Array.h"
#include "Basso_meshing.h"
#include "Basso_ElementType.h"
#include "Basso_stlops.h"

//#define BASSO_IDTYPE unsigned long int


namespace Basso
{

// ****************************************************************************************** //

/**
\brief VBC for meshing simple geometries

The derived classes allow for the generation of various finite element meshes on some
simple geometries like quadralaterals, blocks, toroids, etc. In general the class does not
store the mesh data structures explicitly, but rather provides member functions that 
set them in external data arrays.

This class is templated so that it works with doubles and floats or any other reasonable 
numeric class.
*/
template<class T>
class Basso_MeshPrimative
{
public:
	typedef T value_type;

    /**void constructor*/
	Basso_MeshPrimative() {}

    /**
    Sets the node coordinate and element connectivity matrices.  This member function allows
    for dynamic resizing of the arrays passed to it.
    \param node - The node coordinate array.  The node coords are in column format.
    \param conn - the element connectivity array.  The elements are stored in column format.
    */
	virtual void SetMesh( Basso_Array2D<T> &node, Basso_Array2D<BASSO_IDTYPE> &conn ) const = 0;
	
	/**Returns the type of element used in the domain.*/
	virtual Basso_ElementType MeshElementType() const = 0;
	
	/**Returns the base element used in the mesh.*/
	//virtual void MeshBaseElement( const Basso_ParentElement *ebase ) const { }
	
	/**Returns the spacial dimension of the problem.*/
    virtual int SpacialDimension() const = 0;
	
	/**Returns the spacial dimension of the problem.*/
	virtual int SDIM() const { return SpacialDimension(); }
	
protected:
	

};


}// end of namespace
#endif

