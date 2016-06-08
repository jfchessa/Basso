/*! \file Basso_Line3.h

 Basso ver 1.0	

\author Jack Chessa, jfchessa@utep.edu
\date Wed November 28, 2007

*/

#ifndef _BASSO_LINE3_H_
#define _BASSO_LINE3_H_

#include "Basso_defs.h"
#include "Basso_shape_functions.h"
#include "Basso_ParentElement.h"

namespace Basso {

/**
 \brief Three node 1D parent element

  Class to represent a connonical three node 1D parent element

	Parent coordinate space and geometric nomenclature for Basso_Line3 

\image html line3.jpg
                         
*/

class Basso_Line3 : public Basso_ParentElement {

public:
	
	// constructors
	Basso_Line3() { }
	
	// destructors
	virtual ~Basso_Line3() { }
	
	// member functions
	virtual Basso_ElementType	Type() const { return Basso_LINE3; }
	virtual Basso_ElementShape 	Shape() const { return Basso_LINE; }
	virtual int			NumNodes() const { return 3; }     
	virtual int			NumEdges() const { return 0; }      
	virtual int			NumFaces() const { return 0; }      
	virtual int			Dimension() const { return 1; }
	virtual int			Order()	const { return 2; }
	
	virtual void 		ParentCoord( Basso_Array< Basso_Point > &xi ) 	const;
	virtual Basso_Point   	Centroid() const { return Basso_Point( 0.0 ); }
	virtual void    	NumNodesOnFaces( Basso_iVector &nn_face ) 	const;
	virtual void 		NumNodesOnEdges( Basso_iVector &nn_edge ) 	const;
	virtual void    	NodesOnFace ( int f, Basso_iVector &face_nodeids ) 		const;
	virtual void 		NodesOnEdge ( int e, Basso_iVector &edge_nodeids )			const;
	virtual void    	FaceElementType( Basso_Array<Basso_ElementType> &face_basis ) 	const;
	virtual void 		EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const;
	
	virtual void 		Quadrature( int order, Basso_QuadratureRule &qrule ) const { qrule.SetQuadratureRule( Basso_GAUSS1D_QUAD, order ); }
	
	virtual void 	Na( const Basso_Point &p, Basso_nVector &shapefunct ) 	const 
	    {   
#ifdef ALLOW_DYNAMIC_RESIZE 
            shapefunct.Newsize( NumNodes() ); 
#endif
            shape_line3( shapefunct.Data(), shapefunct.Inc(), p.x() ); 
        }
        
	virtual void 	DNa( const Basso_Point &p, Basso_nMatrix &grad_shape ) 	const 
	    { 
#ifdef ALLOW_DYNAMIC_RESIZE
            grad_shape.Newsize( NumNodes(), Dimension() );
#endif
	        dshape_line3( grad_shape.Data(), grad_shape.LDA(), p.x() ); 
	    }

	virtual void    	Na( Basso_nVector &shapefunct )    const { Na( Centroid(), shapefunct ); }
	virtual void 		DNa( Basso_nMatrix &grad_shape ) const { DNa( Centroid(), grad_shape ); }

private:
	
};

void Basso_Line3::ParentCoord( Basso_Array< Basso_Point > &xi ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( xi.Length()!=NumNodes() )
		xi.Resize(NumNodes());
#endif
	xi[0]=Basso_Point(-1.0);
	xi[1]=Basso_Point(1.0);
}

void 	Basso_Line3::NumNodesOnFaces( Basso_iVector &nn_face ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_face.Length()!=NumEdges() )
		nn_face.Resize( NumEdges() );
#endif
}

void 	Basso_Line3::NumNodesOnEdges( Basso_iVector &nn_edge ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_edge.Length()!=NumEdges() )
		nn_edge.Resize( NumEdges() );
#endif 
}

void Basso_Line3::NodesOnFace ( int e, Basso_iVector &face_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_nodeids.Length()!=0 )
		face_nodeids.Resize(0);
#endif
}

void Basso_Line3::NodesOnEdge ( int e, Basso_iVector &edge_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_nodeids.Length()!=0 )
		edge_nodeids.Resize(0);
#endif
}

void Basso_Line3::FaceElementType( Basso_Array< Basso_ElementType > &face_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_basis.Length()!=NumFaces() )
		face_basis.Resize(NumFaces());
#endif
}

void Basso_Line3::EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_basis.Length()!=NumEdges() )
		edge_basis.Resize(NumEdges());
#endif
}


} // end namespace

#endif




