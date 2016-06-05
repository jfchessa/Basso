/*! \file Basso_Quad8.h
Basso_Quad8
 Basso ver 1.0	

\author Jack Chessa, jfchessa@utep.edu
\date Wed November 28 2007

*/

#ifndef _QUAD8_BASSO_H_
#define _QUAD8_BASSO_H_

#include "Basso_defs.h"
#include "Basso_shape_functions.h"
#include "Basso_ParentElement.h"

namespace Basso {

/**
 \brief Eight node quadrilateral parent element

  Class to represent a connonical eight node quadrilateral parent element

	Parent coordinate space and geometric nomenclature for Basso_Quad8 


\image html quad8.pdf
                         
*/

class Basso_Quad8 : public Basso_ParentElement {

public:
	
	// constructors
	Basso_Quad8() { }
	
	// destructors
	virtual ~Basso_Quad8() { }
	
	// member functions
	virtual Basso_ElementType	Type() const { return Basso_QUAD8; }
	virtual Basso_ElementShape Shape() const { return Basso_QUADRILATERAL; }
	virtual int				NumNodes() const { return 8; }     
	virtual int				NumEdges() const { return 4; }      
	virtual int				NumFaces() const { return 0; }      
	virtual int				Dimension() const { return 2; }
	virtual int				Order()	const { return 3; }
	
	virtual void 		ParentCoord( Basso_Array< Basso_Point > &xi ) 	const;
	virtual Basso_Point      	Centroid() const { return Basso_Point( 0.0, 0.0 ); }
	virtual void    	NumNodesOnFaces( Basso_iVector &nn_face ) 	const;
	virtual void 		NumNodesOnEdges( Basso_iVector &nn_edge ) 	const;
	virtual void    	NodesOnFace ( int f, Basso_iVector &face_nodeids ) 		const;
	virtual void 		NodesOnEdge ( int e, Basso_iVector &edge_nodeids )			const;
	virtual void    	FaceElementType( Basso_Array<Basso_ElementType> &face_basis ) 	const;
	virtual void 		EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const;
	
	virtual void 		Quadrature( int order, Basso_QuadratureRule &qrule ) const { qrule.SetQuadratureRule( Basso_GAUSS2D_QUAD, order ); }
	
	virtual void 	Na( const Basso_Point &p, Basso_nVector &shapefunct ) 	const 
        {   
#ifdef ALLOW_DYNAMIC_RESIZE 
            shapefunct.Newsize( NumNodes() ); 
#endif
            shape_quad8( shapefunct.Data(), shapefunct.Inc(), p.x(), p.y() ); 
        }
        
	virtual void 	DNa( const Basso_Point &p, Basso_nMatrix &grad_shape ) 	const 
	    { 
#ifdef ALLOW_DYNAMIC_RESIZE
            grad_shape.Newsize( NumNodes(), Dimension() );
#endif
	        dshape_quad8( grad_shape.Data(), grad_shape.LDA(), p.x(), p.y() ); 
	    }
	
	virtual void    	Na( Basso_nVector &shapefunct )    const { Na( Centroid(), shapefunct ); }
	virtual void 		DNa( Basso_nMatrix &grad_shape ) const { DNa( Centroid(), grad_shape ); }

private:
	
};

void Basso_Quad8::ParentCoord( Basso_Array< Basso_Point > &xi ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( xi.Length()!=NumNodes() )
		xi.Resize(NumNodes());
#endif
	xi[0]=Basso_Point(-1.0,-1.0);
	xi[1]=Basso_Point( 1.0,-1.0);
	xi[2]=Basso_Point( 1.0, 1.0);
	xi[3]=Basso_Point(-1.0, 1.0);
	xi[4]=Basso_Point( 0.0,-1.0);
	xi[5]=Basso_Point( 1.0, 0.0);
	xi[6]=Basso_Point( 0.0, 1.0);
	xi[7]=Basso_Point(-1.0, 0.0);
}

void 	Basso_Quad8::NumNodesOnFaces( Basso_iVector &nn_face ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_face.Length()!=NumEdges() )
		nn_face.Resize( NumEdges() );
#endif
}

void 	Basso_Quad8::NumNodesOnEdges( Basso_iVector &nn_edge ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_edge.Length()!=NumEdges() )
		nn_edge.Resize( NumEdges() );
#endif
	for ( int i=0; i<NumEdges(); ++i )
		nn_edge[i]=3; 
}

void Basso_Quad8::NodesOnFace ( int e, Basso_iVector &face_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_nodeids.Length()!=0 )
		face_nodeids.Resize(0);
#endif
}

void Basso_Quad8::NodesOnEdge ( int e, Basso_iVector &edge_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_nodeids.Length()!=3 )
		edge_nodeids.Resize(3);
#endif
	
	switch (e) {
	
		case 0:
		edge_nodeids[0]=0; edge_nodeids[1]=1; edge_nodeids[2]=4;
		break;
		
		case 1:
		edge_nodeids[0]=1; edge_nodeids[1]=2; edge_nodeids[2]=5;
		break;
		
		case 2:
		edge_nodeids[0]=2; edge_nodeids[1]=3; edge_nodeids[2]=6;
		break;
		
		case 3:
		edge_nodeids[0]=3; edge_nodeids[1]=4; edge_nodeids[2]=7;
		break;
		
		default:
		Basso_Warning("Basso_Quad8::NodesOnFace","face id out of range");
		break;
		
	}
}

void Basso_Quad8::FaceElementType( Basso_Array< Basso_ElementType > &face_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_basis.Length()!=NumFaces() )
		face_basis.Resize(NumFaces());
#endif
}

void Basso_Quad8::EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_basis.Length()!=NumEdges() )
		edge_basis.Resize(NumEdges());
#endif
	for ( int i=0; i<NumEdges(); ++i )
		edge_basis[i]=Basso_LINE3; 
}


} // end namespace

#endif




