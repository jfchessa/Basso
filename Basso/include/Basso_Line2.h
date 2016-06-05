/*! \file Basso_Line2.h

 Basso ver 1.0	

\author Jack Chessa, jfchessa@utep.edu
\date Wed Apr 6 2007

*/

#ifndef _BASSO_LINE2_H_
#define _BASSO_LINE2_H_

#include "Basso_defs.h"
#include "Basso_shape_functions.h"
#include "Basso_ParentElement.h"

namespace Basso {

/**
 \brief Two node 1D parent element class

  Class to represent a connonical two node 1D parent element

	Parent coordinate space and geometric nomenclature for Basso_Line2 

\image html line2.pdf
                         
*/

class Basso_Line2 : public Basso_ParentElement {

public:
	
	// constructors
	Basso_Line2() { }
	
	// destructors
	virtual ~Basso_Line2() { }
	
	// member functions
	virtual Basso_ElementType	Type() const { return Basso_LINE2; }
	virtual Basso_ElementShape 	Shape() const { return Basso_LINE; }
	virtual int		NumNodes() const { return 2; }     
	virtual int		NumEdges() const { return 0; }      
	virtual int		NumFaces() const { return 0; }      
	virtual int		Dimension() const { return 1; }
	virtual int		Order()	const { return 1; }
	
	virtual void 		ParentCoord( Basso_Array< Basso_Point > &xi ) 	const;
	virtual Basso_Point Centroid() const { return Basso_Point( 0.0 ); }
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
            shape_line2( shapefunct.Data(), shapefunct.Inc(), p.x() ); 
        }
        
	virtual void 	DNa( const Basso_Point &p, Basso_nMatrix &grad_shape ) 	const 
	    { 
#ifdef ALLOW_DYNAMIC_RESIZE
            grad_shape.Newsize( NumNodes(), Dimension() );
#endif
	        dshape_line2( grad_shape.Data(), grad_shape.LDA(), p.x() ); 
	    }

	virtual void    Na( Basso_nVector &shapefunct )	const { Na( Centroid(), shapefunct ); }
	virtual void 	DNa( Basso_nMatrix &grad_shape ) const { DNa( this->Centroid(), grad_shape ); }

private:
	
};

void Basso_Line2::ParentCoord( Basso_Array< Basso_Point > &xi ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( xi.Length()!=NumNodes() )
		xi.Resize(NumNodes());
#endif
	xi[0]=Basso_Point(-1.0);
	xi[1]=Basso_Point(1.0);
}

void 	Basso_Line2::NumNodesOnFaces( Basso_iVector &nn_face ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_face.Length()!=NumEdges() )
		nn_face.Resize( NumEdges() );
#endif
}

void 	Basso_Line2::NumNodesOnEdges( Basso_iVector &nn_edge ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_edge.Length()!=NumEdges() )
		nn_edge.Resize( NumEdges() );
#endif 
}

void Basso_Line2::NodesOnFace ( int e, Basso_iVector &face_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_nodeids.Length()!=0 )
		face_nodeids.Resize(0);
#endif
}

void Basso_Line2::NodesOnEdge ( int e, Basso_iVector &edge_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_nodeids.Length()!=0 )
		edge_nodeids.Resize(0);
#endif
}

void Basso_Line2::FaceElementType( Basso_Array< Basso_ElementType > &face_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_basis.Length()!=NumFaces() )
		face_basis.Resize(NumFaces());
#endif
}

void Basso_Line2::EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_basis.Length()!=NumEdges() )
		edge_basis.Resize(NumEdges());
#endif
}


/*
 void Basso_Line2::Na( const Basso_Point &p, nArray &shapefunct ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( shapefunct.Length()!=NumNodes() )
		shapefunct.Resize(NumNodes());	
#endif
	shapefunct[0]=1-p.x();
	shapefunct[1]=p.x();
}

void Basso_Line2::DNa( const Basso_Point &p, nMatrix &grad_shape ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( mat_nrows(grad_shape)!=NumNodes() || mat_ncols(grad_shape)!=Dimension() )
		Resize( grad_shape, NumNodes(), Dimension() );	
#endif
	grad_shape(0,0)=-1.0;  
	grad_shape(1,0)= 1.0; 
}
*/

} // end namespace

#endif




