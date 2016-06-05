/*! \file Basso_Point1.h

 Basso ver 1.0	

\author Jack Chessa, jfchessa@utep.edu
\date Wed Apr 6 2007

*/

#ifndef _BASSO_POINT1_H_
#define _BASSO_POINT1_H_

#include "Basso_defs.h"
#include "Basso_shape_functions.h"
#include "Basso_ParentElement.h"

namespace Basso {

/**
 \brief Single node 1D parent element class

  Class to represent a connonical one node 1D parent element

	Parent coordinate space and geometric nomenclature for Basso_Point1 

\image html line2.pdf
                         
*/

class Basso_Point1 : public Basso_ParentElement {

public:
	
	// constructors
	Basso_Point1() { }
	
	// destructors
	virtual ~Basso_Point1() { }
	
	// member functions
	virtual Basso_ElementType		Type() const { return Basso_POINT1; }
	virtual Basso_ElementShape 	Shape() const { return Basso_POINT; }
	virtual int		NumNodes() const { return 1; }     
	virtual int		NumEdges() const { return 0; }      
	virtual int		NumFaces() const { return 0; }      
	virtual int		Dimension() const { return 0; }
	virtual int		Order()	const { return 0; }
	
	virtual void 		ParentCoord( Basso_Array< Basso_Point > &xi ) 	const;
	virtual Basso_Point  	Centroid() const { return Basso_Point( 0.0 ); }
	virtual void     	NumNodesOnFaces( Basso_iVector &nn_face ) 	const;
	virtual void 		NumNodesOnEdges( Basso_iVector &nn_edge ) 	const;
	virtual void    	NodesOnFace ( int f, Basso_iVector &face_nodeids ) 		const;
	virtual void 		NodesOnEdge ( int e, Basso_iVector &edge_nodeids )			const;
	virtual void    	FaceElementType( Basso_Array<Basso_ElementType> &face_basis ) 	const;
	virtual void 		EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const;
	
	virtual void 		Quadrature( int order, Basso_QuadratureRule &qrule ) const { qrule.SetQuadratureRule( Basso_GAUSS1D_QUAD, 1 ); }
	
	virtual void 		Na( const Basso_Point &p, Basso_nVector &shapefunct ) const 
		{ 
			if ( p.x()==0.0 )
				shapefunct[0]=1.0;
			else
				shapefunct[0]=0.0;
		}
		
	virtual void 		DNa( const Basso_Point &p, Basso_nMatrix &grad_shape ) 	const 
		{ 
			if ( p.x()==0.0 )
				grad_shape[0][0]=0.0; //BASSO_BIG_NUMBER;
			else
				grad_shape[0][0]=0.0;
		}
	virtual void    Na( Basso_nVector &shapefunct ) 	const { Na( Centroid(), shapefunct ); }
	virtual void 	DNa( Basso_nMatrix &grad_shape ) const { DNa( Centroid(), grad_shape ); }

private:
	
};

void Basso_Point1::ParentCoord( Basso_Array< Basso_Point > &xi ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( xi.Length()!=NumNodes() )
		xi.Resize(NumNodes());
#endif
	xi[0]=Basso_Point(-1.0);
	xi[1]=Basso_Point(1.0);
}

void 	Basso_Point1::NumNodesOnFaces( Basso_iVector &nn_face ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_face.Length()!=NumEdges() )
		nn_face.Resize( NumEdges() );
#endif
}

void 	Basso_Point1::NumNodesOnEdges( Basso_iVector &nn_edge ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_edge.Length()!=NumEdges() )
		nn_edge.Resize( NumEdges() );
#endif 
}

void Basso_Point1::NodesOnFace ( int e, Basso_iVector &face_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_nodeids.Length()!=0 )
		face_nodeids.Resize(0);
#endif
}

void Basso_Point1::NodesOnEdge ( int e, Basso_iVector &edge_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_nodeids.Length()!=0 )
		edge_nodeids.Resize(0);
#endif
}

void Basso_Point1::FaceElementType( Basso_Array< Basso_ElementType > &face_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_basis.Length()!=NumFaces() )
		face_basis.Resize(NumFaces());
#endif
}

void Basso_Point1::EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_basis.Length()!=NumEdges() )
		edge_basis.Resize(NumEdges());
#endif
}

/*
 void Basso_Point1::Na( const Basso_Point &p, Basso_Vector &shapefunct ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( shapefunct.Length()!=NumNodes() )
		shapefunct.Resize(NumNodes());	
#endif
	shapefunct[0]=1-p.x();
	shapefunct[1]=p.x();
}

void Basso_Point1::DNa( const Basso_Point &p, Basso_nMatrix &grad_shape ) const
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




