/*! \file Basso_Pyramid5.h

 Basso ver 1.0	

\author Jack Chessa, jfchessa@utep.edu
\date Wed Apr 6 2007

*/

#ifndef _BASSO_PYRAMID5_H_
#define _BASSO_PYRAMID5_H_

#include "Basso_defs.h"
#include "Basso_shape_functions.h"
#include "Basso_ParentElement.h"

namespace Basso {

/**
 \brief Five node pyramid element

  Class to represent a connonical five node pyramid element

	Parent coordinate space and geometric nomenclature for Basso_Pyramid5 


\image html pyramid5.jpg
                         
*/

class Basso_Pyramid5 : public Basso_ParentElement {

public:
	
	// constructors
	Basso_Pyramid5() { }
	
	// destructors
	virtual ~Basso_Pyramid5() { }
	
	// member functions
	virtual Basso_ElementType	Type() const { return Basso_PYRAMID5; }
	virtual Basso_ElementShape Shape() const { return Basso_PYRAMID; }
	virtual int				NumNodes() const { return 5; }     
	virtual int				NumEdges() const { return 8; }      
	virtual int				NumFaces() const { return 5; }      
	virtual int				Dimension() const { return 3; }
	virtual int				Order()	const { return 1; }
	
	virtual void 		ParentCoord( Basso_Array< Basso_Point > &xi ) 	const;
	virtual Basso_Point     	Centroid() const { return Basso_Point( 0, 0, 0 ); }
	virtual void    	NumNodesOnFaces( Basso_iVector &nn_face ) 	const;
	virtual void 		NumNodesOnEdges( Basso_iVector &nn_edge ) 	const;
	virtual void    	NodesOnFace ( int f, Basso_iVector &face_nodeids ) 		const;
	virtual void 		NodesOnEdge ( int e, Basso_iVector &edge_nodeids )			const;
	virtual void    	FaceElementType( Basso_Array<Basso_ElementType> &face_basis ) 	const;
	virtual void 		EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const;
	
	virtual void 		Quadrature( int order, Basso_QuadratureRule &qrule ) const 
		{ //qrule.SetQuadratureRule( Basso_TETRA_QUAD, order ); 
			Basso_Warning("Basso_Pyramid5::Quadrature","not yet defined"); 
		}
	
	virtual void 		Na( const Basso_Point &p, Basso_nVector &shapefunct ) const;
        
	virtual void 		DNa( const Basso_Point &p, Basso_nMatrix &gshape ) 	const;
	    	
	virtual void    Na( Basso_nVector &shapefunct )	const { Na( Centroid(), shapefunct ); }
	virtual void 	DNa( Basso_nMatrix &grad_shape ) const { DNa( Centroid(), grad_shape ); }

private:
	
};

void Basso_Pyramid5::ParentCoord( Basso_Array< Basso_Point > &xi ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( xi.Length()!=NumNodes() )
		xi.Resize(NumNodes());
#endif
	xi[0]=Basso_Point(-1.0,-1.0,-1.0);
	xi[1]=Basso_Point( 1.0,-1.0,-1.0);
	xi[2]=Basso_Point( 1.0, 1.0,-1.0);
	xi[3]=Basso_Point(-1.0, 1.0,-1.0);
	xi[3]=Basso_Point( 0.0, 0.0, 1.0);
}

void 	Basso_Pyramid5::NumNodesOnFaces( Basso_iVector &nn_face ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_face.Length()!=NumEdges() )
		nn_face.Resize( NumEdges() );
#endif
	nn_face[0]=4;
	for ( int i=1; i<NumFaces(); ++i )
		nn_face[i]=3;
}

void 	Basso_Pyramid5::NumNodesOnEdges( Basso_iVector &nn_edge ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_edge.Length()!=NumEdges() )
		nn_edge.Resize( NumEdges() );
#endif
	for ( int i=0; i<NumEdges(); ++i )
		nn_edge[i]=2; 
}

void Basso_Pyramid5::NodesOnFace ( int f, Basso_iVector &face_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_nodeids.Length()!=3 )
		face_nodeids.Resize(3);
#endif
	switch (f) {
	
		case 0:
		face_nodeids[0]=0; face_nodeids[1]=1; face_nodeids[2]=2; face_nodeids[2]=3;
		break;
		
		case 1:
		face_nodeids[0]=0; face_nodeids[1]=1; face_nodeids[2]=4;
		break;
		
		case 2:
		face_nodeids[0]=1; face_nodeids[1]=2; face_nodeids[2]=4;
		break;
		
		case 3:
		face_nodeids[0]=2; face_nodeids[1]=3; face_nodeids[2]=4;
		break;
		
		case 4:
		face_nodeids[0]=3; face_nodeids[1]=0; face_nodeids[2]=4;
		break;
		
		default:
		Basso_Warning("Basso_Pyramid5::NodesOnFace","face id out of range");
		break;
		
	}
}

void Basso_Pyramid5::NodesOnEdge ( int e, Basso_iVector &edge_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_nodeids.Length()!=2 )
		edge_nodeids.Resize(2);
#endif
	
	switch (e) {
	
		case 0:
		edge_nodeids[0]=0; edge_nodeids[1]=1;
		break;
		
		case 1:
		edge_nodeids[0]=1; edge_nodeids[1]=2;
		break;
		
		case 2:
		edge_nodeids[0]=2; edge_nodeids[1]=3;
		break;
		
		case 3:
		edge_nodeids[0]=3; edge_nodeids[1]=0;
		break;
		
		case 4:
		edge_nodeids[0]=0; edge_nodeids[1]=4;
		break;
		
		case 5:
		edge_nodeids[0]=1; edge_nodeids[1]=4;
		break;
		
		case 6:
		edge_nodeids[0]=2; edge_nodeids[1]=4;
		break;
		
		case 7:
		edge_nodeids[0]=3; edge_nodeids[1]=4;
		break;
		
		default:
		Basso_Warning("Basso_Pyramid5::NodesOnFace","face id out of range");
		break;
		
	}
}

void Basso_Pyramid5::FaceElementType( Basso_Array< Basso_ElementType > &face_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_basis.Length()!=NumFaces() )
		face_basis.Resize(NumFaces());
#endif
    face_basis[0]=Basso_QUAD4;
	for ( int i=1; i<NumEdges(); ++i )
		face_basis[i]=Basso_TRIA3;
}

void Basso_Pyramid5::EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_basis.Length()!=NumEdges() )
		edge_basis.Resize(NumEdges());
#endif
	for ( int i=0; i<NumEdges(); ++i )
		edge_basis[i]=Basso_LINE2; 
}


 void Basso_Pyramid5::Na( const Basso_Point &p, Basso_nVector &shapefunct ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( shapefunct.Length()!=NumNodes() )
		shapefunct.Resize(NumNodes());	
#endif

	shapefunct[0]=0.125*( 1 - 2*p.x() - 2*p.y() - p.z() - 2*p.x()*p.y()*p.z() );
	shapefunct[1]=0.125*( 1 + 2*p.x() - 2*p.y() - p.z() + 2*p.x()*p.y()*p.z() );
	shapefunct[2]=0.125*( 1 + 2*p.x() + 2*p.y() - p.z() - 2*p.x()*p.y()*p.z() );
	shapefunct[3]=0.125*( 1 - 2*p.x() + 2*p.y() - p.z() + 2*p.x()*p.y()*p.z() );
	shapefunct[4]=0.5*(1+p.z());
}

void Basso_Pyramid5::DNa( const Basso_Point &p, Basso_nMatrix &grad_shape ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( grad_shape.NumRows()!=NumNodes() || grad_shape.NumCols()!=Dimension() )
		grad_shape.Resize( NumNodes(), Dimension() );	
#endif
	Basso_Numeric r=p.x(), s=p.y(), t=p.z();
	grad_shape(0,0)=-0.25*(1+s*t); grad_shape(0,1)=-0.25*(1+r*t);  grad_shape(0,2)=-0.125*(1+2*r*s);
	grad_shape(1,0)= 0.25*(1+s*t); grad_shape(1,1)= 0.25*(-1+r*t); grad_shape(1,2)= 0.125*(-1+2*r*s); 
	grad_shape(2,0)= 0.25*(1-s*t); grad_shape(2,1)= 0.25*(1-r*t);  grad_shape(2,2)=-0.125*(1+2*r*s);  
	grad_shape(3,0)=0.25*(-1+s*t); grad_shape(3,1)= 0.25*(1+r*t);  grad_shape(3,2)= 0.125*(-1+2*r*s); 
	grad_shape(4,0)= 0.00;         grad_shape(4,1)= 0.00;          grad_shape(4,2)= 0.5; 
}

} // end namespace

#endif




