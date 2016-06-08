/*! \file Basso_Hexa8.h
Basso_Hexa8
 Basso ver 1.0	

\author Jack Chessa, jfchessa@utep.edu
\date Wed November 28 2007

*/

#ifndef _HEXA8_BASIS_BASSO_H_
#define _HEXA8_BASIS_BASSO_H_

#include "Basso_defs.h"
#include "Basso_shape_functions.h"
#include "Basso_ParentElement.h"

namespace Basso {

/**
 \brief Eight node hexahedral basis

  Class to represent a connonical eight node hexahedral parent element

	Parent coordinate space and geometric nomenclature for Basso_Hexa8 
	
	\image html hexa8.jpg
                         
*/

class Basso_Hexa8 : public Basso_ParentElement {

public:
	
	// constructors
	Basso_Hexa8() { }
	
	// destructors
	virtual ~Basso_Hexa8() { }
	
	// member functions
	virtual Basso_ElementType	Type() const { return Basso_HEXA8; }
	virtual Basso_ElementShape Shape() const { return Basso_HEXAHEDRA; }
	virtual int				NumNodes() const { return 8; }     
	virtual int				NumEdges() const { return 12; }      
	virtual int				NumFaces() const { return 6; }      
	virtual int				Dimension() const { return 3; }
	virtual int				Order()	const { return 2; }
	
	virtual void 		ParentCoord( Basso_Array< Basso_Point > &xi ) 	const;
	virtual Basso_Point      	Centroid() const { return Basso_Point( 0.0, 0.0, 0.0 ); }
	virtual void    	NumNodesOnFaces( Basso_iVector &nn_face ) 	const;
	virtual void 		NumNodesOnEdges( Basso_iVector &nn_edge ) 	const;
	virtual void    	NodesOnFace ( int f, Basso_iVector &face_nodeids ) 		const;
	virtual void 		NodesOnEdge ( int e, Basso_iVector &edge_nodeids )			const;
	virtual void    	FaceElementType( Basso_Array<Basso_ElementType> &face_basis ) 	const;
	virtual void 		EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const;
	
	virtual void 		Quadrature( int order, Basso_QuadratureRule &qrule ) const { qrule.SetQuadratureRule( Basso_GAUSS3D_QUAD, order ); }
	
	virtual void 		Na( const Basso_Point &p, Basso_nVector &shapefunct ) const 
	        {   
#ifdef ALLOW_DYNAMIC_RESIZE 
            shapefunct.Newsize( NumNodes() ); 
#endif
            shape_hexa8(shapefunct.Data(), shapefunct.Inc(), p.x(), p.y(), p.z() ); 
        }
        
	virtual void 		DNa( const Basso_Point &p, Basso_nMatrix &gshape ) 	const
	    { 
#ifdef ALLOW_DYNAMIC_RESIZE
            gshape.Newsize( NumNodes(), Dimension() );
#endif
	        dshape_hexa8( gshape.Data(), gshape.LDA(), p.x(), p.y(), p.z() ); 
	    }	
	 
	virtual void    	Na( Basso_nVector &shapefunct )    const { Na( Centroid(), shapefunct ); }
	virtual void 		DNa( Basso_nMatrix &grad_shape ) const { DNa( Centroid(), grad_shape ); }

private:
	
};

void Basso_Hexa8::ParentCoord( Basso_Array< Basso_Point > &xi ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( xi.Length()!=NumNodes() )
		xi.Resize(NumNodes());
#endif
	xi[0]=Basso_Point(-1.0,-1.0,-1.0);
	xi[1]=Basso_Point( 1.0,-1.0,-1.0);
	xi[2]=Basso_Point( 1.0, 1.0,-1.0);
	xi[3]=Basso_Point(-1.0, 1.0,-1.0);
	xi[4]=Basso_Point(-1.0,-1.0, 1.0);
	xi[5]=Basso_Point( 1.0,-1.0, 1.0);
	xi[6]=Basso_Point( 1.0, 1.0, 1.0);
	xi[7]=Basso_Point(-1.0, 1.0, 1.0);
}

void 	Basso_Hexa8::NumNodesOnFaces( Basso_iVector &nn_face ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_face.Length()!=NumFaces() )
		nn_face.Resize( NumFaces() );
#endif
	for ( int i=0; i<NumFaces(); ++i )
		nn_face[i]=4;
}

void 	Basso_Hexa8::NumNodesOnEdges( Basso_iVector &nn_edge ) 	const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( nn_edge.Length()!=NumEdges() )
		nn_edge.Resize( NumEdges() );
#endif
	for ( int i=0; i<NumEdges(); ++i )
		nn_edge[i]=2; 
}

void Basso_Hexa8::NodesOnFace ( int e, Basso_iVector &face_nodeids ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_nodeids.Length()!=0 )
		face_nodeids.Resize(0);
#endif	
	switch (e) {
	
		case 0:
		face_nodeids[0]=0; face_nodeids[1]=1; face_nodeids[2]=5; face_nodeids[3]=4;
		break;
		
		case 1:
		face_nodeids[0]=1; face_nodeids[1]=2; face_nodeids[2]=6; face_nodeids[3]=5;
		break;
		
		case 2:
		face_nodeids[0]=2; face_nodeids[1]=3; face_nodeids[2]=7; face_nodeids[3]=6;
		break;
		
		case 3:
		face_nodeids[0]=0; face_nodeids[1]=4; face_nodeids[2]=7; face_nodeids[3]=3;
		break;
		
		case 4:
		face_nodeids[0]=0; face_nodeids[1]=3; face_nodeids[2]=2; face_nodeids[3]=1;
		break;
		
		case 5:
		face_nodeids[0]=4; face_nodeids[1]=5; face_nodeids[2]=6; face_nodeids[3]=7;
		break;
		
		default:
		Basso_Warning("Basso_Hexa8::NodesOnFace","face id out of range");
		break;
		
	}
}

void Basso_Hexa8::NodesOnEdge ( int e, Basso_iVector &edge_nodeids ) const
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
		edge_nodeids[0]=4; edge_nodeids[1]=5; 
		break;
	
		case 5:
		edge_nodeids[0]=5; edge_nodeids[1]=6; 
		break;
	
		case 6:
		edge_nodeids[0]=6; edge_nodeids[1]=7;  
		break;
	
		case 7:
		edge_nodeids[0]=7; edge_nodeids[1]=4;  
		break;

		case 8:
		edge_nodeids[0]=0; edge_nodeids[1]=4; 
		break;
		
		case 9:
		edge_nodeids[0]=1; edge_nodeids[1]=5; 
		break;
	
		case 10:
		edge_nodeids[0]=2; edge_nodeids[1]=6; 
		break;
	
		case 11:
		edge_nodeids[0]=3; edge_nodeids[1]=7;  
		break;
		
		default:
		Basso_Warning("Basso_Hexa8::NodesOnEdge","edge id out of range");
		break;
		
	}
}

void Basso_Hexa8::FaceElementType( Basso_Array< Basso_ElementType > &face_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( face_basis.Length()!=NumFaces() )
		face_basis.Resize(NumFaces());
#endif
	for ( int i=0; i<NumFaces(); ++i )
		face_basis[i]=Basso_QUAD4;
}

void Basso_Hexa8::EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const
{
#ifdef ALLOW_DYNAMIC_RESIZE 
	if ( edge_basis.Length()!=NumEdges() )
		edge_basis.Resize(NumEdges());
#endif
	for ( int i=0; i<NumEdges(); ++i )
		edge_basis[i]=Basso_LINE2; 
}


} // end namespace

#endif






