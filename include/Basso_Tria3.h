/*! \file Basso_Tria3.h 

\brief Three node triangle element

  Basso ver 1.0
 
\author Jack Chessa, jfchessa@utep.edu
\date Friday, September 26, 2008

*/


#ifndef _BASSO_TRIA3_H_
#define _BASSO_TRIA3_H_

#include "Basso_defs.h"
#include "Basso_shape_functions.h"
#include "Basso_ParentElement.h"

namespace Basso
{

	/** 
		\brief Basso_Tria3 parent element class
		
		A dervied class to abstract a three node triangular element
		
\image html tria3.pdf		
	
	*/
	class Basso_Tria3 : public Basso_ParentElement
	{
		public:
			
			// constructors
			Basso_Tria3() {}
			
			// destructors
			virtual ~Basso_Tria3() {}

			// member functions
			virtual Basso_ElementType	Type() const { return Basso_TRIA3; }
			virtual Basso_ElementShape Shape() const { return Basso_TRIANGLE; }
			virtual int				NumNodes() const { return 3; }     
			virtual int				NumEdges() const { return 3; }      
			virtual int				NumFaces() const { return 0; }      
			virtual int				Dimension() const { return 2; }
			virtual int				Order()	const { return 1; }

			virtual void 			ParentCoord( Basso_Array< Basso_Point > &xi ) 	const;
			virtual Basso_Point   	Centroid() const { return Basso_Point( 0.33333333333333333, 0.33333333333333333 ); }
			virtual void    		NumNodesOnFaces( Basso_iVector &nn_face ) 	const;
			virtual void 			NumNodesOnEdges( Basso_iVector &nn_edge ) 	const;
			virtual void    		NodesOnFace ( int f, Basso_iVector &face_nodeids ) 		const;
			virtual void 			NodesOnEdge ( int e, Basso_iVector &edge_nodeids )			const;
			virtual void    		FaceElementType( Basso_Array<Basso_ElementType> &face_basis ) 	const;
			virtual void 			EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const;

			virtual void 		Quadrature( int order, Basso_QuadratureRule &qrule ) const { qrule.SetQuadratureRule( Basso_TRIA_QUAD, order ); }
	
	        virtual void 	Na( const Basso_Point &p, Basso_nVector &shapefunct ) 	const 
            {   
#ifdef ALLOW_DYNAMIC_RESIZE 
                shapefunct.Newsize( NumNodes() ); 
#endif
                shape_tria3( shapefunct.Data(), shapefunct.Inc(), p.x(), p.y() ); 
            }
        
	        virtual void 	DNa( const Basso_Point &p, Basso_nMatrix &grad_shape ) 	const 
	        { 
#ifdef ALLOW_DYNAMIC_RESIZE
                grad_shape.Newsize( NumNodes(), Dimension() );
#endif
	            dshape_tria3( grad_shape.Data(), grad_shape.LDA(), p.x(), p.y() ); 
	        }

			virtual void    Na( Basso_nVector &shapefunct )	const { Na( Centroid(), shapefunct ); }
			virtual void 	DNa( Basso_nMatrix &grad_shape ) const { DNa( Centroid(), grad_shape ); }

		protected:
			
			
	};

	void Basso_Tria3::ParentCoord( Basso_Array< Basso_Point > &xi ) const
	{
	#ifdef ALLOW_DYNAMIC_RESIZE 
		if ( xi.Length()!=NumNodes() )
			xi.Resize(NumNodes());
	#endif
		xi[0]=Basso_Point(0.0,0.0);
		xi[1]=Basso_Point(1.0,0.0);
		xi[2]=Basso_Point(0.0,1.0);
	}

	void 	Basso_Tria3::NumNodesOnFaces( Basso_iVector &nn_face ) 	const
	{
	#ifdef ALLOW_DYNAMIC_RESIZE 
		if ( nn_face.Length()!=NumEdges() )
			nn_face.Resize( NumEdges() );
	#endif
	}

	void 	Basso_Tria3::NumNodesOnEdges( Basso_iVector &nn_edge ) 	const
	{
	#ifdef ALLOW_DYNAMIC_RESIZE 
		if ( nn_edge.Length()!=NumEdges() )
			nn_edge.Resize( NumEdges() );
	#endif
		for ( int i=0; i<NumEdges(); ++i )
			nn_edge[i]=2; 
	}

	void Basso_Tria3::NodesOnFace ( int e, Basso_iVector &face_nodeids ) const
	{
	#ifdef ALLOW_DYNAMIC_RESIZE 
		if ( face_nodeids.Length()!=0 )
			face_nodeids.Resize(0);
	#endif
	}

	void Basso_Tria3::NodesOnEdge ( int e, Basso_iVector &edge_nodeids ) const
	{
	#ifdef ALLOW_DYNAMIC_RESIZE 
		if ( edge_nodeids.Length()!=2 )
			edge_nodeids.Resize(2);
	#endif

		switch (e) {

			case 0:
			edge_nodeids[0]=1; edge_nodeids[1]=2;
			break;

			case 1:
			edge_nodeids[0]=2; edge_nodeids[1]=0;
			break;

			case 2:
			edge_nodeids[0]=0; edge_nodeids[1]=1;
			break;

			default:
			Basso_Warning("Basso_Tria3::NodesOnFace","face id out of range");
			break;

		}
	}

	void Basso_Tria3::FaceElementType( Basso_Array< Basso_ElementType > &face_basis ) const
	{
	#ifdef ALLOW_DYNAMIC_RESIZE 
		if ( face_basis.Length()!=NumFaces() )
			face_basis.Resize(NumFaces());
	#endif
	}

	void Basso_Tria3::EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const
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


