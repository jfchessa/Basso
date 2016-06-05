/*! \file ParentElement.h 

\brief Abstract base clase ParentElement

  Basso ver 1.0
 
\author Jack Chessa, jfchessa@utep.edu
\date Monday, September 22, 2008

*/


#ifndef _BASSO_PARENT_ELEMENT_H_
#define _BASSO_PARENT_ELEMENT_H_

#include "Basso_defs.h"
#include "Basso_Array.h"
#include "Basso_iVector.h"
#include "Basso_nVector.h"
#include "Basso_nMatrix.h"
#include "Basso_Point.h"
#include "Basso_QuadratureRule.h"
#include "Basso_ElementType.h"
//#include "Basso_basic_elements.h"


namespace Basso
{

/** 
	\brief Enumerated type to define element shape/geometry
**/
enum Basso_ElementShape
	{	Basso_POINT, Basso_LINE, Basso_TRIANGLE, Basso_QUADRILATERAL, Basso_TETRAHEDRA, Basso_HEXAHEDRA, Basso_PRISM, Basso_PYRAMID };

/** 
	\brief Enumerated type to element types
**/
/**
enum ElementType
		{ Basso_NONE=0, Basso_POINT1, Basso_LINE2,	Basso_TRIA3,	Basso_QUAD4, Basso_TETRA4, Basso_HEXA8, Basso_PRISM6, Basso_PYRAMID5,
			Basso_LINE3, Basso_TRIA6, Basso_QUAD9, Basso_QUAD8, Basso_TETRA10, Basso_HEXA20, Basso_HEXA27, Basso_PYRAMID13, Basso_PYRAMID14, Basso_PRISM15, Basso_PRISM18 };
*/

/** returns the number of nodes in a Basso_ElementType */
int basis_type_numnodes( const Basso_ElementType &bt )
{
			switch ( bt ) 
			{

				case Basso_NONE:
				return 0;
				break;
				
				case Basso_POINT1:
				return 1;
				break;

				case Basso_LINE2:
				return 2;
				break;

				case Basso_TRIA3:
				return 3;
				break;

				case Basso_QUAD4:
				return 4;
				break;

				case Basso_TETRA4:
				return 4;
				break;

				case Basso_HEXA8:
				return 8;
				break;

				case Basso_PRISM6:
				return 6;
				break;

				case Basso_PYRAMID5:
				return 5;
				break;

				case Basso_LINE3:
				return 3;
				break;

				case Basso_TRIA6:
				return 6;
				break;

				case Basso_QUAD9:
				return 9;
				break;

				case Basso_QUAD8:
				return 8;
				break;

				case Basso_TETRA10:
				return 10;
				break;

				case Basso_HEXA20:
				return 20;
				break;

				case Basso_HEXA27:
				return 27;
				break;

				case Basso_PYRAMID13:
				return 13;
				break;
				
				case Basso_PYRAMID14:
				return 14;
				break;
				
				case Basso_PRISM15:
				return 15;
				break;				

				case Basso_PRISM18:
				return 18;
				break;
				
				default:
				Basso_Warning("basis_type_numnodes","unknown Basso_ElementType");
				return 0;

			}
}
		
		std::ostream &operator << ( std::ostream &out, const Basso_ElementType &bt )
		{
			switch ( bt ) {

				case Basso_NONE:
				out << "NONE";
				return out;
				break;
				
				case Basso_POINT1:
				out << "POINT1";
				return out;
				break;

				case Basso_LINE2:
				out << "LINE2";
				return out;
				break;

				case Basso_TRIA3:
				out << "TRIA3";
				return out;
				break;

				case Basso_QUAD4:
				out << "QUAD4";
				return out;
				break;

				case Basso_TETRA4:
				out << "TETRA4";
				return out;
				break;

				case Basso_HEXA8:
				out << "HEXA8";
				return out;
				break;

				case Basso_PRISM6:
				out << "PRISIM6";
				return out;
				break;

				case Basso_PYRAMID5:
				out << "PYRAMID5";
				return out;
				break;

				case Basso_LINE3:
				out << "LINE3";
				return out;
				break;

				case Basso_TRIA6:
				out << "TRIA6";
				return out;
				break;

				case Basso_QUAD9:
				out << "QUAD9";
				return out;
				break;

				case Basso_QUAD8:
				out << "QUAD8";
				return out;
				break;

				case Basso_TETRA10:
				out << "TETRA10";
				return out;
				break;

				case Basso_HEXA20:
				out << "HEXA20";
				return out;
				break;

				case Basso_HEXA27:
				out << "HEXA27";
				return out;
				break;

				case Basso_PYRAMID13:
				out << "PYRAMD13";
				return out;
				break;
				
				case Basso_PYRAMID14:
				out << "PYRAMID14";
				return out;
				break;
				
				case Basso_PRISM15:
				out << "PRISIM15";
				return out;
				break;				

				case Basso_PRISM18:
				out << "PRISM18";
				return out;
				break;
				
				default:
				out << "UNKNOWN";
				Basso_Warning("basis_type_numnodes","unknown Basso_ElementType");
				return out;

			}
			return out;
		}

/** 
	\brief Basso_ParentElement class
	
	A virtual base class that defines a parent (cononical) element.   This class does not hold any
	connectivity data.
	
*/
class Basso_ParentElement 
{
	public:
		
		Basso_ParentElement()  {}
		
		virtual ~Basso_ParentElement() {}
		
		/** Returns the basis type */
		virtual Basso_ElementType	Type() const=0;

		/** Returns the element shape */
		virtual Basso_ElementShape Shape() const=0;

		/** Returns the number of nodes in the element */
		virtual int	NumNodes() const=0;

		/** Returns the number of edges in the element */               
		virtual int	NumEdges() const=0;

		/** Returns the number of faces in the element */                
		virtual int	NumFaces() const=0;

		/** Returns the dimension of the parent element space */                
		virtual int	Dimension() const=0;

		/** Returns the polynomial order of the basis */
		virtual int	Order()	const=0; 

		/** Returns the parent element coordinates */
		virtual void ParentCoord( Basso_Array< Basso_Point > &pts ) const=0;
		//virtual void ParentCoord( Basso_nMatrix &pts ) const=0;

		/** Returns the point in parent element coordinates of the element centroid */
		virtual Basso_Point Centroid() const=0;

		/** Returns the number of nodes on the element faces given in canonical order 
		*/
		virtual void NumNodesOnFaces( Basso_iVector &nn_face ) const=0; 

	  	/** Returns the number of nodes on the element edges given in canonical order 
	  	*/
	  	virtual void NumNodesOnEdges( Basso_iVector &nn_edege ) const=0;

		/** Returns the local node numbers of the nodes on face f 
			\param f The face id
		*/
		virtual void NodesOnFace ( int f, Basso_iVector &face_nodeids ) const=0;

		/** Returns the local node ids of the nodes on the edge e given in cononical order */
		virtual void NodesOnEdge ( int e, Basso_iVector &edge_nodeids ) const=0;

		/** return geometry of each face in canonical order */
		virtual void FaceElementType( Basso_Array<Basso_ElementType> &face_basis ) const=0;

		/** return geometry of each edge in canonical order	*/
		virtual void EdgeElementType( Basso_Array< Basso_ElementType > &edge_basis ) const=0;

		/** Compute the quadrature rule to integrate on the Basis ooa polynomial function of order p */
		virtual void Quadrature( int p, Basso_QuadratureRule &quad ) const=0;

		// NEED TO ADD THIS
		/** returns shape function at the point p, but uses standard pointers to store the I/O 
		\param p The point in the parent coordinate space at which the shape function is to be evaluated
		\param shapefunct On return has the vector field shape function.
		*/
		//virtual void Na( const Basso_Numeric *p, Basso_Numeric *shapefunct ) const=0;
		
		/** returns shape function at the point p 
		\param p The point in the parent coordinate space at which the shape function is to be evaluated
		\param shapefunct On return has the vector field shape function.
		*/
		virtual void Na( const Basso_Point &p, Basso_nVector &shapefunct ) const=0;

		/** returns shape function at the element centroid */
		virtual void Na( Basso_nVector &shapefunct ) const=0;
		
		// NEED TO ADD THIS
        //virtual void Nv( const Basso_Numeric *p, Basso_Numeric *shapefunct, int sdim=-1 ) const;
		
		/** returns shape function for a vector field at the point p 
		\param p The point in the parent coordinate space at which the shape function is to be evaluated
		\param shapefunct On return has the vector field shape function matrix.
		
		    Nv = [ N1  0  0 ;  
		            0 N1  0 ; 
		            0  0 N1 ;
		           N2  0  0 ;  
		            0 N2  0 ; 
		            0  0 N2 ;
		            :  :  : ;
		           Nn  0  0 ;  
		            0 Nn  0 ; 
		            0  0 Nn  ]  for sdim = 3
		*/
        virtual void Nv( const Basso_Point &p, Basso_nMatrix &shapefunct, int sdim=-1 ) const;

		// NEED TO ADD THIS
		//virtual void DNa( const Basso_Numeric *p, Basso_Numeric *grad_shape ) const=0;
		
		/** Compute the gradient of the shape functions w.r.t. the element space at the point p.  The  
		   derivative is oriented column-wise and the shape funciton id row-wise. 
		\param p The point in the parent coordinate space at which the shape function is to be evaluated.	
		
			DNa = [	N1,xi  N1,eta  N1,zeta ;
					N2,xi  N2,eta  N2,zeta ;
					N3,xi  N3,eta  N3,zeta ;
						: 		: 		:	;
					Nn,xi  Nn,eta  Nn,zeta ]   (Fortran matrix storage)
					
					
		*/
		virtual void DNa( const Basso_Point &p, Basso_nMatrix &grad_shape ) const=0;

		/** Compute the gradient of the shape functions w.r.t. the element space at the element centroid */
		virtual void DNa( Basso_nMatrix &grad_shape ) const=0;
		
		
		/** 
		Compute the gradient of the shape functions w.r.t. the real space
		\param p The point in the parent coordinate space at which the shape function is to be evaluated. 
		\param gcoord the global coordinates of the element's nodes in column format (each column contans
		  a nodes coordinats)
		\param grad_shape on return has the gradient of the shape functions
		\param DNScr - scratch space to hold temproray data (dNdxi)
		\param jacScr - scratch space to hold temproray data sdim x sdim (jac matrix) 
		\param sdim spacial dimension of the problem
		The return value is the determinat of the element Jacobian
		*/
		virtual Basso_Numeric DNx( const Basso_Point &p, const Basso_nMatrix &gcoord, Basso_nMatrix &grad_shape,
		 		   Basso_nMatrix &DNScr, Basso_nMatrix &jacScr, int sdim=-1 ) const
		{
			if ( sdim<0 )
				sdim=this->Dimension();
				
#ifdef BASSO_BOUNDS_CHECK
            if ( DNScr.M()<this->NumNodes() || DNScr.N()<this->Dimension() )
            {
#ifdef ALLOW_DYNAMIC_RESIZE 
                DNScr.Newsize( this->NumNodes(), this->Dimension() );
#else
                Basso_Warning("Basso_ParentElement::DNx","incorrect size of DNScr");
                return 0.0;
#endif
            }
            
            if ( jacScr.M()<sdim || jacScr.N()<sdim )
            {
#ifdef ALLOW_DYNAMIC_RESIZE 
                jacScr.Newsize( sdim, sdim );
#else
                Basso_Warning("Basso_ParentElement::DNx","incorrect size of jacScr");
                return 0.0;
#endif
            }
#endif		
            
            grad_shape.Newsize( this->NumNodes(), sdim );
	
			this->DNa( p, DNScr );
			return grad_shapefunc( this->NumNodes(), sdim, this->Dimension(), 
									gcoord.Data(), gcoord.LDA(),
									DNScr.Data(), DNScr.LDA(), 
									jacScr.Data(), jacScr.LDA(), 
									grad_shape.Data(), grad_shape.LDA() );
		}

		/** 
		Compute the gradient of the shape functions w.r.t. the real space
		\param p The point in the parent coordinate space at which the shape function is to be evaluated. 
		\param gcoord the global coordinates of the element's nodes in column format (each column contans
		  a nodes coordinates)
		\param grad_shape on return has the gradient of the shape functions
		\param sdim spacial dimension of the problem
		The return value is the determinat of the element Jacobian.  
		*/
		virtual Basso_Numeric DNx( const Basso_Point &p, const Basso_nMatrix &gcoord, Basso_nMatrix &grad_shape,
		 		   int sdim=-1 ) const
		{
			if ( sdim<0 )
				sdim=this->Dimension();
				
            Basso_nMatrix DNScr( this->NumNodes(), this->Dimension() ), jacScr( sdim, sdim );
                
            return this->DNx( p, gcoord, grad_shape, DNScr, jacScr, sdim );       
		}

		virtual void Print( ostream &out ) const 
		{
			out << "Basso_ParentElement " << this->Type();
		}

		/** 
		Compute the B-matrix
		\param p The point in the parent coordinate space at which the shape function is to be evaluated. 
		\param gcoord the global coordinates of the element's nodes in column format (each column contans
		  a nodes coordinats)
		\param bmat on return has the element B-matrix
		\param sdim spacial dimension of the problem
		The return value is the determinat of the element Jacobian.  NOte this has a static variable and
		is not really thread safe.
		*/		
		virtual Basso_Numeric Bmatrix( const Basso_Point &p, const Basso_nMatrix &gcoord, Basso_nMatrix &bmat,
		 	 int sdim=-1 ) const
		{
			if ( sdim<0 )
				sdim=this->Dimension();
				
            int vdim = voigt_dim(sdim);
            
            bmat.Newsize( vdim, sdim*(this->NumNodes()) );
            
            Basso_nMatrix dNdx( this->NumNodes(), sdim );
            Basso_Numeric detj = DNx( p, gcoord, dNdx, sdim );
            
            form_bmatrix( this->NumNodes(), sdim, dNdx.Data(), dNdx.LDA(), bmat.Data(), bmat.LDA() );
            
            return detj;
            	    
		}
		
		/**
		Compute the determinant of the element Jacobian.
		\param p The point in the parent coordinate space at which the shape function is to be evaluated. 
		\param gcoord the global coordinates of the element's nodes in column format (each column contans
		  a nodes coordinats)
		\param sdim spacial dimension of the problem
		\param DNScr - scratch space to hold temproray data (dNdxi)
		\param jacScr - scratch space to hold temproray data sdim x sdim (jac matrix) 
		The return value is the determinat of the element Jacobian.  NOte this has a static variable and
		is not really thread safe.
		*/
		virtual Basso_Numeric DetJacobian( const Basso_Point &p, const Basso_nMatrix &gcoord, int sdim,
		    Basso_nMatrix &jacScr, Basso_nMatrix &DNScr ) const
		{
			
			this->DNa( p, DNScr );
			
			return element_detjac( this->NumNodes(), sdim, this->Dimension(), 
				gcoord.Data(), gcoord.LDA(), 
				DNScr.Data(), DNScr.LDA(), 
				jacScr.Data(), jacScr.LDA() );
			
		}
		
		/**
		Compute the determinant of the element Jacobian.
		\param p The point in the parent coordinate space at which the shape function is to be evaluated. 
		\param gcoord the global coordinates of the element's nodes in column format (each column contans
		  a nodes coordinats)
		\param sdim spacial dimension of the problem
		The return value is the determinat of the element Jacobian.  
		
		This does have to allocate a matrix for the jacobian and the shape function gradients so there is 
		version where you pass this information to aleviate this issue.
		*/
		virtual Basso_Numeric DetJacobian( const Basso_Point &p, const Basso_nMatrix &gcoord, int sdim=-1 ) const
		{
            if ( sdim<0 ) sdim = this->Dimension();
            
			Basso_nMatrix jmat( sdim, sdim ), dNdxi( this->NumNodes(), this->Dimension() );
			
			this->DNa( p, dNdxi );
			
			return DetJacobian( p, gcoord, sdim, jmat, dNdxi );
			
		}
		
		
		
	protected:
		
};
	
    
    void Basso_ParentElement::Nv( const Basso_Point &p, Basso_nMatrix &shapefunct, int sdim ) const
    {
        if ( sdim<0 ) sdim = this->Dimension();
        Basso_Warning("void Basso_ParentElement::Nv","Not yet implemented, sorry :(");
    }

	std::ostream &operator << ( std::ostream &out, const Basso_ParentElement &elem )
	{
		elem.Print(out);
		return out;
	}
	

}// end namespace
#endif


