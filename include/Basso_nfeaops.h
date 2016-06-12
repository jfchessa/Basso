/*! \file Basso_nfeaop.h

Some finite element operations in Basso_feaops.h, but using nMatrix and nVector classes.

\author Jack Chessa, jfchessa@utep.edu
\date Sunday, June 5, 2016

*/

#ifndef _BASSO_NUMERIC_FEA_OPERATIONS_H_
#define _BASSO_NUMERIC_FEA_OPERATIONS_H_

// std includes
#include <set>

// Basso includes
#include "Basso_nVector.h"
#include "Basso_iVector.h"
#include "Basso_nMatrix.h"
#include "Basso_iMatrix.h"

#include "Basso_feaops.h"

namespace Basso 
{

using namespace std;

/**
	Computes the element nodal coordinate matrix from a global node coordinate array of points
	and a element connectivity vector.  The element coordinate matrix returns the node coordinates
	in row form so that each row returns the x y and z coordinates for the respective node in 
	that element.
	     
	     		   |  x1   x2   x3 ... xn  |
			cmat = |  y1   y2   y3 ... yn  |
	     		   |  z1   z2   z3 ... zn  |
	
	\param node Global node coordinate array (col node format)
	\param econn Element connectivity array
	\param cmat Element coordinate matrix (column format) on return.
**/
void element_coordinates( const Basso_nMatrix &node, const Basso_iVector &econn, Basso_nMatrix &cmat )
{
	element_coordinates( node.Data(), node.M(), node.LDA(), econn.Data(), econn.Length(), cmat.Data(), cmat.LDA() );
}

/**
	Computes the element nodal coordinate matrix from a global node coordinate array of points
	and a element connectivity vector.  The element coordinate matrix returns the node coordinates
	in row form so that each row returns the x y and z coordinates for the respective node in 
	that element.
	     
	     		   |  x1   x2   x3 ... xn  |
			cmat = |  y1   y2   y3 ... yn  |
	     		   |  z1   z2   z3 ... zn  |
	
	\param node Global node coordinate array (col node format)
	\param econn Element connectivity array
	\param nne is the number of nodes in the element
	\param cmat Element coordinate matrix (column format) on return.
**/
void element_coordinates( const Basso_nMatrix &node, const BASSO_IDTYPE *econn, int nne, Basso_nMatrix &cmat )
{
	element_coordinates( node.Data(), node.M(), node.LDA(), econn, nne, cmat.Data(), cmat.LDA() );
}

 

/**
Sets the C matrix for plane stress isotropic material
\param E is Young's modulus
\param nu is Poisson's ratio
**/
void cmat_pstress( Basso_nMatrix &cmat, Basso_Numeric E, Basso_Numeric nu )
{ cmat_pstress( cmat.Data(), cmat.LDA(), E, nu ); }

/**
Sets the C matrix for 2D plane strain isotropic material
\param E is Young's modulus
\param nu is Poisson's ratio
**/
void cmat_pstrain( Basso_nMatrix &cmat, Basso_Numeric E, Basso_Numeric nu )
{ cmat_pstrain( cmat.Data(), cmat.LDA(), E, nu ); }

/**
Sets the C matrix for 2D axisymmetric isotropic material
\param E is Young's modulus
\param nu is Poisson's ratio

the order of the strains are as follows
        [ rr, theta, zz, rz ]

**/
void cmat_axisymmetric( Basso_nMatrix &cmat, Basso_Numeric E, Basso_Numeric nu )
{ cmat_axisymmetric( cmat.Data(), cmat.LDA(), E, nu ); }

/**
Sets the C matrix for 3D isotropic material
\param E is Young's modulus
\param nu is Poisson's ratio
**/
void cmat_3d( Basso_nMatrix &cmat, Basso_Numeric E, Basso_Numeric nu )
{ cmat_3d( cmat.Data(), cmat.LDA(), E, nu ); }


/**
	Computes D <- alpha * B^T . C . B  + beta * D  ( trans='n')
	   ( or D <- alpha * B . C . B^T  + beta * D if trans='t' )
	
	If trans='n' B(mxn), C(mxm) and D(nxn)
	If trans='t' B(mxn), C(nxn) and D(mxm)
*/
void multbtcb( char trans, Basso_Numeric alpha, const Basso_nMatrix &B, 
 	         Basso_nMatrix &C, Basso_Numeric beta, Basso_nMatrix &D	)
{
    multbtcb( 'n', B.M(), B.N(), alpha, B.Data(), B.LDA(), C.Data(), C.LDA(), beta, D.Data(), D.LDA() );  
}

/**
Returns an ordered array of the node ids in a connectivity matris
\param conn the element connectivity matrix
\param gnids on return contains the node ids
*/
template < class InTtYpE >
void get_gnids( const Basso_Array2D<InTtYpE> &conn, Basso_Array<InTtYpE> &gnids )
{
	set<InTtYpE> gnidSet;
	
	const InTtYpE *nptr = conn.Data();
	for ( int i=0; i<conn.Length(); ++i, ++nptr )
		gnidSet.insert(*nptr);
	gnids.Resize(gnidSet.size());
	
	typename set<InTtYpE>::const_iterator sitr;
	InTtYpE *gptr = gnids.Data();
	for ( sitr=gnidSet.begin(); sitr!=gnidSet.end(); ++sitr, ++gptr )
		*gptr = *sitr;
}


} // end of namespace


#endif



