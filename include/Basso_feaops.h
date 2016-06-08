/*! \file Basso_feaop.h
@file
@brief Some basic finite element operations

\author Jack Chessa, jfchessa@utep.edu
\date Sunday, Oct 30, 2011

*/

#ifndef _BASSO_FEA_OPERATIONS_H_
#define _BASSO_FEA_OPERATIONS_H_

#include <math.h>

#include "Basso_BLAS_wrappers.h"
#include "Basso_LAPACK_wrappers.h"
#include "Basso_stlops.h"

namespace Basso 
{

using namespace std;


//  ************* ALL MATRICIES ARE ASSUMED STORED IN FORTRAN COLUMN FORMAT *************  //
/*! 
	The gradient matrices are returned in the following format
	dNa = [ dN1dxi dN2dxi dN3dxi ... dNndxi dN1deta dN2deta dN3deta ... dNndeta ..... ]
	
	so for C matrix storage (row storage) this is 
		dNa = [ dN1dxi   dN2dxi   dN3dxi ...   dNndxi ;
				dN1deta  dN2deta  dN3deta ...  dNndeta;
				dN1dzeta dN2dzeta dN3dzeta ... dNndzeta ]
			
	and in fortran (column) storage this is
		dNa = [ dN1dxi   dN1deta  dN1dzeta ;
				dN2dxi   dN2deta  dN2dzeta ;
				dN3dxi   dN3deta  dN3dzeta ;
		  			:         :        :
				dNndxi   dNndeta  dNndzeta ]
				
				
  So in all the documentation below the matrices are in fortran column storage format
*/


/** Forms a shape function for interpolating vector fields.  This function takes a shape function 
array (length nn) and returns a nn x sdim matrix, Nv, as follows ( for 3D )

Nv = [ Ni 0 0;
       0 Ni 0;
       0 0 Ni; 
      ---------
				];
  
\param N scalar field shape function
\param sdim Spacial dimension of the field
\param Nv on return contains the vector shape function
\param nn optional parameter that dictates the number of nodes. The default is the number of rows in N
**/
/*	void form_vect_shapefunct( const BASSO_NUMERIC *N, int nn, int sdim, BASSO_NUMERIC *Nv, int lda )
	{
		BASSO_NUMERIC *Nvptr=Nv;
		int inc=lda-sdim*nn;
		for ( int s=0; s<sdim; ++s, Nvptr+=inc )
			for ( int i=0; i<sdim*nn; ++i, ++Nvptr )
				*Nvptr = 0.0;
				
		Nvptr=Nv;
		for ( int i=0; i<nn; ++i, Nvptr+=sdim )
			for ( int s=0; s<sdim; ++s )
				*(Nvptr+s*(lda+1)) = N[i];
			
	}
	
*/	

template < class NuMeRiC >
/*!
	Forms the B matrix from a shape function gradient matrix
	
	      | Ni,x  0   0  |
	      | 0  Ni,y   0  |
	Bi =  | 0   0   Ni,z |
	      | 0  Ni,z Ni,y |
		  | Ni,z 0  Ni,x |
		  | Ni,y Ni,x  0 |
	
	\param DNa_x gradient of the shape function w.r.t. the spacial coordinates
	\param B on return contains the B matrix
	\param nn (optional) number of nodes (default is num rows in DNa_x)
	\param sdim (optional) spacial dimension (default is num cols in DNa_x)
**/
void form_bmatrix( int nn, int sdim, const NuMeRiC *DNa_x, int lda, NuMeRiC *B, int ldb )
{

	if ( sdim==1 )
		for ( int i=0; i<nn; ++i )
			B[0+i*ldb]=DNa_x[i];
			
	else if ( sdim==2 ) 
		for ( int i=0; i<nn; ++i ) 
		{
			B[0+2*i*ldb] = DNa_x[i];		B[0+(2*i+1)*ldb] = 0.0;
			B[1+2*i*ldb] = 0.0;  			B[1+(2*i+1)*ldb] = DNa_x[i+lda];
			B[2+2*i*ldb] = DNa_x[i+lda];	B[2+(2*i+1)*ldb] = DNa_x[i];	
		}

	else // sdim==3	
		for ( int i=0; i<nn; ++i ) 
		{
			B[0+3*i*ldb] = DNa_x[i];		B[0+(3*i+1)*ldb] = 0.0;				B[0+(3*i+2)*ldb] = 0.0;
			B[1+3*i*ldb] = 0.0;				B[1+(3*i+1)*ldb] = DNa_x[i+lda];	B[1+(3*i+2)*ldb] = 0.0;
			B[2+3*i*ldb] = 0.0;				B[2+(3*i+1)*ldb] = 0.0;				B[2+(3*i+2)*ldb] = DNa_x[i+2*lda];
			B[3+3*i*ldb] = 0.0;				B[3+(3*i+1)*ldb] = DNa_x[i+2*lda];	B[3+(3*i+2)*ldb] = DNa_x[i+lda];
			B[4+3*i*ldb] = DNa_x[i+2*lda];	B[4+(3*i+1)*ldb] = 0.0;				B[4+(3*i+2)*ldb] = DNa_x[i];
			B[5+3*i*ldb] = DNa_x[i+lda];	B[5+(3*i+1)*ldb] = DNa_x[i];		B[5+(3*i+2)*ldb] = 0.0;
		}
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
	\param nn the number of nodes in the connectivity
	\param cmat Element coordinate matrix (column format) on return.
	\param lda the leading dimension of \param cmat
	\param sdim Spacial dimension 
**/
template< class NuMeRiC >
void element_coordinates( const NuMeRiC *node, const int sdim, const int ldn, 
    const BASSO_IDTYPE *econn, const int nne,
    NuMeRiC *cmat, const int ldc )
{		
    const NuMeRiC *nptr;
    const BASSO_IDTYPE *cptr = econn;
    NuMeRiC *coordptr;
	for ( int n=0; n<nne; ++n, ++cptr )
	{
        nptr = node + ldn*(*cptr) ; 
        coordptr = cmat + n*ldc;
	    for ( int i=0; i<sdim; ++i, ++nptr, ++coordptr )
            *coordptr = *nptr;
	}		
}


/**
Computes the Jacobian matrix for an element
	\param cmat element coordinate matrix. This is in columns format (each node is a column)
	\param dN_xi gradient of element shape functions w.r.t the parent coordinate system
	\param jmat jacobian matrix (return value with dimension sdim x edim)
	\param nn number of nodes, optional (default is number of columns in cmat)
	\param sdim spacial dimension, optional (default is number of rows in cmat)
	\param edim element dimension, optional (default is number of columns in dN_xi)
	\param coordRowFormat (default is false) if set to true then the coordinate 
	        matrix, cmat is in row format (the nodes are in a row)

Note, that the use of the optional parameters allows for the use of matrices that are of 
incompatible dimension as long as nn, sdim and edim are defined.  This is done so that 
the resizing of element matricies can be kept at a minimum.
**/
template < class NuMeRiC >
void element_jacobian( int nn, int sdim, int edim, 
	const NuMeRiC *cmat, int ldcmat, 
	const NuMeRiC *dNxi, int lddNxi, 
	NuMeRiC *jmat, int ldjmat, bool coordRowFormat=false )
{
	if ( coordRowFormat )                                   //  jac = cmat*dNdxi
	    gemm( 't', 'n', sdim, edim, nn, 1.0, cmat, ldcmat, 
		        dNxi, lddNxi, 0.0, jmat, ldjmat );	
	else
	    gemm( 'n', 'n', sdim, edim, nn, 1.0, cmat, ldcmat, 
		        dNxi, lddNxi, 0.0, jmat, ldjmat );	    
}




/**
An alternate verion of element_jacobian that augments the mtrix to be squre.
Not sure this is really correct for 1D elements.

*/
template < class NuMeRiC >
void element_jacobian2( int nn, int sdim, int edim, 
	const NuMeRiC *cmat, int ldcmat, 
	const NuMeRiC *dNxi, int lddNxi, 
	NuMeRiC *jmat, int ldjmat, bool coordRowFormat=false )
{
	if ( coordRowFormat )                                   //  jac = cmat*dNdxi
	    gemm( 't', 'n', sdim, edim, nn, 1.0, cmat, ldcmat, 
		        dNxi, lddNxi, 0.0, jmat, ldjmat );	
	else
	    gemm( 'n', 'n', sdim, edim, nn, 1.0, cmat, ldcmat, 
		        dNxi, lddNxi, 0.0, jmat, ldjmat );	    

	if ( sdim==2 && edim==1 )  // 1d element in 2d space
	{
        jmat[0+ldjmat] = -jmat[1];
		jmat[1+ldjmat] =  jmat[0];
//		NuMeRiC a=1.0/sqrt( jmat[0+ldjmat]*jmat[0+ldjmat] + jmat[1+ldjmat]*jmat[1+ldjmat] );
//		jmat[0+ldjmat] *= a; 
//		jmat[1+ldcmat] *= a;
	}

	
	else if ( edim<sdim )  // 1d or 2d element in 3d space
	{
		if ( sdim==3 && edim==1 )  // 1d element in 3d space
		{  
			jmat[0+ldjmat]=-jmat[0+ldjmat]; 
			jmat[1+ldjmat]= jmat[0];
			jmat[2+ldjmat]=0.0;
			NuMeRiC alpha=sqrt( jmat[0+ldjmat]*jmat[0+ldjmat] + jmat[1+ldjmat]*jmat[1+ldjmat] 
							+ jmat[2+ldjmat]*jmat[2+ldjmat] );
			if ( alpha==0 )
			{
				jmat[0+ldjmat]=jmat[1]; 
				jmat[1+ldjmat]=0.0; 
				jmat[2+ldjmat]=-jmat[0];
			}
			NuMeRiC a=1.0/alpha;
			jmat[0+ldjmat] *= a; 
			jmat[1+ldjmat] *= a; 
			jmat[2+ldjmat] *= a;
			
		}	// else sdim=3 and edim=2
		
		jmat[0+2*ldjmat]=jmat[1]*jmat[2+ldjmat]-jmat[1+ldjmat]*jmat[2];
		jmat[1+2*ldjmat]=jmat[2]*jmat[0+ldjmat]-jmat[2+ldjmat]*jmat[0];
		jmat[2+2*ldjmat]=jmat[0]*jmat[1+ldjmat]-jmat[0+ldjmat]*jmat[1];
		NuMeRiC a=1.0/sqrt( jmat[0+2*ldjmat]*jmat[0+2*ldjmat] + jmat[1+2*ldjmat]*jmat[1+2*ldjmat] 
					+ jmat[2+2*ldjmat]*jmat[2+2*ldjmat] );
		jmat[0+2*ldjmat] *= a; 
		jmat[1+2*ldjmat] *= a; 
		jmat[2+2*ldjmat] *= a;		
	}
	
}


template < class NuMeRiC >
NuMeRiC det2by2( NuMeRiC *A, int lda )
{
	return A[0]*A[1+lda]-A[1]*A[lda];
}

template < class NuMeRiC >
NuMeRiC det3by3( NuMeRiC *A, int lda )
{
    return A[0]*A[lda+1]*A[2*lda+2] - A[0]*A[lda+2]*A[2*lda+1] 
		- A[1]*A[lda]*A[2*lda+2] + A[1]*A[lda+2]*A[2*lda] 
		+ A[2]*A[lda]*A[2*lda+1] - A[2]*A[lda+1]*A[2*lda];
}

/**
Computes the Jacobian (determinant) and Jacobian matrix for an element
	\param cmat element coordinate matrix. This is in columns format (each node is a column)
	\param dN_xi gradient of element shape functions w.r.t the parent coordinate system
	\param jmat jacobian matrix (return value with dimension sdim x edim)
	\param nn number of nodes, optional (default is number of columns in cmat)
	\param sdim spacial dimension, optional (default is number of rows in cmat)
	\param edim element dimension, optional (default is number of columns in dN_xi)
	\param coordRowFormat (default is false) if set to true then the coordinate 
	        matrix, cmat is in row format (the nodes are in a row)
	\return the det(J)

Note, that the use of the optional parameters allows for the use of matrices that are of 
incompatible dimension as long as nn, sdim and edim are defined.  This is done so that 
the resizing of element matricies can be kept at a minimum.
**/
template < class NuMeRiC >
NuMeRiC element_detjac( int nn, int sdim, int edim, 
	const NuMeRiC *cmat, int ldcmat, 
	const NuMeRiC *dNxi, int lddNxi, 
	NuMeRiC *jmat, int ldjmat, bool coordRowFormat=false )
{
	NuMeRiC det;
	element_jacobian( nn, sdim, edim, cmat, ldcmat, 
		dNxi,lddNxi, jmat, ldjmat, coordRowFormat );
				
	// get the determinant.
	if ( sdim==edim )
	{
		if ( sdim==1 )
			det = abs(jmat[0]);
		else if ( sdim==2 )
			det = det2by2( jmat, ldjmat );
		else
			det = det3by3( jmat, ldjmat );
			
		return det;	
	}
	
	if ( edim ==1 )
	{
        det = 0.0;
        for  ( int i=0; i<sdim; ++i )
            det += jmat[i]*jmat[i];
        det = sqrt( det );
        return det;
	}
	
	// edim = 2 and sdim = 3
	NuMeRiC xi1=jmat[0], xi2=jmat[1], xi3=jmat[2];
	NuMeRiC eta1=jmat[ldjmat], eta2=jmat[ldjmat+1], eta3=jmat[ldjmat+2];
	det = sqrt( (xi1*xi1+xi2*xi2+xi3*xi3)*(eta1*eta1+eta2*eta2+eta3*eta3) );
	
	return det;
}

template < class NuMeRiC >
NuMeRiC inv2by2( NuMeRiC *A, int lda )
{
	NuMeRiC a=A[0], b=A[1], c=A[lda], d=A[lda+1];
    NuMeRiC det = a*d-b*c, invDet=1.0/det;
    A[0] 		=  invDet*d;   
    A[1]   		= -invDet*b;   
    A[lda]   	= -invDet*c;
    A[lda+1] 	=  invDet*a;
    return det;
}

template < class NuMeRiC >
NuMeRiC inv3by3( NuMeRiC *A, int lda )
{
    NuMeRiC a0=A[0],  a3=A[lda],    a6=A[2*lda];
    NuMeRiC a1=A[1],  a4=A[lda+1],  a7=A[2*lda+1];
    NuMeRiC a2=A[2],  a5=A[lda+2],  a8=A[2*lda+2];
    NuMeRiC det = a0*a4*a8 - a0*a5*a7 - a1*a3*a8 + a1*a5*a6 + a2*a3*a7 - a2*a4*a6;
    NuMeRiC invDet=1.0/det;
    A[0]=invDet*(a4*a8-a5*a7); A[lda]  =invDet*(a5*a6-a3*a8); A[2*lda]  =invDet*(a3*a7-a4*a6);
    A[1]=invDet*(a2*a7-a1*a8); A[lda+1]=invDet*(a0*a8-a2*a6); A[2*lda+1]=invDet*(a1*a6-a0*a7);
    A[2]=invDet*(a1*a5-a2*a4); A[lda+2]=invDet*(a2*a3-a0*a5); A[2*lda+2]=invDet*(a0*a4-a1*a3);
    return det;
}


/**
Computes the gradient of the shape function with respect to the spacial coordinates.c
	\param cmat element coordinate matrix. This is in columns format (each node is a column)
	\param dN_xi gradient of element shape functions w.r.t the parent coordinate system
	\param jac the jacobian matrix of the element (return value is actually the inverse).
				The dimension of the jacobian matrix must typically be sdim x sdim so that 
				it is big enough to store the inverse
	\param dN_x gradient of element shape functions w.r.t the spacial coordinate system (return value)
	\param nn number of nodes, optional (default is number of columns in cmat)
	\param sdim spacial dimension, optional (default is number of rows in cmat)
	\param edim element dimension, optional (default is number of columns in dN_xi)
	\param coordRowFormat (default is false) if set to true then the coordinate 
	        matrix, cmat is in row format (the nodes are in a row)

The return value is the determinate of the Jacobian at that point

Note, that the use of the optional parameters allows for the use of matrices that are of 
incompatible dimension as long as nn, sdim and edim are defined.  This is done so that 
the resizing of element matrices can be kept at a minimum.

Also, dN_x and dN_xi cannot be the same memory space for this function to work correctly 

**/
template < class NuMeRiC >
NuMeRiC grad_shapefunc( int nn, int sdim, int edim, 
	const NuMeRiC *cmat, int ldcmat, 
	const NuMeRiC *dN_xi, int lddNxi, 
	NuMeRiC *jac, int ldjac, 
	NuMeRiC *dN_x, int lddNx, bool coordRowFormat=false )
{
	
	element_jacobian( nn, sdim, edim, cmat, ldcmat, dN_xi, lddNxi, jac, ldjac, coordRowFormat );
    NuMeRiC det;
   
	if ( sdim == edim )
	{
		if ( sdim==1 )
		{
			det = abs(jac[0]);
			NuMeRiC idet=1.0/jac[0];
			jac[0]=idet;
			for ( int I=0; I<nn; ++I )
				dN_x[I] = dN_xi[I]*idet;
			return det;
		}
		
		else if ( sdim==2 )
		{
			det = inv2by2( jac, ldjac );
		}
		
		else // sdim == 3
		{
			det = inv3by3( jac, ldjac );
		}
	}
	
	else if ( edim == 1 )  // edim=1 and sdim = 2 or 3
	{
        det = 0.0;
        for  ( int i=0; i<sdim; ++i )
        {
            det += jac[i]*jac[i];
            
            if ( jac[i] != 0 )     // this can be made more efficeint.
                jac[i] = 1.0/jac[i];
            
            for ( int I=0; I<nn; ++I ) // this can be made more efficeint.
                dN_x[ i*lddNx + I ] = dN_xi[ I ]*jac[i];
        }
        det = sqrt( det );
        return det;
	}
	
	else // edim==2 and sdim==3
	{
       // Basso_Error("grad_shapefunc","Not yet done for sdim=3 and edim=2");
		NuMeRiC xi1=jac[0], xi2=jac[1], xi3=jac[2];
		NuMeRiC eta1=jac[ldjac], eta2=jac[ldjac+1], eta3=jac[ldjac+2];
		jac[ ldjac*2 ] = xi2*eta3 - xi3*eta2;
		jac[ ldjac*2 + 1 ] = xi3*eta1 - xi1*eta3;
		jac[ ldjac*2 + 2 ] = xi1*eta2 - xi2*eta1;
		NuMeRiC a = 1.0/sqrt( pow(jac[ldjac*2],2) + pow(jac[ldjac*2+1],2) + pow(jac[ldjac*2+2],2) );
		jac[ ldjac*2 ] *= a;
		jac[ ldjac*2 + 1 ] *= a;
		jac[ ldjac*2 + 2 ] *= a;
		inv3by3( jac, ldjac );
		det = sqrt( (xi1*xi1+xi2*xi2+xi3*xi3)*(eta1*eta1+eta2*eta2+eta3*eta3) );
	}
    
    // dN_x = dN_dxi * inv(jac)
    NuMeRiC alpha=1.0, beta=0.0;
    gemm( 'n', 'n', nn, sdim, edim, alpha, dN_xi, lddNxi, jac, ldjac, beta, dN_x, lddNx );

	return det;
}

/**
An alternate version of grad_shapefunc.  It works with the jacobian 
matrix from element_jacobian2

*/
template < class NuMeRiC >
NuMeRiC grad_shapefunc2( int nn, int sdim, int edim, 
	const NuMeRiC *cmat, int ldcmat, 
	const NuMeRiC *dN_xi, int lddNxi, 
	NuMeRiC *jac, int ldjac, 
	NuMeRiC *dN_x, int lddNx, bool coordRowFormat=false )
{
	element_jacobian2( nn, sdim, edim, cmat, ldcmat, dN_xi, lddNxi, jac, ldjac, coordRowFormat );
    NuMeRiC det;
   
    if ( sdim==1 )
    {
        det = jac[0];
        dN_x[0] = dN_xi[0]/det;
        return det; 
    }
    
    else if ( sdim==2 )
        det = inv2by2( jac, ldjac );
    
    else // sdim ==3
        det = inv3by3( jac, ldjac );

    // dN_x = dN_dxi * inv(jac)
    NuMeRiC alpha=1.0, beta=0.0;
    gemm( 'n', 'n', nn, sdim, sdim, alpha, dN_xi, lddNxi, jac, ldjac, beta, dN_x, lddNx );
    
	return det;
}

/*
template < class NuMeRiC >
NuMeRiC grad_shapefunc( int nn, int sdim, int edim, 
	const NuMeRiC *cmat, int ldcmat, 
	const NuMeRiC *dN_xi, int lddNxi, 
	NuMeRiC *jac, int ldjac, 
	NuMeRiC *dN_x, int lddNx, bool coordRowFormat=false )
{
	element_jacobian( nn, sdim, edim, cmat, ldcmat, dN_xi, lddNxi, jac, ldjac, coordRowFormat );
    NuMeRiC det;

   int ipiv[edim];
    int info;
	NuMeRiC b[sdim];  // CAN I GET RID OF b?

	info = getrf( sdim, edim, jac, ldjac, ipiv );

	NuMeRiC det=1.0;
	for ( int i=0; i<edim; ++i )
	{
		det *= jac[i+i*ldjac];
		if ( ipiv[i]!=i+1 )
			det *= -1.0;
	}

	for  ( int i=0; i<nn; ++i ) 
	{
		for ( int j=0; j<sdim; ++j ) 
			b[j]=dN_xi[i+j*lddNxi];
			
		info = getrs( 't', 1, edim, jac, ldjac, ipiv, b, sdim );
		
		for ( int j=0; j<sdim; ++j ) 
			dN_x[i+j*lddNx]=b[j];
	}

	return det;
}
*/


/**
Constructs a scatter vector, sctr, assuming a fixed number of dofs for
  each node.  The mapping of the local to global dofs, ldof is the global 
  dofs is as follows 
 
    gdofi=n*nndof+ldofi
 
  where gdofi is a global dof, ldofi is a local dof and nndof is the number
  of dofs for each node.
 
        \param conn - is the element connectivity
	    \param nne - number of nodes in teh connectivity
		\param nndof - the max number of dofs declared at each node
        \param ldof -  is a vector of the local dofs to map
		\param ndot - the number of dofs in ldof
		\param sctr - on retrn has the appropriate scatter vector
*/
void set_scatter( const int *conn, int nne, int nndof, const int *ldof, int ndof, int *sctr )
{
	const int *cptr=conn, *ldptr;
	int *sptr=sctr;
	int i, s;
	for ( i=0; i<nne; ++i, ++cptr )
	{
		ldptr=ldof;
		for ( s=0; s<ndof; ++s, ++ldptr, ++sptr )
			*sptr = (*cptr)*nndof + *ldptr;
	}
}
/**
Constructs a scatter vector, sctr, assuming a fixed number of dofs for
  each node.  The mapping of the local to global dofs, ldof is the global 
  dofs is as follows 
 
    gdofi=n*nndof+ldofi
 
  where gdofi is a global dof, ldofi is a local dof and nndof is the number
  of dofs for each node.
 
        \param conn - is the element connectivity
	    \param nne - number of nodes in teh connectivity
		\param nndof - the max number of dofs declared at each node
		\param sctr - on retrn has the appropriate scatter vector
*/
void set_scatter( const int *conn, int nne, int nndof, int *sctr )
{
	const int *cptr=conn;
	int *sptr=sctr;
	int i, s;
	for ( i=0; i<nne; ++i, ++cptr )
		for ( s=0; s<nndof; ++s, ++sptr )
			*sptr = (*cptr)*nndof + s;
}

/**
Returns the dimension of the Voigt stress vector for a given
spacial dimension \param sdim.
*/
int voigt_dim( int sdim )
{
    switch ( sdim )
    {
        case 1:
        return 1;
        
        case 2:
        return 3;
        
        case 3:
        return 6;
        
        default:
        Basso_Warning("voigt_dim","incorrect value for sdim (1, 2 or 3)");
        return 6;
    }
}

/**
Sets the C matrix for plane stress isotropic material
\param E is Young's modulus
\param nu is Poisson's ratio
**/
template < class NuMeRiC >
void cmat_pstress( NuMeRiC *cmat, int ldc, NuMeRiC E, NuMeRiC nu )
{
	NuMeRiC c1 = E/(1.0-nu*nu), c2=nu*c1, c3=0.5*(1.0-nu)*c1;
    // C=[c1 c2 0;c2 c1 0; 0 0 c3];
	NuMeRiC *cptr = cmat;
	*(cptr++) = c1; *(cptr++) = c2; *(cptr++) = 0.0;  cptr+=ldc-3;
	*(cptr++) = c2; *(cptr++) = c1; *(cptr++) = 0.0;  cptr+=ldc-3;
	*(cptr++) = 0.0; *(cptr++) = 0.0; *(cptr) = c3;  
}

/**
Sets the C matrix for 2D plane strain isotropic material
\param E is Young's modulus
\param nu is Poisson's ratio
**/
template < class NuMeRiC >
void cmat_pstrain( NuMeRiC *cmat, int ldc, NuMeRiC E, NuMeRiC nu )
{
	NuMeRiC c0=E/(1.-2.*nu)/(1+nu), c1=(1.-nu)*c0, c2=nu*c0, c3=0.5*(1.-2.*nu)*c0;
    // C=[c1 c2 0;c2 c1 0; 0 0 c3];
	NuMeRiC *cptr = cmat;
	*(cptr++) = c1; *(cptr++) = c2; *(cptr++) = 0.0;  cptr+=ldc-3;
	*(cptr++) = c2; *(cptr++) = c1; *(cptr++) = 0.0;  cptr+=ldc-3;
	*(cptr++) = 0.0; *(cptr++) = 0.0; *(cptr) = c3;  
}

/**
Sets the C matrix for 2D axisymmetric isotropic material
\param E is Young's modulus
\param nu is Poisson's ratio

the order of the strains are as follows
        [ rr, theta, zz, rz ]

**/
template < class NuMeRiC >
void cmat_axisymmetric( NuMeRiC *cmat, int ldc, NuMeRiC E, NuMeRiC nu )
{
	NuMeRiC c0=E/(1.-2*nu)/(1.+nu), c1=(1.-nu)*c0, c2=nu*c0, c3=0.5*(1.-2.*nu)*c0;
/*
        C=[ c1 c2 c2 0;
            c2 c1 c2 0; 
            c2 c2 c1 0; 
             0  0  0 c3];
*/
	NuMeRiC *cptr = cmat;
	*(cptr++) = c1; *(cptr++) = c2; *(cptr++) = c2; *(cptr++) = 0.0;  cptr+=ldc-4;
	*(cptr++) = c2; *(cptr++) = c1; *(cptr++) = c2; *(cptr++) = 0.0;  cptr+=ldc-4;
	*(cptr++) = c2; *(cptr++) = c2; *(cptr++) = c1; *(cptr++) = 0.0;  cptr+=ldc-4;
	*(cptr++) = 0.0; *(cptr++) = 0.0; *(cptr++) = 0.0;  *(cptr) = c3;  
}

/**
Sets the C matrix for 3D isotropic material
\param E is Young's modulus
\param nu is Poisson's ratio
**/
template < class NuMeRiC >
void cmat_3d( NuMeRiC *cmat, int ldc, NuMeRiC E, NuMeRiC nu )
{
	NuMeRiC c0=E/(1.-2.*nu)/(1+nu), c1=(1.-nu)*c0, c2=nu*c0, c3=0.5*(1.-2.*nu)*c0;
	NuMeRiC *cptr = cmat;
	*(cptr++) = c1; *(cptr++) = c2; *(cptr++) = c2; *(cptr++) = 0.0; *(cptr++) = 0.0; *(cptr++) = 0.0;  cptr+=ldc-6;
	*(cptr++) = c2; *(cptr++) = c1; *(cptr++) = c2; *(cptr++) = 0.0; *(cptr++) = 0.0; *(cptr++) = 0.0;  cptr+=ldc-6;
	*(cptr++) = c2; *(cptr++) = c2; *(cptr++) = c1; *(cptr++) = 0.0; *(cptr++) = 0.0; *(cptr++) = 0.0;  cptr+=ldc-6;
	*(cptr++) =0.0; *(cptr++) =0.0; *(cptr++) =0.0; *(cptr++) =  c3; *(cptr++) = 0.0; *(cptr++) = 0.0;  cptr+=ldc-6;
	*(cptr++) =0.0; *(cptr++) =0.0; *(cptr++) =0.0; *(cptr++) = 0.0; *(cptr++) =  c3; *(cptr++) = 0.0;  cptr+=ldc-6;
	*(cptr++) =0.0; *(cptr++) =0.0; *(cptr++) =0.0; *(cptr++) = 0.0; *(cptr++) = 0.0; *(cptr) = c3; 
}


/**
Computes the effective values (Mises) for a symmetric tensor (matrix).

	\param s - is a pointer array to a symetric tensor of the form
		[ s11 s22 s33 s12 s23 s31 ]
		[ s11 s22 s12 ]
		
	\param sdim - the spacial dimension of the tensor (1, 2 or 3, the default is 3)

**/
template < class NuMeRiC >
NuMeRiC mises_stress( const NuMeRiC *s, int sdim=3 )
{
	NuMeRiC srt2 = 0.707106781186547;
	NuMeRiC svm=0.0;
	switch (sdim)
	{
		case 1:
		svm=fabs( s[0] );
		break;
			
		case 2:
		svm=srt2*pow( pow(s[0]-s[1],2) + s[1]*s[1] + s[0]*s[0] + 6*s[2]*s[2], 0.5 );
		break;
		
		case 3:
		svm=srt2*pow( pow(s[0]-s[1],2) + pow(s[0]-s[2],2) + pow(s[2]-s[1],2) 
			+ 6*( s[3]*s[3] + s[4]*s[4] + s[5]*s[5] ), 0.5 );
		break;
			
		default:
		Basso_Warning("NuMeRiC mises_stress( const NuMeRiC *, int sdim, NuMeRiC * )",
			"Incorrect sdim value, sdim={1,2,3}");
			
	}
	return svm;	
}

/**
Computes the principal values for a symmetric tensor in Voigt form.

	\param s - is a pointer array to a symetric tensor of the form
		[ s11 s22 s33 s12 s23 s31 ]
		[ s11 s22 s12 ]
		
	\param sdim - the spacial dimension of the tensor (1, 2 or 3, the default is 3)
	\param eval - on return contains the eigen values (principal values)
	\param evct - on return will contain the eigen vector

**/
template < class NuMeRiC >
void principal_value( const NuMeRiC *s, int sdim, NuMeRiC *eval, NuMeRiC *evct=NULL )
{
    if ( sdim == 1 )
    {
        eval[0] = s[0];
        if ( evct != NULL )
            evct[0] = 1.0;
        return;
    }
    
    if ( sdim == 2 )
    {
        NuMeRiC s11=s[0], s22=s[1], s12=s[2];
        NuMeRiC dd = s11*s11-2*s11*s22+4*s12*s12+s22*s22, tr2 = s11+s22;
        if ( dd<0 )
        {    
            Basso_Warning("principal_value","negative discriminant");
            return 0;
        }
        NuMeRiC d = pow(dd,0.5);
        eval[0] = 0.5*(tr2-d);
		eval[1] = 0.5*(tr2+d);
       
        if ( evct != NULL )
        {
            NuMeRiC a1 = 0.5*(s11-s22-d)/s12, a2 = 0.5*(s11-s22+d)/s12;
            NuMeRiC a1nrminv = pow( a1*a1 + 1, -0.5 ), a2nrminv = pow( a2*a2 + 1, -0.5 ); 
            evct[0] = a1*a1nrminv;
            evct[1] = a1nrminv;
            evct[2] = a2*a2nrminv;
            evct[3] = a2nrminv;
        }  
        return;      
    }
    
    if ( sdim == 3 )
    {
		// compute the eigen values and vectors
		int lwork=8;
		bool freeVct=false;
		eval = new NuMeRiC[3];
		
		if ( evct = NULL )
		{
			evct = new NuMeRiC[9];
			freeVct = true;
		}
		NuMeRiC work[lwork];
		
		evct[0]=s[0];  
		evct[1]=s[3]; evct[4]=s[1];  
		evct[2]=s[5]; evct[5]=s[4]; evct[8]=s[2];
		int info = syev( 'v', 'l', 3, evct, 3, eval, work, lwork );
		
		if ( freeVct )
			delete [] evct;
			
		return;      
    }
}


/**
Computes the maximum shear stresses from the principal normal stresses
\param s - contains the principal values (s1 s2 s3)
\param sdim -  the spacial dimension
\param tau - on return contains the max shear stresses (tau12, tau23, tau31)
**/
template < class NuMeRiC >
void max_shear( const NuMeRiC *s, int sdim, NuMeRiC *tau )
{
	switch ( sdim )
	{
		case 1:
		tau[0]=0;
		return;
		
		case 2:
		tau[0]=0.5*fabs(s[2]-s[1]);
		return;
		
		case 3:
		tau[0]=0.5*fabs(s[2]-s[1]);
		tau[1]=0.5*fabs(s[2]-s[3]);
		tau[2]=0.5*fabs(s[3]-s[1]);
		return;
		
		default:
		Basso_Warning("void max_shear( const NuMeRiC *s, int sdim, NuMeRiC *tau )","invalid sdim");
	}
	
	return;
}

/**
Computes the octahedral stress components

\param s -  the stress tensor (symmetric in voigt form)

		[ s11 s22 s33 s12 s23 s31 ]
		[ s11 s22 s12 ]
\param sdim - the spacial dimension
\param oct - on return has the octahedral stress tensor  oct[0]=normal octahedral
oct[1] = octahedral shear stress
*/
template < class NuMeRiC >
void octahedral_stress( const NuMeRiC *s, int sdim, NuMeRiC *oct )
{
	oct[0]=0.0;
	NuMeRiC third = 0.33333333333333333333333;
	for ( int i=0; i<sdim; ++i )
		oct[0] += s[i];
	oct[0] *= third;
	
	switch (sdim)
	{
		case 1:
		oct[1]=oct[0];
		break;
			
		case 2:
		oct[1]=third*pow( pow(s[0]-s[1],2) + s[1]*s[1] + s[0]*s[0] + 6*s[2]*s[2], 0.5 );
		break;
		
		case 3:
		oct[1]=third*pow( pow(s[0]-s[1],2) + pow(s[0]-s[2],2) + pow(s[2]-s[1],2) 
			+ 6*( s[3]*s[3] + s[4]*s[4] + s[5]*s[5] ), 0.5 );
		break;
			
		default:
		Basso_Warning("NuMeRiC octahedral_stress( const NuMeRiC *, int sdim, NuMeRiC * )",
			"Incorrect sdim value, sdim={1,2,3}");
			
	}
}

/**
Computes the deviatoric stress components
The reutnr value is the hydrostatic stress
\param s -  the stress tensor (symmetric in voigt form)

		[ s11 s22 s33 s12 s23 s31 ]
		[ s11 s22 s12 ]
		
\param sdim - the spacial dimension
\param dev - on return has the deviatoric stress tensor(if not supplied 
it is not calculated)
*/
template < class NuMeRiC >
NuMeRiC deviatoric_stress( const NuMeRiC *s, int sdim=3, NuMeRiC *dev=NULL )
{
	NuMeRiC hyd=0.0;
	
	for ( int i=0; i<sdim; ++i )
		hyd += s[i];
	hyd *= 0.3333333333333333;	
	
	if ( dev != NULL )
	{
		switch (sdim)
		{
			case 1:
			dev[0]=dev[0]-hyd;
			break;
			
			case 2:
			dev[0]=s[0]-hyd; dev[1]=s[1]-hyd; dev[2]=s[2];
			break;
			
			case 3:
			dev[0]=s[0]-hyd; dev[1]=s[1]-hyd; dev[2]=s[2]-hyd;
			dev[3]=s[3]; dev[4]=s[4]; dev[5]=s[5];
			break;
			
			default:
			Basso_Warning("NuMeRiC deviatoric_stress( const NuMeRiC *s, int sdim=3, NuMeRiC *dev=NULL )",
				"Incorrect sdim value, sdim={1,2,3}");
			
		}
	}
	
	return hyd;
}


/**
Fills in the upper half of the matrix with the lower half (symmetric)
**/
template < class NuMeRiC >
void fill_upper_half( NuMeRiC *A, int m, int lda )
{
    // fill upper triangle of A with lower triangl (symmetry)
    NuMeRiC *lptr, *uptr;
    for ( int i=0; i<m; ++i )
    {
        lptr = A + i*lda + 1 + i;
        uptr = A + (i+1)*lda + i;
	    for ( int j=i+1; j<m; ++j, ++lptr, uptr+=lda )  
		{
			// A_ij = A_ji or A[ i+j*lda ] = A[ j+i*lda ];
            *uptr = *lptr;
		}
	}
}

/**
	Computes D <- alpha * B^T . C . B  + beta * D  ( trans='n')
	   ( or D <- alpha * B . C . B^T  + beta * D if trans='t' )
	
	If trans='n' B(mxn), C(mxm) and D(nxn)
	If trans='t' B(mxn), C(nxn) and D(mxm)
*/
template < class NuMeRiC >
void multbtcb( char trans, int m, int n, NuMeRiC alpha, const NuMeRiC *B, int ldb, 
 	 const NuMeRiC *C, int ldc, NuMeRiC beta, NuMeRiC *D, int ldd	)
{
    
    if ( beta==0.0 && alpha==0.0 ) return;
    
	int i, j, k, l, dinc, cinc, binc, dnr, cnr;
	const NuMeRiC *Cptr, *B1ptr, *B2ptr;
    NuMeRiC *Dptr;
	bool tr;
	
	if ( trans == 'n' || trans == 'N' )	
	{
		dnr=n;
		cnr=m;
		tr=false;
	}
	else if ( trans == 't' || trans == 'T'  )
	{
		dnr=m;
		cnr=n;
		tr=true;
	}
	else
	{
		cout << "\n*** error in trans parameter in multbtcb ***\n";
		return;
	}
	dinc=ldd-dnr;
	cinc=ldc-cnr;
	binc=ldb-m;
	
	Dptr=D;
	i=0;
	if ( beta == 0.0 ) 
	{ 
		for ( j=0; j<=i; ++j, Dptr+=dinc+j  )  // D -> alpha * op( B^T C B )   JUST LOWER HALF
			for ( i=j; i<dnr; ++i, ++Dptr )
			{
				*Dptr = 0.0; // add in beta*D
				Cptr = C;
				B1ptr = &(B[i*ldb]);
				for ( k=0; k<cnr; ++k, Cptr+=cinc, ++B1ptr )
				{
				    B2ptr = &(B[j*ldb]);
					for ( l=0; l<cnr; ++l, ++Cptr, ++B2ptr )  
					{ 
						// D += alpha*B_ki*C_kl*B_lj  same as alpha*B[k+i*ldb]*C[l+k*ldc]*B[l+j*ldb]
                        *Dptr += alpha*(*B1ptr)*(*Cptr)*(*B2ptr);
					}
				}
			}   
	}
			     
	else if ( beta == 1.0 ) 
	{
		for ( j=0; j<=i; ++j, Dptr+=dinc+j  )  // D -> alpha * op( B^T C B )   JUST LOWER HALF
			for ( i=j; i<dnr; ++i, ++Dptr )
			{
				*Dptr = *Dptr; // add in beta*D
				Cptr = C;
				B1ptr = &(B[i*ldb]);
				for ( k=0; k<cnr; ++k, Cptr+=cinc, ++B1ptr )
				{
				    B2ptr = &(B[j*ldb]);
					for ( l=0; l<cnr; ++l, ++Cptr, ++B2ptr )  
					{ 
						// D += alpha*B_ki*C_kl*B_lj  same as alpha*B[k+i*ldb]*C[l+k*ldc]*B[l+j*ldb]
                        *Dptr += alpha*(*B1ptr)*(*Cptr)*(*B2ptr);
					}
				}
			}
	}
			   
	else
	{
		for ( j=0; j<=i; ++j, Dptr+=dinc+j  )  // D -> alpha * op( B^T C B )   JUST LOWER HALF
			for ( i=j; i<dnr; ++i, ++Dptr )
			{
				*Dptr = beta*( *Dptr ); // add in beta*D
				Cptr = C;
				B1ptr = &(B[i*ldb]);
				for ( k=0; k<cnr; ++k, Cptr+=cinc, ++B1ptr )
				{
				    B2ptr = &(B[j*ldb]);
					for ( l=0; l<cnr; ++l, ++Cptr, ++B2ptr )  
					{ 
						// D += alpha*B_ki*C_kl*B_lj  same as alpha*B[k+i*ldb]*C[l+k*ldc]*B[l+j*ldb]
                        *Dptr += alpha*(*B1ptr)*(*Cptr)*(*B2ptr);
					}
				}
			}
	}
	
    fill_upper_half( D, dnr, ldd );
	
}

/** Same as multbtcb but uses BLAS
	Computes D <- alpha * B^T . C . B  + beta * D  ( trans='n')
	   ( or D <- alpha * B . C . B^T  + beta * D if trans='t' )
	
	If trans='n' B(mxn), C(mxm) and D(nxn)
	If trans='t' B(mxn), C(nxn) and D(mxm)
	
	It seems that the multbtcb is about 25% faster than this routine so I would 
	recommend using that in lieu of multbtcb2
*/
template < class NuMeRiC >
void multbtcb2( char trans, int m, int n, NuMeRiC alpha, const NuMeRiC *B, int ldb, 
 	 const NuMeRiC *C, int ldc, NuMeRiC beta, NuMeRiC *D, int ldd	)
{
	if ( trans == 'n' || trans == 'N'  )
	{
		NuMeRiC A[n*m];
		int lda=n;
		gemm('t','n',n,m,m,1.0,B,ldb,C,ldc,0.0,A,lda); 	 // A -> B^T * C
		gemm('n','n',n,n,m,alpha,A,lda,B,ldb,1.0,D,ldd); // D -> beta*D + alpha * A * B;
	}
	
	else
	{
		Basso_Warning("function multbtcb2","trans = 't' not yet implemented");
	}
/*gemm( const char transa, const char transb, const int m, const int n, const int k, 
					const float alpha, const float *A, const int lda,
 					const float *B, const int ldb, 
					const float beta, float *C, const int ldc )
					*/
}
/*********************************************************************************************************
**
**                             	 Q U A D R A T U R E     F U N C T I O N S
**
**********************************************************************************************************/

/**
Computes the n Gaussian quadrature rule
\param pts - pointer to the array that holds the quadarature points [s1, s2, s3, ...] (must be allocated)
\param wts - pointer to the array of weights [w1 w2 ..] (must be allocated)
\param n - The number of points
\param lda -  The leading dimension of pts
**/
template <class NuMeRiC>
void quadrature_gauss1d( NuMeRiC *pts, NuMeRiC *wts, int npts, int lda=1 )
{
    NuMeRiC *ptPtr=pts, *wtPtr=wts;
    
    switch (npts) 
    {
    
        case 1:
        *ptPtr = 0.0;   *wtPtr = 2.0; 
        break;
        
        case 2:
        *ptPtr      = -0.577350269189626;       *wtPtr      = 1.0; 
        *(ptPtr+=lda)  =  0.577350269189626;       *(++wtPtr)  = 1.0; 
        break;
        
        case 3:
        *ptPtr      = -0.774596669241483;       *wtPtr      = 0.555555555555556; 
        *(ptPtr+=lda)  =  0.000000000000000;       *(++wtPtr)  = 0.888888888888889; 
        *(ptPtr+=lda)  =  0.774596669241483;       *(++wtPtr)  = 0.555555555555556;
        break; 
        
        case 4:
        *ptPtr      = -0.861134311594053;       *wtPtr      = 0.347854845137454; 
        *(ptPtr+=lda)  = -0.339981043584856;       *(++wtPtr)  = 0.652145154862546; 
        *(ptPtr+=lda)  =  0.339981043584856;       *(++wtPtr)  = 0.652145154862546;
        *(ptPtr+=lda)  =  0.861134311594053;       *(++wtPtr)  = 0.347854845137454;
        break; 
        
        case 5:
        *ptPtr      = -0.906179845938664;       *wtPtr      = 0.236726885056189; 
        *(ptPtr+=lda)  = -0.538469310105683;       *(++wtPtr)  = 0.478628670499366; 
        *(ptPtr+=lda)  =  0.000000000000000;       *(++wtPtr)  = 0.568888888888889; 
        *(ptPtr+=lda)  =  0.538469310105683;       *(++wtPtr)  = 0.478628670499366; 
        *(ptPtr+=lda)  =  0.906179845938664;       *(++wtPtr)  = 0.236726885056189;
        break; 
        
        case 6:
        *ptPtr      =   0.932469514203152;   *wtPtr      = 0.171324492379170; 
        *(ptPtr+=lda)  =  -0.932469514203152;   *(++wtPtr)  = 0.171324492379170; 
        *(ptPtr+=lda)  =   0.661209386466265;   *(++wtPtr)  = 0.360761573048139;
        *(ptPtr+=lda)  =  -0.661209386466265;   *(++wtPtr)  = 0.360761573048139;
        *(ptPtr+=lda)  =   0.238619186003152;   *(++wtPtr)  = 0.467913934572691;
        *(ptPtr+=lda)  =  -0.238619186003152;   *(++wtPtr)  = 0.467913934572691;
        break; 
        
        case 7:
		*ptPtr      =   0.949107912342759;   *wtPtr      = 0.129484966168870;
        *(ptPtr+=lda)  =  -0.949107912342759;   *(++wtPtr)  = 0.129484966168870;
        *(ptPtr+=lda)  =   0.741531185599394;   *(++wtPtr)  = 0.279705391489277;
        *(ptPtr+=lda)  =  -0.741531185599394;   *(++wtPtr)  = 0.279705391489277;
		*(ptPtr+=lda)  =   0.405845151377397;   *(++wtPtr)  = 0.381830050505119;
		*(ptPtr+=lda)  =  -0.405845151377397;   *(++wtPtr)  = 0.381830050505119;
		*(ptPtr+=lda)  =   0.000000000000000;   *(++wtPtr)  = 0.417959183673469;
        break; 
        
        case 8:
        *ptPtr      =   0.960289856497536;   *wtPtr      = 0.101228536290376;
		*(ptPtr+=lda)  =  -0.960289856497536;   *(++wtPtr)  = 0.101228536290376;
		*(ptPtr+=lda)  =   0.796666477413627;   *(++wtPtr)  = 0.222381034453374;
		*(ptPtr+=lda)  =  -0.796666477413627;   *(++wtPtr)  = 0.222381034453374;
		*(ptPtr+=lda)  =   0.525532409916329;   *(++wtPtr)  = 0.313706645877887;
		*(ptPtr+=lda)  =  -0.525532409916329;   *(++wtPtr)  = 0.313706645877887;
		*(ptPtr+=lda)  =   0.183434642495650;   *(++wtPtr)  = 0.362683783378362;
		*(ptPtr+=lda)  =  -0.183434642495650;   *(++wtPtr)  = 0.362683783378362;
        break; 
    
    }
    
}

/**
Computes the nxn Gaussian quadrature rule
\param pts - pointer to the array that holds the quadarature points [s1, t1, s2, t2, s3, t3, ...] (must be allocated 2*n*n)
\param wts - pointer to the array of weights [w1 w2 w3 ..] (must be allocated n*n)
\param n - The number of points in one dimension (nxn rule) 
For now this only works for n = { 1, 2, 3, 4 }
\param lda -  The leading dimension of pts (default is 2)
**/
template <class NuMeRiC>
void quadrature_gauss2d( NuMeRiC *pts, NuMeRiC *wts, int n, int lda=2 )
{
    NuMeRiC *ptPtr=pts, *wtPtr=wts;
	int pi=lda-1; 
    
    switch (n) 
	{
			
		case 1:
		*ptPtr = 0.0;      *(++ptPtr) = 0.0;  *wtPtr = 4.0;
		break;

		case 2:
		
        *ptPtr     =  0.57735026918963;  *(++ptPtr) =  0.57735026918963;     *wtPtr = 1.0;
		*(ptPtr+=pi) =  0.57735026918963;  *(++ptPtr) = -0.57735026918963; *(++wtPtr) = 1.0;
		*(ptPtr+=pi) = -0.57735026918963;  *(++ptPtr) =  0.57735026918963; *(++wtPtr) = 1.0;
		*(ptPtr+=pi) = -0.57735026918963;  *(++ptPtr) = -0.57735026918963; *(++wtPtr) = 1.0;
		
//        *ptPtr     =  1;  *(++ptPtr) =  2;     *wtPtr = 1.0;
//		*(ptPtr+=pi) =  3;  *(++ptPtr) = 4; *(++wtPtr) = 1.0;
//		*(ptPtr+=pi) = 5;  *(++ptPtr) =  6; *(++wtPtr) = 1.0;
//		*(ptPtr+=pi) = 7;  *(++ptPtr) = 8; *(++wtPtr) = 1.0;
		
		break;

		case 3:
        *ptPtr = 0.77459666924148;  *(++ptPtr) =  0.77459666924148;     *wtPtr = 0.30864197530864;
        *(ptPtr+=pi) =  0.77459666924148;  *(++ptPtr) =  -0.77459666924148; *(++wtPtr) = 0.30864197530864;
        *(ptPtr+=pi) =  0.77459666924148;  *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.49382716049383;
        *(ptPtr+=pi) =  -0.77459666924148;  *(++ptPtr) =  0.77459666924148; *(++wtPtr) = 0.30864197530864;
        *(ptPtr+=pi) =  -0.77459666924148;  *(++ptPtr) =  -0.77459666924148; *(++wtPtr) = 0.30864197530864;
        *(ptPtr+=pi) =  -0.77459666924148;  *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.49382716049383;
        *(ptPtr+=pi) =  0.00000000000000;  *(++ptPtr) =  0.77459666924148; *(++wtPtr) = 0.49382716049383;
        *(ptPtr+=pi) =  0.00000000000000;  *(++ptPtr) =  -0.77459666924148; *(++wtPtr) = 0.49382716049383;
        *(ptPtr+=pi) =  0.00000000000000;  *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.79012345679012;
		break;

		case 4:
		*ptPtr = 0.86113431159405;  *(++ptPtr) =  0.86113431159405;     *wtPtr = 0.12100299328560;
		*(ptPtr+=pi) =  0.86113431159405;  *(++ptPtr) = -0.86113431159405; *(++wtPtr) = 0.12100299328560;
		*(ptPtr+=pi) =  0.86113431159405;  *(++ptPtr) = 0.33998104358486; *(++wtPtr) = 0.22685185185185;
		*(ptPtr+=pi) =  0.86113431159405;  *(++ptPtr) = -0.33998104358486; *(++wtPtr) = 0.22685185185185;
		*(ptPtr+=pi) =  -0.86113431159405;  *(++ptPtr) = 0.86113431159405; *(++wtPtr) = 0.12100299328560;
		*(ptPtr+=pi) =  -0.86113431159405;  *(++ptPtr) = -0.86113431159405; *(++wtPtr) = 0.12100299328560;
		*(ptPtr+=pi) =  -0.86113431159405;  *(++ptPtr) = 0.33998104358486; *(++wtPtr) = 0.22685185185185;
		*(ptPtr+=pi) =  -0.86113431159405;  *(++ptPtr) = -0.33998104358486; *(++wtPtr) = 0.22685185185185;
		*(ptPtr+=pi) =  0.33998104358486;  *(++ptPtr) = 0.86113431159405; *(++wtPtr) = 0.22685185185185;
		*(ptPtr+=pi) =  0.33998104358486;  *(++ptPtr) = -0.86113431159405; *(++wtPtr) = 0.22685185185185;
		*(ptPtr+=pi) =  0.33998104358486;  *(++ptPtr) = 0.33998104358486; *(++wtPtr) = 0.42529330301069;
		*(ptPtr+=pi) =  0.33998104358486;  *(++ptPtr) = -0.33998104358486; *(++wtPtr) = 0.42529330301069;
		*(ptPtr+=pi) =  -0.33998104358486;  *(++ptPtr) = 0.86113431159405; *(++wtPtr) = 0.22685185185185;
		*(ptPtr+=pi) =  -0.33998104358486;  *(++ptPtr) = -0.86113431159405; *(++wtPtr) = 0.22685185185185;
		*(ptPtr+=pi) =  -0.33998104358486;  *(++ptPtr) = 0.33998104358486; *(++wtPtr) = 0.42529330301069;
		*(ptPtr+=pi) =  -0.33998104358486;  *(++ptPtr) = -0.33998104358486; *(++wtPtr) = 0.42529330301069;
		break;
			
		default:
		Basso_Error("quadrature_gauss2d","too many points");
		break;
			
	}
    
}

/**
Computes the nxnxn Gaussian quadrature rule
\param pts - pointer to the array that holds the quadarature points [s1, t1, r1, s2, t2, r2, ...] (must be allocated 3*n*n)
\param wts - pointer to the array of weights [w1 w2 w3 ..] (must be allocated n*n)
\param n - The number of points in one dimension (nxn rule) 
\param lda -  The leading dimension of pts (default is 3)
For now this only works for n = { 1, 2, 3 }
**/
template <class NuMeRiC>
void quadrature_gauss3d( NuMeRiC *pts, NuMeRiC *wts, int n, int lda=3 )
{
    NuMeRiC *ptPtr=pts, *wtPtr=wts;
	int pi=lda-2;
	
    switch (n) 
	{
		case 1:
		*ptPtr = 0.0; *(++ptPtr) = 0.0; *(++ptPtr) = 0.0; *wtPtr = 8.0;
		break;

		case 2:
		*ptPtr =      0.57735026918963; *(++ptPtr) =  0.57735026918963; *(++ptPtr) =  0.57735026918963; *wtPtr     = 1.0;
		*(ptPtr+=pi) =  0.57735026918963; *(++ptPtr) =  0.57735026918963; *(++ptPtr) = -0.57735026918963; *(++wtPtr) = 1.0;
		*(ptPtr+=pi) =  0.57735026918963; *(++ptPtr) = -0.57735026918963; *(++ptPtr) =  0.57735026918963; *(++wtPtr) = 1.0;
		*(ptPtr+=pi) =  0.57735026918963; *(++ptPtr) = -0.57735026918963; *(++ptPtr) = -0.57735026918963; *(++wtPtr) = 1.0;
		*(ptPtr+=pi) = -0.57735026918963; *(++ptPtr) =  0.57735026918963; *(++ptPtr) =  0.57735026918963; *(++wtPtr) = 1.0;
		*(ptPtr+=pi) = -0.57735026918963; *(++ptPtr) =  0.57735026918963; *(++ptPtr) = -0.57735026918963; *(++wtPtr) = 1.0;
		*(ptPtr+=pi) = -0.57735026918963; *(++ptPtr) = -0.57735026918963; *(++ptPtr) =  0.57735026918963; *(++wtPtr) = 1.0;
		*(ptPtr+=pi) = -0.57735026918963; *(++ptPtr) = -0.57735026918963; *(++ptPtr) = -0.57735026918963; *(++wtPtr) = 1.0;
		break;

		case 3:
		*ptPtr =      0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++ptPtr) =  0.77459666924148; *wtPtr     = 0.17146776406036;
		*(ptPtr+=pi) =  0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++wtPtr) = 0.17146776406036;
		*(ptPtr+=pi) =  0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.27434842249657;
		*(ptPtr+=pi) =  0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++wtPtr) = 0.17146776406036;
		*(ptPtr+=pi) =  0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++wtPtr) = 0.17146776406036;
		*(ptPtr+=pi) =  0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.27434842249657;
		*(ptPtr+=pi) =  0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++ptPtr) =  0.77459666924148; *(++wtPtr) = 0.27434842249657;
		*(ptPtr+=pi) =  0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++ptPtr) = -0.77459666924148; *(++wtPtr) = 0.27434842249657;
		*(ptPtr+=pi) =  0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.43895747599451;
		*(ptPtr+=pi) = -0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++wtPtr) = 0.17146776406036;
	    *(ptPtr+=pi) = -0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++wtPtr) = 0.17146776406036;
	    *(ptPtr+=pi) = -0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.27434842249657;
	    *(ptPtr+=pi) = -0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++wtPtr) = 0.17146776406036;
	    *(ptPtr+=pi) = -0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++wtPtr) = 0.17146776406036;
	    *(ptPtr+=pi) = -0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.27434842249657;
	    *(ptPtr+=pi) = -0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++ptPtr) =  0.77459666924148; *(++wtPtr) = 0.27434842249657;
	    *(ptPtr+=pi) = -0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++ptPtr) = -0.77459666924148; *(++wtPtr) = 0.27434842249657;
	    *(ptPtr+=pi) = -0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.43895747599451;
	    *(ptPtr+=pi) =  0.00000000000000; *(++ptPtr) =  0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++wtPtr) = 0.27434842249657;
	    *(ptPtr+=pi) =  0.00000000000000; *(++ptPtr) =  0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++wtPtr) = 0.27434842249657;
	    *(ptPtr+=pi) =  0.00000000000000; *(++ptPtr) =  0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.43895747599451;
		*(ptPtr+=pi) =  0.00000000000000; *(++ptPtr) = -0.77459666924148; *(++ptPtr) =  0.77459666924148; *(++wtPtr) = 0.27434842249657;
		*(ptPtr+=pi) =  0.00000000000000; *(++ptPtr) = -0.77459666924148; *(++ptPtr) = -0.77459666924148; *(++wtPtr) = 0.27434842249657;
		*(ptPtr+=pi) =  0.00000000000000; *(++ptPtr) = -0.77459666924148; *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.43895747599451;
		*(ptPtr+=pi) =  0.00000000000000; *(++ptPtr) =  0.00000000000000; *(++ptPtr) =  0.77459666924148; *(++wtPtr) = 0.43895747599451;
		*(ptPtr+=pi) =  0.00000000000000; *(++ptPtr) =  0.00000000000000; *(++ptPtr) = -0.77459666924148; *(++wtPtr) = 0.43895747599451;
		*(ptPtr+=pi) =  0.00000000000000; *(++ptPtr) =  0.00000000000000; *(++ptPtr) =  0.00000000000000; *(++wtPtr) = 0.70233196159122;
		break;
			
		default:
		Basso_Error("quadrature_gauss3d","too many points");
		break;
			
	}
    
}


/** 
	computes 2D tensorial product quadrature rules from 1D rules.
	\param q1 input rule
	\param q2 second input rule
	\param q3 output tensorial rule 
*/
template <class NuMeRiC>
void tensorial_quadrature_rule( const NuMeRiC *pts1, NuMeRiC *wts1, int npts1,
        const NuMeRiC *pts2, NuMeRiC *wts2, int npts2,
        NuMeRiC *pts3, NuMeRiC *wts3, int ldwts3 )
{

    const NuMeRiC *p1=pts1, *w1=wts1, *p2, *w2;
    NuMeRiC *p3=pts3, *w3=wts3;

	int n=0, i, j;
    for ( i=0; i<npts1; ++i ) 
	{
        p2=pts2; w2=wts2;
        for ( j=0; j<npts2; ++j ) 
		{
            *p3 = *p1;  ++p3;              // xi
            *p3 = *p2;  p3 += ldwts3-2;    // eta
            
            *w3 = (*w1)*(*w2);  ++w3;
           
            ++w2; ++p2;
		}
        ++w1; ++p1;
 	}
}
/** 
	computes 3D tensorial product quadrature rules from 1D rules.
	\param q1 input rule
	\param q2 second input rule
	\param q3 third input rule
	\param q4 output tensorial rule 
*/
template <class NuMeRiC>
void tensorial_quadrature_rule( const NuMeRiC *pts1, NuMeRiC *wts1, int npts1,
        const NuMeRiC *pts2, NuMeRiC *wts2, int npts2,
		const NuMeRiC *pts3, NuMeRiC *wts3, int npts3,
        NuMeRiC *pts4, NuMeRiC *wts4, int ldwts4 )
{

    const NuMeRiC *p1=pts1, *w1=wts1, *p2, *w2, *p3=pts3, *w3=wts3;
    NuMeRiC *p=pts3, *w=wts3;

	int n=0, i, j, k;
    for ( i=0; i<npts1; ++i ) 
	{
        p2=pts2; w2=wts2;
        for ( j=0; j<npts2; ++j ) 
		{
        	p3=pts3; w3=wts3;
        	for ( k=0; k<npts3; ++k ) 
			{
            	*p = *p1;  ++p;              // xi
            	*p = *p2;  ++p;              // eta
            	*p = *p3;  p += ldwts4-3;    // zeta
            
            	*w = (*w1)*(*w2)*(*w3);  ++w;
           
            	++w3; ++p3;
			}
       		++w2; ++p2;
		}
       	++w1; ++p1;
 	}
}

/** 
	computes 2D or 3D tensorial product quadrature rules from 1D rules.
	\param q1 input rule
	\param q output tensorial rule 
*/
template <class NuMeRiC>
void tensorial_quadrature_rule( const NuMeRiC *pts1, NuMeRiC *wts1, int npts1,
        NuMeRiC *pts, NuMeRiC *wts, int ldp, int sdim )
{
	if ( sdim==1 )
	{
		const NuMeRiC *p1=pts1, *w1=wts1;
		NuMeRiC *p=pts, *w=wts;
		for ( int i=0; i<npts1; ++i )
		{
			*p = *p1; p+=ldp-1; ++p1;
			*w = *w1; ++w; ++w1;
		}
	}
	else if ( sdim==2 )
		tensorial_quadrature_rule( pts1, wts1, npts1, pts1, wts1, npts1,
				pts, wts, ldp );
	
	else if ( sdim==3 )
		tensorial_quadrature_rule( pts1, wts1, npts1, pts1, wts1, npts1,
		 	pts1, wts1, npts1, pts, wts, ldp );

}


/** 
Computes quadrature rules for triangular domains 
	npts = 1, 3, 4, 5, 6, 7, 9, 12, and 13
\param pts - pointer to the array that holds the quadarature points [s1, t1, .. s2, t2, .., ...] (must be allocated 3*n*n)
\param wts - pointer to the array of weights [w1 w2 w3 ..] (must be allocated n*n)
\param npts - The number of pointsin the rule
\param ldp - leading dimension of pts (defalut is 2)
The return value is the order of the integration
*/
template <class NuMeRiC>	
int quadrature_tria( NuMeRiC *pts, NuMeRiC *wts, int npts, int ldp=2  )  
{
	NuMeRiC *ptPtr=pts, *wtPtr=wts;
	int ord, pi=ldp-1;
	switch (npts)
	{
		case 1:
		ord=1;
		*ptPtr = 0.333333333333; *(++ptPtr) = 0.333333333333; *wtPtr = 0.5; 
		break;

		case 3:
		ord=2;
		*ptPtr = 0.1666666666667; *(++ptPtr) = 0.1666666666667; *wtPtr     = 0.1666666666667; ptPtr += pi;
		*ptPtr = 0.6666666666667; *(++ptPtr) = 0.1666666666667; *(++wtPtr) = 0.1666666666667; ptPtr += pi;
		*ptPtr = 0.1666666666667; *(++ptPtr) = 0.6666666666667; *(++wtPtr) = 0.1666666666667;
		break;

		case 4:
		ord=3;
		*ptPtr = 0.333333333333333; *(++ptPtr) = 0.333333333333333; *wtPtr     = -0.281250000000; ptPtr += pi;
		*ptPtr = 0.600000000000000; *(++ptPtr) = 0.200000000000000; *(++wtPtr) = 0.2604166666667; ptPtr += pi;
		*ptPtr = 0.200000000000000; *(++ptPtr) = 0.600000000000000; *(++wtPtr) = 0.2604166666667; ptPtr += pi;
		*ptPtr = 0.200000000000000; *(++ptPtr) = 0.200000000000000; *(++wtPtr) = 0.2604166666667;
		break;

		case 6:
		ord= 4;
		*ptPtr = 0.816847572980459; *(++ptPtr) = 0.091576213509770; *wtPtr     = 0.054975871827616; ptPtr += pi;
		*ptPtr = 0.091576213509770; *(++ptPtr) = 0.816847572980459; *(++wtPtr) = 0.054975871827616; ptPtr += pi;
		*ptPtr = 0.091576213509770; *(++ptPtr) = 0.091576213509770; *(++wtPtr) = 0.054975871827616; ptPtr += pi;
		*ptPtr = 0.108103018168070; *(++ptPtr) = 0.445948490915965; *(++wtPtr) = 0.111690794839005; ptPtr += pi;
		*ptPtr = 0.445948490915965; *(++ptPtr) = 0.108103018168070; *(++wtPtr) = 0.111690794839005; ptPtr += pi;
		*ptPtr = 0.445948490915965; *(++ptPtr) = 0.445948490915965; *(++wtPtr) = 0.111690794839005;
		break;

		case 7:
		ord=5;
		*ptPtr = 0.1012865073235; *(++ptPtr) = 0.1012865073235; *wtPtr     = 0.06296959027240; ptPtr += pi;
		*ptPtr = 0.7974269853531; *(++ptPtr) = 0.1012865073235; *(++wtPtr) = 0.06296959027240; ptPtr += pi;
		*ptPtr = 0.1012865073235; *(++ptPtr) = 0.7974269853531; *(++wtPtr) = 0.06296959027240; ptPtr += pi;
		*ptPtr = 0.4701420641051; *(++ptPtr) = 0.0597158717898; *(++wtPtr) = 0.06619707639425; ptPtr += pi;
		*ptPtr = 0.4701420641051; *(++ptPtr) = 0.4701420641051; *(++wtPtr) = 0.06619707639425; ptPtr += pi;
		*ptPtr = 0.0597158717898; *(++ptPtr) = 0.4701420641051; *(++wtPtr) = 0.06619707639425; ptPtr += pi;
		*ptPtr = 0.3333333333333; *(++ptPtr) = 0.3333333333333; *(++wtPtr) = 0.11250000000000;
		break;

		case 9:
		ord=5;
		*ptPtr = 0.124949503233232; *(++ptPtr) = 0.437525248383384; *wtPtr     = 0.102975252380443; ptPtr += pi;
		*ptPtr = 0.437525248383384; *(++ptPtr) = 0.124949503233232; *(++wtPtr) = 0.102975252380443; ptPtr += pi;
		*ptPtr = 0.437525248383384; *(++ptPtr) = 0.437525248383384; *(++wtPtr) = 0.102975252380443; ptPtr += pi;
		*ptPtr = 0.797112651860071; *(++ptPtr) = 0.165409927389841; *(++wtPtr) = 0.031845707143112; ptPtr += pi;
		*ptPtr = 0.797112651860071; *(++ptPtr) = 0.037477420750088; *(++wtPtr) = 0.031845707143112; ptPtr += pi;
		*ptPtr = 0.165409927389841; *(++ptPtr) = 0.797112651860071; *(++wtPtr) = 0.031845707143112; ptPtr += pi;
		*ptPtr = 0.165409927389841; *(++ptPtr) = 0.037477420750088; *(++wtPtr) = 0.031845707143112; ptPtr += pi;
		*ptPtr = 0.037477420750088; *(++ptPtr) = 0.797112651860071; *(++wtPtr) = 0.031845707143112; ptPtr += pi;
		*ptPtr = 0.037477420750088; *(++ptPtr) = 0.165409927389841; *(++wtPtr) = 0.031845707143112;
		break;

		case 12:
		ord=6;
		*ptPtr = 0.873821971016996; *(++ptPtr) = 0.063089014491502; *wtPtr     = 0.025422453185103; ptPtr += pi;
		*ptPtr = 0.063089014491502; *(++ptPtr) = 0.873821971016996; *(++wtPtr) = 0.025422453185103; ptPtr += pi;
		*ptPtr = 0.063089014491502; *(++ptPtr) = 0.063089014491502; *(++wtPtr) = 0.025422453185103; ptPtr += pi;
		*ptPtr = 0.501426509658179; *(++ptPtr) = 0.249286745170910; *(++wtPtr) = 0.058393137863189; ptPtr += pi;
		*ptPtr = 0.249286745170910; *(++ptPtr) = 0.501426509658179; *(++wtPtr) = 0.058393137863189; ptPtr += pi;
		*ptPtr = 0.249286745170910; *(++ptPtr) = 0.249286745170910; *(++wtPtr) = 0.058393137863189; ptPtr += pi;
		
		*ptPtr = 0.636502499121399; *(++ptPtr) = 0.310352451033785; *(++wtPtr) = 0.041290537809187; ptPtr += pi;
		*ptPtr = 0.636502499121399; *(++ptPtr) = 0.053145049844816; *(++wtPtr) = 0.041290537809187; ptPtr += pi;
		*ptPtr = 0.310352451033785; *(++ptPtr) = 0.636502499121399; *(++wtPtr) = 0.041290537809187; ptPtr += pi;
		*ptPtr = 0.310352451033785; *(++ptPtr) = 0.053145049844816; *(++wtPtr) = 0.041290537809187; ptPtr += pi;
		*ptPtr = 0.053145049844816; *(++ptPtr) = 0.636502499121399; *(++wtPtr) = 0.041290537809187; ptPtr += pi;
		*ptPtr = 0.053145049844816; *(++ptPtr) = 0.310352451033785; *(++wtPtr) = 0.041290537809187;
		break;

		case 13:
		ord=7;
		
		*ptPtr = 0.0651301029022; *(++ptPtr) = 0.0651301029022; *wtPtr     = 0.02667361780440; ptPtr += pi;
		*ptPtr = 0.8697397941956; *(++ptPtr) = 0.0651301029022; *(++wtPtr) = 0.02667361780440; ptPtr += pi;
		*ptPtr = 0.0651301029022; *(++ptPtr) = 0.8697397941956; *(++wtPtr) = 0.02667361780440; ptPtr += pi;		
		*ptPtr = 0.3128654960049; *(++ptPtr) = 0.0486903154253; *(++wtPtr) = 0.03855688044515; ptPtr += pi;
		*ptPtr = 0.6384441885698; *(++ptPtr) = 0.3128654960049; *(++wtPtr) = 0.03855688044515; ptPtr += pi;
		*ptPtr = 0.0486903154253; *(++ptPtr) = 0.6384441885698; *(++wtPtr) = 0.03855688044515; ptPtr += pi;
		
		*ptPtr = 0.6384441885698; *(++ptPtr) = 0.0486903154253; *(++wtPtr) = 0.03855688044515; ptPtr += pi;
		*ptPtr = 0.3128654960049; *(++ptPtr) = 0.6384441885698; *(++wtPtr) = 0.03855688044515; ptPtr += pi;
		*ptPtr = 0.0486903154253; *(++ptPtr) = 0.3128654960049; *(++wtPtr) = 0.03855688044515; ptPtr += pi;
		
		*ptPtr = 0.2603459660790; *(++ptPtr) = 0.2603459660790; *(++wtPtr) = 0.08780762881660; ptPtr += pi;
		*ptPtr = 0.4793080678419; *(++ptPtr) = 0.2603459660790; *(++wtPtr) = 0.08780762881660; ptPtr += pi;
		*ptPtr = 0.2603459660790; *(++ptPtr) = 0.4793080678419; *(++wtPtr) = 0.08780762881660; ptPtr += pi;
		
		*ptPtr = 0.3333333333333; *(++ptPtr) = 0.3333333333333; *(++wtPtr) = -0.07478502223385;
		break;

		default:
		ord=0;
		break;

	}
	return ord;
}

/**
	Computes quadrature rules tor tetrahedra
	
	npts = 1, 4, or 5
	\param pts - pointer to the array that holds the quadarature points [s1, t1, r1, .. s2, t2, r2, .., ...] (must be allocated 3*n*n)
    \param wts - pointer to the array of weights [w1 w2 w3 ..] (must be allocated n*n)
    \param npts - The number of pointsin the rule
    \param ldp - leading dimension of pts (defalut is 3)
    The return value is the order of the integration
*/	
template <class NuMeRiC>
int quadrature_tetra( NuMeRiC *pts, NuMeRiC *wts, int npts, int ldp=3 ) 
{
	NuMeRiC *ptPtr=pts, *wtPtr=wts;
	int ord, pi=ldp-2;
	switch (npts)
	{
		case 1:
		ord=1;
		*ptPtr = 0.25000000000000; *(++ptPtr) = 0.25000000000000; *(++ptPtr) = 0.25000000000000; *wtPtr = 0.16666666666667; 
		break;

		case 4:
		ord=2;
		*ptPtr = 0.58541020000000; *(++ptPtr) = 0.13819660000000; *(++ptPtr) = 0.13819660000000; *wtPtr     = 0.04166666666667; ptPtr += pi;
		*ptPtr = 0.13819660000000; *(++ptPtr) = 0.58541020000000; *(++ptPtr) = 0.13819660000000; *(++wtPtr) = 0.04166666666667; ptPtr += pi;
		*ptPtr = 0.13819660000000; *(++ptPtr) = 0.13819660000000; *(++ptPtr) = 0.58541020000000; *(++wtPtr) = 0.04166666666667; ptPtr += pi;
		*ptPtr = 0.13819660000000; *(++ptPtr) = 0.13819660000000; *(++ptPtr) = 0.13819660000000; *(++wtPtr) = 0.04166666666667;
		break;

		case 5:
		ord=3;
		*ptPtr = 0.25000000000000; *(++ptPtr) = 0.25000000000000; *(++ptPtr) = 0.25000000000000; *wtPtr     =-0.13333333333333; ptPtr += pi;
		*ptPtr = 0.50000000000000; *(++ptPtr) = 0.16666666666667; *(++ptPtr) = 0.16666666666667; *(++wtPtr) = 0.07500000000000; ptPtr += pi;
		*ptPtr = 0.16666666666667; *(++ptPtr) = 0.50000000000000; *(++ptPtr) = 0.16666666666667; *(++wtPtr) = 0.07500000000000; ptPtr += pi;
		*ptPtr = 0.16666666666667; *(++ptPtr) = 0.16666666666667; *(++ptPtr) = 0.50000000000000; *(++wtPtr) = 0.07500000000000; ptPtr += pi;
		*ptPtr = 0.16666666666667; *(++ptPtr) = 0.16666666666667; *(++ptPtr) = 0.16666666666667; *(++wtPtr) = 0.07500000000000;
		break;

		default:
		ord=0;
		break;

	}
	return ord;
}

// ******************************************************************************************
//       S O M E      B A S I C    E L E M E N T    S T I F F N E S S     M A T R I X     
//                   C O M P U T A T I O N     F U N C T I O N S 

/**
Computes the stiffness matrix for a plane stress three node triangular 
element with linear isotropic elastic material
\param coord - The node coordiante matrix [x1 y1... x2 y2 .. x3 y3 ..]
\parm ldc - The leading dimension of the coordinate matrix. Must be at least 2.
\param E - YOung's modulus
\param nu - Poission's ratio
\param thk - the thickness of the plate
\param ke - on return has the stiffness matrix
\param ldk - the leading dimension of the stiffness matrix
*/
template <class NuMeRiC>
void kmat_t3( const NuMeRiC *coord, const int &ldc, 
const NuMeRiC &E, const NuMeRiC &nu, const NuMeRiC &thk, NuMeRiC *ke, const int &ldk )
{
	NuMeRiC x1=coord[0],     y1=coord[1],
	        x2=coord[ldc],   y2=coord[1+ldc],
			x3=coord[2*ldc], y3=coord[1+2*ldc];
	NuMeRiC detJ = x1*y2-x1*y3-x2*y1+x2*y3+x3*y1-x3*y2;
	NuMeRiC EE = E/(1-nu*nu), nm1=0.5*(nu-1.0), np1=0.5*(nu+1.0), alpha=0.5*EE*thk/detJ;
	NuMeRiC y23=y2-y3, x32=x3-x2, y31=y3-y1, x13=x1-x3, y12=y1-y2, x21=x2-x1;
		
	NuMeRiC *kptr = ke;	
	*kptr = alpha*(y23*y23-nm1*x32*x32);  	++kptr; 		// k(0,0)
	*kptr = alpha*(np1*x32*y23); 			++kptr; 		// k(1,0)
	*kptr = alpha*(y31*y23-nm1*x13*x32);	++kptr; 		// k(2,0)
	*kptr = alpha*(nu*x13*y23-nm1*x32*y31);	++kptr; 		// k(3,0)
	*kptr = alpha*(y12*y23-nm1*x21*x32);	++kptr; 		// k(4,0)
	*kptr = alpha*(nu*x21*y23-nm1*x32*y12);	kptr += ldk-5; 	// k(5,0)
	
	*kptr = ke[1]; 							++kptr; 		// k(0,1)
	*kptr = alpha*(x32*x32-nm1*y23*y23); 	++kptr; 		// k(1,1)
	*kptr = alpha*(nu*x32*y31-nm1*x13*y23);	++kptr; 		// k(2,1)
	*kptr = alpha*(x13*x32-nm1*y31*y23);	++kptr; 		// k(3,1)
	*kptr = alpha*(nu*x32*y12-nm1*x21*y23);	++kptr; 		// k(4,1)
	*kptr = alpha*(x21*x32-nm1*y23*y12);	kptr += ldk-5; 	// k(5,1)
	
	*kptr = ke[2]; 							++kptr; 		// k(0,2)
	*kptr = ke[ldk+2]; 						++kptr; 		// k(1,2)
	*kptr = alpha*(y31*y31-nm1*x13*x13);	++kptr; 		// k(2,2)
	*kptr = alpha*(np1*x13*y31);			++kptr; 		// k(3,2)
	*kptr = alpha*(y12*y31-nm1*x21*x13);	++kptr; 		// k(4,2)
	*kptr = alpha*(nu*x21*y31-nm1*x13*y12);	kptr += ldk-5; 	// k(5,2)
	
	*kptr = ke[3]; 							++kptr; 		// k(0,3)
	*kptr = ke[ldk+3]; 						++kptr; 		// k(1,3)
	*kptr = ke[2*ldk+3];					++kptr; 		// k(2,3)
	*kptr = alpha*(x13*x13-nm1*y31*y31);	++kptr; 		// k(3,3)
	*kptr = alpha*(nu*x13*y12-nm1*x21*y31);	++kptr; 		// k(4,3)
	*kptr = alpha*(x21*x13-nm1*y31*y12);	kptr += ldk-5; 	// k(5,3)
	
	*kptr = ke[4]; 							++kptr; 		// k(0,4)
	*kptr = ke[ldk+4]; 						++kptr; 		// k(1,4)
	*kptr = ke[2*ldk+4];					++kptr; 		// k(2,4)
	*kptr = ke[3*ldk+4];					++kptr; 		// k(3,4)
	*kptr = alpha*(y12*y12-nm1*x21*x21);	++kptr; 		// k(4,4)
	*kptr = alpha*(np1*x21*y12);			kptr += ldk-5; 	// k(5,4)
	
	*kptr = ke[5]; 							++kptr; 		// k(0,5)
	*kptr = ke[ldk+5]; 						++kptr; 		// k(1,5)
	*kptr = ke[2*ldk+5];					++kptr; 		// k(2,5)
	*kptr = ke[3*ldk+5];					++kptr; 		// k(3,5)
	*kptr = ke[4*ldk+5];					++kptr; 		// k(4,5)
	*kptr = alpha*(x21*x21-nm1*y12*y12);	kptr += ldk-5; 	// k(5,5)
}

/**
Computes the stiffness matrix for a plane stress three node triangular 
element with linear isotropic elastic material
\param iblock - the number of elements passed
\param coord - The node coordiante matrix. Each column holds the nodal coordinates for that
		element.   [ x1 y1 z1 x2 y2 z2 x3 y3 z3; x1 y1 z1 x2 y2 z2 x3 y3 z3; .... ]^T
\parm ldc - The leading dimension of the coordinate matrix. Must be at least 9.
\param E - Young's modulus
\param nu - Poission's ratio
\param thk - the thickness of the plate
\param ke - on return has the lower half of each element stiffness matrix in a column
		[ k11 k21 k31 k41 k51 k61 k22 k32 k42 k52 k62 k33 k43 k53 k63 k44 k54 k64 k55 k65 k66 ]
\param ldk - the leading dimension of the stiffness matrix must be at least 21
*/
template <class NuMeRiC>
void kmat_t3( const int &nblock, const NuMeRiC *coord, const int &ldc, 
const NuMeRiC &E, const NuMeRiC &nu, const NuMeRiC &thk, NuMeRiC *ke, const int &ldk )
{
	NuMeRiC x1, y1, x2, y2, x3, y3, detJ;
	NuMeRiC EE = E/(1-nu*nu), nm1=0.5*(nu-1.0), np1=0.5*(nu+1.0), alpha=0.5*EE*thk/detJ;
	NuMeRiC y23, x32, y31, x13, y12, x21;
	
	NuMeRiC *kptr = ke;	
	for ( int k=0; k<nblock; ++k )
	{
		x1=coord[0+k*ldc];	y1=coord[1+k*ldc];
		x2=coord[4+k*ldc];	y2=coord[5+k*ldc];
		x3=coord[7+k*ldc];	y3=coord[8+k*ldc];
		detJ = x1*y2-x1*y3-x2*y1+x2*y3+x3*y1-x3*y2;
		y23=y2-y3; x32=x3-x2; y31=y3-y1; x13=x1-x3; y12=y1-y2; x21=x2-x1;
		
		*kptr = alpha*(y23*y23-nm1*x32*x32);  	++kptr; 		// k(0,0)
		*kptr = alpha*(np1*x32*y23); 			++kptr; 		// k(1,0)
		*kptr = alpha*(y31*y23-nm1*x13*x32);	++kptr; 		// k(2,0)
		*kptr = alpha*(nu*x13*y23-nm1*x32*y31);	++kptr; 		// k(3,0)
		*kptr = alpha*(y12*y23-nm1*x21*x32);	++kptr; 		// k(4,0)
		*kptr = alpha*(nu*x21*y23-nm1*x32*y12);	++kptr; 		// k(5,0)
	
		*kptr = alpha*(x32*x32-nm1*y23*y23); 	++kptr; 		// k(1,1)
		*kptr = alpha*(nu*x32*y31-nm1*x13*y23);	++kptr; 		// k(2,1)
		*kptr = alpha*(x13*x32-nm1*y31*y23);	++kptr; 		// k(3,1)
		*kptr = alpha*(nu*x32*y12-nm1*x21*y23);	++kptr; 		// k(4,1)
		*kptr = alpha*(x21*x32-nm1*y23*y12);	++kptr; 		// k(5,1)
	
		*kptr = alpha*(y31*y31-nm1*x13*x13);	++kptr; 		// k(2,2)
		*kptr = alpha*(np1*x13*y31);			++kptr; 		// k(3,2)
		*kptr = alpha*(y12*y31-nm1*x21*x13);	++kptr; 		// k(4,2)
		*kptr = alpha*(nu*x21*y31-nm1*x13*y12);	++kptr; 		// k(5,2)
	
		*kptr = alpha*(x13*x13-nm1*y31*y31);	++kptr; 		// k(3,3)
		*kptr = alpha*(nu*x13*y12-nm1*x21*y31);	++kptr; 		// k(4,3)
		*kptr = alpha*(x21*x13-nm1*y31*y12);	++kptr; 		// k(5,3)
	
		*kptr = alpha*(y12*y12-nm1*x21*x21);	++kptr; 		// k(4,4)
		*kptr = alpha*(np1*x21*y12);			++kptr; 		// k(5,4)
	
		*kptr = alpha*(x21*x21-nm1*y12*y12);	kptr += ldk-20; // k(5,5)
	}
}

} // end of namespace


#endif



