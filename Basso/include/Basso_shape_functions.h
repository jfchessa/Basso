/*! \file Basso_shape_functions.h


	\brief Functions for computing element shape functions and derivatives

	\author Jack Chessa, jfchessa@utep.edu
	\date November 23, 2007

*/

#ifndef _BASSO_SHAPE_FUNCTIONS_H_
#define _BASSO_SHAPE_FUNCTIONS_H_

namespace Basso
{

#ifndef ONETHIRD
#define ONETHIRD 0.33333333333333333333333333333333
#endif

#ifndef ONESIXTH
#define ONESIXTH 0.16666666666666666666666666666667
#endif

using namespace std;
	
	
	/****************************************************************************************************** 
	 ***                                                                                                ***
	 ***              Functions to various element shape functions and parent space gradients           ***
	 ***                                                                                                ***
	 ******************************************************************************************************/
/** 
	The gradient matrices are returned in the follwoing format
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
*/


/** 
	Computes the shape function for a two node line element
		\param Na the array where the shape function values are returned
		\param inc increment in the Na vector where the shape function values will be places (default=1)
		\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
		
		\image html line2.pdf
		
*/
template < class NuMeRiC >
int shape_line2( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0 )
{	
	NuMeRiC *ptr=Na;
	*ptr 		= 0.5*(1.0-xi);  ptr+=inc;
	*ptr 		= 0.5*(1.0+xi);
	return 2;
}

/**
Computes the gradient of the shape function for a two node line element w.r.t to the parent
coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
		\image html line2.pdf
*/
template < class NuMeRiC >
int dshape_line2( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0 )
{
	NuMeRiC *ptr=dNa;
	*(ptr++)	= -0.5;  	 
	*(ptr)	 	=  0.5; 
	return 2;
}

/** 
	Computes the shape function for a three node line element
		\param Na the array where the shape function values are returned
		\param inc increment in the Na vector where the shape function values will be places (default=1)
		\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]

		\image html line3.pdf		
*/
template < class NuMeRiC >
int shape_line3( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0 )
{	
	Na[0]=0.5*(xi*xi-xi);
	Na[inc]=0.5*(xi*xi+xi);	
	Na[2*inc]=1.0-xi*xi;
	return 3;
}

/**
Computes the gradient of the shape function for a three node line element w.r.t to the parent
coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
	\param edim (optional) on returns contains the dimension of the element in the parament space
	
		\image html line3.pdf

*/
template < class NuMeRiC >
int dshape_line3( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0   )
{	
	dNa[0]=xi-0.5; 
	dNa[1]=xi+0.5; 
	dNa[2]=-2*xi;
	return 3;
}

/** 
	Computes the shape function for a three node triangle element
		\param Na the array where the shape function values are returned
		\param inc increment in the Na vector where the shape function values will be places (default=1)
		\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[0,1]
		\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[0,1-xi]

		\image html tria3.pdf
		
*/
template < class NuMeRiC >
int shape_tria3( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD )
{
	NuMeRiC *ptr=Na;
	*ptr 	= 1.0-xi-eta;  	ptr+=inc; 
	*ptr 	= xi;  			ptr+=inc; 	 
	*ptr 	= eta;
	return 3;
}

/**
Computes the gradient of the shape function for a three node triangle element w.r.t to the parent
coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[0,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[0,1-xi]
	\param edim (optional) on returns contains the dimension of the element in the parament space

		\image html tria3.pdf
*/
template < class NuMeRiC >
int dshape_tria3( NuMeRiC *dNa, int lda, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD )
{	
	NuMeRiC *ptr=dNa;
	*(ptr++)	= -1.0;  	*(ptr++)	=  1.0; 	*(ptr++) 	=  0.0;
    ptr += lda-3;
	*(ptr++)	= -1.0;  	*(ptr++)	=  0.0; 	*(ptr)   	=  1.0;
	return 3;
}

/**
		Computes the gradient of the shape functions w.r.t. to the spacial coordinates:
			dNa = [ dN1dx dN2dx dN3dx dN1dy dN2dy dN3dy ]
		
		The return value is the area of the triangle	
			\image html tria3.pdf
*/
template< class NuMeRiC >
NuMeRiC gradshape_tria3( NuMeRiC *dNa, int lda, NuMeRiC *x1, NuMeRiC *x2, NuMeRiC *x3, int sdim=2 )
{
	NuMeRiC *ptr=dNa;
	NuMeRiC invarea2;
 	if ( sdim==2 )
	{
		invarea2 = 1.0 / ( x2[0]*x3[1]-x2[1]*x3[0] + x3[0]*x1[1]-x3[1]*x1[0]
				+ x1[0]*x2[1]-x1[1]*x2[0] );
		
		*(ptr++) = invarea2*( x2[1] - x3[1] );  	
		*(ptr++) = invarea2*( x3[1] - x1[1] );  	
		*(ptr++) = invarea2*( x1[1] - x2[1] );  
		ptr += lda-3;	
		*(ptr++) = invarea2*( x3[0] - x2[0] );
		*(ptr++) = invarea2*( x1[0] - x3[0] );
		*(ptr)   = invarea2*( x2[0] - x1[0] );
	}
	else if ( sdim==3 )
	{
		cout << "gradshape_tria3 with sdim=3 not yet implemented\n";
	}
	return 0.5/invarea2;
}


/** 
	Computes the shape function for a four node triangle element (cubic bubble)
		\param Na the array where the shape function values are returned
		\param inc increment in the Na vector where the shape function values will be places (default=1)
		\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[0,1]
		\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[0,1-xi]

		\image html tria4.pdf

*/
template < class NuMeRiC >
int shape_tria4( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD )
{
	Na[0]=1-3*(xi+eta)+4*xi*eta+2*(xi*xi+eta*eta);	
	Na[inc]=xi*(2*xi-1);
	Na[2*inc]=eta*(2*eta-1);
	Na[3*inc]=4*xi*(1-xi-eta);
	return 4;
}

/**
Computes the gradient of the shape function for a four node triangle element (cubic bubble) 
w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[0,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[0,1-xi]
	\param edim (optional) on returns contains the dimension of the element in the parament space

		\image html tria4.pdf

*/
template < class NuMeRiC >
int dshape_tria4( NuMeRiC *dNa, int lda, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD  )
{
	dNa[0]=-1-3*eta; dNa[0+lda]=-1-3*xi;
	dNa[1]=1-3*eta;  dNa[1+lda]=-3*xi;
	dNa[2]=-3*eta;   dNa[2+lda]=1-3*xi;
	dNa[3]=9*eta;    dNa[3+lda]= 9*xi;
	return 4;
}

/** 
	Computes the shape function for a six node triangle element
		\param Na the array where the shape function values are returned
		\param inc increment in the Na vector where the shape function values will be places (default=1)
		\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[0,1]
		\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[0,1-xi]

		\image html tria6.pdf
	
*/
template < class NuMeRiC >
int shape_tria6( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD )
{	
	Na[0]=1-3*(xi+eta)+4*xi*eta+2*(xi*xi+eta*eta);	
	Na[inc]=xi*(2*xi-1);
	Na[2*inc]=eta*(2*eta-1);
	Na[3*inc]=4*xi*(1-xi-eta);	
	Na[4*inc]=4*xi*eta;
	Na[5*inc]=4*eta*(1-xi-eta);
	return 6;
}

/**
Computes the gradient of the shape function for a six node triangle element 
w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[0,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[0,1-xi]
	\param edim (optional) on returns contains the dimension of the element in the parament space

		\image html tria6.pdf
*/
template < class NuMeRiC >
int dshape_tria6( NuMeRiC *dNa, int lda, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD  )
{	
	dNa[0]=4*(xi+eta)-3;   dNa[0+lda]=4*(xi+eta)-3;
	dNa[1]=4*xi-1;         dNa[1+lda]=0.0; 
	dNa[2]=0.0;            dNa[2+lda]=4*eta-1;
	dNa[3]=4*(1-eta-2*xi); dNa[3+lda]=-4*xi;
	dNa[4]=4*eta;          dNa[4+lda]=4*xi;
	dNa[5]=-4*eta;         dNa[5+lda]=4*(1-xi-2*eta);
	return 6;
}

/** 
	Computes the shape function for a four node quadrilateral element
		\param Na the array where the shape function values are returned
		\param inc increment in the Na vector where the shape function values will be places (default=1)
		\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
		\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[-1,1]


		\image html quad4.pdf	
*/
template < class NuMeRiC >
int shape_quad4( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0 )
{	
	Na[0]=0.25*(1-xi)*(1-eta);	
	Na[inc]=0.25*(1+xi)*(1-eta);
	Na[2*inc]=0.25*(1+xi)*(1+eta);
	Na[3*inc]=0.25*(1-xi)*(1+eta);
	return 4;
}

/**
	Computes the gradient of the shape function for a four node quadrilateral element
		w.r.t to the parent coordinate system
		\param dNa the array where the shape function gradient values are returned
		\param lda the leading dimension of dNa
		\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
		\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[-1,1]
		\param edim (optional) on returns contains the dimension of the element in the parament space

	\image html quad4.pdf
		
*/
template < class NuMeRiC >
int dshape_quad4( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0 )
{	
	dNa[0]=-0.25*(1-eta);	dNa[0+lda]=-0.25*(1-xi);
	dNa[1]= 0.25*(1-eta);	dNa[1+lda]=-0.25*(1+xi);
	dNa[2]= 0.25*(1+eta);	dNa[2+lda]= 0.25*(1+xi);
  	dNa[3]=-0.25*(1+eta);	dNa[3+lda]= 0.25*(1-xi);
	return 4;
}

/** 
		Computes the shape function for a eight node quadrilateral element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[-1,1]

		\image html quad8.pdf	
*/
template < class NuMeRiC >
int shape_quad8( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0 )
{
	NuMeRiC r=0.5*(xi+1.0), s=0.5*(eta+1.0);

	Na[0] =   ( r - 1.0 ) * ( s - 1.0 ) * ( 1.0 - 2.0 * r - 2.0 * s );
	Na[inc] =   r * ( s - 1.0 ) * ( 1.0 - 2.0 * r + 2.0 * s );
	Na[2*inc] =   r * s * ( 2.0 * r + 2.0 * s - 3.0 );
	Na[3*inc] =   ( r - 1.0 ) * s * ( 2.0 * r - 2.0 * s + 1.0 );
	Na[4*inc] =   4.0 * r * ( r - 1.0 ) * ( s - 1.0 );
	Na[5*inc] = - 4.0 * r * s * ( s - 1.0 );
	Na[6*inc] = - 4.0 * r * ( r - 1.0 ) * s;
	Na[7*inc] =   4.0 * ( r - 1.0 ) * s * ( s - 1.0 );
	return 8;
}

	/**
	Computes the gradient of the shape function for a eight node quadrilateral element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[-1,1]
	\param edim (optional) on returns contains the dimension of the element in the parament space

		\image html quad8.pdf
	*/
template < class NuMeRiC >
	int dshape_quad8( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0 )
	{
		
		NuMeRiC r=0.5*(xi+1.0), s=0.5*(eta+1.0);

		dNa[0] = ( s - 1.0  ) * ( - 4.0  * r - 2.0  * s + 3.0  );
		dNa[1] = ( s - 1.0  ) * ( - 4.0  * r + 2.0  * s + 1.0  );
		dNa[2] =   s         * (   4.0  * r + 2.0  * s - 3.0  );
		dNa[3] =   s         * (   4.0  * r - 2.0  * s - 1.0  );
		dNa[4] =   4.0  * ( 2.0  * r - 1.0  )     * ( s - 1.0  );
		dNa[5] = - 4.0  *                     s * ( s - 1.0  );
		dNa[6] = - 4.0  * ( 2.0  * r - 1.0  ) * s;
		dNa[7] =   4.0  *                     s * ( s - 1.0  );

		dNa[0+lda] = ( r - 1.0  ) * ( - 4.0  * s - 2.0  * r + 3.0  );
		dNa[1+lda] =   r *       (   4.0  * s - 2.0  * r - 1.0  );
		dNa[2+lda] =   r *       (   4.0  * s + 2.0  * r - 3.0  );
		dNa[3+lda] = ( r - 1.0  ) * ( - 4.0  * s + 2.0  * r + 1.0  );
		dNa[4+lda] =   4.0  * r * ( r - 1.0  );
		dNa[5+lda] = - 4.0  * r               * ( 2.0  * s - 1.0  );
		dNa[6+lda] = - 4.0  * r * ( r - 1.0  );
		dNa[7+lda] =   4.0  *     ( r - 1.0  ) * ( 2.0  * s - 1.0  );

		return 8;
	}

	/** 
		Computes the shape function for a nine node quadrilateral element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[-1,1]

		\image html quad9.pdf
	
	*/
template < class NuMeRiC >
	int shape_quad9( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0 )
	{
	
      	Na[0]=0.25*xi*eta*(xi-1)*(eta-1);
      	Na[inc]=0.25*xi*eta*(xi+1)*(eta-1);
      	Na[2*inc]=0.25*xi*eta*(xi+1)*(eta+1);
      	Na[3*inc]=0.25*xi*eta*(xi-1)*(eta+1);
      	Na[4*inc]=-0.5*eta*(xi+1)*(xi-1)*(eta-1);
      	Na[5*inc]=-0.5*xi*(xi+1)*(eta+1)*(eta-1);
      	Na[6*inc]=-0.5*eta*(xi+1)*(xi-1)*(eta+1);
      	Na[7*inc]=-0.5*xi*(xi-1)*(eta+1)*(eta-1);
      	Na[8*inc]=(xi+1)*(xi-1)*(eta+1)*(eta-1);
		return 9;
		
	}

	/**
	Computes the gradient of the shape function for a nine node quadrilateral element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[-1,1]
	\param edim (optional) on returns contains the dimension of the element in the parament space

		\image html quad9.pdf
		
	*/
template < class NuMeRiC >
	int dshape_quad9( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0 )
	{

		dNa[0]=0.25*eta*(2*xi-1)*(eta-1);		dNa[0+lda]=0.25*xi*(xi-1)*(2*eta-1);
	    dNa[1]=0.25*eta*(2*xi+1)*(eta-1);		dNa[1+lda]=0.25*xi*(xi+1)*(2*eta-1);
	    dNa[2]=0.25*eta*(2*xi+1)*(eta+1);		dNa[2+lda]=0.25*xi*(xi+1)*(2*eta+1);
	    dNa[3]=0.25*eta*(2*xi-1)*(eta+1);		dNa[3+lda]=0.25*xi*(xi-1)*(2*eta+1);
	    dNa[4]=-xi*eta*(eta-1);	  			dNa[4+lda]=-0.5*(xi+1)*(xi-1)*(2*eta-1);
	    dNa[5]=-0.5*(2*xi+1)*(eta+1)*(eta-1);	dNa[5+lda]=-xi*eta*(xi+1);
	    dNa[6]=-xi*eta*(eta+1);				dNa[6+lda]=-0.5*(xi+1)*(xi-1)*(2*eta+1);
	    dNa[7]=-0.5*(2*xi-1)*(eta+1)*(eta-1);	dNa[7+lda]=-xi*eta*(xi-1);
		dNa[8]=2*xi*(eta*eta-1);				dNa[8+lda]=2*eta*(xi*xi-1);	
			return 9;	
		
	}
	
	
	/** 
		Computes the shape function for a twelve node quadrilateral element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[-1,1]

		\image html quad12.pdf

	*/
template < class NuMeRiC >
	int shape_quad12( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0 )
	{
		
		NuMeRiC r=0.5*(xi+1.0), s=0.5*(eta+1.0);
		NuMeRiC a = 0.0, b = 1.0/3.0, c = 2.0/3.0, d = 1.0; 
	    NuMeRiC corner = 9.0  * ( ( 2.0  * r - 1.0  ) * ( 2.0  * r - 1.0  ) 
			+ ( 2.0  * s - 1.0  ) * ( 2.0  * s - 1.0  ) ) - 10.0;

		Na[0] =     0.125   * ( r - d ) * ( s - d ) * corner;
		Na[inc] =   - 0.125   * ( r - a ) * ( s - d ) * corner;
		Na[2*inc] =    0.125   * ( r - a ) * ( s - a ) * corner;
		Na[3*inc] =   - 0.125   * ( r - d ) * ( s - a ) * corner;
		Na[4*inc] =  - 13.5     * ( r - a ) * ( r - c ) * ( r - d ) * ( s - d );
		Na[5*inc] =    13.5     * ( r - a ) * ( r - b ) * ( r - d ) * ( s - d );
		Na[6*inc] =    13.5     * ( r - a ) * ( s - a ) * ( s - c ) * ( s - d );
		Na[7*inc] =  - 13.5     * ( r - a ) * ( s - a ) * ( s - b ) * ( s - d );
		Na[8*inc] = - 13.5     * ( r - a ) * ( r - b ) * ( r - d ) * ( s - a );
		Na[9*inc] =   13.5     * ( r - a ) * ( r - c ) * ( r - d ) * ( s - a );
		Na[10*inc] =    13.5     * ( r - d ) * ( s - a ) * ( s - b ) * ( s - d );
		Na[11*inc] =  - 13.5     * ( r - d ) * ( s - a ) * ( s - c ) * ( s - d );
		return 12;

		
	}

	/**
	Computes the gradient of the shape function for a twelve node quadrilateral element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at xi=[-1,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at eta=[-1,1]
	\param edim (optional) on returns contains the dimension of the element in the parament space


	\image html quad12.pdf

	*/
template < class NuMeRiC >
	int dshape_quad12( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0 )
	{

		NuMeRiC r=0.5*(xi+1.0), s=0.5*(eta+1.0);
		NuMeRiC a = 0.0, b = 1.0/3.0, c = 2.0/3.0, d = 1.0; 
		NuMeRiC corner = 9.0  * ( ( 2.0  * r - 1.0  ) * ( 2.0  * r - 1.0  ) 
				+ ( 2.0  * s - 1.0  ) * ( 2.0  * s - 1.0  ) ) - 10.0;
				
		NuMeRiC dcdr = 36.0  * ( 2.0  * r - 1.0  );
		dNa[0] =  0.125 * ( s - d ) * ( ( r - d ) * dcdr + corner );
	    dNa[4] =  - 13.5  * ( s - d ) * ( 3.0  * r * r 
			- 2.0  * ( a + c + d ) * r + a * c + c * d + d * a ); 
	    dNa[5] =    13.5  * ( s - d ) * ( 3.0  * r * r 
			- 2.0  * ( a + b + d ) * r + a * b + b * d + d * a );
		dNa[1] = - 0.125  * ( s - d ) * ( ( r - a ) * dcdr + corner );
		dNa[11] = - 13.5  * ( s - a ) * ( s - c ) * ( s - d ); 
		dNa[6] =   13.5  * ( s - a ) * ( s - c ) * ( s - d );
		dNa[10] =   13.5  * ( s - a ) * ( s - b ) * ( s - d );
		dNa[7] = - 13.5  * ( s - a ) * ( s - b ) * ( s - d );
		dNa[3] = - 0.125  * ( s - a ) * ( ( r - d ) * dcdr + corner );
	    dNa[9] =   13.5  * ( s - a ) * ( 3.0  * r * r 
			- 2.0  * ( a + c + d ) * r + a * c + c * d + d * a ); 
	    dNa[8] = - 13.5  * ( s - a ) * ( 3.0  * r * r 
			- 2.0  * ( a + b + d ) * r + a * b + b * d + d * a );
		dNa[2] = 0.125  * ( s - a ) * ( ( r - a ) * dcdr + corner );

		NuMeRiC dcds = 36.0  * ( 2.0  * s - 1.0  );
		dNa[0+lda] =  0.125  * ( r - d ) * ( corner + ( s - d ) * dcds );
		dNa[4+lda] =  - 13.5  * ( r - a ) * ( r - c ) * ( r - d ) ;
		dNa[5+lda] =  13.5  * ( r - a ) * ( r - b ) * ( r - d );
		dNa[1+lda] = - 0.125   * ( r - a ) * ( corner + ( s - d ) * dcds );
	    dNa[11+lda] =  - 13.5  * ( r - d ) * ( 3.0  * s * s 
		- 2.0  * ( a + c + d ) * s + a * c + c * d + d * a );
	    dNa[6+lda] =  13.5  * ( r - a ) * ( 3.0  * s * s 
			- 2.0  * ( a + c + d ) * s + a * c + c * d + d * a );
	    dNa[10+lda] =  13.5  * ( r - d ) * ( 3.0  * s * s 
			- 2.0  * ( a + b + d ) * s + a * b + b * d + d * a );
	    dNa[7+lda] =  - 13.5  * ( r - a ) * ( 3.0  * s * s 
			- 2.0  * ( a + b + d ) * s + a * b + b * d + d * a );
		dNa[3+lda] =  - 0.125  * ( r - d ) * ( corner + ( s - a ) * dcds );
		dNa[9+lda] = 13.5  * ( r - a ) * ( r - c ) * ( r - d ) ;
		dNa[8+lda] = - 13.5  * ( r - a ) * ( r - b ) * ( r - d ) ;
		dNa[2+lda] = 0.125  * ( r - a ) * ( corner + ( s - a ) * dcds );

		return 12;		
		
	}
	
	
	/** 
		Computes the shape function for a four node tetrahedal element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at xi
			\param eta the point in the parent coordinate space where the shape functions are evalueated at eta
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at


			\image html tetra4.pdf
	
	*/
template < class NuMeRiC >
	int shape_tetra4( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=ONESIXTH, const NuMeRiC &eta=ONESIXTH, 
		const NuMeRiC &zeta=ONESIXTH )
	{
	
		Na[0]=1-xi-eta-zeta;
		Na[inc]=xi;
		Na[2*inc]=eta;
		Na[3*inc]=zeta;
		return 4;
	}

	/**
	Computes the gradient of the shape function for a four node tetrahedal element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at  
	\param eta the point in the parent coordinate space where the shape functions are evalueated at  
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at  
	\param edim (optional) on returns contains the dimension of the element in the parament space

	\image html tetra4.pdf

	*/
	template < class NuMeRiC >
	int dshape_tetra4( NuMeRiC *dNa, int lda, const NuMeRiC &xi=ONESIXTH, const NuMeRiC &eta=ONESIXTH, 
		const NuMeRiC &zeta=ONESIXTH )
	{

		dNa[0]=-1.0;	dNa[0+lda]=-1.0;	dNa[0+2*lda]=-1.0;
		dNa[1]= 1.0;	dNa[1+lda]= 0.0;	dNa[1+2*lda]= 0.0;
		dNa[2]= 0.0;	dNa[2+lda]= 1.0;	dNa[2+2*lda]= 0.0;
		dNa[3]= 0.0;	dNa[3+lda]= 0.0;	dNa[3+2*lda]= 1.0;

		return 4;
		
	}

	/** 
		Computes the shape function for a ten node tetrahedal element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at 
			\param eta the point in the parent coordinate space where the shape functions are evalueated at 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at


			\image html tetra10.pdf
	
	*/
	template < class NuMeRiC >
	int shape_tetra10( NuMeRiC *Na, int inc, const NuMeRiC &xi=ONESIXTH, const NuMeRiC &eta=ONESIXTH, 
		const NuMeRiC &zeta=ONESIXTH )
	{

		NuMeRiC l2=xi, l3=eta, l4=zeta, l1=1-xi-eta-zeta;

		Na[0] = l1*(2*l1-1); 
		Na[inc] = l2*(2*l2-1); 
		Na[2*inc] = l3*(2*l3-1); 
		Na[3*inc] = l4*(2*l4-1); 
		Na[4*inc] = 4*l1*l2;
		Na[5*inc] = 4*l3*l2;
		Na[6*inc] = 4*l1*l3;
		Na[7*inc] = 4*l1*l4;
		Na[8*inc] = 4*l4*l2;
		Na[9*inc]= 4*l3*l4;
		return 10;
	}

	/**
	Computes the gradient of the shape function for a ten node tetrahedal element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at  
	\param eta the point in the parent coordinate space where the shape functions are evalueated at  
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at  
	\param edim (optional) on returns contains the dimension of the element in the parament space


	\image html tetra10.pdf

	*/
	template < class NuMeRiC >
	int dshape_tetra10( NuMeRiC *dNa, int lda, const NuMeRiC &xi=ONESIXTH, const NuMeRiC &eta=ONESIXTH, 
		const NuMeRiC &zeta=ONESIXTH )
	{

		NuMeRiC l2=xi, l3=eta, l4=zeta, l1=1-xi-eta-zeta;
		
		dNa[0]=-(4*l1-1); dNa[0+lda]=-(4*l1-1); dNa[0+2*lda]=-(4*l1-1);
		dNa[1]= 4*l2-1; 	dNa[1+lda]= 0.0; 		dNa[1+2*lda]= 0.0;
		dNa[2]= 0.0;		dNa[2+lda]= 4*l3-1;	dNa[2+2*lda]= 0.0;
		dNa[3]= 0.0;		dNa[3+lda]= 0.0;		dNa[3+2*lda]= 4*l4-1;
		dNa[4]=4*(l1-l2);	dNa[4+lda]=-4*l2;		dNa[4+2*lda]=-4*l2;
		dNa[5]= 4*l3;		dNa[5+lda]= 4*l2;		dNa[5+2*lda]= 0.0;
		dNa[6]=-4*l3;		dNa[6+lda]=4*(l1-l3);	dNa[6+2*lda]=-4*l3;
		dNa[7]=-4*l4;		dNa[7+lda]=-4*l4;		dNa[7+2*lda]=4*(l1-l4);
		dNa[8]= 4*l4;		dNa[8+lda]= 0.0;		dNa[8+2*lda]= 4*l2;   
		dNa[9]= 0.0;		dNa[9+lda]= 4*l4;		dNa[9+2*lda]= 4*l3;

		return 10;
		
	}

	/** 
		Computes the shape function for a eight node hexahedral element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[-1,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[-1,1] 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1]


			\image html hexa8.pdf
	
	*/
	template < class NuMeRiC >
	int shape_hexa8( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.0 )
	{

		Na[0] = 0.125*(1-xi)*(1-eta)*(1-zeta);
		Na[inc] = 0.125*(1+xi)*(1-eta)*(1-zeta);
		Na[2*inc] = 0.125*(1+xi)*(1+eta)*(1-zeta);
		Na[3*inc] = 0.125*(1-xi)*(1+eta)*(1-zeta);
		Na[4*inc] = 0.125*(1-xi)*(1-eta)*(1+zeta);
		Na[5*inc] = 0.125*(1+xi)*(1-eta)*(1+zeta);
		Na[6*inc] = 0.125*(1+xi)*(1+eta)*(1+zeta);
		Na[7*inc] = 0.125*(1-xi)*(1+eta)*(1+zeta);
		return 8;	
	}

	/**
	Computes the gradient of the shape function for a eight node hexahedral element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evaluated at  
	\param eta the point in the parent coordinate space where the shape functions are evaluated at  
	\param zeta the point in the parent coordinate space where the shape functions are evaluated at  
	\param edim (optional) on returns contains the dimension of the element in the parent space

			\image html hexa8.pdf

	*/
	template < class NuMeRiC >
	int dshape_hexa8( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.0 )
	{
		
		dNa[0]=0.125*(-1+eta+zeta-eta*zeta);	dNa[0+lda]=0.125*(-1+xi+zeta-xi*zeta);	dNa[0+2*lda]=0.125*(-1+xi+eta-xi*eta);
		dNa[1]=0.125*(1-eta-zeta+eta*zeta);	    dNa[1+lda]=0.125*(-1-xi+zeta+xi*zeta);	dNa[1+2*lda]=0.125*(-1-xi+eta+xi*eta);
		dNa[2]=0.125*(1+eta-zeta-eta*zeta);		dNa[2+lda]=0.125*(1+xi-zeta-xi*zeta);	dNa[2+2*lda]=0.125*(-1-xi-eta-xi*eta);	
		dNa[3]=0.125*(-1-eta+zeta+eta*zeta);	dNa[3+lda]=0.125*(1-xi-zeta+xi*zeta);	dNa[3+2*lda]=0.125*(-1+xi-eta+xi*eta);
		dNa[4]=0.125*(-1+eta-zeta+eta*zeta);	dNa[4+lda]=0.125*(-1+xi-zeta+xi*zeta);	dNa[4+2*lda]=0.125*(1-xi-eta+xi*eta);
		dNa[5]=0.125*(1-eta+zeta-eta*zeta);		dNa[5+lda]=0.125*(-1-xi-zeta-xi*zeta);	dNa[5+2*lda]=0.125*(1+xi-eta-xi*eta);
		dNa[6]=0.125*(1+eta+zeta+eta*zeta);		dNa[6+lda]=0.125*(1+xi+zeta+xi*zeta);	dNa[6+2*lda]=0.125*(1+xi+eta+xi*eta);
		dNa[7]=0.125*(-1-eta-zeta-eta*zeta);	dNa[7+lda]=0.125*(1-xi+zeta-xi*zeta);	dNa[7+2*lda]=0.125*(1-xi+eta-xi*eta);

		return 8;
		
	}

	/** 
		Computes the shape function for a twenty node hexahedral element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[-1,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[-1,1] 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1]

			\image html hexa20.pdf
	
	*/
	template < class NuMeRiC >
	int shape_hexa20( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.0 )
	{

		Na[0] = 0.125*(1-xi)*(1-eta)*(1-zeta)*(-xi-eta-zeta-2);
		Na[inc] = 0.25*(1-xi)*(1-eta)*(1-zeta*zeta);
		Na[2*inc] = 0.125*(1-xi)*(1-eta)*(1+zeta)*(-xi-eta+zeta-2);
		Na[3*inc] = 0.25*(1-xi*xi)*(1-eta)*(1+zeta);
		Na[4*inc] = 0.125*(1+xi)*(1-eta)*(1+zeta)*(xi-eta+zeta-2);
		Na[5*inc] = 0.25*(1+xi)*(1-eta)*(1-zeta*zeta);
		Na[6*inc] = 0.125*(1+xi)*(1-eta)*(1-zeta)*(xi-eta-zeta-2);
		Na[7*inc] = 0.25*(1-xi*xi)*(1-eta)*(1-zeta);
		Na[8*inc] = 0.25*(1-xi)*(1-eta*eta)*(1-zeta);
		Na[9*inc] = 0.25*(1-xi)*(1-eta*eta)*(1+zeta);
		Na[10*inc] = 0.25*(1+xi)*(1-eta*eta)*(1+zeta);
		Na[11*inc] = 0.25*(1+xi)*(1-eta*eta)*(1-zeta);
		Na[12*inc] = 0.125*(1-xi)*(1+eta)*(1-zeta)*(-xi+eta-zeta-2);
		Na[13*inc] = 0.25*(1-xi)*(1+eta)*(1-zeta*zeta);
		Na[14*inc] = 0.125*(1-xi)*(1+eta)*(1+zeta)*(-xi+eta+zeta-2);
		Na[15*inc] = 0.25*(1-xi*xi)*(1+eta)*(1+zeta);
		Na[16*inc] = 0.125*(1+xi)*(1+eta)*(1+zeta)*(+xi+eta+zeta-2);
		Na[17*inc] = 0.25*(1+xi)*(1+eta)*(1-zeta*zeta);
		Na[18*inc] = 0.125*(1+xi)*(1+eta)*(1-zeta)*(xi+eta-zeta-2);
		Na[19*inc] = 0.25*(1-xi*xi)*(1+eta)*(1-zeta);
		return 20;
	}

	/**
	Computes the gradient of the shape function for a twenty node hexahedral element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at  
	\param eta the point in the parent coordinate space where the shape functions are evalueated at  
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at  
	\param edim (optional) on returns contains the dimension of the element in the parament space

			\image html hexa20.pdf

	*/
	template < class NuMeRiC >
	int dshape_hexa20( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.0 )
	{
		
		dNa[0]=-0.125*(1-eta)*(1-zeta)*(-zeta-xi-eta-2)-0.125*(1-eta)*(1-xi)*(1-zeta); 
		dNa[0+lda]=-0.125*(1-xi)*(1-zeta)*(-zeta-xi-eta-2)-0.125*(1-eta)*(1-xi)*(1-zeta); 
		dNa[0+2*lda]=-0.125*(1-eta)*(1-xi)*(-zeta-xi-eta-2) -0.125*(1-eta)*(1-xi)*(1-zeta);
		dNa[1]=-0.25*(1-eta)*(1-zeta*zeta);
		dNa[1+lda]=-0.25*(1-xi)*(1-zeta*zeta);
		dNa[1+2*lda]=-0.5*(1-eta)*(1-xi)*zeta; 
		dNa[2]=-0.125*(1-eta)*(zeta+1)*(zeta-xi-eta-2)-0.125*(1-eta)*(1-xi)*(zeta+1); 
		dNa[2+lda]=-0.125*(1-xi)*(zeta+1)*(zeta-xi-eta-2)-0.125*(1-eta)*(1-xi)*(zeta+1); 
		dNa[2+2*lda]= 0.125*(1-eta)*(1-xi)*(zeta-xi-eta-2)+0.125*(1-eta)*(1-xi)*(zeta+1); 
		dNa[3]=-0.5*(1-eta)*xi*(zeta+1);
		dNa[3+lda]=-0.25*(1-xi*xi )*(zeta+1);
		dNa[3+2*lda]= 0.25*(1-eta)*(1-xi*xi); 
		dNa[4]= 0.125*(1-eta)*(zeta+1)*(zeta+xi-eta-2)+0.125*(1-eta)*(xi+1)*(zeta+1); 
		dNa[4+lda]= -0.125*(xi+1)*(zeta+1)*(zeta+xi-eta-2)-0.125*(1-eta)*(xi+1)*(zeta+1); 
		dNa[4+2*lda]= 0.125*(1-eta)*(xi+1)*(zeta+xi-eta-2)+0.125*(1-eta)*(xi+1)*(zeta+1);
		dNa[5]= 0.25*(1-eta)*(1-zeta*zeta);
		dNa[5+lda]=-0.25*(xi+1)*(1-zeta*zeta);
		dNa[5+2*lda]=-0.5*(1-eta)*(xi+1)*zeta; 
		dNa[6]= 0.125*(1-eta)*(1-zeta)*(-zeta+xi-eta-2)+0.125*(1-eta)*(xi+1)*(1-zeta); 
		dNa[6+lda]=-0.125*(xi+1)*(1-zeta)*(-zeta+xi-eta-2)-0.125*(1-eta)*(xi+1)*(1-zeta); 
		dNa[6+2*lda]=-0.125*(1-eta)*(xi+1)*(-zeta+xi-eta-2)-0.125*(1-eta)*(xi+1)*(1-zeta); 
		dNa[7]=-0.5*(1-eta)*xi*(1-zeta);
		dNa[7+lda]=-0.25*(1-xi*xi)*(1-zeta);
		dNa[7+2*lda]=-0.25*(1-eta)*(1-xi*xi); 
		dNa[8]=-0.25*(1-eta )*(1-zeta); 
		dNa[8+lda]=-0.5*eta*(1-xi)*(1-zeta); 
		dNa[8+2*lda]=-0.25*(1-eta*eta )*(1-xi) ;
		dNa[9]=-0.25*(1-eta*eta)*(zeta+1); 
		dNa[9+lda]=-0.5*eta*(1-xi)*(zeta+1);
		dNa[9+2*lda]= 0.25*(1-eta*eta)*(1-xi);
		dNa[10]= 0.25*(1-eta*eta)*(zeta+1); 
		dNa[10+lda]=-0.5*eta*(xi+1)*(zeta+1);
		dNa[10+2*lda]= 0.25*(1-eta*eta)*(xi+1);
		dNa[11]= 0.25*(1-eta*eta)*(1-zeta); 
		dNa[11+lda]=-0.5*eta*(xi+1)*(1-zeta); 
		dNa[11+2*lda]=-0.25*(1-eta*eta)*(xi+1) ;
		dNa[12]=-0.125*(eta+1)*(1-zeta)*(-zeta-xi+eta-2)-0.125*(eta+1)*(1-xi)*(1-zeta); 
		dNa[12+lda]= 0.125*(1-xi)*(1-zeta)*(-zeta-xi+eta-2)+0.125*(eta+1)*(1-xi)*(1-zeta); 
		dNa[12+2*lda]=-0.125*(eta+1)*(1-xi)*(-zeta-xi+eta-2)-0.125*(eta+1)*(1-xi)*(1-zeta);
		dNa[13]=-0.25*(eta+1)*(1-zeta*zeta);
		dNa[13+lda]= 0.25*(1-xi)*(1-zeta*zeta); 
		dNa[13+2*lda]=-0.5*(eta+1)*(1-xi)*zeta; 
		dNa[14]=-0.125*(eta+1)*(zeta+1)*(zeta-xi+eta-2)-0.125*(eta+1)*(1-xi)*(zeta+1); 
		dNa[14+lda]= 0.125*(1-xi)*(zeta+1)*(zeta-xi+eta-2)+0.125*(eta+1)*(1-xi)*(zeta+1); 
		dNa[14+2*lda]= 0.125*(eta+1)*(1-xi)*(zeta-xi+eta-2)+0.125*(eta+1)*(1-xi)*(zeta+1); 
		dNa[15]=-0.5*(eta+1)*xi*(zeta+1);
		dNa[15+lda]= 0.25*(1-xi*xi)*(zeta+1);
		dNa[15+2*lda]= 0.25*(eta+1)*(1-xi*xi); 
		dNa[16]= 0.125*(eta+1)*(zeta+1)*(zeta+xi+eta-2)+0.125*(eta+1)*(xi+1)*(zeta+1);
		dNa[16+lda]= 0.125*(xi+1)*(zeta+1)*(zeta+xi+eta-2)+0.125*(eta+1)*(xi+1)*(zeta+1); 
		dNa[16+2*lda]= 0.125*(eta+1)*(xi+1)*(zeta+xi+eta-2)+0.125*(eta+1)*(xi+1)*(zeta+1);	
		dNa[17]= 0.25*(eta+1)*(1-zeta);
		dNa[17+lda]= 0.25*(xi+1)*(1-zeta*zeta); 
		dNa[17+2*lda]=-0.5*(eta+1)*(xi+1)*zeta; 
		dNa[18]= 0.125*(eta+1)*(1-zeta)*(-zeta+xi+eta-2)+0.125*(eta+1)*(xi+1)*(1-zeta);
		dNa[18+lda]= 0.125*(xi+1)*(1-zeta)*(-zeta+xi+eta-2)+0.125*(eta+1)*(xi+1)*(1-zeta);
		dNa[18+2*lda]=-0.125*(eta+1)*(xi+1)*(-zeta+xi+eta-2)-0.125*(eta+1)*(xi+1)*(1-zeta);
		dNa[19]=-0.5*(eta+1)*xi*(1-zeta);
		dNa[19+lda]= 0.25*(1-xi*xi)*(1-zeta);
		dNa[19+2*lda]=-0.25*(eta+1)*(1-xi*xi);

		return 20;
		
	}

	/** 
		Computes the shape function for a twenty seven node hexahedral element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[-1,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[-1,1] 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1]

			\image html hexa27.pdf
	
	*/
	template < class NuMeRiC >
	int shape_hexa27( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.0 )
	{
		cout << "\nshape_hexa27" << " not yet implemented\n";
		return 27;
	}

	/**
	Computes the gradient of the shape function for a twenty seven node hexahedral element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at  
	\param eta the point in the parent coordinate space where the shape functions are evalueated at  
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at 
	\param edim (optional) on returns contains the dimension of the element in the parament space 

			\image html hexa27.pdf

	*/
	template < class NuMeRiC >
	int dshape_hexa27( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.0 )
	{
	    cout << "\ndshape_hexa27" << " not yet implemented\n";
	    return 27;	
	}
	
	/** 
		Computes the shape function for a six node prism element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[0,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[0,1-xi] 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1]

			\image html prism6.pdf
	
	*/
	template < class NuMeRiC >
	int shape_prism6( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD, 
		const NuMeRiC &zeta=0.0 )
	{

		Na[0]=(1-xi-eta)*(1-zeta);
		Na[inc]=xi*(1-zeta);
		Na[2*inc]=eta*(1-zeta);
		Na[3*inc]=(1-xi-eta)*zeta;
		Na[4*inc]=xi*zeta;
		Na[5*inc]=eta*zeta;
		return 6;
	}

	/**
	Computes the gradient of the shape function for a six node prism element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[0,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[0,1-xi] 
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1]  
	\param edim (optional) on returns contains the dimension of the element in the parament space

			\image html prism6.pdf

	*/
	template < class NuMeRiC >
	int dshape_prism6( NuMeRiC *dNa, int lda, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD, 
		const NuMeRiC &zeta=0.0 )
	{

		dNa[0]= zeta-1.0; 	dNa[0+lda]= zeta-1.0;		dNa[0+2*lda]= xi+eta-1.0;
		dNa[1]= 1.0-zeta;		dNa[1+lda]= 0.0;			dNa[1+2*lda]= -xi;
		dNa[2]= 0.0;			dNa[2+lda]= 1.0-zeta;		dNa[2+2*lda]= -eta ;
		dNa[3]= -zeta;		dNa[3+lda]= -zeta;		dNa[3+2*lda]= 1-xi-eta; 
		dNa[4]= zeta;			dNa[4+lda]= 0.0;			dNa[4+2*lda]= xi ;
		dNa[5]= 0.0;			dNa[5+lda]= zeta;			dNa[5+2*lda]= eta;
		//3;
		return 6;
	}
	
	/** 
		Computes the shape function for a twelve node prism element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[0,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[0,1-xi] 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1]

			\image html prism12.pdf
	
	*/
	template < class NuMeRiC >
	int shape_prism12( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD, 
		const NuMeRiC &zeta=0.0 )
	{	
        cout << "\nshape_prism12" << " not yet implemented\n";
		return 12;
	}

	/**
	Computes the gradient of the shape function for a twelve node prism element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[0,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[0,1-xi] 
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1] 
	\param edim (optional) on returns contains the dimension of the element in the parament space 

			\image html prism12.pdf

	*/
	template < class NuMeRiC >
	int dshape_prism12( NuMeRiC *dNa, int lda, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD, 
		const NuMeRiC &zeta=0.0  )
	{
        cout << "\ndshape_prism12" << " not yet implemented\n";
	
		return 12;
	}
	
	/** 
		Computes the shape function for a fiveteen node prism element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[0,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[0,1-xi] 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1]

			\image html prism15.pdf	
	*/
	template < class NuMeRiC >
	int shape_prism15( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD, 
		const NuMeRiC &zeta=0.0 )
	{
        cout << "\nshape_prism15" << " not yet implemented\n";
	
	    return 15;	
	}

	/**
	Computes the gradient of the shape function for a eighteen node prism element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[0,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[0,1-xi] 
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1] 
	\param edim (optional) on returns contains the dimension of the element in the parament space
	*/
	template < class NuMeRiC >
	int dshape_prism15( NuMeRiC *dNa, int lda, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD, 
		const NuMeRiC &zeta=0.0 )
	{
	
        cout << "\ndshape_prism15" << " not yet implemented\n";

	    return 15;	
	}
	
	/** 
		Computes the shape function for a eghteen node prism element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[0,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[0,1-xi] 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1]

			\image html prism18.pdf	
	*/
	template < class NuMeRiC >
	int shape_prism18( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD, 
		const NuMeRiC &zeta=0.0 )
	{
	
		Na[0] = 2*(1-3*(xi+eta)+4*xi*eta+2*(xi*xi+eta*eta))*(1-zeta)*(0.5-zeta);
		Na[inc] = 2*xi*(2*xi-1)*(1-zeta)*(0.5-zeta);
		Na[2*inc] = 2*eta*(2*eta-1)*(1-zeta)*(0.5-zeta);
		Na[3*inc] = 2*(1-3*(xi+eta)+4*xi*eta+2*(xi*xi+eta*eta))*zeta*(0.5-zeta);
		Na[4*inc] = 2*xi*(2*xi-1)*zeta*(0.5-zeta);
		Na[5*inc] = 2*eta*(2*eta-1)*zeta*(0.5-zeta);
		Na[6*inc] = 8*xi*(1-xi-eta)*(1-zeta)*(0.5-zeta);
		Na[7*inc] = 8*xi*eta*(1-zeta)*(0.5-zeta);
		Na[8*inc] = 8*eta*(1-xi-eta)*(1-zeta)*(0.5-zeta);    
		Na[9*inc] = 4*(1-3*(xi+eta)+4*xi*eta+2*(xi*xi+eta*eta))*zeta*(1-zeta);
		Na[10*inc] = 4*xi*(2*xi-1)*zeta*(1-zeta); 	   
		Na[11*inc] = 4*eta*(2*eta-1)*zeta*(1-zeta);
		Na[12*inc] = 16*xi*(1-xi-eta)*zeta*(1-zeta);
		Na[13*inc] = 16*xi*eta*zeta*(1-zeta);
		Na[14*inc] = 16*eta*(1-xi-eta)*zeta*(1-zeta);  
		Na[15*inc] = 8*xi*(1-xi-eta)*zeta*(0.5-zeta);
		Na[16*inc] = 8*xi*eta*zeta*(0.5-zeta);
		Na[17*inc] = 8*eta*(1-xi-eta)*zeta*(0.5-zeta);
		return 18;
	}

	/**
	Computes the gradient of the shape function for a eighteen node prism element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[0,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[0,1-xi] 
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[-1,1] 
	\param edim (optional) on returns contains the dimension of the element in the parament space

			\image html prism18.pdf	

	*/
	template < class NuMeRiC >
	int dshape_prism18( NuMeRiC *dNa, int lda, const NuMeRiC &xi=ONETHIRD, const NuMeRiC &eta=ONETHIRD, 
		const NuMeRiC &zeta=0.0 )
	{

		dNa[ 0 ]= 2 *(4 *xi + 4* eta - 3)*(0.5 - zeta)*(1 - zeta); 
		dNa[ 0 +lda]= 2 *(4 *xi + 4 *eta - 3)*(0.5 - zeta)*(1 - zeta);
	  dNa[ 0 +2*lda]= -2*(2*(xi*xi+eta*eta)-3*(xi+eta)+4* eta*xi+1)*(1-zeta)-2*(2*(xi*xi+eta*eta) 
			-3*(xi+eta+4*eta*xi+1))*(0.5-zeta) ;
		dNa[ 1 ]= 2 *(2* xi - 1)*(0.5 - zeta)*(1 - zeta) + 4 *xi* (0.5 - zeta)*(1 - zeta) ;
		dNa[ 1 +lda]= 0.0 ;
		dNa[ 1 +2*lda]= - 2 *xi* (2 *xi - 1)*(1 - zeta) - 2 *xi* (2* xi - 1)*(0.5 - zeta) ;
		dNa[ 2 ]= 0.0 ;
		dNa[ 2 +lda]= 2* (2 *eta - 1)*(0.5 - zeta)*(1 - zeta)+ 4 *eta *(0.5 - zeta)*(1 - zeta) ;
		dNa[ 2 +2*lda]= - 2 *eta *(2* eta - 1)*(1 - zeta) - 2* eta *(2*eta - 1)*(0.5 - zeta) ;
		dNa[ 3 ]= 2 *(4* xi + 4* eta - 3)*(0.5 - zeta)* zeta ;
		dNa[ 3 +lda]= 2 *(4* xi + 4 *eta - 3)*(0.5 - zeta)* zeta ;
	    dNa[ 3 +2*lda]= 2*(2*(xi*xi+eta*eta)-3*(xi+eta)+4*eta*xi+1)*(0.5-zeta)-2*(2*(xi*xi+eta*eta) 
			-3*(xi+eta)+4*eta*xi+1)*zeta ;
		dNa[ 4 ]= 2*(2* xi - 1)*(0.5 - zeta)* zeta + 4* xi* (0.5 - zeta)*zeta ;
		dNa[ 4 +lda]= 0.0 ;
		dNa[ 4 +2*lda]= 2* xi* (2* xi - 1)*(0.5 - zeta) - 2 *xi* (2 *xi - 1)* zeta ;
		dNa[ 5 ]= 0.0 ;
		dNa[ 5 +lda]= 2 *(2*eta - 1)*(0.5 - zeta)* zeta + 4 *eta* (0.5 - zeta)* zeta ;
		dNa[ 5 +2*lda]= 2* eta *(2 *eta - 1)*(0.5 - zeta) - 2* eta* (2 *eta - 1) *zeta ;
		dNa[ 6 ]= 8 *(- xi - eta + 1)*(0.5 - zeta)*(1 - zeta)- 8 *xi *(0.5 - zeta)*(1 - zeta) ;
		dNa[ 6 +lda]= - 8 *xi *(0.5 - zeta)*(1 - zeta) ;
		dNa[ 6 +2*lda]= - 8* (- xi - eta + 1)* xi* (1 - zeta) - 8 *(- xi - eta + 1) *xi* (0.5 - zeta); 
		dNa[ 7 ]= 8* eta *(0.5 - zeta)*(1 - zeta) ;
		dNa[ 7 +lda]= 8* xi *(0.5 - zeta)*(1 - zeta) ;
		dNa[ 7 +2*lda]= - 8 *eta* xi* (1 - zeta) - 8 *eta* xi* (0.5 - zeta) ;
		dNa[ 8 ]= - 8 *eta* (0.5 - zeta)*(1 - zeta) ;
		dNa[ 8 +lda]= 8 *(- xi - eta + 1)*(0.5 - zeta)*(1 - zeta)  - 8 *eta *(0.5 - zeta)*(1 - zeta);
		dNa[ 8 +2*lda]= - 8 *eta* (- xi - eta + 1)*(1 - zeta) - 8* eta* (- xi - eta + 1)*(0.5 - zeta) ;
		dNa[ 9 ]= 4* (4* xi + 4* eta - 3)*(1 - zeta) *zeta ;
		dNa[ 9 +lda]= 4* (4 *xi + 4 *eta - 3)*(1 - zeta) *zeta ;
	    dNa[ 9 +2*lda]= 4* (2* (xi*xi  + eta*eta ) - 3* (xi + eta) + 4 *eta *xi + 1)*(1 - zeta) 
			- 4 *(2* (xi*xi  + eta*eta ) - 3 *(xi + eta) + 4* eta *xi + 1)* zeta ;
		dNa[ 10 ]= 4* (2* xi - 1)*(1 - zeta)* zeta + 8 *xi *(1 - zeta) *zeta ;
		dNa[ 10 +lda]= 0.0 ;
		dNa[ 10 +2*lda]= 4* xi *(2* xi - 1)*(1 - zeta) - 4* xi* (2* xi - 1)* zeta ;
		dNa[ 11 ]= 0.0 ;
		dNa[ 11 +lda]= 4 *(2 *eta - 1)*(1 - zeta)* zeta + 8 *eta* (1 - zeta) *zeta ;
		dNa[ 11 +2*lda]= 4* eta *(2 *eta - 1)*(1 - zeta) - 4 *eta* (2* eta - 1) *zeta ;
		dNa[ 12 ]= 16* (- xi - eta + 1)*(1 - zeta) *zeta - 16* xi* (1 - zeta)* zeta ;
		dNa[ 12 +lda]= - 16 *xi* (1 - zeta) *zeta ;
		dNa[ 12 +2*lda]= 16* (- xi - eta + 1)* xi* (1 - zeta) - 16* (- xi - eta + 1)* xi *zeta ;
		dNa[ 13 ]= 16* eta* (1 - zeta) *zeta ;
		dNa[ 13 +lda]= 16 *xi* (1 - zeta)* zeta ;
		dNa[ 13 +2*lda]= 16 *eta* xi* (1 - zeta) - 16* eta* xi *zeta ;
		dNa[ 14 ]= - 16* eta* (1 - zeta)* zeta ;
		dNa[ 14 +lda]= 16 *(- xi - eta + 1)*(1 - zeta)* zeta - 16 *eta *(1 - zeta)* zeta ;
		dNa[ 14 +2*lda]= 16 *eta* (- xi - eta + 1)*(1 - zeta) - 16* eta* (- xi - eta + 1) *zeta ; 
		dNa[ 15 ]= 8* (- xi - eta + 1)*(0.5 - zeta) *zeta - 8* xi* (0.5 - zeta) *zeta ;
		dNa[ 15 +lda]= - 8 *xi* (0.5 - zeta) *zeta ;
		dNa[ 15 +2*lda]= 8 *(- xi - eta + 1)* xi *(0.5 - zeta) - 8* (- xi - eta + 1) *xi *zeta ;
		dNa[ 16 ]= 8 *eta *(0.5 - zeta) *zeta ;
		dNa[ 16 +lda]= 8 *xi* (0.5 - zeta) *zeta ;
		dNa[ 16 +2*lda]= 8 *eta* xi* (0.5 - zeta) - 8* eta* xi* zeta ;
		dNa[ 17 ]= - 8 *eta* (0.5 - zeta) *zeta ;
		dNa[ 17 +lda]= 8* (- xi - eta + 1)*(0.5 - zeta) *zeta - 8* eta* (0.5 - zeta) *zeta ;
		dNa[ 17 +2*lda]= 8 *eta *(- xi - eta + 1)*(0.5 - zeta) - 8 *eta* (- xi - eta + 1) *zeta ;
		//3;
		return 18;
	}
	
	
	
	/** 
		Computes the shape function for a five node pyramid element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[-1,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[-1,1] 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[0,1]


			\image html pyramid5.pdf	
	
	*/
	template < class NuMeRiC >
	int shape_pyramid5( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.125 )
	{
	    cout << "\nshape_pyramid5" << " not yet implemented\n";
		return 5;
	}

	/**
	Computes the gradient of the shape function for a five node pyramid element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[-1,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[-1,1] 
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[0,1]
	\param edim (optional) on returns contains the dimension of the element in the parament space
	
		\image html pyramid5.pdf
	
	*/
	template < class NuMeRiC >
	int dshape_pyramid5( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.125 )
	{

        cout << "\ndshape_pyramid5" << " not yet implemented\n";
		return 5;	
	}
	
	/** 
		Computes the shape function for a thirteen node pyramid element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[-1,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[-1,1] 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[0,1]

				\image html pyramid13.pdf
	
	*/
	template < class NuMeRiC >
	int shape_pyramid13( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.125 )
	{
        cout << "\nshape_pyramid13" << " not yet implemented\n";
		return 13;
	}

	/**
	Computes the gradient of the shape function for a thirteen node pyramid element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[-1,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[-1,1] 
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[0,1]
	\param edim (optional) on returns contains the dimension of the element in the parament space

		\image html pyramid13.pdf
	*/
	template < class NuMeRiC >
	int dshape_pyramid13( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.125 )
	{
        cout << "\ndshape_pyramid13" << " not yet implemented\n";
	
		return 13;
	}
	
	/** 
		Computes the shape function for a fourteen node pyramid element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[-1,1]
			\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[-1,1] 
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[0,1]
	
		\image html pyramid14.pdf	
	
	*/
	template < class NuMeRiC >
	int shape_pyramid14( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.125 )
	{
        cout << "\nshape_pyramid14" << " not yet implemented\n";
		return 14;
	}

	/**
	Computes the gradient of the shape function for a fourteen node pyramid element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at, xi=[-1,1]
	\param eta the point in the parent coordinate space where the shape functions are evalueated at, eta=[-1,1] 
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at, zeta=[0,1]
	\param edim (optional) on returns contains the dimension of the element in the parament space
	
		\image html pyramid14.pdf
	
	*/
	template < class NuMeRiC >
	int dshape_pyramid14( NuMeRiC *dNa, int lda, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.125 )
	{
	
        cout << "\ndshape_pyramid14" << " not yet implemented\n";
		
		//3;	
		return 14;
	}	
	
	/** 
		Computes the shape function for a one node point element
			\param Na the array where the shape function values are returned
			\param inc increment in the Na vector where the shape function values will be places (default=1)
			\param xi the point in the parent coordinate space where the shape functions are evalueated at 
			\param eta the point in the parent coordinate space where the shape functions are evalueated at  
			\param zeta the point in the parent coordinate space where the shape functions are evalueated at 

	\image html point1.pdf
	
	*/
	template < class NuMeRiC >
	int shape_point1( NuMeRiC *Na, int inc=1, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.0 )
	{


		Na[0]=1.0;
		return 1;
	}

	/**
	Computes the gradient of the shape function for a one node point element
	w.r.t to the parent coordinate system
	\param dNa the array where the shape function gradient values are returned
	\param lda the leading dimension of dNa
	\param xi the point in the parent coordinate space where the shape functions are evalueated at 
	\param eta the point in the parent coordinate space where the shape functions are evalueated at 
	\param zeta the point in the parent coordinate space where the shape functions are evalueated at 
	\param edim (optional) on returns contains the dimension of the element in the parament space
	*/
	template < class NuMeRiC >
	int dshape_point1( NuMeRiC *dNa, int ldadNa, const NuMeRiC &xi=0.0, const NuMeRiC &eta=0.0, const NuMeRiC &zeta=0.125 )
	{

		dNa[0]=0.0; 
		//1;
		return 1;
	
	}

} // end namespace
#endif



