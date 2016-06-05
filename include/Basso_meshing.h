/*! \file meshing.h

Some basic mesh generation functions

\author Jack Chessa, jfchessa@utep.edu
\date Sunday, Oct 30, 2011

*/

#ifndef _BASSO_MESHING_H_
#define _BASSO_MESHING_H_

#include "Basso_defs.h"

#include <math.h>
#include <set>
#include <map>

namespace Basso
{

using namespace std;

enum MeshBiasType { BiasNONE, BiasPOWER, BiasGEOMETRIC, BiasBELL };

//--------------------------------------------------------------------------
/**
Computes a 1D array of nodes between the points
    
    (x1) ----------- (x2)
    
with \param nx points along the u-direction  

\param x1 - coordinate of first point 
\param x2 - coordinate of the second point
\param nx - number of points
the return value is a pointer to the next node ot be added in pts
*/
template <class VAL_TYPEx, class VAL_TYPE>
VAL_TYPE *node_array_1d_none( int nx, VAL_TYPEx x1, VAL_TYPEx x2, VAL_TYPE *pts )
{
    if ( nx<=2 )
    {
        pts[0]=x1; pts[1]=x2;
        return pts+2;
    }
    VAL_TYPE *pt=pts;
    VAL_TYPE h=(x2-x1)/(nx-1);
    for ( int i=0; i<nx; ++i, ++pt )
    {
        *pt = x1 + i*h;
    }
	return pt;
}

/**
Computes a 1D array of nodes between the points with a geometric bias
    
    (x1) ----------- (x2)
    
with \param nx points along the u-direction  

\param x1 - coordinate of first point 
\param x2 - coordinate of the second point
\param nx - number of points
\param bias = 'GEOMETRIC' spacing is such that the size of the spacing at x2 is
               1/b of the spacing at x1.
The return value is a pointer to the next node to be added in pts.
*/
template <class VAL_TYPEx, class VAL_TYPE>
VAL_TYPE *node_array_1d_gometric( BASSO_IDTYPE n, VAL_TYPEx x1, VAL_TYPEx x2, VAL_TYPE *pts, VAL_TYPEx bias )
{
    if ( n<=2 )
    {
        pts[0]=x1; pts[1]=x2;
        return pts+2;
    }
	VAL_TYPE expn=1.0/(n-2.0);
    VAL_TYPE rv=pow(bias,expn), d=1.0;
    pts[0]=0.0;
    for ( BASSO_IDTYPE i=1; i<n; ++i )
    {
        pts[i] = pts[i-1]+d;
        d=d/rv;
    }
    d=1.0/pts[n-1];
    for ( BASSO_IDTYPE i=0; i<n; ++i )
    {
        pts[i] = pts[i]*d;
        pts[i] = (1.0-pts[i])*x1+pts[i]*x2;
    }
	return pts+n;
}

/**
Computes a 1D array of nodes between the points with a exponential bias.

    
    (x1) ----------- (x2)
    
with \param nx points along the u-direction  

\param x1 - coordinate of first point 
\param x2 - coordinate of the second point
\param nx - number of points
\param bias = 'spacing is x^b (where b is the \param bias)
The return value is a pointer to the next node to be added in pts.
*/
template <class VAL_TYPEx, class VAL_TYPE>
VAL_TYPE *node_array_1d_power( BASSO_IDTYPE n, VAL_TYPEx x1, VAL_TYPEx x2, VAL_TYPE *pts, VAL_TYPEx bias )
{
    node_array_1d_none( n, 0.0, 1.0, pts );
    for ( BASSO_IDTYPE i=0; i<n; ++i ) 
    {
        pts[i] = pow(pts[i],bias);
        pts[i] = (1.0-pts[i])*x1+pts[i]*x2;
    }
    return pts+n;
}

/**
Computes a 1D array of nodes between the points with a bell shaped bias.

    
    (x1) ----------- (x2)
    
with \param nx points along the u-direction  

\param x1 - coordinate of first point 
\param x2 - coordinate of the second point
\param nx - number of points
\param bias = spacing is such that the size of the spacing  in the middle is closer
                    1/2*(tanh(b*(s-1/2))+1) where s is a parametric coord between -1 and 1
The return value is a pointer to the next node to be added in pts.
*/
template <class VAL_TYPEx, class VAL_TYPE>
VAL_TYPE *node_array_1d_bellcurve( BASSO_IDTYPE n, VAL_TYPEx x1, VAL_TYPEx x2, VAL_TYPE *pts, VAL_TYPEx bias )
{
    node_array_1d_none( n, 0.0, 1.0, pts );
    for ( BASSO_IDTYPE i=0; i<n; ++i ) 
        pts[i] = 0.5*(tanh(bias*(pts[i]-0.5))+1.0);
    for ( BASSO_IDTYPE i=0; i<n; ++i ) 
        pts[i] = pts[i]-pts[0];
    VAL_TYPE f=1.0/pts[n-1];
    for ( BASSO_IDTYPE i=0; i<n; ++i ) 
        pts[i] = pts[i]*f;
    for ( BASSO_IDTYPE i=0; i<n; ++i ) 
        pts[i] = (1.0-pts[i])*x1+pts[i]*x2;
	
	return pts+n;
}


/**
Computes a 1D array of nodes between the points
    
    (x1) ----------- (x2)
    
with \param nx points along the u-direction  The node spacing can be biased using the bias1.

\param x1 - coordinate of first point 
\param x2 - coordinate of the second point
\param bias, \param biasfact -  the bias type and value (default is none)

\param bias = 'POWER' spacing is x^b (where b is the \param biasfactor)
\param bias = 'GEOMETRIC' spacing is such that the size of the spacing at x2 is
               1/b of the spacing at x1.
\param bias = 'BELL' spacing is such that the size of the spacing  in the middle is closer
                    1/2*(tanh(b*(s-1/2))+1) where s is a parametric coord between -1 and 1
The return value is a pointer to the next node to be added in pts.
*/
template <class VAL_TYPEx, class VAL_TYPE>
VAL_TYPE *node_array_1d( BASSO_IDTYPE nx, VAL_TYPEx x1, VAL_TYPEx x2, VAL_TYPE *pts, 
                MeshBiasType bias=BiasNONE, VAL_TYPEx biasfact=1.0 )
{
    switch (bias)
    {
        
        case BiasGEOMETRIC:
        return node_array_1d_gometric( nx, x1, x2, pts, biasfact );
 
        
        case BiasPOWER:
        return node_array_1d_power( nx, x1, x2, pts, biasfact );
   
        
        case BiasBELL:
        return node_array_1d_bellcurve( nx, x1, x2, pts, biasfact );

        
        default:   // No Bias
        return node_array_1d_none( nx, x1, x2, pts );
    }
}

/**
Computes a 2D array of nodes (in the z=0 plane) between the points
    
        (x4,y4) ----------- (x3,y3)
          |                    |
         |                    |
        |                    |
       |                    |
    (x1,y1) ----------- (x2,y2)
    
with n1 points along the u-direction (p1-p2) and n2 points in the v-direction
(p1-p4).  The node spacing can be biased using the bias1 and bias2.

The node coordinates are returned in \param pts as follows

pts = [  p1x, p1y, 0, 0, ..., pt2x, pt2y, 0 ,0, ... ]

where the 0s pad if \param ldp is greater than 2.

\param x1, \param y1 - x and y coodinates of first point 
\param x2, \param y2 - x and y coodinates of the second point
\param x3, \param y3 - x and y coodinates of the third point
\param x4, \param y4 - x and y coodinates of the fourth point
\param n1, \param n2 - number of nodes in the u and v directions respectively
\param pts - the array to hold the node coordinates
\param ldp - the leading dimension of \param pts (>=2)
\param bias1, \param b1 -  the biast type and value for the u-direction (default is none)
\param bias2, \param b2 -  the biast type and value for the v-direction (default is none)
The return value is a pointer to the next node to be added in pts.

**/
template <class VAL_TYPEx, class VAL_TYPE>
VAL_TYPE *node_array_2d( VAL_TYPEx x1, VAL_TYPEx y1, 
                VAL_TYPEx x2, VAL_TYPEx y2, 
                VAL_TYPEx x3, VAL_TYPEx y3, 
                VAL_TYPEx x4, VAL_TYPEx y4,
                BASSO_IDTYPE n1, BASSO_IDTYPE n2, VAL_TYPE *pts, int ldp=2,
                MeshBiasType bias1=BiasNONE, VAL_TYPE bf1=1.0, 
                MeshBiasType bias2=BiasNONE, VAL_TYPE bf2=1.0 )
{
   // BASSO_IDTYPE nn=n1*n2;
    VAL_TYPE xis[n1], etas[n2], one=1.0;
    node_array_1d( n1, -one, one, xis,  bias1, bf1 );
    node_array_1d( n2, -one, one, etas, bias2, bf2 );
    
    VAL_TYPE *ptr=pts, xpt, ypt, Ns1, Ns2, Ns3, Ns4, xi, eta;
    for ( BASSO_IDTYPE j=0; j<n2; ++j )
    {
        eta=etas[j];
        for ( BASSO_IDTYPE i=0; i<n1; ++i )
        {
            xi=xis[i];
            Ns1 = 0.25*(1.0-xi)*(1.0-eta);
            Ns2 = 0.25*(1.0+xi)*(1.0-eta);
            Ns3 = 0.25*(1.0+xi)*(1.0+eta);
            Ns4 = 0.25*(1.0-xi)*(1.0+eta);
            xpt = Ns1*x1 + Ns2*x2 + Ns3*x3 + Ns4*x4;
            ypt = Ns1*y1 + Ns2*y2 + Ns3*y3 + Ns4*y4;
            *ptr=xpt; ++ptr; 
            *ptr=ypt; ptr += ldp-1; 
        }
    }
	return ptr;
}
//**********************************************************
template <class VAL_TYPEx, class VAL_TYPE>
VAL_TYPE *node_array_2d( VAL_TYPEx x1, VAL_TYPEx y1, VAL_TYPEx z1, 
                VAL_TYPEx x2, VAL_TYPEx y2, VAL_TYPEx z2, 
                VAL_TYPEx x3, VAL_TYPEx y3, VAL_TYPEx z3, 
                VAL_TYPEx x4, VAL_TYPEx y4, VAL_TYPEx z4,
                BASSO_IDTYPE n1, BASSO_IDTYPE n2, VAL_TYPE *pts, int ldp=2,
                MeshBiasType bias1=BiasNONE, VAL_TYPEx bf1=1.0, 
                MeshBiasType bias2=BiasNONE, VAL_TYPEx bf2=1.0 )
{
    //BASSO_IDTYPE nn=n1*n2;
    VAL_TYPEx xis[n1], etas[n2], one=1.0, fourth=0.25;
    node_array_1d( n1, -one, one, xis,  bias1, bf1 );
    node_array_1d( n2, -one, one, etas, bias2, bf2 );

    VAL_TYPE *ptr=pts, xpt, ypt, zpt, Ns1, Ns2, Ns3, Ns4, xi, eta;
    for ( BASSO_IDTYPE j=0; j<n2; ++j )
    {
        eta=etas[j];
        for ( BASSO_IDTYPE i=0; i<n1; ++i )
        {
            xi=xis[i];
            Ns1 = fourth*(one-xi)*(one-eta);
            Ns2 = fourth*(one+xi)*(one-eta);
            Ns3 = fourth*(one+xi)*(one+eta);
            Ns4 = fourth*(one-xi)*(one+eta);
            xpt = Ns1*x1 + Ns2*x2 + Ns3*x3 + Ns4*x4;
            ypt = Ns1*y1 + Ns2*y2 + Ns3*y3 + Ns4*y4;
            zpt = Ns1*z1 + Ns2*z2 + Ns3*z3 + Ns4*z4;
			//cout << "(" << xi << "," << eta << "):\t\t" << xpt << " " << ypt << " " << zpt << endl;
            *ptr=xpt; ++ptr; 
            *ptr=ypt; ++ptr; 
            *ptr=zpt; ptr += ldp-2; 
        }
    }
	return ptr;
}

/**
Computes a 2D triangular array of nodes (in the z=0 plane) between the points

                 (x3, y3)
                 /    \
               /       \
             /          \
           /             \
         /                \    
       /                   \
    (x1,y1) ------------(x2,y2)
    
with n1 points along the u-direction (p1-p2) and n2 points in the v-direction
(p1-p3).  The node spacing can be biased using the bias1 and bias2.

The node coordinates are returned in \param pts as follows

pts = [  p1x, p1y, 0, 0, ..., pt2x, pt2y, 0 ,0, ... ]

where the 0s pad if \param ldp is greater than 2.

\param x1, \param y1 - x and y coodinates of first point 
\param x2, \param y2 - x and y coodinates of the second point
\param x3, \param y3 - x and y coodinates of the third point
\param n1, \param n2 - number of nodes in the u and v directions respectively
\param pts - the array to hold the node coordinates
\param ldp - the leading dimension of \param pts (>=2)
\param bias1, \param b1 -  the biast type and value for the u-direction (default is none)
\param bias2, \param b2 -  the biast type and value for the v-direction (default is none)
The return value is a pointer to the next node to be added in pts.

**/
template <class VAL_TYPEx, class VAL_TYPE>
VAL_TYPE *node_array_2d( VAL_TYPEx x1, VAL_TYPEx y1, 
                VAL_TYPEx x2, VAL_TYPEx y2, 
                VAL_TYPEx x3, VAL_TYPEx y3,
                BASSO_IDTYPE n1, BASSO_IDTYPE n2, VAL_TYPE *pts, int ldp=2,
                MeshBiasType bias1=BiasNONE, VAL_TYPEx bf1=1.0, 
                MeshBiasType bias2=BiasNONE, VAL_TYPEx bf2=1.0 )
{
    BASSO_IDTYPE nn=n1*n2;
	VAL_TYPE xis[n1], etas[n2], one=1.0;;
    node_array_1d( n1, -one, one, xis,  bias1, bf1 );
    node_array_1d( n2, -one, one, etas, bias2, bf2 );
    
    VAL_TYPE *ptr=pts, xpt, ypt, Ns1, Ns2, Ns3, xi, eta;
    for ( BASSO_IDTYPE j=0; j<n2; ++j )
    {
        eta=etas[j];
        for ( BASSO_IDTYPE i=0; i<n1; ++i )
        {
            xi=xis[i];
    
            Ns1 = 1.0-xi-eta;
            Ns2 = xi;
            Ns3 = eta;
            xpt = Ns1*x1 + Ns2*x2 + Ns3*x3;
            ypt = Ns1*y1 + Ns2*y2 + Ns3*y3;
            *ptr=xpt; ++ptr; 
            *ptr=ypt; ptr += ldp-1; 
        
        }
    }
	return ptr;
}

/**
Computes a 3D array of nodes with n1 points along the u-direction (p1-p2), 
n2 points in the v-direction (p1-p4) and n3 points in the w-direction (p1-p5).
The node spacing can be biased using the bias1, bias2 and bias3.

The node coordinates are returned in \param pts as follows

pts = [  p1x, p1y, p1z, 0, 0, ..., pt2x, pt2y, pt2z, 0 ,0, ... ]

where the 0s pad if \param ldp is greater than 3.

\param x1, \param y1, \param z1 - x, y and z coodinates of first point 
\param x2, \param y2, \param z2 - x, y and z coodinates of the second point
\param x3, \param y3, \param z3 - x, y and z coodinates of the third point
\param x4, \param y4, \param z4 - x, y and z coodinates of the fourth point
\param n1, \param n2, \param n3- number of nodes in the u, v and w directions respectively
\param pts - the array to hold the node coordinates
\param ldp - the leading dimension of \param pts (>=3)
\param bias1, \param b1 -  the biast type and value for the u-direction (default is none)
\param bias2, \param b2 -  the biast type and value for the v-direction (default is none)
The return value is a pointer to the next node to be added in pts.

**/
template <class VAL_TYPEx, class VAL_TYPE>
VAL_TYPE *node_array_3d( VAL_TYPEx x1, VAL_TYPEx y1, VAL_TYPEx z1, 
                VAL_TYPEx x2, VAL_TYPEx y2, VAL_TYPEx z2, 
                VAL_TYPEx x3, VAL_TYPEx y3, VAL_TYPEx z3, 
                VAL_TYPEx x4, VAL_TYPEx y4, VAL_TYPEx z4, 
                VAL_TYPEx x5, VAL_TYPEx y5, VAL_TYPEx z5, 
                VAL_TYPEx x6, VAL_TYPEx y6, VAL_TYPEx z6, 
                VAL_TYPEx x7, VAL_TYPEx y7, VAL_TYPEx z7, 
                VAL_TYPEx x8, VAL_TYPEx y8, VAL_TYPEx z8,
                BASSO_IDTYPE n1, BASSO_IDTYPE n2, BASSO_IDTYPE n3, VAL_TYPE *pts, int ldp=3,
                MeshBiasType bias1=BiasNONE, VAL_TYPEx bf1=1.0, 
                MeshBiasType bias2=BiasNONE, VAL_TYPEx bf2=1.0, 
                MeshBiasType bias3=BiasNONE, VAL_TYPEx bf3=1.0 )
{
    BASSO_IDTYPE nn=n1*n2*n3;
    VAL_TYPE xis[n1], etas[n2], zetas[n3];
	VAL_TYPE one = 1.0;
    node_array_1d( n1, -one, one, xis,   bias1, bf1 );
    node_array_1d( n2, -one, one, etas,  bias2, bf2 );
    node_array_1d( n3, -one, one, zetas, bias3, bf3 );
    
    VAL_TYPE *ptr=pts, xpt, ypt, zpt, 
        Ns1, Ns2, Ns3, Ns4, Ns5, Ns6, Ns7, Ns8, xi, eta, zeta;
    for ( BASSO_IDTYPE k=0; k<n3; ++k )
    {
        zeta=zetas[k];
        for ( BASSO_IDTYPE j=0; j<n2; ++j )
        {
            eta=etas[j];
            for ( BASSO_IDTYPE i=0; i<n1; ++i )
            {
                xi=xis[i];
                Ns1 = 0.125*(1.0-xi)*(1.0-eta)*(1.0-zeta);
                Ns2 = 0.125*(1.0+xi)*(1.0-eta)*(1.0-zeta);
                Ns3 = 0.125*(1.0+xi)*(1.0+eta)*(1.0-zeta);
                Ns4 = 0.125*(1.0-xi)*(1.0+eta)*(1.0-zeta);
                Ns5 = 0.125*(1.0-xi)*(1.0-eta)*(1.0+zeta);
                Ns6 = 0.125*(1.0+xi)*(1.0-eta)*(1.0+zeta);
                Ns7 = 0.125*(1.0+xi)*(1.0+eta)*(1.0+zeta);
                Ns8 = 0.125*(1.0-xi)*(1.0+eta)*(1.0+zeta);
                xpt = Ns1*x1 + Ns2*x2 + Ns3*x3 + Ns4*x4 + Ns5*x5 + Ns6*x6 + Ns7*x7 + Ns8*x8;
                ypt = Ns1*y1 + Ns2*y2 + Ns3*y3 + Ns4*y4 + Ns5*y5 + Ns6*y6 + Ns7*y7 + Ns8*y8;
                zpt = Ns1*z1 + Ns2*z2 + Ns3*z3 + Ns4*z4 + Ns5*z5 + Ns6*z6 + Ns7*z7 + Ns8*z8;
                *ptr=xpt; ++ptr; 
                *ptr=ypt; ++ptr; 
                *ptr=zpt; ptr += ldp-2; 
            }
        }
    }
	return ptr;
    
}

//--------------------------------------------------------------------------

/**
Uses a connectivity pattern, \param connPtrn, to generate an element connectivity
in 1D.  

\param connPtrn - the connectivity of the first element (on return this will have
                  the connectivity of the connectivity of on past the last element)
\param nne - the number of nodes in \param connPtrn
\param num_u - the number of elements to be constructed in the u-direction
\param conn - on return will have the element connectivity
\param ldc - the leading dimension of conn (the default is nne)
\param inc_u - the amount that the \param connPtrn will be incremented each step
               in the u-direction.
The return value is a pointer to the next element that would be inserted in conn.
      
The connectivity ir returnd in \param conn as follows 

    conn = [ e1n1, e1n2, ..., e1nne, 0, 0, ..., e2n1, e2n2, ..., e2nne, 0, 0, ...,
              nen1, nen2, ..., nenne, 0, 0, ...  ]
      
where einj is the jth node of element i.  The zeros are padding if ldc>nne.              
**/
BASSO_IDTYPE *gen_conn_1d( BASSO_IDTYPE *connPtrn, int nne, BASSO_IDTYPE num_u, BASSO_IDTYPE *conn, int ldc=-1, int inc_u=1 )
{
    if ( ldc<nne ) ldc=nne;
    
    BASSO_IDTYPE *nPtr=conn;
    for ( BASSO_IDTYPE c=0; c<num_u; ++c )
    {
        for ( int n=0; n<nne; ++n )
        {
            *nPtr = connPtrn[n];
            ++nPtr; 
            connPtrn[n] += inc_u;
        }
        nPtr += ldc-nne;
    }
	return nPtr;
}
/**
Uses a connectivity pattern, \param connPtrn, to generate an element connectivity
in 2D.  

\param connPtrn - the connectivity of the first element (on return this will have
                  the connectivity of the connectivity of on past the last element)
\param nne - the number of nodes in \param connPtrn
\param num_u - the number of elements to be constructed in the u-direction
\param num_v - the number of elements to be constructed in the v-direction
\param conn - on return will have the element connectivity
\param ldc - the leading dimension of conn (the default is nne)
\param inc_u - the amount that the \param connPtrn will be incremented each step
               in the u-direction.
\param inc_v - the amount that the \param connPtrn will be incremented each step
               in the v-direction.

The return value is a pointer to the next element that would be inserted in conn.

The connectivity ir returnd in \param conn as follows 

    conn = [ e1n1, e1n2, ..., e1nne, 0, 0, ..., e2n1, e2n2, ..., e2nne, 0, 0, ...,
              nen1, nen2, ..., nenne, 0, 0, ...  ]
      
where einj is the jth node of element i.  The zeros are padding if ldc>nne.              
**/
BASSO_IDTYPE *gen_conn_2d( BASSO_IDTYPE *connPtrn, int nne, BASSO_IDTYPE num_u, BASSO_IDTYPE num_v, BASSO_IDTYPE *conn, int ldc=-1,
        int inc_u=1, int inc_v=2 )
{
    if ( ldc<nne ) ldc=nne;
    
    BASSO_IDTYPE *nPtr=conn;
    for ( BASSO_IDTYPE r=0; r<num_v; ++r )
    {
        for ( BASSO_IDTYPE c=0; c<num_u; ++c )
        {
            for ( int n=0; n<nne; ++n )
            {
                *nPtr = connPtrn[n];
                ++nPtr; 
                connPtrn[n] += inc_u;
            }
            nPtr += ldc-nne;
        }
        for ( int n=0; n<nne; ++n )
            connPtrn[n] += inc_v-inc_u;
    }
	return nPtr;
}

/**
Uses a connectivity pattern, \param connPtrn, to generate an element connectivity
in 3D.  

\param connPtrn - the connectivity of the first element (on return this will have
                  the connectivity of the connectivity of on past the last element)
\param nne - the number of nodes in \param connPtrn
\param num_u - the number of elements to be constructed in the u-direction
\param num_v - the number of elements to be constructed in the v-direction
\param num_w - the number of elements to be constructed in the w-direction
\param conn - on return will have the element connectivity
\param ldc - the leading dimension of conn (the default is nne)
\param inc_u - the amount that the \param connPtrn will be incremented each step
               in the u-direction.
\param inc_v - the amount that the \param connPtrn will be incremented each step
               in the v-direction.
\param inc_w - the amount that the \param connPtrn will be incremented each step
               in the w-direction.
			   
The return value is a pointer to the next element that would be inserted in conn.
               
The connectivity ir returnd in \param conn as follows 

    conn = [ e1n1, e1n2, ..., e1nne, 0, 0, ..., e2n1, e2n2, ..., e2nne, 0, 0, ...,
              nen1, nen2, ..., nenne, 0, 0, ...  ]
      
where einj is the jth node of element i.  The zeros are padding if ldc>nne.              
**/
BASSO_IDTYPE *gen_conn_3d( BASSO_IDTYPE *connPtrn, int nne, BASSO_IDTYPE num_u, BASSO_IDTYPE num_v, BASSO_IDTYPE num_w,
    BASSO_IDTYPE *conn, int ldc=-1, int inc_u=1, int inc_v=2, int inc_w=-999999999 )
{
    if ( ldc<nne ) ldc=nne;
	if ( inc_w == -999999999 ) inc_w = num_u + 3;
        
    BASSO_IDTYPE *nPtr=conn;
    
    for ( BASSO_IDTYPE l=0; l<num_w; ++l )
    {
        for ( BASSO_IDTYPE r=0; r<num_v; ++r )
        {   
            for ( BASSO_IDTYPE c=0; c<num_u; ++c )
            {
                for ( int n=0; n<nne; ++n )
                {
                    *nPtr = connPtrn[n];
                    ++nPtr; 
                    connPtrn[n] += inc_u;
                }
                nPtr += ldc-nne;
            }
            for ( int n=0; n<nne; ++n )
                connPtrn[n] += inc_v-inc_u;
        }
        for ( int n=0; n<nne; ++n )
            connPtrn[n] += inc_w-inc_v;
    } 
	return nPtr;	
}


//  --------------------   R E N U M B E R I N G     F U N C T I O N S   -------------------- //

/** 

 Renumbers the nodes so that they are zero (default) offset and continous.  Nodes
 that are not referenced in element are removed.  For this to work correctly
 it is assumed that the nids are in assending order.

	\param element - element connectivity matrix.  On return will have the renumbered 
	                connectivity.
	\param ne - number of elements in  element
	\param nne - number of nodes in each element in   element
	\param lde - the leading dimension of  element (pads each element connectivity)
	\param nid - a vector of the node ids (this is unchanged)
	\param nn - the number of nodes 
	\param node - the array of node coordinates.  On return this will have the eliminated 
	            nodes removed and the other node coordinates shifted to fill in these 
	            removed nodes.
	\param sdim - the spacial dimension of node (sdim=2 (x,y) sdim=3 (x,y,z)) 
	\param ldn - the leading dimension of node (pads the dimesion of each node)
	\param nodeMap - on return will have a map from the old node numbers to the new 
	                offset node numbering with the removed nodes. 
	\param offset - the offset to be used on the node numbering (the default is zero).  If this
	                is set to something other than zero then the connectivity and the map will 
	                be affected (i.e. they will not "point" directly to the location in node)

    teh return value is the number of nodes returned in node that are used in element.
**/
template <class VAL_TYPE>
BASSO_IDTYPE remove_unused_nodes( BASSO_IDTYPE *element, BASSO_IDTYPE ne, int nne, int lde, 
	const BASSO_IDTYPE *nid, BASSO_IDTYPE nn, VAL_TYPE *node, int sdim, int ldn, map<BASSO_IDTYPE,BASSO_IDTYPE> &nodeMap, int offset=0 )
{
	// set up the nodeMap with the nodes in element
	BASSO_IDTYPE *nptr=element;
	for ( BASSO_IDTYPE e=0; e<ne; ++e, nptr += lde-nne )
		for ( int i=0; i<nne; ++i, ++nptr )
			nodeMap[ *nptr ] = -1;
			
	// go back and renumber
	map<BASSO_IDTYPE,BASSO_IDTYPE>::iterator mitr;
	BASSO_IDTYPE n=offset-1;
	for ( mitr=nodeMap.begin(); mitr!=nodeMap.end(); ++mitr )
		mitr->second = ++n;
	BASSO_IDTYPE newNumNode = n+1;
		
	// now remove nodes that are not used
	//BASSO_IDTYPE o=0;
	for ( n=0; n<nn; ++n )
	{
		mitr = nodeMap.find( nid[n] );
		if ( mitr != nodeMap.end() ) // this node is used
		{
			if ( mitr->second != n ) // we need to move this node
			{
				BASSO_IDTYPE nnew = mitr->second;
				node[ nnew*ldn ] = node[ n*ldn ];
				node[ nnew*ldn + 1 ] = node[ n*ldn + 1 ];
				node[ nnew*ldn + 2 ] = node[ n*ldn + 2 ];
			}
		}
	}
	
	// renumber the connectivity
	nptr=element;
	for ( BASSO_IDTYPE e=0; e<ne; ++e, nptr += lde-nne )
		for ( int i=0; i<nne; ++i, ++nptr )
			*nptr = nodeMap[ *nptr ];
		
	return newNumNode;
}

/**
Constructs an invers node map from a node id array.  
\param nid - a vector of the node ids (this is unchanged)
\param nn - the number of nodes 
\param invMap - on return will have a map from the old node numbers to the new 
	               offset node numbering with the removed nodes. 
	
so invMap[ old nod id ] => directly to the row of the node coordinate matrix.
*/
void inv_nodeid_map( const BASSO_IDTYPE *nid, BASSO_IDTYPE nn, map<BASSO_IDTYPE,BASSO_IDTYPE> &invMap, int offset=0 )
{
	const BASSO_IDTYPE *idptr = nid;
	for ( BASSO_IDTYPE i=0; i<nn; ++i, ++idptr )
		invMap[ *idptr ] = i + offset;
}



/**
Renumbers the connectivity so that it points to the row of the node coordinate matrix directly

	\param element - element connectivity matrix.  On return will have the renumbered 
	                connectivity.
	\param ne - number of elements in   element
	\param nne - number of nodes in each element in  element
	\param lde - the leading dimension of  lement (pads each element connectivity)
	\param nid - a vector of the node ids (this is unchanged)
	\param nn - the number of nodes 
	\param offset - the offset to be used on the node numbering (the default is zero).  If this
	                is set to something other than zero then the connectivity and the map will 
	                be affected (i.e. they will not "point" directly to the location in node)


*/
void renumber_connectivity( BASSO_IDTYPE *element, BASSO_IDTYPE ne, int nne, int lde, const BASSO_IDTYPE *nid, BASSO_IDTYPE nn, int offset=0  )
{
	map<BASSO_IDTYPE,BASSO_IDTYPE> invMap;
	inv_nodeid_map( nid, nn, invMap, offset );
	
	// renumber the connectivity
	BASSO_IDTYPE *nptr=element;
	for ( BASSO_IDTYPE e=0; e<ne; ++e, nptr += lde-nne )
		for ( int i=0; i<nne; ++i, ++nptr )
			*nptr = invMap[ *nptr ];
		
}

/**
Renumbers the connectivity so that it points to the row of the node coordinate matrix directly

	\param element - element connectivity matrix.  On return will have the renumbered 
	                connectivity.
	\param ne - total number of entries in element to be renumbered
	\param nid - a vector of the node ids (this is unchanged)
	\param nn - the number of nodes 
	\param offset - the offset to be used on the node numbering (the default is zero).  If this
	                is set to something other than zero then the connectivity and the map will 
	                be affected (i.e. they will not "point" directly to the location in node)


*/
void renumber_connectivity( BASSO_IDTYPE *element, BASSO_IDTYPE ne, const BASSO_IDTYPE *nid, BASSO_IDTYPE nn, int offset=0  )
{
	map<BASSO_IDTYPE,BASSO_IDTYPE> invMap;
	inv_nodeid_map( nid, nn, invMap, offset );
	
	// renumber the connectivity
	BASSO_IDTYPE *nptr=element;
	for ( BASSO_IDTYPE e=0; e<ne; ++e, ++nptr )
		*nptr = invMap[ *nptr ];
		
}


}
#endif

