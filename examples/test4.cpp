#include <iostream>

#include "Basso_defs.h"

#include "Basso_Point1.h"
#include "Basso_Line2.h"
#include "Basso_Line3.h"
#include "Basso_Tria3.h"
#include "Basso_Tria6.h"
#include "Basso_Quad4.h"
#include "Basso_Quad8.h"
#include "Basso_Tetra4.h"
#include "Basso_Tetra10.h"
#include "Basso_Hexa8.h"

#include "Basso_iMatrix.h"
#include "Basso_nfeaops.h"

using namespace std;
using namespace Basso;

int main(int argc, char* argv[])
{
	cout << "Test of the Element classes in Basso\n";
	
	Basso_Quad4 quad4;
	
	Basso_nMatrix nodes(2,6);
	nodes[0][0] = 0.0; nodes[0][1] = 0.0; 
	nodes[1][0] = 1.0; nodes[1][1] = 0.0; 
	nodes[2][0] = 2.0; nodes[2][1] = 0.0; 
	nodes[3][0] = 0.0; nodes[3][1] = 1.0; 
	nodes[4][0] = 1.0; nodes[4][1] = 1.0; 
	nodes[5][0] = 2.0; nodes[5][1] = 1.0; 
	
	cout << nodes << "\n";
	
	Basso_iMatrix conn(4,2);
	conn[0][0]=0; conn[0][1]=1; conn[0][2]=4; conn[0][3]=3;
	conn[1][0]=1; conn[1][1]=2; conn[1][2]=5; conn[1][3]=4;
	
	cout << conn << "\n";
	
	Basso_nMatrix ecoord(2,4), bmat(3,8), kmat(8,8), cmat(3,3);
	Basso_Numeric E=10e6, nu=.3, jac;
	Basso_QuadratureRule qrule;
	quad4.Quadrature(2,qrule);
	cmat_pstrain( cmat, E, nu );
	cout << cmat << "\n";
	for ( int e=0; e<conn.N(); ++e )
	{
		element_coordinates( nodes, conn[e], conn.M(), ecoord );
		cout << "\n-------------- e=" << e << " --------------\n" << ecoord << "\n";
		kmat.Zero();
		Basso_QuadratureRule::ConstIterator qitr;
		for ( qitr=qrule.Begin(); qitr!=qrule.End(); ++qitr )
		{
			jac = quad4.Bmatrix( qitr->Pt(), ecoord, bmat );
			cout << bmat << "\n";
			multbtcb( 'n', jac*(qitr->Wt()), bmat, cmat, 1.0, kmat );
		}
		cout << " kmat:\n" << kmat << "\n";
	}
	
	return 0;
}