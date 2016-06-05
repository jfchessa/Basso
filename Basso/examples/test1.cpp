#include <iostream>

#include "Basso_defs.h"
#include "Basso_shape_functions.h"
#include "Basso_feaops.h"

using namespace std;
using namespace Basso;

int main(int argc, char* argv[])
{
	cout << "Find an element jacobian matrix for CST element\n";
	
	double nodeCoord[] = { 0.0, 0.0,    3.0, 0.0,   0.5, 2.0 };
	double dNa[6], jmat[4];
	
	dshape_tria3( dNa, 3 );
	
	double *nptr=dNa;
	for ( int i=0; i<2; ++i )
	{
		for ( int I=0; I<3; ++I, ++nptr )
			cout << *nptr << " ";
		cout << "\n";
	}
	
	element_jacobian( 3, 2, 2, nodeCoord, 2, dNa, 3, jmat, 2 );
				
	nptr=jmat;
	for ( int i=0; i<2; ++i )
	{
		for ( int j=0; j<2; ++j, ++nptr )
			cout << *nptr << " ";
		cout << "\n";
	}		 
	
	return 0;
}