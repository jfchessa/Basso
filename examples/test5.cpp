#include <iostream>

#include "Basso_Numeric.h"
#include "Basso_Point.h"
#include "Basso_MeshBlockHexa8.h"
#include "Basso_iMatrix.h"
#include "Basso_nMatrix.h"

using namespace std;
using namespace Basso;

int main(int argc, char* argv[])
{
	cout << "Test of the meshing classes in Basso\n";
	
	Basso_MeshBlockHexa8<Basso_Numeric> mesh;
	
	Basso_Point p0( 0.0, 0.0, 0.0 ), p1( 10.0, 5.0, 3.0 );
	mesh.SetCorners(p0,p1);
	mesh.SetNumNodes(3,2,2);
	
	int nn = mesh.NumNodes(), ne=mesh.NumElements();
	Basso_nMatrix nodes(3,nn);
	Basso_iMatrix conn(8,ne);
	
	mesh.SetMesh(nodes,conn);
	
	cout << nodes << "\n" << conn << "\n";
	
	return 0;
}