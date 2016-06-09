#include <iostream>
#include <sstream>

#include "Basso_Numeric.h"
#include "Basso_Point.h"
#include "Basso_PMeshBlockHexa8.h"
#include "Basso_iMatrix.h"
#include "Basso_nMatrix.h"

#include "Basso_write_ensight.h"

using namespace std;
using namespace Basso;

int main(int argc, char* argv[])
{
	cout << "Test of the meshing classes in Basso and writing to ensight\n";
	
	Basso_PMeshBlockHexa8<Basso_Numeric> mesh;
	
	Basso_Point p0( 0.0, 0.0, 0.0 ), p1( 10.0, 5.0, 3.0 );
	mesh.SetCorners(p0,p1);
	mesh.SetElementSize(0.5);
	mesh.SetNumPartitions(4,2,2);
	
	Basso_nMatrix nodes;
	Basso_iMatrix conn;
	Basso_Array<BASSO_IDTYPE> gnids;
	for ( int p=0; p<mesh.NumPartitions(); ++p )
	{
		
		mesh.SetPID(p);
		int nn = mesh.NumNodes(), ne=mesh.NumElements();
		nodes.Resize(3,nn);
		conn.Resize(8,ne);
		gnids.Resize(nn);
		mesh.SetMesh(nodes,conn);
		mesh.GlobalNodeIDs(gnids);
		
		cout << gnids << "\n";
		
		ostringstream ss;
		ss << p;
		string filename = "test6_"+ss.str();
		ensight_ascii_fegeometry( filename+".geom", nodes.Data(), nodes.N(), nodes.M(), nodes.LDA(),
        conn.Data(), conn.N(), conn.M(), conn.LDA(), "hexa8" );
		ensight_case( filename, filename+".geom" ); 
	}
	
	
	return 0;
}