// std includes
#include <iostream>
#include <sstream>


// basso includes
#include "Basso_defs.h"
#include "Basso_Numeric.h"
#include "Basso_Point.h"
#include "Basso_MeshBlockHexa8.h"
#include "Basso_iMatrix.h"
#include "Basso_nMatrix.h"
#include "Basso_ParentElement.h"
#include "Basso_write_ensight.h"
#include "Basso_Hexa8.h"
#include "Basso_Quad4.h"
#include "Basso_nfeaops.h"
#include "Basso_DOFMap.h"

#include "Basso_PointBCSet.h"
#include "Basso_FaceBCSet.h"

using namespace std;
using namespace Basso;


/*
*/
int main(int argc, char* argv[])
{
	cout << "Set up a mesh and a dofmap\n";
	
	// Define the finite element mesh
	Basso_Point p0( 0.0, 0.0, 0.0 ), p1( 10.0, 10.0, 10.0 );
	Basso_MeshBlockHexa8<Basso_Numeric> mesh(p0,p1);
	mesh.SetElementSize(5.0);
	
	int nn = mesh.NumNodes(), ne=mesh.NumElements();
	Basso_nMatrix nodes(3,nn);
	Basso_iMatrix conn(8,ne);
	Basso_Array<BASSO_IDTYPE> gnids(nn);
		
	mesh.SetMesh(nodes,conn);
	
	// set up the global dof map
	Basso_DOFMap GDOFMap(nn,3);
	
	// compute the rhs
	Basso_nVector rhs(3*nn); 
	Basso_Array2D<BASSO_IDTYPE> faceConn;
	mesh.SideSegs(faceConn,Basso_MeshBlockHexa8<BASSO_IDTYPE>::kPos1Face);
	
	Basso_Quad4 quad4;
	Basso_nVector ff=string{"10.0 0.0 0.0"};
	Basso_iVector ldof=string("0 1 2");
	Basso_FaceBCSet trac(faceConn,&quad4,ff,ldof); 
	trac.AddForce(nodes,GDOFMap,rhs);
	cout << faceConn << "\n" << rhs << "\n";
	

	// Get the fixed dofs
	Basso_PointBCSet spcs;
	set<BASSO_IDTYPE> nfix;
	mesh.SideNodes(nfix,Basso_MeshBlockHexa8<BASSO_IDTYPE>::kNeg1Face);
	spcs.AddPointBCs(nfix,0);
	nfix.clear();
	mesh.SideNodes(nfix,Basso_MeshBlockHexa8<BASSO_IDTYPE>::kNeg2Face);
	spcs.AddPointBCs(nfix,1);
	nfix.clear();
	mesh.SideNodes(nfix,Basso_MeshBlockHexa8<BASSO_IDTYPE>::kNeg3Face);
	spcs.AddPointBCs(nfix,2);
	cout << spcs < "\n";
	
	return 0;
}