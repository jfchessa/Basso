// std includes
#include <iostream>
#include <sstream>

// Trilinos includes
#include "Epetra_Version.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

// basso includes
#include "Basso_defs.h"
#include "Basso_Numeric.h"
#include "Basso_Point.h"
#include "Basso_PMeshBlockHexa8.h"
#include "Basso_iMatrix.h"
#include "Basso_nMatrix.h"
#include "Basso_ParentElement.h"
#include "Basso_Hexa8.h"
#include "Basso_Quad4.h"
#include "Basso_nfeaops.h"

#include "Basso_PointBCSet.h"
#include "Basso_FaceBCSet.h"


// tbasso includes
#include "TBasso_DOFMap.h"
#include "TBasso_FECrsMatrix.h"
#include "TBasso_FEVector.h"
#include "TBasso_FESolver.h"
#include "TBasso_PVTK.h"

using namespace std;
using namespace Basso;


/*
  This example requires 8 processors
*/
int main(int argc, char* argv[])
{
	cout << "Set up a mesh and a dofmap\n";

#ifdef HAVE_MPI
	MPI_Init(&argc,&argv);
	Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
	Epetra_SerialComm comm;
#endif
	
	// Define the finite element mesh
	Basso_Point p0( 0.0, 0.0, 0.0 ), p1( 10.0, 5.0, 2.5 );
	Basso_PMeshBlockHexa8<Basso_Numeric> mesh(p0,p1,comm.MyPID());
	mesh.SetElementSize(0.1);
	mesh.SetNumPartitions(2,1,1);
	
	if ( comm.NumProc() != mesh.NumPartitions() )
	{
		cout << "Requires " << mesh.NumPartitions() << " processes\n";
		return 1;
	}
	int nn = mesh.NumNodes(), ne=mesh.NumElements();
	Basso_nMatrix nodes(3,nn);
	Basso_iMatrix conn(8,ne);
	Basso_Array<BASSO_IDTYPE> gnids(nn);
		
	mesh.SetMesh(nodes,conn);
	mesh.GlobalNodeIDs(gnids);
	
	// set up the global dof map
	TBasso_DOFMap GDOFMap(comm,gnids,3);
	
	// Compute the global stiffness matrix
 	TBasso_FECrsMatrix Kmat(GDOFMap,conn);

	Basso_Hexa8 hexa8;
	Basso_nMatrix ecoord(3,8), bmat(6,3*8), ke(3*8,3*8), cmat(6,6);
	Basso_Array<BASSO_IDTYPE> sctr(3*8);
	Basso_Numeric E=10e6, nu=.3, jac;
	Basso_QuadratureRule qrule;
	hexa8.Quadrature(2,qrule);
	cmat_3d( cmat, E, nu );
	for ( int e=0; e<conn.NumCols(); ++e )
	{
		element_coordinates( nodes, conn[e], conn.NumRows(), ecoord );
		ke.Zero();
		GDOFMap.SetScatter(conn[e],conn.NumRows(),sctr);
		Basso_QuadratureRule::ConstIterator qitr;
		for ( qitr=qrule.Begin(); qitr!=qrule.End(); ++qitr )
		{
			jac = hexa8.Bmatrix( qitr->Pt(), ecoord, bmat );
			multbtcb( 'n', jac*(qitr->Wt()), bmat, cmat, 1.0, ke );
		}
		Kmat.Scatter(ke,sctr);
	}
	Kmat.FillComplete();
	
	// compute the rhs
	TBasso_FEVector rhs(GDOFMap); 
	Basso_Array2D<BASSO_IDTYPE> faceConn;
	mesh.SideSegs(faceConn,Basso_MeshBlockHexa8<BASSO_IDTYPE>::kNeg1Face);
	Basso_Quad4 quad4;
	Basso_nVector ff=string{"10.0 0.0 0.0"};
	Basso_iVector ldof=string("0 1 2");
	Basso_FaceBCSet trac(faceConn,&quad4,ff,ldof); 
	trac.AddForce(nodes,GDOFMap,rhs);
	
	rhs.FillComplete();
	
	// Get the fixed dofs
	set<BASSO_IDTYPE> nfix;
	mesh.SideNodes(nfix,Basso_MeshBlockHexa8<BASSO_IDTYPE>::kPos1Face);
	Basso_PointBCSet spcs;
	spcs.AddPointBCs(nfix,0);
	spcs.AddPointBCs(nfix,1);
	spcs.AddPointBCs(nfix,2);
	
	// solve the system
	TBasso_FEVector disp(GDOFMap); 
	TBasso_FESolver solver(Kmat,disp,rhs,GDOFMap,spcs);
	
	solver.Solve();
	
	TBasso_LocalVector dlocal( GDOFMap );
	dlocal.PullValues(disp);
	
	// write the results
	TBasso_VTKPUnstructuredGrid results(comm,"test1");
	results.SetMesh(nodes,conn,&hexa8);
	Basso_VTKFloatDataArray ddata( "displacement", dlocal.Data(), nn, 3 );
	results.AddPointData(&ddata);
	results.WriteFile();
	
#ifdef HAVE_MPI
	MPI_Finalize() ; 
#endif
	
	return 0;
}