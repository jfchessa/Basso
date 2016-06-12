// std includes
#include <iostream>
#include <sstream>

// mpi includes
#include "mpi.h"

// Trilinos includes
#include "Epetra_Version.h"
#ifdef HAVE_MPI
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
#include "Basso_write_ensight.h"
#include "Basso_Hexa8.h"
#include "Basso_nfeaops.h"

// tbasso includes
#include "TBasso_DOFMap.h"
#include "TBasso_FECrsMatrix.h"
//#include "TBasso_FEVector.h"

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
	
	if ( comm.NumProc() != 8 )
	{
		cout << "Requires 8 cpus\n";
		return 1;
	}
	
	// Define the finite element mesh
	Basso_Point p0( 0.0, 0.0, 0.0 ), p1( 10.0, 10.0, 10.0 );
	Basso_PMeshBlockHexa8<Basso_Numeric> mesh(p0,p1,comm.MyPID());
	mesh.SetElementSize(1.0);
	mesh.SetNumPartitions(2,2,2);
	
	int nn = mesh.NumNodes(), ne=mesh.NumElements();
	Basso_nMatrix nodes(3,nn);
	Basso_iMatrix conn(8,ne);
	Basso_Array<BASSO_IDTYPE> gnids(nn);
		
	mesh.SetMesh(nodes,conn);
	mesh.GlobalNodeIDs(gnids);
	
	// set up the global dof map
	TBasso_DOFMap GDOFMap(gnids,3);
	GDOFMap.FixDOFMap(comm);
	
	// Compute the global stiffness matrix
 	TBasso_FECrsMatrix Kmat(GDOFMap,conn);
	Basso_Array<BASSO_IDTYPE> sctr(3*8);
	for ( int e=0; e<conn.NumCols(); ++e )
	{
		GDOFMap.SetScatter(conn[e],8,sctr);
		Kmat.InitialFill(sctr);
	}
	Basso_Hexa8 hexa8;
	Basso_nMatrix ecoord(3,8), bmat(6,3*8), ke(3*8,3*8), cmat(6,6);
	Basso_Numeric E=10e6, nu=.3, jac;
	Basso_QuadratureRule qrule;
	hexa8.Quadrature(2,qrule);
	cmat_3d( cmat, E, nu );
	for ( int e=0; e<conn.N(); ++e )
	{
		element_coordinates( nodes, conn[e], conn.M(), ecoord );
		ke.Zero();
		GDOFMap.SetScatter(conn[e],conn.M(),sctr);
		Basso_QuadratureRule::ConstIterator qitr;
		for ( qitr=qrule.Begin(); qitr!=qrule.End(); ++qitr )
		{
			jac = hexa8.Bmatrix( qitr->Pt(), ecoord, bmat );
			cout << bmat << "\n";
			multbtcb( 'n', jac*(qitr->Wt()), bmat, cmat, 1.0, ke );
		}
		Kmat.ScatterMatrix(ke,sctr);
	}
	Kmat.FillComplete();
	
	// compute the rhs
	//TBasso_FEVector rhs(GDOFMap); 
	
	// solve the system
	
	// write the mesh
	ostringstream ss;
	ss << comm.MyPID();
	string filename = "test1_"+ss.str();
	ensight_ascii_fegeometry( filename+".geom", nodes.Data(), nodes.N(), nodes.M(), nodes.LDA(),
        conn.Data(), conn.N(), conn.M(), conn.LDA(), "hexa8" );
	ensight_case( filename, filename+".geom" ); 
	
#ifdef HAVE_MPI
	MPI_Finalize() ; 
#endif
	
	return 0;
}