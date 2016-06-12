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

// tbasso includes
#include "TBasso_DOFMap.h"
#include "TBasso_FECrsMatrix.h"
//#include "TBasso_FEVector.h"

using namespace std;
using namespace Basso;

/*
class TBasso_FELinearOperator
{
public:
	TBasso_FELinearOperator( const TBasso_DOFMap &gdof ) { gdof_ptr_ = &gdof; }
	
	void StiffnessMatrix( TBasso_FECrsMatrix &Kmat,  const Basso_nMatrix &node, 
	       const Basso_iMatrix &conn, const Basso_ParentElement *element, 
		   const Basso_nMatrix &cmat, Basso_Numeric a=1.0 ) const;

protected:
	const TBasso_DOFMap *gdof_ptr_;
	
};

void TBasso_FELinearOperator::StiffnessMatrix( TBasso_FECrsMatrix &Kmat,  const Basso_nMatrix &node, 
		const Basso_iMatrix &conn, const Basso_ParentElement *element, const Basso_nMatrix &cmat,
		Basso_Numeric a ) const
{
}			
*/


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
	//TBasso_FELinearOperator feop(GDOFMap);
	Basso_Hexa8 hexa8;
	Basso_nMatrix cmat(6,6);
	//feop.StiffnessMatrix( Kmat, nodes, conn, &hexa8, cmat );
	//Kmat.FillComplete();
	
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