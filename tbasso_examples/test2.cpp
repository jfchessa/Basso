// std includes
#include <iostream>
#include <set>

// mpi includes
#include "mpi.h"

// Trilinos includes
#include "Epetra_Version.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// basso includes
#include "Basso_defs.h"
#include "Basso_Numeric.h"
#include "Basso_iMatrix.h"
#include "Basso_nMatrix.h"
#include "Basso_nfeaops.h"

// tbasso includes
#include "TBasso_DOFMap.h"

using namespace std;
using namespace Basso;


/*
  This example requires 2 processors
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
	
	if ( comm.NumProc() != 2 )
	{
		cout << "Requires 2 cpus\n";
		return 1;
	}

	/*
	
	This is a simple example of a domain decomposition for standard finite elements using
	Trilinos.  This is outlined in the ppt file "Basso and Trilnions Data Structures for FEA Data"
	
	This file should be run with exactly two threads to work correctly.
	
	25 ------- 26  ------- 27  ------- 28  ------- 29
	|          |           |           |           |
	|     5    |     6     |     8     |     9     |
	|          |           |           |           |
	20 ------- 21  ------- 22  ------- 23  ------- 24
	|          |           |           |           |
	|     0    |     1     |     2     |    3      |
	|          |           |           |           |
	10 ------- 11  ------- 12  ------- 13  ------- 14
	
	PID: 0 ownes elements 0, 1, 5 and 6 and nodes 10, 11, 12, 21, 22, 27, 25, and 26
	PID: 1 ownes elements 2, 3, 8 and 9 and nodes 13, 14, 23, 24, 27, 28 and 29
	
	All nodes have two dofs 
	
	Symmetry bc are envoked on the left and bottome edges and a 100 unit load is uniformly 
	distributed on the right edge (in the normal direction).
	
	----------------------------------------------
	Jack Chessa
	Associate Professor of Mechancial Engineering
	The University of Texs at El Paso
	jchessa@utep.edu
	
	*/
	
	Basso_iMatrix myConn(4,4);
	switch ( comm.MyPID() ) 
	{
	case 0:
		myConn = "10 11 20 21;11 12 21 22;21 22 26 27;20 21 25 26";
		break;
	
	case 1:
		myConn = "12 13 22 23;13 14 23 24;23 24 28 29;22 23 27 28";
		break;
	}
	
	cout << myConn << "\n";
	
	Basso_Array<int> mygnids;
	get_gnids( myConn, mygnids );
	
	cout << mygnids << "\n";
	
	TBasso_DOFMap gdofMap(mygnids,2);
	gdofMap.FixDOFMap(comm);
	
	cout << gdofMap << "\n";
	
#ifdef HAVE_MPI
	MPI_Finalize() ; 
#endif
	
	return 0;
}