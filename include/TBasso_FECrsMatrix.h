/*! \file TBasso_FECrsMatrix.h

	\brief Defines the class TBasso_FECrsMatrix

This class requires Trilinos Epetra	
	
\author Jack Chessa, jfchessa@utep.edu
\date Thursday, June 8, 2016

*/

#ifndef _TRILINOS_BASSO_FE_CRS_MATRIX_H_
#define _TRILINOS_BASSO_FE_CRS_MATRIX_H_

// std includes

// trilinos includes
#include "Epetra_FECrsMatrix.h"

// basso includes
#include "Basso_defs.h"

namespace Basso
{
/**
 \brief Compressed row sparse matrix
 This is a wrapper class for the Epetra_FECrsMatrix class.  You have to first set the sparsity of the 
 matrix using InitialFill.  Then fix the sparsity with the FixSparsity member function.  At that point 
 the ScatterMatrix function can be useed to add local finite element matricies.  When the fill is done
 then sue teh FillComplete function.
*/
class TBasso_FECrsMatrix 
{
public:
	/**Constructor
	\param map The global dof map
	\param nnzr The approximate number of non-zeros in a row.  The default is 81.
	*/
	TBasso_FECrsMatrix( const TBasso_DOFMap &map, int nnzr=81 ) 
		: crs_( Copy, map.OwnedGDOFMap(), nnzr ) { }
		
	/**Constructor that also does thie initial zero fill in
	\param map The global dof map
	\param conn the local conenctivity
	*/
	TBasso_FECrsMatrix( const TBasso_DOFMap &map,  const Basso_iMatrix &conn ) 
		: crs_( Copy, map.OwnedGDOFMap(), 81*map.NumDofPerNode() ) { InitialFill(map,conn); }
		
	/**Sets the initial sparsity entites to zero.  This needs to be done before adding and
	element stiffness matricies.  This will set the values in the FECrsMatrix to zero corresponding
	to the square submatrix defined by the rows and colums in sctr.
	\param sctr the global dofs to be filled with zeros
	*/	
	void InitialFill( const Basso_Array<BASSO_IDTYPE> &sctr );
	
	/**Sets the initial sparsity entites to zero for a set of elements given by an element connectivity matrix. 
	This needs to be done before adding and element stiffness matricies.  This will set the values in the 
	FECrsMatrix to zero corresponding global dofs of the connectivity matrix with the dofs in the global dof map gdof
	\param gdof the global dof map to define the global rows and columns to be filled
	\param conn the connectivity matrix
	*/	
	void InitialFill( const TBasso_DOFMap &gdof, const Basso_iMatrix &conn );
	
	/**Scatters a square matrix, ke, to the rows and columns given by sctr.  These dofs must be initialized
	first with an InitialFill call.
	\param ke the square matrix to be added into the sparse matrix
	\param sctr the global dofs (row and coulumn) to add the matrix into.
	*/
	void Scatter( const Basso_nMatrix &ke, const Basso_Array<BASSO_IDTYPE> &sctr );
	
	/**To be called when teh matrix fill in is complete.  This will move all the off processor
	values onto the processor.  This needs to be done before you atempt to solve the system.
	This is simple a call to the Epetra_FECrsMatrix GlobalAssemble member function, so you can see
	that documention for more details.
	*/
	void FillComplete() { crs_.GlobalAssemble(); }
	
	/**Returns a non-constant reference to the underlying Epetra_FECrsMatrix class.*/
	Epetra_FECrsMatrix &EpetraObject()                 { return crs_; }
	const Epetra_FECrsMatrix &EpetraObject()     const { return crs_; }
	
	/**Prints the matrix.  This function os overloaded with the << operator.*/
	virtual void Print( std::ostream &out=BASSO_STDOUT ) const;
	
	/** Write to an ASCII matlab file */
	int WriteMatlab( const string &filename ) const;

protected:
	Epetra_FECrsMatrix crs_;
	
};


int TBasso_FECrsMatrix::WriteMatlab( const string &filename ) const
{	
	std::ofstream outfile;
	
    outfile.open( filename.c_str(), ios::out );
    if ( !outfile )
        return 1;

	outfile << "NOt yet implemented\n";
	
/* 	const Epetra_Comm &lcomm = crs_.Map().Comm();
	int mnnzr = crs_.MaxNumEntries();
	double values[mnnzr];
	int indices[mnnzr];
	for ( int p=0; p<lcomm.NumProc(); ++p )
	{
		if ( lcomm.MyPID() == p )
			for ( int i=0; i<crs_.NumMyRows(); ++i )
			{
				int I = crs_.GRID(i), nnz;
				crs_.ExtractMyRowCopy( i, mnnzr, nnz, values, indices);
				for ( int j=0; j<nnz; ++ j )
					outfile << setw(10) << I << setw(10) << indices[j] << setw(15) << values[j] << "\n";
			}
		lcomm.Barrier();		
	}
    lcomm.Barrier(); */
	
	outfile.close();
	
	return 0;
}

void TBasso_FECrsMatrix::Print( std::ostream &out ) const
{
	out << "TBasso_FECrsMatrix wrapper for Epetra_FECrsMatrix" << crs_;
}

void TBasso_FECrsMatrix::InitialFill( const Basso_Array<BASSO_IDTYPE> &sctr ) 
{
	Basso_nVector zeros(sctr.Length()*sctr.Length());
	crs_.InsertGlobalValues(sctr.Length(),sctr.Data(),zeros.Data());
}

void TBasso_FECrsMatrix::InitialFill( const TBasso_DOFMap &gdof, const Basso_iMatrix &conn )
{
	int nn = conn.NumRows(), ndofpn=gdof.NumDofPerNode();
	int nz = (nn*ndofpn);
	Basso_nVector zeros(nz*nz);
	Basso_Array<BASSO_IDTYPE> sctr(nz);
	for ( int e=0; e<conn.NumCols(); ++e )
	{
		gdof.SetScatter( conn[e], conn.NumRows(), sctr );
		//cout << sctr << " TBasso_FECrsMatrix.h:111\n";
		crs_.InsertGlobalValues(sctr.Length(),sctr.Data(),zeros.Data());
	}
}

void TBasso_FECrsMatrix::Scatter( const Basso_nMatrix &ke, const Basso_Array<BASSO_IDTYPE> &sctr ) 
{ 
	crs_.SumIntoGlobalValues(sctr.Length(),sctr.Data(),ke.Data());
}

std::ostream &operator << ( std::ostream &out, const TBasso_FECrsMatrix &A )
{
	A.Print( out );
	return out;
}

} // end of namspace
#endif