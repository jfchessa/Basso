/*! \file TBasso_FEVector.h

	\brief Defines the class TBasso_FEVector

This class requires Trilinos Epetra	
	
\author Jack Chessa, jfchessa@utep.edu
\date Thursday, June 8, 2016

*/

#ifndef _TRILINOS_BASSO_FE_VECTOR_H_
#define _TRILINOS_BASSO_FE_VECTOR_H_

// std includes

// trilinos includes
#include "Epetra_FEVector.h"

// basso includes
#include "Basso_defs.h"

namespace Basso
{
/**
 \brief Vector for fe type processing
 This is a wrapper class for the Epetra_Vector class. 
*/

class  TBasso_FEVector 
{
public:
	TBasso_FEVector( const TBasso_DOFMap &map ) : fev_( map.OwnedGDOFMap() ) {}
	
	/**Scatters a local vector, fe, to the global locations given by sctr.  
	\param fe the local vector to be added into the global vector.
	\param sctr the global dofs to add the vector into.
	*/
	void Scatter( const Basso_nVector &fe, const Basso_Array<BASSO_IDTYPE> &sctr )
	      { fev_.SumIntoGlobalValues( sctr.Length(), sctr.Data(), fe.Data() ); }
	
	/**Gather any overlapping/shared data into the non-overlapping partitioning defined by the Map that was 
	passed to this vector at construction time. Data imported from other processors is stored on the owning 
	processor with a "sumInto" or accumulate operation. This is a collective method – every processor must 
	enter it before any will complete it.
	*/
	void FillComplete() { fev_.GlobalAssemble(); }
	
	/**Prints the vector.  This function is overloaded with the << operator.*/
	void Print( std::ostream &out=BASSO_STDOUT ) const;
	
	Epetra_FEVector &EpetraObject()             { return fev_; }
	const Epetra_FEVector &EpetraObject()     const { return fev_; }
	
	/** Write to an ASCII matlab file */
	int WriteMatlab( const string &filename ) const;
	
protected:
    Epetra_FEVector fev_; 
};
	
void TBasso_FEVector::Print( std::ostream &out ) const
{
	out << "TBasso_FEVector wrapper for Epetra_FEVector" << fev_;
}

int TBasso_FEVector::WriteMatlab( const string &filename ) const
{	
	std::ofstream outfile;
	
    outfile.open( filename.c_str(), ios::out );
    if ( !outfile )
        return 1;

	outfile << "NOt yet implemented\n";
/* 	for ( indx_type i=0; i<rowDim; ++i )
	{
		for ( indx_type j=0; j<colDim; ++j )
    		outfile << (*this)[j][i] << " ";
		outfile << "\n";
	} */
   
	outfile.close();
	
	return 0;
}

std::ostream &operator << ( std::ostream &out, const TBasso_FEVector &A )
{
	A.Print( out );
	return out;
}

} // end of namspace
#endif