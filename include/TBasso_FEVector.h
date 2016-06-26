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
#include "Epetra_Vector.h"
#include "Epetra_Import.h"

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
	
	TBasso_FEVector( const Epetra_BlockMap &map ) : fev_( map ) {}
	
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

	outfile << "NOt yet implemented\n";  // Use EpetraExt
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

/**
 *  \brief Class to contain local data
 * This class is a wrapper for an Epetra_Vector based on a shared GDOF map.  This vector
 * can "pull" the data from a FEVector based on an "owned" dof map so that the local shared
 * values are accessable. This class constructs a static Epetra_Import at insatiation.
 */
 class TBasso_LocalVector 
 {
	 public:
		 /**
		 * Constructor
		 \param map the gdof map used to get the shared and owned gdof maps
		 */
		TBasso_LocalVector(  const TBasso_DOFMap &map  ) 
			: fl(map.SharedGDOFMap()), imprtr(map.SharedGDOFMap(),map.OwnedGDOFMap()) {  }
		
		 /**
		 * Constructor
		 \param sharedmap the shared gdof map 
		 \param ownedmap  the owned gdof map
		 */	
		TBasso_LocalVector(  const Epetra_BlockMap &sharedmap, const Epetra_BlockMap &ownedmap ) 
			: fl(sharedmap), imprtr(sharedmap,ownedmap) { }
		
		/**\return the number of local values in the vector*/
		//BASSO_IDTYPE NumLocalValues( ) const { return fl.MyLength(); }
		BASSO_IDTYPE Length( ) const         { return fl.MyLength(); }
		
		/** Returns a constant pointer to the local vector data
		\return constant pointer to the local vector data.
		*/
		const Basso_Numeric *Data( ) const 
		{ 
			Basso_Numeric *dptr; 
			fl.ExtractView( &dptr );
			return dptr;
		}
		
		/** Gathers the local values according to the gather ldof vector
		\param gthr the array of local indices (ldofs) to pull gather the values
		\param vals the array where the values are places (no size chacking)
		*/
		void Gather( const Basso_Array<BASSO_IDTYPE> &gthr, Basso_nVector &vals ) const;
		/** pulls the data from a vector based on the owned gdof map to the local shared gdof
		values
		\param v the vecto whcih to pull down values from.  Must use the same owned gdof map as the
		    map used in teh constructor.
		*/
		void PullValues(const TBasso_FEVector &v)
		{ fl.Import(v.EpetraObject(),imprtr,Insert); }
		
		/**Prints the vector.  This function is overloaded with the << operator.*/
		void Print( std::ostream &out=BASSO_STDOUT ) const;
		
	protected:
		Epetra_Vector fl;
		Epetra_Import imprtr;
 };
 
void TBasso_LocalVector::Gather( const Basso_Array<BASSO_IDTYPE> &gthr, Basso_nVector &vals ) const	
{
	const BASSO_IDTYPE *gptr = gthr.Data();
	for ( int i=0; i<gthr.Length(); ++i, ++gptr )
		vals[i] = fl[*gptr];
}

void TBasso_LocalVector::Print( std::ostream &out ) const
{
	out << "TBasso_LocalkVector wrapper for Epetra_Vector" << fl;
}

std::ostream &operator << ( std::ostream &out, const TBasso_LocalVector &A )
{
	A.Print( out );
	return out;
}

} // end of namspace
#endif