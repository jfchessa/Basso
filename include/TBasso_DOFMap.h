/*! \file TBasso_DOFMap.h

	\brief Defines the class TBasso_DOFMap

This class requires Trilinos Epetra	
	
\author Jack Chessa, jfchessa@utep.edu
\date Thursday, June 8, 2016

*/

#ifndef _TRILINOS_BASSO_DOF_MAP_H_
#define _TRILINOS_BASSO_DOF_MAP_H_

// std includes
#include <map>
#include <iomanip>

// trilinos includes
#include "Epetra_Map.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

// basso includes
#include "Basso_defs.h"
#include "Basso_Array.h"
#include "Basso_DOFMap.h"

namespace Basso
{

	
/**
 \brief DOF map with a fixed number of active dofs per node.
	
	This class requires Trilinos Epetra	
*/
class TBasso_DOFMap : public Basso_DOFMap
{
public:
	/**
	Constructor
	\param gnids - an array of the global node ids that are supported by the processor.  This array needs to 
	map onto the rows of the local node coordinate matrix
	\param nldof the number of dofs active on each node.  Default value is 6.
	*/
	TBasso_DOFMap( const Basso_Array<BASSO_IDTYPE> &gnids, int nldof=6 ) 
				: gdof_map_(nldof,gnids.Length())
		{ fixed_=false; SetGNIDs(gnids); owned_gdof_map_=NULL; }

	
	/**
	Constructor that fully sets up the GDOFMap.  If calls FixDOFMap.
	\param Comm -  The Epetra mpi commuication operator
	\param gnids - an array of the global node ids that are supported by the processor.  This array needs to 
	map onto the rows of the local node coordinate matrix
	\param nldof the number of dofs active on each node.  Default value is 6.
	*/	
	TBasso_DOFMap( const Epetra_Comm &Comm, const Basso_Array<BASSO_IDTYPE> &gnids, int nldof=6 ) 
				: gdof_map_(nldof,gnids.Length())
		{ fixed_=false; SetGNIDs(gnids); owned_gdof_map_=NULL; FixDOFMap(Comm); }
		
	/** Deconstructor */
	~TBasso_DOFMap() { if ( owned_gdof_map_!=NULL ) delete owned_gdof_map_; }
		
	/** Sets the global node ids that are supported by the procesor.
	\param gnids - an array of the global node ids that are supported by the processor.  This array needs to 
	map onto the rows of the local node coordinate matrix
	*/	
	void SetGNIDs( const Basso_Array<BASSO_IDTYPE> &gnids ) { gnids_ = gnids; }
	
	/**
	Initialized the global dof map.  This takes the global node ids supported on each processor
	and determines a unique global id 
	\param Comm - The mpi commuication operator
	*/
	void FixDOFMap( const Epetra_Comm &Comm );
	
	/**Returns the number of active dofs for each node.*/
	int NumDofPerNode() const { return gdof_map_.NumRows(); }
	
	/**Returns the number of nodes supported on the local processor.*/
	int NumMyNodes()    const { return gnids_.Length(); }
	
	/**Initializes a scatter vector for the gdof map.  There is no bounds checking to expidite 
	the member function.  The form of the scatter vector is as follows:
		sctr = [ node0_dof0, node0_dof1, node0_dof2 ..., node1_dof0, node1_dof2, ... ] 
	\param econn a vector of the local node ids for which to compute the scatter vector. Do not use global node ids here.
	\param sctr on return has the scatter vector to the global dofs.  It should be of length econn.Length()*NumDofPerNode(). 
	*/
	void SetScatter( const Basso_Array<BASSO_IDTYPE> &econn, Basso_Array<BASSO_IDTYPE> &sctr ) const;
	
	/**Initializes a scatter vector for the gdof map.  There is no bounds checking to expidite 
	the member function.  The form of the scatter vector is as follows:
		sctr = [ node0_dof0, node0_dof1, node0_dof2 ..., node1_dof0, node1_dof2, ... ] 
	\param econn a pointer to an array of the local node ids for which to compute the scatter vector. Do not use global node ids here.
	\param nne the length of the array in econn
	\param sctr on return has the scatter vector to the global dofs.  It should be of length econn.Length()*NumDofPerNode(). 
	*/
	void SetScatter( const BASSO_IDTYPE *econn, int nne, Basso_Array<BASSO_IDTYPE> &sctr ) const;
	
	
	/**Initializes a scatter vector for the gdof map using only particular dofs.  There is no bounds checking to expidite 
	the member function.  The form of the scatter vector is as follows:
		sctr = [ node0_dof0, node0_dof1, node0_dof2 ..., node1_dof0, node1_dof2, ... ] 
	\param econn a pointer to an array of the local node ids for which to compute the scatter vector. Do not use global node ids here.
	\param nne the length of the array in econn
	\param sctr on return has the scatter vector to the global dofs.  It should be of length econn.Length()*NumDofPerNode(). 
	\param ldofs an array that has the local dofs to be used in forming the scatter vector.  Again there is no error checking here.
	*/
	void SetScatter( const BASSO_IDTYPE *econn, int nne, Basso_Array<BASSO_IDTYPE> &sctr, const Basso_Array<int>  &ldofs ) const;
	
	/**Returns the global dof given a lodal node id and a local dof
	 \param lnid the local node id
	\param ldof the local dof  for the lnid node to return the global dof.
	\return the globaldof corresponding to the ldof dof of node lnid.	 
	 */
	BASSO_IDTYPE GDOF( BASSO_IDTYPE lnid, int ldof ) const 
                { if ( !fixed_) return 0; return gdof_map_[lnid][ldof]; }
		
	/** Returns a pointer to a one to one gdof map that can be used to construct
	Epetra fintie element objects (Epetra_FECsrMatrix, etc. )
	\return a constant pointer to the gdof map
	*/
	const Epetra_Map &OwnedGDOFMap() const { return *owned_gdof_map_; }
	
	/**Sets the dofs in the map 
	\param ldofs on return has the local dofs active in the map
	*/
	void LDOFs( Basso_Array<int> & ldofs ) const;
	
	/**Returns the number of total global dofs defined by the map*/
	BASSO_IDTYPE MaxAllGDOF() const 
	{ 
		if ( owned_gdof_map_ == NULL )
			return 0;
		return owned_gdof_map_->MaxAllGID64(); 
	}
	/**Returns the number of total global dofs defined by the map*/
	BASSO_IDTYPE NumTotalGDOFs() const 
	{ 
		if ( owned_gdof_map_ == NULL )
			return 0;
		return owned_gdof_map_->NumGlobalElements64(); 
	}
	
	/**Returns the number of global dofs on the processor*/
	BASSO_IDTYPE NumMyGDOFs() const 
	{ 
		if ( owned_gdof_map_ == NULL )
			return 0;
		return owned_gdof_map_->NumMyElements(); 
	}
	
	/**Prints the map.  Also is used for  << overlaod.*/
	virtual void Print( std::ostream &out=BASSO_STDOUT ) const;
	
protected:
	
	Basso_Array2D<BASSO_IDTYPE> gdof_map_;
	Basso_Array<BASSO_IDTYPE> gnids_;
	bool fixed_;
	Epetra_Map *owned_gdof_map_;
};

void TBasso_DOFMap::LDOFs( Basso_Array<int> & ldofs ) const
{
	ldofs.Resize( NumDofPerNode() );
	for ( int s=0; s<ldofs.Length(); ++ s )
		ldofs[s] = s;
}

void TBasso_DOFMap::Print( std::ostream &out ) const 
{
 	if ( !fixed_ ) 
	{
		out << "TBasso_DOFMap, not yet set\n";
		return;
	}
	
	int nn=NumMyNodes();
	if ( nn>1000 ) nn=1000;
	
	out << "TBasso_DOFMap, PID=" << owned_gdof_map_->Comm().MyPID() <<
	       "\n      lnid      gnid      ldof      gdof\n"; 
	for ( int i=0; i<nn; ++i )
		for ( int s=0; s<NumDofPerNode(); ++s )
			out << setw(10) << i << setw(10) << gnids_[i] << setw(10) << s << setw(10) << gdof_map_[i][s] << "\n"; 
}

void TBasso_DOFMap::SetScatter( const Basso_Array<BASSO_IDTYPE> &econn, Basso_Array<BASSO_IDTYPE> &sctr ) const
{
	SetScatter(econn.Data(), econn.Length(), sctr);		
}


void TBasso_DOFMap::SetScatter( const BASSO_IDTYPE *econn, int nne, Basso_Array<BASSO_IDTYPE> &sctr ) const
{
	const BASSO_IDTYPE *nptr = econn;
	BASSO_IDTYPE *cptr = sctr.Data(), I;
	int s;
	for ( I=0; I<nne; ++I, ++nptr )
		for ( s=0; s<NumDofPerNode(); ++s, ++cptr ) 
			*cptr = gdof_map_[*nptr][s];		
}

void TBasso_DOFMap::SetScatter( const BASSO_IDTYPE *econn, int nne, Basso_Array<BASSO_IDTYPE> &sctr, const Basso_Array<int>  &ldofs ) const
{
	const BASSO_IDTYPE *nptr = econn;
	BASSO_IDTYPE *cptr = sctr.Data(), I;
	int s;
	for ( I=0; I<nne; ++I, ++nptr )
	{
		const int *sptr=ldofs.Data();
		for ( s=0; s<ldofs.Length(); ++s, ++cptr, ++sptr ) 
			*cptr = gdof_map_[*nptr][*sptr];
	}		
}

void TBasso_DOFMap::FixDOFMap( const Epetra_Comm &Comm )
{
   /* ================= Compute the owned nodes ids on that processor ================= 
   /*	
	 1) setup a non-subjective Epetra_Map containing the nodes ids that
		  are in the element connectivity of this processor
	*/	 
	
	Epetra_Map imap0(-1,gnids_.Length(),gnids_.Data(),0,Comm);      // temp
	Epetra_IntVector supportNID(imap0);                             // temp
	supportNID.PutValue(Comm.MyPID());  
    
	/*
    2) define an Epetra_IntVector that will contain the owning processor id (PID)
	    for each node.  This map can be arbitrailly set up by Trilinos, but will be
	    subjective	
    */
    int numGlobalNodes=imap0.MaxAllGID()+1;
    Epetra_Map imap1(numGlobalNodes,0,Comm);   // temp
	Epetra_IntVector ownedPID(imap1);          // temp
    Epetra_Export exptr(imap0,imap1);          // temp
    ownedPID.Export(supportNID,exptr,Insert);	
	
	/*
	3) Import the data from the shared map to the subjective IntVector
	*/
    Epetra_Import imptr(imap0,imap1);          // temp
    supportNID.Import(ownedPID,imptr,Insert);
    
	int numMyNodes=0;
    for ( int i=0; i<supportNID.Map().NumMyElements(); ++i )
        if ( supportNID[i] == Comm.MyPID() )
            ++numMyNodes;    
		
	Epetra_IntSerialDenseVector myNIDs(numMyNodes);           // myNIDs; This vector has the gnids
    int j=0; // ** set myNIDs **                              // of the nodes "owned "by this processor
    for ( int i=0; i<supportNID.Map().NumMyElements(); ++i )
        if ( supportNID[i] == Comm.MyPID() )
            myNIDs[j++] = imap0.GID(i);
		
	/* ================= Compute the owned GDOF ids ================= */
	int numSharedNodes = gnids_.Length();	
 	Epetra_IntSerialDenseVector nodeNDOF(numSharedNodes);
	int numMyGDOFs = numMyNodes*NumDofPerNode();
	Epetra_IntSerialDenseVector myGDOFs(numMyGDOFs);  // GDOF of ownded nodes
	int gdof = 0;
	for ( int pid=0; pid<Comm.NumProc(); ++pid )
	{
		if ( Comm.MyPID()==pid )
		{
			int ii=0;
			for ( int n=0; n<numMyNodes; ++n )
				for ( int s=0; s<NumDofPerNode(); ++s, ++ii, ++gdof )
					myGDOFs[ii] = gdof;
		}
		Comm.Broadcast( &gdof, 1, pid );
	}
	
	owned_gdof_map_ = new Epetra_Map(-1,numMyGDOFs,myGDOFs.Values(),0,Comm);
	
	/* ============ Compute the GDOF for the shared nids ============= */ 
	Epetra_Map ownedNodeMap(-1,numMyNodes,myNIDs.Values(),0,Comm);  
	Epetra_Map sharedNodeMap(-1,numSharedNodes,gnids_.Data(),0,Comm);
	Epetra_IntVector startMyGDOF(ownedNodeMap), startSharedGDOF(sharedNodeMap);   
	for  ( int i=0; i<numMyNodes; ++i )
		startMyGDOF[i]=myGDOFs[NumDofPerNode()*i];
	
	Epetra_Import imprt(sharedNodeMap,ownedNodeMap);  
	
	startSharedGDOF.Import(startMyGDOF,imprt,Insert);
	
	for  ( int i=0; i<numSharedNodes; ++i )
		for ( int s=0; s<NumDofPerNode(); ++s )
			gdof_map_[i][s] = startSharedGDOF[i]+s;
		
	// done! so set as fixed
	fixed_ = true;

}

std::ostream &operator << ( std::ostream &out, const TBasso_DOFMap &A )
{
	A.Print( out );
	return out;
}

} // end of namspace
#endif