
#ifndef  _TRILINOS_BASSO_PARALLEL_VTK_H_
#define  _TRILINOS_BASSO_PARALLEL_VTK_H_

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <list>

// Trilinos includes
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// Basso includes
#include "Basso_defs.h"
#include "Basso_VTK.h"

using namespace std;

namespace Basso
{


/**
 *  \brief Serial VTK XML based output for unstrictured grids
 */
class TBasso_VTKPUnstructuredGrid : public Basso_VTKUnstructuredGrid
{
	public:
	/**
	 *  constructor
	 *  \param p the process id from the mpi communicator
	 *  \param jobname The name of the resulting .vtu file  (jobname.vtu)
	 */
		TBasso_VTKPUnstructuredGrid( const Epetra_Comm &Comm, const string &jobname="results" ) 
			: comm(Comm), Basso_VTKUnstructuredGrid(jobname,Comm.MyPID()), pjob(jobname) { }
	/**
	 *  Prints the vtk file to the ostream
	 *  \param out the output stream to print to	 
	 */
		virtual void PrintPVTK( std::ostream &out=BASSO_STDOUT  ) const;
	/**
	 *  Writes the full VTK file 
	 */
		virtual void WriteFile( ) const;
	/**
	 *  Returns the filename with the extension	
	 */
		virtual string FileName( ) const { return pjob+FileExtension(); }
		
		virtual string FileName( int pid ) const 
		{ 
			string pfname = pjob;
			pfname += IDTag(pid);
			pfname += (Basso_VTKUnstructuredGrid(*this).FileExtension());
			return pfname;
		}
		
	protected:
	/**
	 *  the extension of the file
	 */
		virtual string FileExtension( ) const { return ".pvtu"; }
		
	protected:
		const Epetra_Comm &comm;
		string pjob;

};

void TBasso_VTKPUnstructuredGrid::WriteFile( ) const 
{
	Basso_VTKUnstructuredGrid(*this).WriteFile();  // write the serial .vtu file
	
	if ( comm.MyPID() != 0 ) return;
	
    ofstream outfile;
	outfile.open( FileName().c_str(), ios::out );
    if ( !outfile )
	{
		Basso_Warning("Basso_VTKPUnstructuredGrid::WriteFile","cannot open file");
        return;
	}
	
	outfile << "<?xml version=\"1.0\"?>\n";
	PrintPVTK(outfile);
	outfile.close();
}

void TBasso_VTKPUnstructuredGrid::PrintPVTK( std::ostream &out ) const
{
	
	if ( comm.MyPID() != 0 ) return;
	
	out << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >\n"
		<< "<PUnstructuredGrid GhostLevel=\"0\">\n" ;
		
	out << "<PPoints>\n";
	this->points.PrintInfo(out);
	out << "</PPoints>\n";
	
	out << "<PCells>\n";
	this->cells.PrintInfo(out);
	out << "</PCells>\n";
	
	if ( !CellDataEmpty( ) )
	{
		out << "<PCellData>\n";
		this->celldata.PrintInfo(out);
		out << "</PCellData>\n";
	}
	
	if ( !PointDataEmpty( ) )
	{
		out << "<PPointData>\n";
		this->pointdata.PrintInfo(out);
		out << "</PPointData>\n";
	}
	
	for ( int pid=0; pid<comm.NumProc(); ++pid )
		out << "<Piece Source=\""+FileName(pid)+"\"/>\n";
		
	out << "</PUnstructuredGrid>\n</VTKFile>";
}

std::ostream &operator << ( std::ostream &out, const TBasso_VTKPUnstructuredGrid &A )
{
	A.PrintPVTK( out );
	A.Print( out );
	return out;
}

	
} // end namespace
#endif

