/*! \file ensight.h

Functions to write Ensight data files

\author Jack Chessa, jfchessa@utep.edu
\date Sunday, Oct 30, 2011

*/

#ifndef _BASSO_WRITE_ENSIGHT_H_
#define _BASSO_WRITE_ENSIGHT_H_

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <list>

#include "Basso_defs.h"

namespace Basso
{

#define FLOAT_PREC_ENSIGHT 5
#define FLOAT_FWID_ENSIGHT 12
#define INT_FWID_ENSIGHT 10
#define STR_FWID_ENSIGHT 80

using namespace std;

/**
performs a basic ofstream open with some checking
*/
bool ensight_open_file( const string &filename, ofstream &outfile )
{
    outfile.open( filename.c_str(), ios::out );
    if ( outfile )
        return false;
    return true;
}

/** opens an existing geometry file and adds an element connectivity
Only works for ASCII files
*/
int ensight_add_connectivity( const string &filename, const BASSO_IDTYPE *conn, BASSO_IDTYPE ne, int nne, int ldc, 
        const string &etype, int conn_offset=1 )
{
    ofstream outfile;
    outfile.open( filename.c_str(), ios::app );
    if ( !outfile )
        return 1;
  
    // write the element connectivity
    outfile << etype << "\n" << setw(10) << ne << "\n";
    for ( BASSO_IDTYPE e=0; e<ne; ++e )
    {
        for ( int i=0; i<nne; ++i )
            outfile << setw(10) << conn[e*ldc+i] + conn_offset;
        outfile << "\n";
    }
    outfile.close();
    
    return 0;
    
}

/**
Writes a geometry file for a rectangular structured grid with uniform spacing
\param filename -  the file that the data is to be written to (include the extension)
\param nx - the number of points in the x-axis direction
\param ny - the number of points in the y-axis direction
\param nz - the number of points in the z-axis direction
\param x0 - The x-coordinate of the origin (starting point of the grid)
\param y0 - The y-coordinate of the origin (starting point of the grid)
\param z0 - The z-coordinate of the origin (starting point of the grid)
\param deltx - The increment in x in the grid
\param delty - The increment in y in the grid
\param deltz - The increment in z in the grid
*/
template < class NuMeRiC >
int ensight_ascii_uniform( const string &filename, BASSO_IDTYPE nx, BASSO_IDTYPE ny, BASSO_IDTYPE nz,
	NuMeRiC x0, NuMeRiC y0, NuMeRiC z0, NuMeRiC deltx, NuMeRiC delty, NuMeRiC deltz, int partid=1 )
{
    ofstream outfile;
    if ( ensight_open_file(filename,outfile) )
        return 1;
  
    outfile << "Ensight Gold Geometry File"
            << "\nStructured Grid Geometry"
            << "\nnode id off"
            << "\nelement id off"
            << "\npart\n" << setw(10) << partid << "\na uniform rectilinear finite difference grid";

    // write the grid parameters     
    outfile << "\nblock uniform\n" 
			<< setw(10) << nx << setw(10) << ny << setw(10) << nz << "\n";
			
    outfile.precision(5);
    outfile << setw(12) << std::scientific << x0 << endl
			<< setw(12) << std::scientific << y0 << endl
			<< setw(12) << std::scientific << z0 << endl
			<< setw(12) << std::scientific << deltx << endl
			<< setw(12) << std::scientific << delty << endl
			<< setw(12) << std::scientific << deltz << endl;
	return 0;
}

/**
Writes finite element type geometry
\param filename -  the file that the data is to be written to (include the extension)
\param node -  the node coordinates [x1 y1 z1.. x2 y2 z2 ..]
\param nn - the number of nodes to be written
\param sdim - the number of dimensions in \param node (i.e.  sdim=2 only x and y are given)
\param ldn - the leading dimension of \param node
\param conn -  the element connectivity matrix [ e1n1 e1n2 ... e2n1 e2n2 .... ]
\param ne - the number of elements
\param nne - the number of nodes in an elemetn connectivity
\param ldc - the leading dimension of \param conn
\param etype - the appropriate element type for conn.  This is a required by ensight - 
	point, bar2, bar3, tria3, tria6, quad4, quad8, tetra4, tetra10, pyramid5, pyramid13,
	hexa8 hexa20, penta6, penta15
*/
template < class NuMeRiC >
int ensight_ascii_fegeometry( const string &filename, const NuMeRiC *node, BASSO_IDTYPE nn, int sdim, int ldn,
        const BASSO_IDTYPE *conn, BASSO_IDTYPE ne, int nne, int ldc, const string &etype, int partid=1, int conn_offset=1 )
{
    ofstream outfile;
    if ( ensight_open_file(filename,outfile) )
        return 1;
  
    outfile << "Ensight Gold Geometry File"
            << "\nFinite Element Geometry"
            << "\nnode id off"
            << "\nelement id off"
            << "\npart\n" << setw(10) << partid << "\na finite element part";
            
    // write the node coordiantes        
    outfile << "\ncoordinates\n" << setw(10) << nn << "\n";
    outfile.precision(5);
    for ( int s=0; s<sdim; ++s )
        for ( BASSO_IDTYPE i=0; i<nn; ++i )
            outfile << setw(12) << std::scientific << node[ ldn*i + s ] << endl;
    for ( BASSO_IDTYPE i=0; i<(3-sdim)*nn; ++i )
        outfile << setw(12) << std::scientific << 0.0 << endl;
    
    // write the element connectivity
    outfile.close();
    ensight_add_connectivity( filename, conn, ne, nne, ldc, etype, conn_offset );
/*    
    outfile << etype << "\n" << setw(10) << ne << "\n";
    for ( BASSO_IDTYPE e=0; e<ne; ++e )
    {
        for ( int i=0; i<nne; ++i )
            outfile << setw(10) << conn[e*ldc+i] + 1;
        outfile << "\n";
    }
    outfile.close();
*/    
    return 0;
    
}

/**
 Writes scalar grid data to an ensight gold fromat.  If etype is given then
 it assumes per element data else per node data.

 etype Line2, Line3, Tria3, Tria6, Quad4, Quad8, 
       Tetra4, Tetra10, Hexa8, Hexa20

\param filename - the name of the file to be written
\param data - the pointer to the data array
\param nd - the number of data points (nodes,elements)
\param etype - If per element data is given this is the element type (the default is none)
\param partid - an optional part id (default is 1)

 etype Line2, Line3, Tria3, Tria6, Quad4, Quad8, 
       Tetra4, Tetra10, Hexa8, Hexa20

*/
 
template < class NuMeRiC >
int ensight_ascii_field( const string &filename, const NuMeRiC *data, BASSO_IDTYPE nd, 
            const string &etype="none", int partid=1 )
{
    
    ofstream outfile;
    if ( ensight_open_file(filename,outfile) )
        return 1;
        
     outfile << "Scalar Data"
            << "\npart\n" << setw(10) << partid << "\n";
    
    if ( etype=="none" )  // nodal data     
        outfile << "coordinates\n";

    else  // element data
        outfile << etype << "\n";
         
    const NuMeRiC *dptr=data; 
	outfile.precision(5);
    for ( BASSO_IDTYPE i=0; i<nd; ++i, ++dptr )
        outfile << setw(12) << std::scientific << *dptr << "\n";
    
    return 0;
}

/**
 Writes general data to an ensight gold fromat. This is typically used by other
function calls. 
\param filename - the name of the file to be written
\param description - a string the is put in the file as a description
\param data - the pointer to the data array
\param nd - the number of data points (nodes,elements)
\param sdim - the spacial dimension of the first rank of the field
\param vdim - the spacial dimension of the second rank of the field
\param ldd - the leading dimension of the data array
\param etype - If per element data is given this is the element type (the default is none)
\param partid - an optional part id (default is 1)
 If etype is given then it assumes per element data else per node data.

 etype Line2, Line3, Tria3, Tria6, Quad4, Quad8, 
       Tetra4, Tetra10, Hexa8, Hexa20

*/
template < class NuMeRiC >
int ensight_ascii_genfield( const string &filename, const string &description,
            const NuMeRiC *data, BASSO_IDTYPE nd, int sdim, int vdim, int ldd, 
            const string &etype="none", int partid=1 )
{    
    ofstream outfile;
    if ( ensight_open_file(filename,outfile) )
        return 1;
        
     outfile << description
            << "\npart\n" << setw(10) << partid << "\n";
    
    if ( etype=="none" )  // nodal data     
        outfile << "coordinates\n";

    else  // element data
        outfile << etype << "\n";
         
    //const NuMeRiC *dptr=data; 
	outfile.precision(5);
    for ( int s=0; s<sdim; ++s )
        for ( BASSO_IDTYPE i=0; i<nd; ++i )
            outfile << setw(12) << std::scientific << data[ i*ldd + s ] << "\n";
    
    for ( BASSO_IDTYPE s=0; s<(vdim-sdim)*nd; ++s )
        outfile << setw(12) << 0.0 << "\n";
        
    return 0;
}

/**
Returns a string equal to \param filename with a timestamp added onto the end.
\param filename - the base filename on which the timestamp will be added
\param nstamp - the width of the timestamp field
\param step - the time step
Example:
add_timestamp( "datafile", 4, 32 )
will return the following string "datafile0032"
*/
string add_timestamp( const string &filename, int nstamp, int step )
{
    stringstream ss (stringstream::in | stringstream::out);
    ss << right;
    ss.width(nstamp);
    ss << setfill('0');
    ss << step;
    return filename + ss.str();
}

/**
 Writes vector grid data to an ensight gold fromat.  If etype is given then
 it assumes per element data else per node data.


\param filename - the name of the file to be written
\param data - the pointer to the data array
\param nd - the number of data points (nodes,elements)
\param sdim - the spacial dimension of the vector field
\param ldd - the leading dimension of the data array
\param etype - If per element data is given this is the element type (the default is none)
\param partid - an optional part id (default is 1)

 etype Line2, Line3, Tria3, Tria6, Quad4, Quad8, 
       Tetra4, Tetra10, Hexa8, Hexa20

 
*/
template < class NuMeRiC >
int ensight_ascii_field( const string &filename, const NuMeRiC *data, BASSO_IDTYPE nd, int sdim, int ldd, 
            const string &etype="none", int partid=1 )
{    
    return ensight_ascii_genfield( filename, "Vector Data", data, nd, sdim, 3, ldd, etype, partid );
}

template < class NuMeRiC >
int ensight_ascii_tnsfield( const string &filename, const NuMeRiC *data, BASSO_IDTYPE nd, int sdim, int ldd, 
            const string &etype="none", int partid=1 )
{
    return ensight_ascii_genfield( filename, "Tensor Data", data, nd, sdim, 6, ldd, etype, partid );
}

void write_case_variables( ofstream &outfile, const string &jobname, int ntpad,
    const list<string> &scalarvar, 
    const list<string> &vectorvar, 
    const list<string> &tensorvar, 
    const list<string> &scalarcell, 
    const list<string> &vectorcell, 
    const list<string> &tensorcell )
{
    if ( scalarvar.size()==0 && 
         vectorvar.size()==0 &&
         tensorvar.size()==0 &&
         scalarcell.size()==0 && 
         vectorcell.size()==0 &&
         tensorcell.size()==0 )  // all empty
    return;
    
    outfile << "VARIABLE\n";
    
    string onback="";
    if ( ntpad>0 )
        onback=string(ntpad,'*');
    
    list<string>::const_iterator vitr;
    for ( vitr=scalarvar.begin(); vitr!=scalarvar.end(); ++vitr )
        outfile << "scalar per node:  " << jobname+"."+(*vitr)+onback << "\n";
    
    for ( vitr=vectorvar.begin(); vitr!=vectorvar.end(); ++vitr )
        outfile << "vector per node:  " << jobname+"."+(*vitr)+onback << "\n"; 
    
    for ( vitr=tensorvar.begin(); vitr!=tensorvar.end(); ++vitr )
        outfile << "tensor symm per node:  " << jobname+"."+(*vitr)+onback << "\n";
         
    for ( vitr=scalarcell.begin(); vitr!=scalarcell.end(); ++vitr )
        outfile << "constant per element:  " << jobname+"."+(*vitr)+onback << "\n";
    
    for ( vitr=vectorcell.begin(); vitr!=vectorcell.end(); ++vitr )
        outfile << "vector per element:  " << jobname+"."+(*vitr)+onback << "\n"; 
    
    for ( vitr=tensorcell.begin(); vitr!=tensorcell.end(); ++vitr )
        outfile << "tensor symm per element:  " << jobname+"."+(*vitr)+onback << "\n";     
}

/**
Writes an Ensight case file for non-time dependent data. 
\param jobname - the string with the jobname.  The case fill will be jobname.case, so one should not put the file extension on the string
\param geomfile - the string withthe geometry file (needs the file extension)
\param scalarvar - a list<string> that contains the names of the per node scalar variable data
\param vectorvar - a list<string> that contains the names of the per node vector variable data
\param tensorvar - a list<string> that contains the names of the per node tensor variable data
\param scalarcell - a list<string> that contains the names of the per element scalar variable data
\param vectorcell - a list<string> that contains the names of the per element vector variable data
\param tensorcell - a list<string> that contains the names of the per element tensor variable data
*/
int ensight_case( const string &jobname, const string &geomfile, 
    const list<string> &scalarvar=list<string>(), 
    const list<string> &vectorvar=list<string>(), 
    const list<string> &tensorvar=list<string>(), 
    const list<string> &scalarcell=list<string>(), 
    const list<string> &vectorcell=list<string>(), 
    const list<string> &tensorcell=list<string>() )
{
    ofstream outfile;
    if ( ensight_open_file(jobname+".case",outfile) )
        return 1;
        
    outfile << "FORMAT  \ntype:     ensight gold\n";
    outfile << "GEOMETRY\nmodel:    " << geomfile << "\n";
    write_case_variables(outfile,jobname,0,
            scalarvar,vectorvar,tensorvar,scalarcell,vectorcell,tensorcell);
      
    return 0;              
}

/**
Writes an Ensight case file for time dependent data. 
\param jobname - the string with the jobname.  The case fill will be jobname.case, so one should not put the file extension on the string
\param geomfile - the string withthe geometry file (needs the file extension)
\param times - a pointer to the array that holds the time 
\param nsteps - the length of times.
\param scalarvar - a list<string> that contains the names of the per node scalar variable data
\param vectorvar - a list<string> that contains the names of the per node vector variable data
\param tensorvar - a list<string> that contains the names of the per node tensor variable data
\param scalarcell - a list<string> that contains the names of the per element scalar variable data
\param vectorcell - a list<string> that contains the names of the per element vector variable data
\param tensorcell - a list<string> that contains the names of the per element tensor variable data
*/
template < class NuMeRiC >
int ensight_case( const string &jobname, const string &geomfile, const NuMeRiC *times, int nsteps,
    const list<string> &scalarvar=list<string>(), 
    const list<string> &vectorvar=list<string>(), 
    const list<string> &tensorvar=list<string>(), 
    const list<string> &scalarcell=list<string>(), 
    const list<string> &vectorcell=list<string>(), 
    const list<string> &tensorcell=list<string>() )
{
    ofstream outfile;
    if ( ensight_open_file(jobname+".case",outfile) )
        return 1;
        
    outfile << "FORMAT  \ntype:     ensight gold\n";
    outfile << "GEOMETRY\nmodel:    " << geomfile << "\n";
    write_case_variables(outfile,jobname,5,
        scalarvar,vectorvar,tensorvar,scalarcell,vectorcell,tensorcell);
    
    outfile << "TIMES\ntime set: 1\nnumber of steps: " << nsteps
        << "\nfilename start number: 0\nfilename increment: 1\ntime values:\n";
    int c=0;
    for ( int i=0; i<nsteps; ++i, ++c ) 
    {
        if ( c == 6 )
        {
            outfile << "\n";
            c=0;
        }
        outfile << setw(12) << std::scientific << times[i] << " ";
    }
    return 0;              
}


} // end namespace
#endif

