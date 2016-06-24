
#ifndef  _BASSO_VTK_H_
#define  _BASSO_VTK_H_

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <list>

#include "Basso_defs.h"
#include "Basso_ParentElement.h"

using namespace std;

namespace Basso
{

/**
 *  \brief class for a VTK data array
 *  
 *  Abstract base class for various VTK data array classes
 */
class Basso_VTKDataArray
{
	public:
		/** constructor
		 *  \param name the name of the data
		 *  \param np the number of points in the data array.  One point may have multiple components.  The 
		 *  length of the data array is np*nc
		 *  \param nc the number of components at each point.
		 */
		Basso_VTKDataArray( const string &name="none", BASSO_IDTYPE np=0, int nc=1 ) 
			: dataname(name), numpts(np), numcomp(nc)  {  } 
			
		/**
		 *  Sets the basic information of the data array
		 *  \param name the name of the data array
		 *  \param np the number of points in the data array.  One point may have multiple components.  The 
		 *  length of the data array is np*nc
		 *  \param nc the number of components at each point.
		 */
		void SetDataInfo(const string &name, BASSO_IDTYPE np, int nc)
		{
			dataname = name;
			numpts = np;
			numcomp = nc;
		}
			
		/**
		 *  Prints the data arry in the vtk xml format to an output stream
		 *  \param out - the output stream to write the array to.		 
		 */
		virtual void Print( std::ostream &out=BASSO_STDOUT  ) const = 0;
		
		/**
		 *  The number of points in the data array
		 *  \return the number of points in the data array
		 */
		BASSO_IDTYPE NumPoints( ) const { return numpts; }
		
		/**
		 *  Returns the total length of the data array.  Length = Num_points * num_components
		 *  \return the length of the array.
		 */
		BASSO_IDTYPE Length( ) const { return numpts*numcomp; }
		
		/**
		 *  Returns the number of components for each data point in the data array.
		 *  \return the number of compoennts for eacy ppoint.
		 */
		int NumComponents( ) const { return numcomp; }
		
		/**
		 *  \return the name of the data array.
		 */
		const string &Name( ) const { return dataname; }
		
	
		virtual string DataType( ) const = 0; 
		
		/**
		 *  \return The format of the data array (ascii, binary, appended)
		 */
		virtual string FormatType( ) const { return "ascii"; } 
		
		/**
		 *  \return A string with the number of componets
		 */
		virtual string NumberOfComponents( ) const 
		{
			if ( numcomp<2 )
				return "";
			ostringstream ss;
			ss << numcomp;
			return " NumberOfComponents=\""+ss.str()+"\"";
		}		
		
	protected:
		string dataname;
		BASSO_IDTYPE numpts;
		int  numcomp;
		
};

/**
 *  \brief Float VTK data array
 */
class Basso_VTKFloatDataArray : public Basso_VTKDataArray
{
	public: 
		/** constructor
		 *  \param name the name of the data
		 *  \param data a constant pointer to the data.  This is a shallow copy.
		 *  \param np the number of points in the data array.  One point may have multiple components.  The 
		 *  length of the data array is np*nc
		 *  \param nc the number of components at each point.
		 */
		Basso_VTKFloatDataArray( const string &name="none", const Basso_Numeric *data=NULL,  BASSO_IDTYPE np=0, int nc=1 ) 
				:  Basso_VTKDataArray(name,np,nc), dptr(data)  { }
		
		/**
		 *  Sets the basic information of the data array
		 *  \param name the name of the data array
		 *  \param np the number of points in the data array.  One point may have multiple components.  The 
		 *  length of the data array is np*nc
		 *  \param nc the number of components at each point.
		 */	
		void SetDataArray(const string &name, const Basso_Numeric *data,  BASSO_IDTYPE np, int nc=1 )
		{
			SetDataInfo(name,np,nc);
			dptr = data;
		}
		
		/**
		 *  Returns a string of the variable type.
		 */
		virtual string DataType( ) const 
		{ 
			int nbit = 8*sizeof(Basso_Numeric); 
			ostringstream ss;
			ss << nbit;
			return "Float"+ss.str(); 
		} 
		
		/**
		 *  Prints the data arry in the vtk xml format to an output stream
		 *  \param out - the output stream to write the array to.		 
		 */
		virtual void Print( std::ostream &out=BASSO_STDOUT  ) const;
		
	protected:
		const Basso_Numeric *dptr;
	
};


void Basso_VTKFloatDataArray::Print( std::ostream &out ) const 
{
		
	out << "<DataArray type=\""+DataType()+"\" Name=\""+Name()+"\" format=\""+FormatType()+"\" " 
			<< NumberOfComponents()+">\n";
	
	const Basso_Numeric *ptr = dptr;
	int ls=80, i;
	for ( i=0; i<Length(); ++i, ++ptr )
	{
		out << *ptr << " ";
		if ( (i+1)%ls == 0 )
			out << "\n";
	}	
	if ( i%ls != 0 ) out << "\n";
	out << "</DataArray>\n"; 
}

std::ostream &operator << ( std::ostream &out, const Basso_VTKFloatDataArray &A )
{
	A.Print( out );
	return out;
}

/**
 *  \brief integer (BASSO_IDTYPE) VTK data array
 */
class Basso_VTKIntDataArray : public Basso_VTKDataArray
{
	public: 
		/** constructor
		 *  \param name the name of the data
		 *  \param data a constant pointer to the data.  This is a shallow copy.
		 *  \param np the number of points in the data array.  One point may have multiple components.  The 
		 *  length of the data array is np*nc
		 *  \param nc the number of components at each point.
		 */
		Basso_VTKIntDataArray( const string &name="none", const BASSO_IDTYPE *data=NULL,  BASSO_IDTYPE np=0, int nc=1 ) 
				:  Basso_VTKDataArray(name,np,nc), dptr(data)  { }
		/**
		 *  Sets the basic information of the data array
		 *  \param name the name of the data array
		 *  \param np the number of points in the data array.  One point may have multiple components.  The 
		 *  length of the data array is np*nc
		 *  \param nc the number of components at each point.
		 */		
		void SetDataArray(const string &name, const BASSO_IDTYPE *data,  BASSO_IDTYPE np, int nc=1 )
		{
			SetDataInfo(name,np,nc);
			dptr = data;
		}		
		/**
		 *  Returns a string of the variable type.
		 */		
		virtual string DataType( ) const 
		{ 
			int nbit = 8*sizeof(BASSO_IDTYPE); 
			ostringstream ss;
			ss << nbit;
			return "Int"+ss.str(); 
		} 
		
		/**
		 *  Prints the data arry in the vtk xml format to an output stream
		 *  \param out - the output stream to write the array to.		 
		 */
		virtual void Print( std::ostream &out=BASSO_STDOUT  ) const;
		
	protected:
		const BASSO_IDTYPE *dptr;
	
};
void Basso_VTKIntDataArray::Print( std::ostream &out ) const 
{
		
	out << "<DataArray type=\""+DataType()+"\" Name=\""+Name()+"\" format=\""+FormatType()+"\" " 
			<< NumberOfComponents()+">\n";
	
	const BASSO_IDTYPE *ptr = dptr;
	int ls=80, i;
	for ( i=0; i<Length(); ++i, ++ptr )
	{
		out << *ptr << " ";
		if ( (i+1)%ls == 0 )
			out << "\n";
	}	
	if ( i%ls != 0 ) out << "\n";
	out << "</DataArray>\n"; 
} 

std::ostream &operator << ( std::ostream &out, const Basso_VTKIntDataArray &A )
{
	A.Print( out );
	return out;
}

/**
 *  \brief unsigned int VTK data array
 */
class Basso_VTKUIntDataArray : public Basso_VTKDataArray
{
	public: 
		/** constructor
		 *  \param name the name of the data
		 *  \param data a constant pointer to the data.  This is a shallow copy.
		 *  \param np the number of points in the data array.  One point may have multiple components.  The 
		 *  length of the data array is np*nc
		 *  \param nc the number of components at each point.
		 */
		Basso_VTKUIntDataArray( const string &name="none", const unsigned int *data=NULL,  BASSO_IDTYPE np=0, int nc=1 ) 
				:  Basso_VTKDataArray(name,np,nc), dptr(data)  { }
				
		/**
		 *  Sets the basic information of the data array
		 *  \param name the name of the data array
		 *  \param np the number of points in the data array.  One point may have multiple components.  The 
		 *  length of the data array is np*nc
		 *  \param nc the number of components at each point.
		 */		
		void SetDataArray(const string &name, const unsigned int *data,  BASSO_IDTYPE np, int nc=1 )
		{
			SetDataInfo(name,np,nc);
			dptr = data;
		}
		
		/**
		 *  Returns a string of the variable type.
		 */
		virtual string DataType( ) const 
		{ 
			int nbit = 8*sizeof(unsigned int); 
			ostringstream ss;
			ss << nbit;
			return "UInt"+ss.str(); 
		} 
		/**
		 *  Prints the data arry in the vtk xml format to an output stream
		 *  \param out - the output stream to write the array to.		 
		 */
		virtual void Print( std::ostream &out=BASSO_STDOUT  ) const;
		
	protected:
		const unsigned int *dptr;
		//const uint8_t *dptr;
	
};
void Basso_VTKUIntDataArray::Print( std::ostream &out ) const 
{
		
	out << "<DataArray type=\""+DataType()+"\" Name=\""+Name()+"\" format=\""+FormatType()+"\" " 
			<< NumberOfComponents()+">\n";
	
	const unsigned int *ptr = dptr;
	int ls=80, i;
	for ( i=0; i<Length(); ++i, ++ptr )
	{
		out << *ptr << " ";
		if ( (i+1)%ls == 0 )
			out << "\n";
	}	
	if ( i%ls != 0 ) out << "\n";
	out << "</DataArray>\n"; 
} 

std::ostream &operator << ( std::ostream &out, const Basso_VTKUIntDataArray &A )
{
	A.Print( out );
	return out;
}

/**
 *  \brief Abstracts the VTK Points data
 *  
 */
class Basso_VTKPoints 
{
	public:
		/** constructor
		 *  \param ncoords a constant pointer to the data.  This is a shallow copy and assumes the nodes are in (x1 y1 z1, x2 y2 z2 ..) order
		 *  \param np the number of points (nodes) in the data array. 
		 *  \param nc the number of components at each point (default is 3).
		 */
		 Basso_VTKPoints ( const Basso_Numeric *ncoords=NULL, BASSO_IDTYPE np=0, int nc=3 )
			: coords("node coordinates",ncoords,np,nc) { }
		
		/** Sets the data 
		 *  \param ncoords a constant pointer to the data.  This is a shallow copy and assumes the nodes are in (x1 y1 z1, x2 y2 z2 ..) order
		 *  \param np the number of points (nodes) in the data array. 
		 *  \param nc the number of components at each point (default is 3).
		 */
		void SetPoints( const Basso_Numeric *ncoords, BASSO_IDTYPE np, int nc=3  )
			{ coords.SetDataArray("node coordinates",ncoords,np,nc); }
	
		/**
		 *  Prints the Points data in the vtk xml format to an output stream
		 *  \param out - the output stream to write the array to.		 
		 */
		virtual void Print( std::ostream &out=BASSO_STDOUT  ) const;
		
		/**\return teh number of nodes (points) in the class*/
		BASSO_IDTYPE NumPoints( ) const { return coords.NumPoints(); }
		
	protected:
		Basso_VTKFloatDataArray coords;
};

void Basso_VTKPoints::Print( std::ostream &out  ) const
{
	out << "<Points>\n" 
		<< coords 
	    << "</Points>\n";
}

std::ostream &operator << ( std::ostream &out, const Basso_VTKPoints &A )
{
	A.Print( out );
	return out;
}

/**
 *  \brief Abstracts the VTK Cell data
 *  
 */
class Basso_VTKCells 
{
	public:
	
		/** 
		 *  void constructor
		 */
		 Basso_VTKCells (  )  { }
			
		/** constructor tha tinitializes from a block of elements of the same type.
		 *  \param conn a constant pointer to the element conectivity data.  This is a shallow copy and 
		 *     assumes the connectivities are in (n1 n2 n3 .. n1 n2 n3 ...) order
		 *  \param etype the elemetn topology type
		 *  \param np is the totla length of the array conn. (Not the number of elements)
		 */
		 Basso_VTKCells ( const BASSO_IDTYPE *conn, const Basso_ParentElement *etype,
			BASSO_IDTYPE np )
			: connectivity("connectivity",conn,np) { SetElementData(etype); }
	
		/** initializes from a block of elements of the same type.
		 *  \param conn a constant pointer to the element conectivity data.  This is a shallow copy and 
		 *     assumes the connectivities are in (n1 n2 n3 .. n1 n2 n3 ...) order
		 *  \param etype the elemetn topology type
		 *  \param np is the totla length of the array conn. (Not the number of elements)
		 */
		void SetCells( const BASSO_IDTYPE *conn, const Basso_ParentElement *etype, BASSO_IDTYPE np  )
			{ connectivity.SetDataArray("connectivity",conn,np); SetElementData(etype); }
		/**
		 *  Prints the Cell data in the vtk xml format to an output stream
		 *  \param out - the output stream to write the array to.		 
		 */
		virtual void Print( std::ostream &out=BASSO_STDOUT  ) const;
		
		/**
		 *  \return The number of cells/elements in the class.
		 */
		BASSO_IDTYPE NumCells( ) const { return types.NumPoints(); }
		
	protected:
		void SetElementData( const Basso_ParentElement *etype );
		unsigned int ElementTypeID( const Basso_ElementType &et ) const;
		
	protected:
		Basso_VTKIntDataArray connectivity;
		Basso_VTKIntDataArray offsets;
		Basso_VTKUIntDataArray types;
		
		Basso_Array<BASSO_IDTYPE> offdata;
		Basso_Array<unsigned int> typedata;
};

void Basso_VTKCells::SetElementData( const Basso_ParentElement *etype )
{
	
	int nne=etype->NumNodes(), eid=ElementTypeID(etype->Type());
	BASSO_IDTYPE ne = connectivity.Length()/nne;
	
	offdata.Resize(ne);
	typedata.Resize(ne);
	
	BASSO_IDTYPE off=nne;
	for ( int e=0; e<ne; ++e, off+=nne )
	{
		typedata[e] = eid;
		offdata[e] = off;
	}
	
	offsets.SetDataArray("offsets",offdata.Data(),offdata.Length());
	types.SetDataArray("types",typedata.Data(),offdata.Length());
}

void Basso_VTKCells::Print( std::ostream &out  ) const
{
	out << "<Cells>\n" 
		<< connectivity << offsets << types 
	    << "</Cells>\n";
}		

unsigned int Basso_VTKCells::ElementTypeID( const Basso_ElementType &et ) const
{
	switch (et)
	{
		case Basso_LINE2:
			return 3; 
		case Basso_TRIA3:
			return 5;  	
		case Basso_QUAD4 :
			return 9; 	
		case Basso_TETRA4 :
			return 10; 	
		case Basso_HEXA8 :
			return 12; 	
		case Basso_PRISM6 :
			return 13; 	
		case Basso_PYRAMID5 :
			return 14; 	
		case Basso_LINE3 :
			return 21; 	
		case Basso_TRIA6 :
			return 22; 	
		case Basso_QUAD8 :
			return 23; 	
		case Basso_TETRA10 :
			return 24; 	
		case Basso_HEXA20:
			return 25;  
		case Basso_POINT1 :
			return 1; 
		default:
			Basso_Warning("Basso_VTKCells::ElementTypeID","Element type not supported");
			return 1; 
	}
}

std::ostream &operator << ( std::ostream &out, const Basso_VTKCells &A )
{
	A.Print( out );
	return out;
}

/**
 *  \brief Abstracts the data for a VTK PointData
 *  Basically a list of Basso_VTKDataArray.  Currently only supports scalar and vector data types.
 */
class Basso_VTKPointData
{
	public:
		/**
		 *  void constructor.
		 */
		Basso_VTKPointData( ) { }
		
		/**
		 *  Adds a data array to the point data set.
		 *  \param data a pointer to the data array.  This is a shallow copy.
		 */
		virtual void AddDataArray( const Basso_VTKDataArray *data ) { datalist.push_back(data); } 
		/**
		 *  Prints the Point data array in the vtk xml format to an output stream
		 *  \param out - the output stream to write the array to.		 
		 */
		virtual void Print( std::ostream &out=BASSO_STDOUT  ) const;
		/**
		 *  \return The number of data arrays in the the class.
		 */
		virtual int NumDataArrays( ) const { return datalist.size(); }
		
		/**
		 *  \return a string with the names of the Scalars and Vectors data variables.
		 */
		virtual string NamesString( ) const;
		
	protected:
		virtual string MyNameIs() const { return "PointData"; }
		
	protected:
		list< const Basso_VTKDataArray * > datalist;
};

string Basso_VTKPointData::NamesString( ) const
{
	set<string> scalars, vectors;
	for ( list< const Basso_VTKDataArray * >::const_iterator itr=datalist.begin();
				itr!=datalist.end(); ++itr )
		if ( (*itr)->NumComponents() == 1 )
			scalars.insert((*itr)->Name());
		else
			vectors.insert((*itr)->Name());
		
	string dnames="";
	
	if ( !scalars.empty() )
	{		
		dnames += " Scalars=\" ";
		for ( set<string>::const_iterator itr=scalars.begin(); itr!=scalars.end(); ++itr )
			dnames += (*itr)+" ";
		dnames += "\"";
	}
	
	if ( !vectors.empty() )
	{		
		dnames += " Vectors=\" ";
		for ( set<string>::const_iterator itr=vectors.begin(); itr!=vectors.end(); ++itr )
			dnames += (*itr)+" ";
		dnames += "\"";
	}
	
	return dnames;
}
		
void Basso_VTKPointData::Print( std::ostream &out  ) const
{
	out << "<"+MyNameIs() << NamesString() << ">\n"; 
	for ( list<const Basso_VTKDataArray*>::const_iterator itr=datalist.begin();
				itr!=datalist.end(); ++itr )
		(*itr)->Print(out);		
	out << "</"+MyNameIs()+">\n";
}		

std::ostream &operator << ( std::ostream &out, const Basso_VTKPointData &A )
{
	A.Print( out );
	return out;
}

/**
 *  \brief Abstracts the data for a VTK CellData
 *  Basically a list of Basso_VTKDataArray.  Currently only supports scalar and vector data types.
 */
class Basso_VTKCellData : public Basso_VTKPointData
{
	public:
		/**
		 *  void constructor.
		 */
		Basso_VTKCellData( ) : Basso_VTKPointData( ) { }
		
	protected:
		virtual string MyNameIs() const { return "CellData"; }
		
};

std::ostream &operator << ( std::ostream &out, const Basso_VTKCellData &A )
{
	A.Print( out );
	return out;
}

/**
 *  \brief
 */
class Basso_VTKFile
{
	public:
		Basso_VTKFile ( const string &jobname ) : job(jobname) { }
	
		virtual void Print( std::ostream &out=BASSO_STDOUT  ) const=0;
		
		virtual void AddPointData( const Basso_VTKDataArray *data ) 
				{ pointdata.AddDataArray(data); } 
		virtual void AddCellData( const Basso_VTKDataArray *data ) 
				{ celldata.AddDataArray(data); } 
	
		virtual void WriteFile( ) const;
		virtual string FileName( ) const { return job+FileExtension(); }
		
	protected:
		virtual string FileExtension( ) const = 0;
		virtual void PrintData( std::ostream &out=BASSO_STDOUT ) const;
		
	protected:
		Basso_VTKCellData celldata;
		Basso_VTKPointData pointdata;
		
		string job;
};

void Basso_VTKFile::PrintData( std::ostream &out ) const
{
	pointdata.Print(out);
	celldata.Print(out);
}

void Basso_VTKFile::WriteFile( ) const 
{
    ofstream outfile;
	outfile.open( FileName().c_str(), ios::out );
    if ( !outfile )
	{
		Basso_Warning("Basso_VTKFile::WriteFile","cannot open file");
        return;
	}
	
	outfile << "<?xml version=\"1.0\"?>\n";
	Print(outfile);
    
	outfile.close();
}

/**
 *  \brief
 */
class Basso_VTKUnstructuredGrid : public Basso_VTKFile
{
	public:
		Basso_VTKUnstructuredGrid( const string &jobname="results" ) : Basso_VTKFile(jobname) { }
		
		void SetMesh( const Basso_nMatrix &nodes, const Basso_Array2D<BASSO_IDTYPE> &conn,
					const Basso_ParentElement *etype )
		{
			points.SetPoints(nodes.Data(),nodes.NumCols(),nodes.NumRows());
			cells.SetCells(conn.Data(),etype,conn.Length());
		}
	
		virtual void Print( std::ostream &out=BASSO_STDOUT  ) const;
	
	protected:
		virtual string FileExtension( ) const { return ".vtu"; }
	
	protected:
	
		Basso_VTKPoints points;
		Basso_VTKCells  cells;

};

void Basso_VTKUnstructuredGrid::Print( std::ostream &out ) const
{
	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >\n"
		<< "<UnstructuredGrid>\n" 
		<< "<Piece NumberOfPoints=\"" << points.NumPoints() 
					<< "\" NumberOfCells=\"" << cells.NumCells() << "\">\n";  
		
	points.Print(out);
	cells.Print(out);
	PrintData(out);
	
	out << "</Piece>\n</UnstructuredGrid>\n</VTKFile>";
}

std::ostream &operator << ( std::ostream &out, const Basso_VTKUnstructuredGrid &A )
{
	A.Print( out );
	return out;
}

/**
 *  \brief Serial VTK XML based output for unstrictured grids
 */
class Basso_VTKPUnstructuredGrid 
{
	public:
	/**
	 *  constructor
	 *  \param jobname The name of the resulting .vtu file  (jobname.vtu)
	 */
		Basso_VTKPUnstructuredGrid(const string &jobname="results" ) : job(jobname) { }
	/**
	 *  Prints the vtk file to the ostream
	 *  \param out the output stream to print to	 
	 */
		virtual void Print( std::ostream &out=BASSO_STDOUT  ) const;
	/**
	 *  Writes the full VTK file 
	 */
		virtual void WriteFile( ) const;
	/**
	 *  Returns the filename with the extension	
	 */
		virtual string FileName( ) const { return job+FileExtension(); }
		
	protected:
	/**
	 *  the extension of the file
	 */
		virtual string FileExtension( ) const { return ".pvtu"; }
		
	protected:
		list < Basso_VTKUnstructuredGrid > grids;
		string job;
};

void Basso_VTKPUnstructuredGrid::WriteFile( ) const 
{
    ofstream outfile;
	outfile.open( FileName().c_str(), ios::out );
    if ( !outfile )
	{
		Basso_Warning("Basso_VTKPUnstructuredGrid::WriteFile","cannot open file");
        return;
	}
	
	outfile << "<?xml version=\"1.0\"?>\n";
	Print(outfile);
	outfile.close();
}

void Basso_VTKPUnstructuredGrid::Print( std::ostream &out ) const
{
	out << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >\n"
		<< "<PUnstructuredGrid GhostLevel=\"0\">\n" ;
		
	/* 	<PPointData>
		</PPointData> 
		<PCellData>
		</PCellData> 
		<PPoints>
		</PPoints> 
		<Piece Source=”unstructuredGrid0.vtu”/>  */
		
	out << "</PUnstructuredGrid>\n</VTKFile>";
}

std::ostream &operator << ( std::ostream &out, const Basso_VTKPUnstructuredGrid &A )
{
	A.Print( out );
	return out;
}

	
} // end namespace
#endif

