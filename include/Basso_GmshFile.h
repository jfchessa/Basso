/*! \file Basso_GmshFilet.h

	\brief Defines a class to process gmsh files

	
\author Jack Chessa, jfchessa@utep.edu
\date Thursday, June 8, 2016

*/

#ifndef _BASSO_GMSH_FILE_H_
#define _BASSO_GMSH_FILE_H_

// std includes

// trilinos includes

// basso includes
#include "Basso_defs.h"
#include "Basso_ElementList.h"

namespace Basso
{

class Basso_GmshFilter
{
	public:
		Basso_GmshFilter( const string &line )
		{ Set(line); }
			
		Basso_GmshFilter( unsigned int t=0, unsigned int ge=0, unsigned int ph=0, unsigned int pr=0 )
			{ Set(t,ge,ph,pr); }
		
		void Set(unsigned int t=0, unsigned int pa=0, unsigned int ph=0, unsigned int pr=0)
			{ type=t; geomid=pa; physid=ph; procid=pr; }
			
		void Set( const string &line )
		{ 
			int ntag, eid, nprocid=0, itag=0;
			procid=0;
			istringstream ss(line);
			ss >> eid >> type >> ntag >> physid >> geomid;	
			itag=2;
			if ( ntag > itag )
				ss >> nprocid >> procid;
		}
		
		void Clear( ) { type=0; physid=0; geomid=0; procid=0; }
			
		bool Match( unsigned int t, unsigned int pa, unsigned int ph, unsigned int pr ) const
		{
			if ( type!=0   & type!=t ) return false;
			if ( geomid!=0 & geomid!=pa ) return false;
			if ( physid!=0 & physid!=ph ) return false;
			if ( procid!=0 & procid!=pr ) return false;
			return true;
		}
		
		/* bool Match( const string &line ) const
		{
			int ntags, npart, p, procid, eid, geomid;
			istringstream ss(line);
		//	ss >> eid >> type >> ntags >> physid >> geomid >> npart;
			if ( npart == 0 )
				procid = 0;
		//	ss >> procid;
			return Match(type,physid,geomid,procid);
		} */
		
		bool Match( const Basso_GmshFilter &f ) const
		{ return Match(f.type,f.geomid,f.physid,f.procid); }
		
		void Print( std::ostream &out=BASSO_STDOUT ) const;
		
		bool operator = ( const Basso_GmshFilter &b ) const
			{ return Match(b); }
			
		bool operator < ( const Basso_GmshFilter &b ) const 
		{
			if ( type < b.type ) return true;
			if ( type > b.type ) return false;
			
			if ( geomid < b.geomid ) return true;
			if ( geomid > b.geomid ) return false;
			
			if ( physid < b.physid ) return true;
			if ( physid > b.physid ) return false;
			
			if ( procid < b.procid ) return true;
			if ( procid > b.procid ) return false;
			
			return false;	
		}
		
		unsigned int type;
		unsigned int geomid;
		unsigned int physid;
		unsigned int procid;
};

void Basso_GmshFilter::Print( std::ostream &out ) const
{
	out << "Basso_GmshFilter:";
	if ( type ) out << " type=" << type;
	if ( geomid ) out << " gid=" << geomid;
	if ( physid ) out << " pid=" << physid;
	if ( procid ) out << " partition=" << procid;
	out << " ";
}

std::ostream &operator << ( std::ostream &out, const Basso_GmshFilter &A )
{
	A.Print( out );
	return out;
}

enum GmshElementType 
{  
    Gmsh_NoElem=0, Gmsh_Line2=1, Gmsh_Tria3, Gmsh_Quad4, Gmsh_Tetra4, Gmsh_Hexa8, 
    Gmsh_Prism6, Gmsh_Pyramid5, Gmsh_Line3, Gmsh_Tria6, Gmsh_Quad9, 
    Gmsh_Tetra10, Gmsh_Hexa27, Gmsh_Prism18, Gmsh_Pyramid14, Gmsh_Point1, 
    Gmsh_Quad8, Gmsh_Hexa20, Gmsh_Prism15, Gmsh_Pyramid13, Gmsh_Tria9, 
    Gmsh_Tria10, Gmsh_Tria12inc, Gmsh_Tria15, Gmsh_Tria15inc, Gmsh_Tria21,
    Gmsh_Line4, Gmsh_Line5, Gmsh_Line6, Gmsh_Tetra20, Gmsh_Tetra35, 
    Gmsh_Tetra56, Gmsh_Hexa64=92, Gmsh_Hexa125=93 
};

Basso_ElementType convert_gmsh_etype( GmshElementType etype )
{
	switch (etype)
	{
		case Gmsh_Line2:        return Basso_LINE2;
		case Gmsh_Tria3:        return Basso_TRIA3; 
		case Gmsh_Quad4:        return Basso_QUAD4; 
		case Gmsh_Tetra4:       return Basso_TETRA4; 
		case Gmsh_Hexa8:        return Basso_HEXA8; 
		case Gmsh_Prism6:       return Basso_PRISM6; 
		case Gmsh_Pyramid5:     return Basso_PYRAMID5; 
		case Gmsh_Line3:        return Basso_LINE3; 
		case Gmsh_Tria6:        return Basso_TRIA6; 
		case Gmsh_Quad9:        return Basso_QUAD9; 
		case Gmsh_Tetra10:      return Basso_TETRA10; 
		case Gmsh_Hexa27:       return Basso_HEXA27; 
		case Gmsh_Prism18:      return Basso_PRISM18; 
		case Gmsh_Pyramid14:    return Basso_PYRAMID14; 
		case Gmsh_Point1:       return Basso_POINT1; 
		case Gmsh_Quad8:        return Basso_QUAD8; 
		case Gmsh_Hexa20:       return Basso_HEXA20; 
		case Gmsh_Prism15:      return Basso_PRISM15; 
		case Gmsh_Pyramid13:    return Basso_PYRAMID13; 
		case Gmsh_Tria9:        return Basso_TRIA9; 
		case Gmsh_Tria10:       return Basso_TRIA10; 
		case Gmsh_Tria12inc:    return Basso_TRIA12INC; 
		case Gmsh_Tria15:       return Basso_TRIA15; 
		case Gmsh_Tria15inc:    return Basso_TRIA15INC; 
		case Gmsh_Tria21:       return Basso_TRIA21;
		case Gmsh_Line4:        return Basso_LINE4; 
		case Gmsh_Line5:        return Basso_LINE5; 
		case Gmsh_Line6:        return Basso_LINE6; 
		case Gmsh_Tetra20:      return Basso_TETRA20; 
		case Gmsh_Tetra35:      return Basso_TETRA35; 
		case Gmsh_Tetra56:      return Basso_TETRA56;
		default:                return Basso_NONE;
	}
}

unsigned int gmsh_num_elem_nodes( GmshElementType etype ) 
{
	switch (etype)
	{
		case Gmsh_Line2:        return 2;
		case Gmsh_Tria3:        return 3; 
		case Gmsh_Quad4:        return 4; 
		case Gmsh_Tetra4:       return 4; 
		case Gmsh_Hexa8:        return 8; 
		case Gmsh_Prism6:       return 6; 
		case Gmsh_Pyramid5:     return 5; 
		case Gmsh_Line3:        return 3; 
		case Gmsh_Tria6:        return 6; 
		case Gmsh_Quad9:        return 9; 
		case Gmsh_Tetra10:      return 10; 
		case Gmsh_Hexa27:       return 27; 
		case Gmsh_Prism18:      return 18; 
		case Gmsh_Pyramid14:    return 14; 
		case Gmsh_Point1:       return 1; 
		case Gmsh_Quad8:        return 8; 
		case Gmsh_Hexa20:       return 20; 
		case Gmsh_Prism15:      return 15; 
		case Gmsh_Pyramid13:    return 13; 
		case Gmsh_Tria9:        return 9; 
		case Gmsh_Tria10:       return 10; 
		case Gmsh_Tria12inc:    return 12; 
		case Gmsh_Tria15:       return 15; 
		case Gmsh_Tria15inc:    return 15; 
		case Gmsh_Tria21:       return 21;
		case Gmsh_Line4:        return 4; 
		case Gmsh_Line5:        return 5; 
		case Gmsh_Line6:        return 6; 
		case Gmsh_Tetra20:      return 20; 
		case Gmsh_Tetra35:      return 35; 
		case Gmsh_Tetra56:      return 56; 
		case Gmsh_Hexa64:       return 64; 
		case Gmsh_Hexa125:      return 125; 
		default:                return 0;
	}
}
		
class Basso_GmshElement
{
	public:
	
		Basso_GmshElement( const string &line )
		{ Set(line); }
		
		void Set( const string &line )
		{
			int ntag, id, nprocid=0, itag=0;
			nghostid=0;
			procid=0;
			//ndata=0;
			istringstream ss(line);
			ss >> eid >> id >> ntag >> physid >> geomid;	
			type=static_cast<GmshElementType>(id);
			itag=2;
			if ( ntag > itag )
			{
				ss >> nprocid;
				for  ( int i=0; i<nprocid; ++i )
				{
					ss >> id;
					if ( id > 0 )
						procid = id;
					else
					{
						ghostid[nghostid] = -id;
						nghostid++;
					}
				}
				itag += nprocid+1;
			}
			/* for ( int i=0; i<ntag-itag; ++i )
			{
				cout << "WOHW\n";
				ss >> data[ndata];
				++ndata;
			}  */
			nne = gmsh_num_elem_nodes(type);
			for  ( int i=0; i<nne; ++i )
				ss >> conn[i];
		}
		
		void Print( std::ostream &out=BASSO_STDOUT ) const;
		
	public:
		GmshElementType type;
		unsigned int eid;
		unsigned int geomid;
		unsigned int physid;
		unsigned int procid;
		unsigned int nghostid;
		unsigned int ghostid[10];
		unsigned int nne;
		unsigned int conn[20];
		/* unsigned int ndata;
		int data[5]; */
};

void Basso_GmshElement::Print( std::ostream &out ) const
{
	out << "Basso_GmshElement: eid=" << eid
		<< ", type=" << static_cast<int>(type)
		<< ", physical id=" << physid
		<< ", geometry id=" << geomid
		<< ", partition id=" << procid;
		
	if ( nghostid )
	{
		out << ", ghost partitions = { ";
		for ( int i=0; i<nghostid; ++i )
		{
			out << ghostid[i];
			if ( i!=nghostid-1 )
				out << ", ";
		}
		out << " }";
	}
	
	out << "\n\tconnectivity = { ";
	for ( int i=0; i<nne; ++i )
	{
		out << conn[i];
		if ( i!=nne-1 )
			out << ", ";
	}
	out << " }\n";
	
	/* if ( ndata )
	{
		out << "\tother data = { ";
		for ( int i=0; i<ndata; ++i )
		{
			out << data[i];
			if ( i!=ndata-1 )
				out << ", ";
		}
		out << " }\n";
	} */
}

std::ostream &operator << ( std::ostream &out, const Basso_GmshElement &A )
{
	A.Print( out );
	return out;
}
/**
 \brief Class to process gmsh files

*/
class Basso_GmshFile
{

	public:
		Basso_GmshFile( const string &file ) : gmshfile(file), partid(0) { }
		
		const string &Filename( ) const { return gmshfile; }
		
		void ReadPhysicalIDs( set<int> &pids ) const;
		
		int NumPartitions( ) const;
		
		void SetElementPhysicalIDs( const set<int> &pids ) { physids = pids; }
		
		void SetPartition( int pid ) { partid = pid; }
		
		//void SetNodeSetPhysicalIDs( set<int> &pids );
		
		// void ReadNodeSets ( ) const;
		
		void ReadNodes( Basso_nMatrix &nodes, Basso_Array<BASSO_IDTYPE> &nids, 
				const set<BASSO_IDTYPE> &n2read=set<BASSO_IDTYPE>() ) const;
				
		void ReadElementBlocks( Basso_ElementBlockList &elem ) const;
		
	//protected:
		bool ProcessElement( Basso_GmshFilter &efilt, int etype=0 ) const
		{
			efilt.geomid = 0;		// ignore geometry id
			
			// does the element type pass
			if ( etype != 0 )  // then check
				if ( etype != efilt.type )
					return false;
			
			// does the physid pass
			if ( physids.size() > 0 )  // then check
				if ( physids.count(efilt.physid) == 0 ) // physical id not in set
					return false;
			
			// does the procid pass
			if ( partid > 0 ) // then check
			{
				if ( partid != efilt.procid )
					return false;
			}
			else
				efilt.procid = 0;
			
			return true;
		}
	
	protected:
	
		string gmshfile;
		
		set< int > physids;
		set< int > nodesetids;
		int partid;
		
};

		
void Basso_GmshFile::ReadNodes( Basso_nMatrix &nodes, Basso_Array<BASSO_IDTYPE> &nids, 
				 const set<BASSO_IDTYPE> &n2read ) const
{
	bool readall=false;
	if ( n2read.size( ) == 0 )
		readall=true;
	
	fstream infile;
    infile.open( gmshfile.c_str(), ios::in );
    if ( !infile )
	{
		Basso_Warning("Basso_GmshFile::ReadNodes","could not open file");
        return;
	}
	
    string line;
	BASSO_IDTYPE lnn, nid, numnode=0;
    while ( !infile.eof() )
    {
        getline( infile, line );
        if ( line.compare(0,5,"$Node")==0 )
        {
			getline( infile, line );
            lnn = atoi(line.c_str()); 
            for ( BASSO_IDTYPE n=0; n<lnn; ++n )
            {
				getline( infile, line );
				istringstream ss(line);
				ss >> nid;	
				if ( readall | n2read.count(nid)>0 )
					++numnode;
            }
        }
    }

	nodes.Resize(3,numnode);
	nids.Resize(numnode);
	
	infile.clear();
	infile.seekg( 0, ios::beg );
	
	BASSO_IDTYPE lnid=0;
	Basso_Numeric x, y, z;
    while ( !infile.eof() )
    {
        getline( infile, line );
        if ( line.compare(0,5,"$Node")==0 )
        {
			getline( infile, line );
            lnn = atoi(line.c_str()); 
            for ( int n=0; n<lnn; ++n )
            {
				getline( infile, line );
				istringstream ss(line);
				ss >> nid >> x >> y >>  z;	
				if ( readall |  n2read.count(nid)>0 )
				{
					nodes[lnid][0]=x;
					nodes[lnid][1]=y;
					nodes[lnid][2]=z;
					nids[lnid]=nid;
					++lnid;
				}
            }
        }
    }
	
	infile.close();	
}

int  Basso_GmshFile::NumPartitions( ) const
{
	set<int> pids;
	fstream infile;
    infile.open( gmshfile.c_str(), ios::in );
    if ( !infile )
	{
		Basso_Warning("Basso_GmshFile::ReadPhysicalIDs","could not open file");
        return 0;
	}
	
    string line;
	int lne;
	Basso_GmshFilter efilt;
    while ( !infile.eof() )
    {
        getline( infile, line );
        if ( line.compare(0,8,"$Element")==0 )
        {
			getline( infile, line );
            lne = atoi(line.c_str()); 
            for ( int n=0; n<lne; ++n )
            {
				getline( infile, line );
				efilt.Set(line);
				pids.insert(efilt.procid);
            }
        }
    }
	infile.close();
	return pids.size();
}

void Basso_GmshFile::ReadPhysicalIDs( set<int> &pids ) const
{
	fstream infile;
    infile.open( gmshfile.c_str(), ios::in );
    if ( !infile )
	{
		Basso_Warning("Basso_GmshFile::ReadPhysicalIDs","could not open file");
        return;
	}
	
    string line;
	int lne;
	Basso_GmshFilter efilt;
    while ( !infile.eof() )
    {
        getline( infile, line );
        if ( line.compare(0,8,"$Element")==0 )
        {
			getline( infile, line );
            lne = atoi(line.c_str()); 
            for ( int n=0; n<lne; ++n )
            {
				getline( infile, line );
				efilt.Set(line);
				if ( efilt.procid == partid | partid==0 )  // check that element is on the partition
					pids.insert(efilt.physid);
            }
        }
    }	
	
	infile.close();
}

void Basso_GmshFile::ReadElementBlocks( Basso_ElementBlockList &elems ) const
{	
	fstream infile;
	cout << "opening " << gmshfile.c_str() << "\n";
    infile.open( gmshfile.c_str(), ios::in );
    if ( !infile.is_open() )
	{
		Basso_Warning("Basso_GmshFile::ReadElementBlocks","could not open file");
        return;
	}
	
    string line;
	int lne;
	map<Basso_GmshFilter,unsigned int> blocks;
	Basso_GmshFilter efilt;
    while ( !infile.eof() )
    {
        getline( infile, line );
        if ( line.compare(0,8,"$Element")==0 )
        {
			getline( infile, line );
            lne = atoi(line.c_str()); 
            for ( int n=0; n<lne; ++n )
            {
				getline( infile, line );
				efilt.Set(line);
				if ( ProcessElement(efilt) )
					blocks[efilt] += 1;
            }
        }
    }
	cout << "element blocks " << blocks << "\n";
	
	// initialize the element blocks
	map<Basso_GmshFilter,unsigned int>::const_iterator bitr;
	for ( bitr=blocks.begin(); bitr!=blocks.end(); ++bitr )
		elems.AddElementBlock( bitr->first.type, bitr->second, bitr->first.physid ); 
	
	// read connectivities
	Basso_ElementBlockList::iterator itr;
	bitr = blocks.begin();
	for ( itr=elems.Begin(); itr!=elems.End(); ++itr, ++bitr )
	{
		Basso_ElementType type= static_cast<Basso_ElementType>(bitr->first.type);
		
		BASSO_IDTYPE *connptr = itr->Data();
		BASSO_IDTYPE *eidptr = itr->DataEID();
		int nne = itr->NNE();
		
		infile.clear();
		infile.seekg( 0, ios::beg );		
		while ( !infile.eof() )
		{
			getline( infile, line );
			if ( line.compare(0,8,"$Element")==0 )
			{
				getline( infile, line );
				for ( int n=0; n<lne; ++n )
				{
					getline( infile, line );    		// DOES NOT GO HERE WTH!?!
					
					efilt.Set(line);					// set filter 
					efilt.geomid = 0;  				// ignore geometry id
					
					if ( (bitr->first).Match(efilt) ) 	// does it pass the filter test? 
					{
						Basso_GmshElement newelem( line );
						*eidptr = newelem.eid;
						++eidptr;
						for ( int i=0; i<nne; ++i, ++connptr )
							*connptr = newelem.conn[i];
					}
						
				}
			}
		}
	}
	
	infile.close();
}

} // end of namspace
#endif