/*! \file Basso_ElementList.h

	\brief Defines the class Basso_ElemetList

	
\author Jack Chessa, jfchessa@utep.edu
\date Thursday, June 8, 2016

*/

#ifndef  _BASSO_ELEMENT_LIST_H_
#define  _BASSO_ELEMENT_LIST_H_

// std includes
#include <list>
#include <map>

// trilinos includes

// basso includes
#include "Basso_defs.h"

#include "Basso_Point1.h"
#include "Basso_Line2.h"
#include "Basso_Line3.h"
#include "Basso_Tria3.h"
#include "Basso_Tria6.h"
#include "Basso_Quad4.h"
#include "Basso_Quad8.h"
#include "Basso_Tetra4.h"
#include "Basso_Tetra10.h"
#include "Basso_Hexa8.h"
#include "Basso_Pyramid5.h"

using namespace std;

namespace Basso
{

/**
 * \brief Block of similar elements 
 
 This holds a block of similar elements with the element connectivity, a reference to
 the element formulation and global element ids.
 */
class Basso_ElementBlock : public Basso_Array2D<BASSO_IDTYPE>
{

	public:
	/** Constructor
		\param ne - The number of elements
		\param nne - The number of nodes in each element
		\param pid - An identifying part id number.
	*/
		Basso_ElementBlock( BASSO_IDTYPE ne, int nne, int pid=0 ) 
			: Basso_Array2D<BASSO_IDTYPE>(nne,ne), partid(pid), geids_(ne), localnum(false) 
			{ elemptr_ = NULL; SetGIDs2LIDs( ); }
		
		/** Constructor
		\param type - The type of element in the block
		\param ne - The number of elements in the block
		\param pid - An identifying part id number.
		*/
		Basso_ElementBlock( Basso_ElementType type, BASSO_IDTYPE ne, int pid=0 ) 
			: partid(pid), localnum(false) 
		{ 
			elemptr_ = NULL;
			SetBasis(type);
			this->Resize( elemptr_->NumNodes(), ne );
			geids_.Resize(ne);
			SetGIDs2LIDs( );
		}
		
		/** 
		 *  Copy constructors (deep copy)
		 */
		Basso_ElementBlock( const Basso_ElementBlock &b )  
		{ 
			Copy(b); 
		}	
		
		/** destructor */	
		~Basso_ElementBlock( ) 
		{ 
			if (elemptr_!= NULL)  
			{
				delete elemptr_;
				elemptr_=NULL;
			}
		}
		
		/** copy (This is a deep copy so beware) */
		void Copy( const Basso_ElementBlock &b ) 
		{
			if ( &b == this ) return;
			elemptr_ = NULL;
			partid = b.partid;
			geids_ = b.geids_;
			SetBasis( b.elemptr_->Type() );	// new pointer allocated
		    Basso_Array2D<BASSO_IDTYPE> *baseptr = this;
			baseptr->Copy(b);
		}
		
		
		/** assignemnet operator*/
		Basso_ElementBlock &operator = ( const Basso_ElementBlock &b ) { Copy(b); return *this; }
		Basso_ElementBlock &operator = ( Basso_ElementBlock &b ) { Copy(b); return *this; }
		
		/**
		\return a constant pointer to the element basis
		*/
		const Basso_ParentElement *ElementBasis( ) const { return elemptr_; }
		
		/**
		\return true if the element basis is set otherwise returns false.
		*/
		bool BasisSet( ) const { return (elemptr_==NULL ? false : true); }
		
		/** Pointer to the element connectivity data
		*/
		BASSO_IDTYPE *DataEID( ) { return geids_.Data(); }
		
		/** Constant pointer to the element connectivity data
		*/
		const BASSO_IDTYPE *DataEID( ) const { return geids_.Data(); }
		
		/** \return the number of elements in the block
		*/
		BASSO_IDTYPE NumElements ( ) const { return NumCols(); }
		/** \return the number of nodes in each element of the block
		*/
		BASSO_IDTYPE NNE ( ) const { return NumRows(); }
		
		/**\return the part id for the block */
		int PartID( ) const { return partid; }
		
		/** Adds all the node ids that support the elements in the block.  If the elements have been 
		renumberd to local node id values this will not work correctly.  The set is not initially
		cleared.
		\param gnidset, on return has the global node ids that support the elements in the block.
		*/
		virtual BASSO_IDTYPE AddSupportNIDs( set<BASSO_IDTYPE> &gnidset ) const
		{
			if ( localnum )
			{
				Basso_Warning("Basso_ElementBlock::AddSupportNIDs","elements have been renumbered to be local");
				return 0;
			}
			const BASSO_IDTYPE *nptr =  this->Data();
			for ( BASSO_IDTYPE i=0; i<this->Length(); ++i, ++nptr )
				gnidset.insert(*nptr);
			return gnidset.size();
		}
		
		/**\return true if the element conectivity has been set to be local node ids, else false.*/
		bool NumberingLocal( ) const { return localnum; )
		
		/**Renumbers the connectivity to local node ids
		\param gnids an array of global node ids.  The order will set up the maping. It is advised
		to use assending order, but this is not required.
		*/
		void RenumberLocal( const Basso_Array<BASSO_IDTYPE> &gnids );	
		
		/**Renumbers the connectivity to local node ids
		\param g2lmap a stl map of the global nid -> local nids/
		*/
		void RenumberLocal( const map<BASSO_IDTYPE,BASSO_IDTYPE> &g2lmap );
		
	void Print( std::ostream &out=BASSO_STDOUT ) const;
	
	protected:
		void SetBasis( Basso_ElementType type );
		void SetGIDs2LIDs( )
		{
			for ( int e=0; e<NumElements(); ++e )
				geids_[e]=e;
		}
		
	protected:
		Basso_ParentElement *elemptr_;
		Basso_Array<BASSO_IDTYPE> geids_;
		int partid;
		bool localnum;
		
};

void Basso_ElementBlock::RenumberLocal( const Basso_Array<BASSO_IDTYPE> &gnids )
{
	if ( localnum ) return;
	map<BASSO_IDTYPE,BASSO_IDTYPE> g2lmap;
	for ( BASSO_IDTYPE lnid=0; lnid<gnids.Length(); ++lnid )
		g2lmap[ gnids[lnid] ] = lnid;
	 localnum=true;
}

void Basso_ElementBlock::RenumberLocal( const map<BASSO_IDTYPE,BASSO_IDTYPE> &g2lmap )
{
	if ( localnum ) return;
	BASSO_IDTYPE *cptr = Data();
	map<BASSO_IDTYPE,BASSO_IDTYPE>::const_iterator mitr;
	for ( BASSO_IDTYPE e=0; e<NumElements(); ++e )
		for ( int i=0; i<NNE(); ++i, ++cptr )
		{
			mitr = g2lmap.find(*cptr);
			if ( mitr!=g2lmap.end() )
				*cptr = mitr->second;
			else
				Basso_Warning("Basso_ElementBlock::RenumberLocal","Map not surjective");
		}
		
	 localnum=true;
}
		
void Basso_ElementBlock::Print( std::ostream &out ) const 
{
	out << "Basso_ElementBlock\n"
	    << "\tnumber of elements : " << NumElements()
	    << "\n\telement type       : " << *elemptr_  
		<< "\n\tpart id            : " << partid;
		
	out << "\nleid   geid   connectivity ";
	if ( localnum ) 
		out << "(local node numbering)\n";
	else
		out << "(global node numbering)\n";
	for ( int e=0; e<NumElements(); ++e )
	{
		cout << setw(7) << e << setw(7) << geids_[e];
		for ( int i=0; i<NNE(); ++i )
			cout << setw(7) << (*this)[e][i];
		cout << "\n";
	}
}

std::ostream &operator << ( std::ostream &out, const Basso_ElementBlock &A )
{
	A.Print( out );
	return out;
}
	
void Basso_ElementBlock::SetBasis( Basso_ElementType type )
{
	if (elemptr_!= NULL) delete elemptr_;
	
	switch ( type )
	{
		case Basso_POINT1:
			elemptr_ = new Basso_Point1( );
			break;
		case Basso_LINE2:
			elemptr_ = new Basso_Line2( );
			break;
		case Basso_TRIA3:
			elemptr_ = new Basso_Tria3( );
			break;
		case Basso_QUAD4:
			elemptr_ = new Basso_Quad4( );
			break;
		case Basso_TETRA4:
			elemptr_ = new Basso_Tetra4( );
			break;
		case Basso_HEXA8:
			elemptr_ = new Basso_Hexa8( );
			break;
		case  Basso_LINE3:
			elemptr_ = new Basso_Line3( );
			break;
		case Basso_TRIA6:
			elemptr_ = new Basso_Tria6( );
			break;
		case Basso_TETRA10:
			elemptr_ = new Basso_Tetra10( );
			break;
		case Basso_PYRAMID5:
			elemptr_ = new Basso_Pyramid5( );
			break;
		case Basso_QUAD8:
			elemptr_ = new Basso_Quad8( );
			break;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////


/**
 *  \brief iterator for Basso_ElementBlockList
 */
 class Basso_ElementBlockListIterator 
 {
	public:
	
		Basso_ElementBlockListIterator(  ) { }
				
		Basso_ElementBlockListIterator( list< Basso_ElementBlock * >::iterator &a )
				: litr_(a) { }
				
		Basso_ElementBlockListIterator &operator = ( const Basso_ElementBlockListIterator &a )
			{ litr_ = a.litr_; return *this; }

		Basso_ElementBlockListIterator &operator++ ( ) { ++litr_; return *this; }
		Basso_ElementBlockListIterator &operator-- ( ) { --litr_; return *this; }
	
		Basso_ElementBlock &operator* ( )  { return **litr_; }
		Basso_ElementBlock *operator-> ( ) { return *litr_; }
		
		bool operator != ( const Basso_ElementBlockListIterator &itr )
			{ return ( litr_ == itr.litr_  ? false : true ); }
	
		bool operator == ( const Basso_ElementBlockListIterator &itr ) 
			{ return ( litr_ == itr.litr_  ? true : false ); }
	
	protected:
		list< Basso_ElementBlock * >::iterator litr_;
		
 };
 
/**
 *  \brief constant iterator for Basso_ElementBlockList
 */
 class Basso_ElementBlockListConstIterator 
 {
	public:
	
		Basso_ElementBlockListConstIterator(  ) { }
				
		Basso_ElementBlockListConstIterator( list< Basso_ElementBlock * >::const_iterator &a )
				: litr_(a) { }
				
		Basso_ElementBlockListConstIterator &operator = ( const Basso_ElementBlockListConstIterator &a )
			{ litr_ = a.litr_; return *this; }

		Basso_ElementBlockListConstIterator &operator++ ( ) { ++litr_; return *this; }
		Basso_ElementBlockListConstIterator &operator-- ( ) { --litr_; return *this; }
	
		const Basso_ElementBlock &operator* ( )  { return **litr_; }
		const Basso_ElementBlock *operator-> ( ) { return *litr_; }
		
		bool operator != ( const Basso_ElementBlockListConstIterator &itr )
			{ return ( litr_ == itr.litr_  ? false : true ); }
	
		bool operator == ( const Basso_ElementBlockListConstIterator &itr ) 
			{ return ( litr_ == itr.litr_  ? true : false ); }
	
	protected:
		list< Basso_ElementBlock * >::const_iterator litr_;
		
 };
 
/**
 * \brief List of element blocks 
 */
class Basso_ElementBlockList
{
	public:
		typedef Basso_ElementBlockListIterator iterator;
		typedef Basso_ElementBlockListConstIterator const_iterator;
	
		/** void constructor */
		Basso_ElementBlockList( ) { }
	
		/** Deconstructor*/
		~Basso_ElementBlockList( ) 
		{ 
			list< Basso_ElementBlock *>::iterator itr;
			for ( itr=elemblocks_.begin(); itr!=elemblocks_.end(); ++itr )
				delete *itr;
		}
		
		/** Adds a block of elements to the list
		\param type the element type
		\param ne the number of elements
		\param an identifying part id number (default=0)
		*/
		void AddElementBlock( Basso_ElementType type, BASSO_IDTYPE ne, int pid=0 ) 
		{ 
			elemblocks_.push_back( new Basso_ElementBlock(type,ne,pid) ); 
		}
		
		/** Adds a block of elements to the list (this will probably be removed)
		\param type the element type id
		\param ne the number of elements
		\param an identifying part id number (default=0)
		*/
		void AddElementBlock( unsigned int type, BASSO_IDTYPE ne, int pid=0 ) 
		{ 
			Basso_ElementType et = static_cast<Basso_ElementType>(type);
			AddElementBlock( et, ne, pid ); 
		}
		
		//void RenumberLocal( const Basso_Array<BASSO_IDTYPE> &gnids );
		
		/** 
		 *  Iterators
		 */
		 /**\return an iterator to the begining of the list of element blocks */
		iterator Begin( ) 
		{
			list< Basso_ElementBlock * >::iterator itr = elemblocks_.begin();
			return Basso_ElementBlockListIterator( itr ); 
		}
		
		 /**\return a constant iterator to the begining of the list of element blocks */
		const_iterator Begin( ) 
		{
			list< Basso_ElementBlock * >::const_iterator itr = elemblocks_.begin();=
			return Basso_ElementBlockListConstIterator( itr ); 
		}
		
		 /**\return an iterator to just past the end of the list of element blocks */
		iterator End( ) 
		{
			list< Basso_ElementBlock * >::iterator itr = elemblocks_.end();
			return Basso_ElementBlockListIterator( itr ); 
		}
		/**\return a constant iterator to just past the end of the list of element blocks */
		const_iterator End( ) 
		{
			list< Basso_ElementBlock * >::const_iterator itr = elemblocks_.end();
			return Basso_ElementBlockListConstIterator( itr ); 
		}
		
		/** Adds all the node ids that support the elements in the block.  If the elements have been 
		renumberd to local node id values this will not work correctly.  The set is not initially
		cleared.
		\param gnidset, on return has the global node ids that support the elements in the block.
		*/
		BASSO_IDTYPE AddSupportNIDs( set<BASSO_IDTYPE> &gnidset ) const
		{
			list< Basso_ElementBlock* >::const_iterator itr;
			for ( itr=elemblocks_.begin(); itr!=elemblocks_.end(); ++itr )
				(*itr)->AddSupportNIDs(gnidset);
			return gnidset.size();
		}
		
		/**Renumbers the connectivity to local node ids.  This assumes that the only local nodes are those
		in support of the elements in this block list.
		*/
		void RenumberLocal( );

		/**\return the number of blocks in the block list*/
		int NumBlocks( ) const { return elemblocks_.size(); }
		
		void Print( std::ostream &out=BASSO_STDOUT ) const;
		
	protected:
		list< Basso_ElementBlock * > elemblocks_;
};	

void Basso_ElementBlockList::RenumberLocal(  )
{
	set<BASSO_IDTYPE> gnidset;
	AddSupportNIDs(gnidset);
	
 	map<BASSO_IDTYPE,BASSO_IDTYPE> g2lmap;
	set<BASSO_IDTYPE>::const_iterator sitr;
	BASSO_IDTYPE lnid=0;
	for ( sitr=gnidset.begin(); sitr!=gnidset.end(); ++sitr, ++lnid )
		g2lmap[*sitr]=lnid;
	gnidset.clear();
	
	list< Basso_ElementBlock * >::iterator bitr;
	for ( bitr=elemblocks_.begin(); bitr!=elemblocks_.end(); ++bitr )
		(*bitr)->RenumberLocal(g2lmap);
	
}

void Basso_ElementBlockList::Print( std::ostream &out ) const 
{
	out << "Basso_ElementBlockList : "
	    << "number of blocks=" << NumBlocks() << "\n";
	
	list< Basso_ElementBlock *>::const_iterator itr;	
	int b=0;
	for ( itr=elemblocks_.begin(); itr!=elemblocks_.end(); ++itr, ++b )
		cout << "Block " << b << ": " << **itr;
}

std::ostream &operator << ( std::ostream &out, const Basso_ElementBlockList &A )
{
	A.Print( out );
	return out;
}
 
} // end of namspace
#endif