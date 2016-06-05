#ifndef _BASSO_STLOPS_H_
#define _BASSO_STLOPS_H_

#include "Basso_defs.h"

#include <queue>
#include <list>
#include <set>
#include <map>

namespace Basso
{

using namespace std;

template<class StLTyPe>
void PrintArray1D( const StLTyPe *l, int n, ostream &os=cout )
{
	const StLTyPe *ptr = l;
	for ( int i=0; i<n; ++i, ++ptr )
		os << *ptr << " ";
}

template<class StLTyPe>
void PrintArray2D( const StLTyPe *l, int m, int n, int ldl=-1, ostream &os=cout )
{
    if ( ldl<0 ) ldl=n;
    
	const StLTyPe *ptr = l;
	for ( int i=0; i<m; ++i, ptr += ldl-n )
	{
		for ( int j=0; j<n; ++j, ++ptr )
			os << setw(10) << *ptr << " ";
		os << endl;
	}
}

template<class StLTyPe>
void PrintArray( const StLTyPe *l, int n, ostream &os=cout )  { PrintArray1D(l,n,os); }
	
template<class StLTyPe>
void PrintArray( const StLTyPe *l, int m, int n, int ldl, ostream &os=cout )  { PrintArray2D(l,m,n,ldl,os); }

template<class StLTyPe>
void PrintFortranArray2D( const StLTyPe *l, int m, int n, int ldl=-1, ostream &os=cout )
{
    if ( ldl<0 ) ldl=m;
    
	for ( int i=0; i<m; ++i )
	{
		for ( int j=0; j<n; ++j )
			os << setw(10) << l[i+j*ldl] << " ";
		os << endl;
	}
}

template<class StLTyPe>
void PrintFortranArray( const StLTyPe *l, int m, int n, int ldl=-1, ostream &os=cout )  { PrintFortranArray2D(l,m,n,ldl,os); }

template<class StLTyPe>
void PrintSTL1( const StLTyPe &l,  ostream &os=cout  )
{
	typename StLTyPe::const_iterator itr, eItr=l.end();
    --eItr;
	os << "{ ";
	for ( itr=l.begin(); itr!=eItr; ++itr )
		os << *itr << ", ";
	os << *itr << " }";
}

template<class LiStTyPe>
ostream& operator << (ostream& os, const list<LiStTyPe> &s)
{
    os << "list:";
    PrintSTL1(s,os);
	return os;
}

template<class SeTTyPe>
ostream& operator << (ostream& os, const set<SeTTyPe> &s)
{
    os << "set:";
	PrintSTL1(s,os);
	return os;
}

template<class VeCtTyPe>
ostream& operator << (ostream& os, const vector<VeCtTyPe> &s)
{
    os << "vector:";
	PrintSTL1(s,os);
	return os;
}

template<class QuETyPe>
ostream& operator << (ostream& os, const queue<QuETyPe> &s)
{
    os << "queue:";
	PrintSTL1(s,os);
	return os;
}

template<class TyPe1, class TyPe2>
ostream& operator << (ostream& os, const map<TyPe1,TyPe2> &l)
{
	typename map<TyPe1,TyPe2>::const_iterator itr=l.begin();
	os << "map:{ ";
	os << itr->first << "->" << itr->second << " ";
	for ( ++itr; itr!=l.end(); ++itr )
		os << ", " << itr->first << "->" << itr->second;
	os << "}";
	return os;
}

template<class StLTyPe>
void set2array( const set<StLTyPe> &s, StLTyPe *a )
{
    typename set< StLTyPe >::const_iterator sitr;
    StLTyPe *aptr=a;
    for ( sitr=s.begin(); sitr!=s.end(); ++sitr, ++aptr )
        *aptr = *sitr; 
}

/** Inserts the n elements in an array a into a set */
template<class StLTyPe>
void array2set( const StLTyPe *a, int n, set<StLTyPe> &s  )
{
    const StLTyPe *aptr=a;
	for ( int i=0; i<n; ++i, ++aptr )
		s.insert( *aptr ); 
}

template<class StLTyPe>
void list2array( const list<StLTyPe> &s, StLTyPe *a )
{
    typename std::list<StLTyPe>::const_iterator sitr;
    StLTyPe *aptr=a;
    for ( sitr=s.begin(); sitr!=s.end(); ++sitr, ++aptr )
        *aptr = *sitr;
}

template<class StLTyPe>
void vector2array( const vector<StLTyPe> &s, StLTyPe *a )
{
    typename std::vector<StLTyPe>::const_iterator sitr;
    StLTyPe *aptr=a;
    for ( sitr=s.begin(); sitr!=s.end(); ++sitr, ++aptr )
        *aptr = *sitr;
}

template< class T1, class T2 >
void map2array( const map<T1,T2> &s, T1 *a, T2 *b )
{
    typename std::map<T1,T2>::const_iterator sitr;
    T1 *aptr=a;
    T2 *bptr=b;
    for ( sitr=s.begin(); sitr!=s.end(); ++sitr, ++aptr, ++bptr )
    {
        *aptr = sitr->first;
        *bptr = sitr->second;
    }
}

} // end namespace
#endif

