/*! \file Basso_defs.h 

\brief basic Basso definitions 

Basso ver 1.0

	\author Jack Chessa, jfchessa@utep.edu
	\date Wed Apr 6 2007

*/

#ifndef _BASSO_DEFS_H_
#define _BASSO_DEFS_H_

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdlib.h>

namespace Basso 
{

using namespace std;


/** A define for standard output stream to be used in Basso */
#ifndef BASSO_STDOUT
#define BASSO_STDOUT std::cout
#endif

/** A define for error and warning output stream to be used in Basso */
#ifndef BASSO_ERROUT
#define BASSO_ERROUT std::cout
#endif

/** Defines the precision for numerics (float or double) */
#ifndef BASSO_NUMERIC_TYPE
#define BASSO_NUMERIC_TYPE double
#endif

/** defines the type used for node ids, element ids, dof indicies etc any
long unsigned counting type of index 
This is the type to hold array dimensions for large arrays as well 
example on my Mac
int (-2,147,483,647 to 2,147,483,647)
unsigned int (0 to 4,294,967,295) 
int long (-6e18 to 6e18)
unsigned int long (0 to 18e8)
long long  <- for trilinos
*/
#ifndef BASSO_IDTYPE
#define BASSO_IDTYPE int
#endif


#ifndef BASSO_ARRAY_INDEX
#define BASSO_ARRAY_INDEX unsigned int
#endif

/** value to be assigned to all unallocated data structure pointsrs */
#ifndef BASSO_NULL_PTR
#define BASSO_NULL_PTR NULL
#endif

/** 
*/


/** Explicitly check the bounds when assessing arrays */
#define BASSO_BOUNDS_CHECK
#define ALLOW_DYNAMIC_RESIZE 

#define BASSO_MAX_NODES_PER_ELEMENT 64

/** \brief Handles simple error messages */
class Basso_Error
{
private:

public:
	Basso_Error( string source, string mesg )
		{ BASSO_ERROUT << "ERROR: in " << source << " " << mesg << "\n";exit(1);}
};

/** \brief Handles simple warning messages */
class Basso_Warning
{
public:
	Basso_Warning( string source, string mesg )
		{ BASSO_ERROUT << "WARNING: in " << source << " " << mesg << "\n"; }
};


/**
Trims the leading and trailing spaces from a std::string 
\param pString - the string to be trimmed
\param pWhitespace - the whitespace to be trimmed,  the default is " "
*/
const std::string string_trim(const std::string& pString,
                       const std::string& pWhitespace = " \t")
{
    const size_t beginStr = pString.find_first_not_of(pWhitespace);
    if (beginStr == std::string::npos)
    {
        // no content
        return "";
    }

    const size_t endStr = pString.find_last_not_of(pWhitespace);
    const size_t range = endStr - beginStr + 1;

    return pString.substr(beginStr, range);
}


} // end of namespace

#endif


