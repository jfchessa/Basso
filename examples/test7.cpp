// std includes
#include <iostream>
#include <set>

// basso includes
#include "Basso_defs.h"
#include "Basso_Numeric.h"
#include "Basso_Array.h"
#include "Basso_PointBCSet.h"
#include "Basso_FaceBCSet.h"

// tbasso includes

using namespace std;
using namespace Basso;

int main(int argc, char* argv[])
{
	cout << "Define some SPCs\n";
	
	Basso_Array<BASSO_IDTYPE> nids=string("0 2 4 6 10 12 14 16");
	Basso_PointBCSet spcset;
	spcset.AddPointBCs(nids,0,1.0);
	spcset.AddPointBCs(nids,1,0.0);
	cout << spcset << "\n";
	
	Basso_PointBCSet::ConstIterator sitr;
	for ( sitr=spcset.Begin(); sitr!=spcset.End(); ++sitr )
		cout << "fix node " << sitr->LNID() << " dof " << sitr->LDOF() << " to have a value of " << sitr->Value() << "\n";
	
	return 0;
}