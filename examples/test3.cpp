#include <iostream>

#include "Basso_defs.h"
#include "Basso_QuadratureRule.h"
//#include "Basso_ParentElement.h"

using namespace std;
using namespace Basso;

int main(int argc, char* argv[])
{
	cout << "Test of the Element classes in Basso\n";
	
	Basso_QuadratureRule qrule(Basso_GAUSS1D_QUAD,3);
	
	Basso_QuadratureRule::ConstIterator qitr;
	for ( qitr=qrule.Begin(); qitr!=qrule.End(); ++qitr )
			cout << qitr->Pt() << qitr->Wt() << "\n";
	
	return 0;
}