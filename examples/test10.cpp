#include <iostream>

#include "Basso_defs.h"
#include "Basso_GmshFile.h"

#include <list>
#include <set>
#include <map>
 
using namespace std;
using namespace Basso;

int main(int argc, char* argv[])
{
	
	Basso_GmshFile gmsh("plate.msh");
	gmsh.SetPartition(3);
	
	set<int> physids;
	gmsh.ReadPhysicalIDs(physids);
	cout << "Physical ids on the parttion " << physids << "\n";
	gmsh.SetElementPhysicalIDs(physids);
	
	Basso_ElementBlockList elements;
	gmsh.ReadElementBlocks(elements); 
	
	cout << elements;
	
	set<BASSO_IDTYPE> snids;
	elements.AddSupportNIDs(snids);
	cout << "Support nodes\n" << snids << "\n";
	
	Basso_nMatrix nodes;
	Basso_Array<BASSO_IDTYPE> gnids;
	
	gmsh.ReadNodes(nodes,gnids,snids);
	cout << nodes << "\n";
	
	return 0;
}