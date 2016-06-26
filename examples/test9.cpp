#include <iostream>

#include "Basso_defs.h"

#include "Basso_VTK.h"
#include "Basso_Quad4.h"


#include <list>
 
using namespace std;
using namespace Basso;

int main(int argc, char* argv[])
{
	
	Basso_nMatrix coords("0.0 1.0 2.0 0.0 1.0 2.0;0.0 0.0 0.0 1.0 1.0 1.0;0.0 0.0 0.0 0.0 0.0 0.0");
	Basso_Array2D<BASSO_IDTYPE> conn("0 1; 1 2; 4 5; 3 4");
	Basso_Quad4 quad4;
	
	Basso_nVector temperature("100 110 120 120 140 200");

	Basso_VTKUnstructuredGrid results;
	results.SetMesh(coords,conn,&quad4);
	Basso_VTKFloatDataArray tdata("temperature",temperature.Data(),temperature.Length());
	results.AddPointData(&tdata);
	
	cout << results << "\n";
	
	results.WriteFile();
	
	
	return 0;
}