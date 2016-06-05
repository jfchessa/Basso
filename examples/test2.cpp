#include <iostream>

#include "Basso_defs.h"
#include "Basso_Point.h"

using namespace std;
using namespace Basso;

int main(int argc, char* argv[])
{
	cout << "Test of the Element classes in Basso\n";
	
	Basso_Numeric xref[3] = { 3.1, 4.5, -2.0 };
	
	Basso_Point p1(1,3,-4);
	cout << p1 << "\n"; 
	
	Basso_Point p2(xref);
	cout << p2 << "\n"; 
	
	xref[0] = -3.14156;
	cout << p2 << "\n"; 
	
	p1 = p2;
	
	cout << p1 << "\n";
	p1.x(1) = 6.789;
	
	cout << p1 << "\n" << p2 << "\n";
	
	
	return 0;
}