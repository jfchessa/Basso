#ifndef _BASSO_ELEMENT_TYPE_H_
#define _BASSO_ELEMENT_TYPE_H_

namespace Basso
{



enum Basso_ElementType 
{   
    Basso_NONE=0, Basso_LINE2=1, Basso_TRIA3, Basso_QUAD4, Basso_TETRA4, Basso_HEXA8, Basso_PRISM6, Basso_PYRAMID5, 
    Basso_LINE3, Basso_TRIA6, Basso_QUAD9, Basso_TETRA10, Basso_HEXA27, Basso_PRISM18, Basso_PYRAMID14,
    Basso_POINT1, Basso_QUAD8, Basso_HEXA20, Basso_PRISM15, Basso_PYRAMID13, 
    Basso_TRIA9, Basso_TRIA10, Basso_TRIA12INC, Basso_TRIA15, Basso_TRIA15INC, Basso_TRIA21,
    Basso_LINE4, Basso_LINE5, Basso_LINE6,
    Basso_TETRA20, Basso_TETRA35, Basso_TETRA56 
};


/*
int (*)( NuMeRiC *Na, const NuMeRiC xi )

void shape_function_pointer( int &etype )
{
    switch (etype)
    {
        case 1:
        break;
        
          
    }
}
*/

} // end namespace

#endif


