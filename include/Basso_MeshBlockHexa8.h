#ifndef _BASSO_MESH_BLOCK_HEXA8_H_
#define _BASSO_MESH_BLOCK_HEXA8_H_

#include "Basso_MeshPrimative3D.h"

namespace Basso
{
	
// ****************************************************************************************** //
/**
\brief Hexa8 mesh on a 3D block region.
*/
template<class T>
class Basso_MeshBlockHexa8 : public Basso_MeshPrimative3D<T>
{
public:
    
	typedef T value_type;
    enum FaceID { kNeg2Face, kPos1Face, kPos2Face, kNeg1Face, kNeg3Face, kPos3Face };
    enum EdgeID { kEdge01, kEdge12, kEdge23, kEdge30, kEdge45, kEdge56, kEdge67, kEdge70, 
                    kEdge04, kEdge15, kEdge26, kEdge37  };
    enum VertexID { kVertex0, kVertex1, kVertex2, kVertex3, kVertex4, kVertex5, kVertex6, kVertex7 };
    
	/** Define a block in 3D by its verticies */
	
	Basso_MeshBlockHexa8( );
		
	Basso_MeshBlockHexa8( const T *x1, const T *x2, const T *x3, const T *x4,
		const T *x5, const T *x6, const T *x7, const T *x8 );
		
	Basso_MeshBlockHexa8( const T *x1, const T *x7 );
	
	void SetCorners( const T *x1, const T *x2, const T *x3, const T *x4,
		const T *x5, const T *x6, const T *x7, const T *x8 ); 
		
	void SetCorners( const T *x1, const T *x7 ); 

	virtual void SetElementSize( T he ); 
	virtual void SetNumNodes( BASSO_IDTYPE n1, BASSO_IDTYPE n2, BASSO_IDTYPE n3 ) { nn1=n1; nn2=n2; nn3=n3; }

	virtual BASSO_IDTYPE NumNodes() const { return nn1*nn2*nn3; }
	virtual int NumNodesPerElement() const { return 8; }
	virtual BASSO_IDTYPE NumElements() const { return (nn1-1)*(nn2-1)*(nn3-1); }
	
	virtual void SetMesh( Basso_Array2D<T> &node, Basso_Array2D<BASSO_IDTYPE> &conn ) const;
	virtual Basso_ElementType MeshElementType() const { return Basso_HEXA8; }
	
	virtual void SideSegs( Basso_Array2D<BASSO_IDTYPE> &conn, int id ) const;
	virtual void SideNodes( set<BASSO_IDTYPE> &nids, int id ) const;
	virtual Basso_ElementType SideElementType() const { return Basso_QUAD4; }
	
	virtual void EdgeSegs( Basso_Array2D<BASSO_IDTYPE> &conn, int id ) const;
	virtual void EdgeNodes( set<BASSO_IDTYPE> &nids, int id ) const;
	virtual Basso_ElementType EdgeElementType() const { return Basso_LINE2; }
	
    virtual BASSO_IDTYPE VertexNodeID( int id ) const;

protected:
	void DefaultParameters( )
	{
		bias1 = bias2 = bias3 = BiasNONE; 
		bf1 = bf2 = bf3 = 1.0;
		nn1 = nn2 = nn3 = 2; 
		order = 1;
	}
	
protected:
	T v1[3], v2[3], v3[3], v4[3], v5[3], v6[3], v7[3], v8[3];
	MeshBiasType bias1, bias2, bias3;
	T bf1, bf2, bf3; 
	BASSO_IDTYPE nn1, nn2, nn3;
	int order;
};


template<class T>
void Basso_MeshBlockHexa8<T>::SetElementSize( T he ) 
{ 
	std::cout << "Basso_MeshBlockHexa8<T>::SetElementSize( T he ) not yet impleementd";
}
 
template<class T>
void Basso_MeshBlockHexa8<T>::SetCorners( const T *x1, const T *x2, const T *x3, const T *x4,
		const T *x5, const T *x6, const T *x7, const T *x8 )
{
	v1[0]=x1[0]; v1[1]=x1[1]; v1[2]=x1[2];
	v2[0]=x2[0]; v2[1]=x2[1]; v2[2]=x2[2];
	v3[0]=x3[0]; v3[1]=x3[1]; v3[2]=x3[2];
	v4[0]=x4[0]; v4[1]=x4[1]; v4[2]=x4[2];
	v5[0]=x5[0]; v5[1]=x5[1]; v5[2]=x5[2];
	v6[0]=x6[0]; v6[1]=x6[1]; v6[2]=x6[2];
	v7[0]=x7[0]; v7[1]=x7[1]; v7[2]=x7[2];
	v8[0]=x8[0]; v8[1]=x8[1]; v8[2]=x8[2];	
}

template<class T>
void Basso_MeshBlockHexa8<T>::SetCorners( const T *x1, const T *x7 )
{
	v1[0]=x1[0]; v1[1]=x1[1]; v1[2]=x1[2];
	v2[0]=x7[0]; v2[1]=x1[1]; v2[2]=x1[2];
	v3[0]=x7[0]; v3[1]=x7[1]; v3[2]=x1[2];
	v4[0]=x1[0]; v4[1]=x7[1]; v4[2]=x1[2];	
	v5[0]=x1[0]; v5[1]=x1[1]; v5[2]=x7[2];
	v6[0]=x7[0]; v6[1]=x1[1]; v6[2]=x7[2];
	v7[0]=x7[0]; v7[1]=x7[1]; v7[2]=x7[2];
	v8[0]=x1[0]; v8[1]=x7[1]; v8[2]=x7[2];	
}

template<class T>
Basso_MeshBlockHexa8<T>::Basso_MeshBlockHexa8( ) : Basso_MeshPrimative3D<T>() 
{
	v1[0]=0.0; v1[1]=0.0; v1[2]=0.0;
	v2[0]=0.0; v2[1]=0.0; v2[2]=0.0;
	v3[0]=0.0; v3[1]=0.0; v3[2]=0.0;
	v4[0]=0.0; v4[1]=0.0; v4[2]=0.0;
	v5[0]=0.0; v5[1]=0.0; v5[2]=0.0;
	v6[0]=0.0; v6[1]=0.0; v6[2]=0.0;
	v7[0]=0.0; v7[1]=0.0; v7[2]=0.0;
	v8[0]=0.0; v8[1]=0.0; v8[2]=0.0;
	
	DefaultParameters( );
}
 

template<class T>
Basso_MeshBlockHexa8<T>::Basso_MeshBlockHexa8( const T *x1, const T *x2, const T *x3, const T *x4,
	const T *x5, const T *x6, const T *x7, const T *x8 ) : Basso_MeshPrimative3D<T>() 
{
	SetCorners(x1,x2,x3,x4,x5,x6,x7,x8);
	
	DefaultParameters( );
}
 
template<class T>
Basso_MeshBlockHexa8<T>::Basso_MeshBlockHexa8( const T *x1, const T *x7 ) : Basso_MeshPrimative3D<T>() 
{
	SetCorners(x1,x7);
	
	DefaultParameters( );
}

template<class T>
void Basso_MeshBlockHexa8<T>::SetMesh( Basso_Array2D<T> &node, Basso_Array2D<BASSO_IDTYPE> &conn ) const
{
	node.Resize( 3, NumNodes() );
	node_array_3d( v1[0], v1[1], v1[2],   v2[0], v2[1], v2[2], 
				v3[0], v3[1], v3[2],   v4[0], v4[1], v4[2],
				v4[0], v5[1], v5[2],   v6[0], v6[1], v6[2],
				v7[0], v7[1], v7[2],   v8[0], v8[1], v8[2],
                nn1, nn2, nn3, node.Data(), node.LDA(),
				bias1, bf1, bias2, bf2, bias3, bf3 );
				
	conn.Resize( 8, NumElements() );
	BASSO_IDTYPE cptrn[8];
	cptrn[0]=0;         cptrn[1]=1;         cptrn[2]=nn1+1;         cptrn[3]=nn1;
	cptrn[4]=0+nn1*nn2; cptrn[5]=1+nn1*nn2; cptrn[6]=nn1+1+nn1*nn2; cptrn[7]=nn1+nn1*nn2;
	gen_conn_3d( cptrn, 8, nn1-1, nn2-1, nn3-1, conn.Data(), conn.LDA() );
}

template<class T>
void Basso_MeshBlockHexa8<T>::SideSegs( Basso_Array2D<BASSO_IDTYPE> &conn, int id )  const
{
	BASSO_IDTYPE cptrn[4];
	
	switch ( id )
	{
		case 0:  // -2 normal face 
		
		conn.Resize(4, (nn1-1)*(nn3-1) );
		cptrn[0]=0; cptrn[1]=1; cptrn[2]=nn1*nn2+1; cptrn[3]=nn1*nn2;	
		gen_conn_2d( cptrn, 4, nn1-1, nn3-1, 
			conn.Data(), conn.LDA(), 1, nn1*(nn2-1)+2 );
		break;
		
		case 1:  // +1 normal face 
		
		conn.Resize(4, (nn2-1)*(nn3-1) );
		cptrn[0]=nn1-1; cptrn[1]=cptrn[0]+nn1; 
		cptrn[2]=cptrn[1]+nn1*nn2; cptrn[3]=cptrn[2]-nn1;	
		gen_conn_2d( cptrn, 4, nn2-1, nn3-1, 
			conn.Data(), conn.LDA(), nn1, 2*nn1 );
		break;
		
		case 2:  // +2 normal face 
		
		conn.Resize(4, (nn1-1)*(nn3-1) );
		cptrn[0]=nn1*(nn2-1)+1; cptrn[1]=nn1*(nn2-1); 
		cptrn[2]=nn1*(nn2-1)+nn1*nn2; cptrn[3]=cptrn[2]+1;	
		gen_conn_2d( cptrn, 4, nn1-1, nn3-1, 
			conn.Data(), conn.LDA(), 1, nn1*(nn2-1)+2 );
		break;
		
		case 3:  // -1 normal face 
		
		conn.Resize(4, (nn2-1)*(nn3-1) );
		cptrn[0]=nn1; cptrn[1]=0; 
		cptrn[2]=nn1*nn2; cptrn[3]=cptrn[2]+nn1;	
		gen_conn_2d( cptrn, 4, nn2-1, nn3-1, 
			conn.Data(), conn.LDA(), nn1, 2*nn1 );
		break;
		
		case 4:  // -3 normal face 
		
		conn.Resize(4, (nn2-1)*(nn1-1) );
		cptrn[0]=1; cptrn[1]=0; 
		cptrn[2]=nn1; cptrn[3]=nn1+1;	
		gen_conn_2d( cptrn, 4, nn1-1, nn2-1, 
			conn.Data(), conn.LDA(), 1, 2 );
		break;
		
		case 5:  // +3 normal face 
		
		conn.Resize(4, (nn2-1)*(nn1-1) );
		cptrn[0]=(nn3-1)*nn1*nn2; cptrn[1]=cptrn[0]+1; 
		cptrn[2]=cptrn[1]+nn1; cptrn[3]=cptrn[2]-1;	
		gen_conn_2d( cptrn, 4, nn1-1, nn2-1, 
			conn.Data(), conn.LDA(), 1, 2 );
		break;
		
		default:
		cout << "Basso_MeshBlockHexa8<T>::SideSegs unknown side id=" << id << "\n";
		
	} 
}

template<class T>
void Basso_MeshBlockHexa8<T>::SideNodes( set<BASSO_IDTYPE> &nidset, int id ) const
{
	Basso_Array2D<BASSO_IDTYPE> sideSeg;
	SideSegs( sideSeg, id );
	for ( BASSO_IDTYPE e=0; e<sideSeg.N(); ++e )
		for ( BASSO_IDTYPE i=0; i<sideSeg.M(); ++i )
 			nidset.insert( sideSeg(i,e) );
}

template<class T>
void Basso_MeshBlockHexa8<T>::EdgeSegs( Basso_Array2D<BASSO_IDTYPE> &conn, int id ) const
{
	BASSO_IDTYPE cptrn[2];
	
	switch ( id )
	{
		case 0:  // edge connecting vertex 0 to 1
		
		conn.Resize(2, nn1-1 );
		cptrn[0]=0; cptrn[1]=1;	
		gen_conn_1d( cptrn, 2, nn1-1, conn.Data(), conn.LDA(), 1 );
		break;
					
		case 1:  // edge connecting vertex 1 to 2
		
		conn.Resize(2, nn2-1 );
		cptrn[0]=nn1-1; cptrn[1]=cptrn[0]+nn1;	
		gen_conn_1d( cptrn, 2, nn2-1, conn.Data(), conn.LDA(), nn1 );
		break;
							
		case 2:  // edge connecting vertex 2 to 3
		
		conn.Resize(2, nn1-1 );
		cptrn[0]=nn1*nn2-1; cptrn[1]=cptrn[0]-1;	
		gen_conn_1d( cptrn, 2, nn1-1, conn.Data(), conn.LDA(), -1 );
		break;
			
		case 3:  // edge connecting vertex 3 to 0
		
		conn.Resize(2, nn2-1 );
		cptrn[0]=(nn2-1)*nn1; cptrn[1]=cptrn[0]-nn1;	
		gen_conn_1d( cptrn, 2, nn2-1, conn.Data(), conn.LDA(), -nn1 );
		break;
				
		case 4:  // edge connecting vertex 4 to 5
		
		conn.Resize(2, nn1-1 );
		cptrn[0]=(nn3-1)*nn1*nn2; cptrn[1]=cptrn[0]+1;	
		gen_conn_1d( cptrn, 2, nn1-1, conn.Data(), conn.LDA(), 1 );
		break;
					
		case 5:  // edge connecting vertex 5 to 6
		
		conn.Resize(2, nn2-1 );
		cptrn[0]=(nn3-1)*nn1*nn2+(nn1-1); cptrn[1]=cptrn[0]+nn1;	
		gen_conn_1d( cptrn, 2, nn2-1, conn.Data(), conn.LDA(), nn1 );
		break;
							
		case 6:  // edge connecting vertex 6 to 7
		
		conn.Resize(2, nn1-1 );
		cptrn[0]=nn1*nn2*nn3-1; cptrn[1]=cptrn[0]-1;	
		gen_conn_1d( cptrn, 2, nn1-1, conn.Data(), conn.LDA(), -1 );
		break;
			
		case 7:  // edge connecting vertex 7 to 4
		
		conn.Resize(2, nn2-1 );
		cptrn[0]=nn1*(nn2*nn3-1); cptrn[1]=cptrn[0]-nn1;	
		gen_conn_1d( cptrn, 2, nn2-1, conn.Data(), conn.LDA(), -nn1 );
		break;
						
		default:
		std:cout << "Basso_MeshBlockHexa8<T>::EdgeSegs unknown edge id=" << id << "\n";
		
	} 	
}

template<class T>
void Basso_MeshBlockHexa8<T>::EdgeNodes( set<BASSO_IDTYPE> &nidset, int id ) const
{
	Basso_Array2D<BASSO_IDTYPE> edgeSeg;
	EdgeSegs( edgeSeg, id );
	for ( BASSO_IDTYPE e=0; e<edgeSeg.N(); ++e )
		for ( BASSO_IDTYPE i=0; i<edgeSeg.M(); ++i )
 			nidset.insert( edgeSeg(i,e) );
}

template<class T>
BASSO_IDTYPE Basso_MeshBlockHexa8<T>::VertexNodeID( int id ) const
{
	switch ( id )
	{
	    case 0:
        return 0;
        
        case 1:
        return nn1-1;
        
        case 2:
        return nn1*nn2-1;
        
        case 3:
        return nn1*(nn2-1);
        
        case 4:
        return (nn3-1)*nn1*nn2;
        
        case 5:
        return (nn3-1)*nn1*nn2 + (nn1-1);
        
        case 6:
        return nn1*nn2*nn3-1;
        
        case 7:
        return nn1*(nn2*nn3-1);
	    
	   	default:
		std:cout << "Basso_MeshBlockHexa8<T>::VertexNodeID unknown vertex id=" << id << "\n";
    };
}
/****************************************************************/

} // end namespace
#endif

