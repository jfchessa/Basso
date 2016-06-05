#ifndef _BASSO_MESH_QUADRILATERAL_H_
#define _BASSO_MESH_QUADRILATERAL_H_

#include "Basso_MeshPrimative2D.h"

namespace Basso
{
	

// ****************************************************************************************** //
/** 
\brief VBC for meshes on quadrilaterals
This is a base classe for generating meshes on quadrilateral regions.  
The quadrilateral is defined by four corner nodes, and may be a plane in 3D (I think)
*/
template<class T>
class Basso_MeshQuadrilateral : public Basso_MeshPrimative2D<T>
{
public:
	typedef T value_type;

	/**
	An enumerated type to id teh four edges of the quadrilateral
	*/
    enum EdgeID { kRightEdge=0, kTopEdge, kLeftEdge, kBottomEdge  };
                    
	Basso_MeshQuadrilateral() : Basso_MeshPrimative2D<T>() { DefaultParameters( ); }
	
	/**
	Constructor
	\param x1 - the (x,y,z) coordinate of the first corner
	\param x2 - the (x,y,z) coordinate of the second corner
	\param x3 - the (x,y,z) coordinate of the third corner
	\param x4 - the (x,y,z) coordinate of the fourth corner
	*/
	Basso_MeshQuadrilateral( const T *x1, const T *x2, const T *x3, const T *x4 )
	    {  SetCorners(x1,x2,x3,x4); DefaultParameters( ); }
	
	/**
	Constructor that assumes a rectanglular region between conrer 1 and 3.
	\param x1 - the (x,y,z) coordinate of the first corner
	\param x3 - the (x,y,z) coordinate of the third corner
	*/		
	Basso_MeshQuadrilateral( const T *x1, const T *x3 ) {  SetCorners(x1,x3); DefaultParameters( ); }
	
	/**
	Sets the corner points.
	\param x1 - the (x,y,z) coordinate of the first corner
	\param x2 - the (x,y,z) coordinate of the second corner
	\param x3 - the (x,y,z) coordinate of the third corner
	\param x4 - the (x,y,z) coordinate of the fourth corner
	*/	
	void SetCorners( const T *x1, const T *x2, const T *x3, const T *x4 ); 
	
	/**
	Sets the corner pointsusing only two opposign conrers.
	This ssumes a rectanglular region between conrer 1 and 3.
	\param x1 - the (x,y,z) coordinate of the first corner
	\param x3 - the (x,y,z) coordinate of the third corner
	*/	
	void SetCorners( const T *x1, const T *x3 ); 
	
	/**Sets the number of nodes in the 1 and 2 directions */
	virtual void SetNumNodes( BASSO_IDTYPE n1, BASSO_IDTYPE n2 ) { nn1=n1; nn2=n2; }
	/** returns the number of nodes in the mesh */
    virtual BASSO_IDTYPE NumNodes() const { return nn1*nn2; }
    /** Returns the number of nodes in the 1-direction */
    virtual BASSO_IDTYPE NumNodes1() const { return nn1; }
    /** Returns the number of nodes in the 2-direction */
    virtual BASSO_IDTYPE NumNodes2() const { return nn2; }
	
	virtual void EdgeSegs( Basso_Array2D<BASSO_IDTYPE> &conn, int id ) const;
	virtual void EdgeNodes( set<BASSO_IDTYPE> &nids, int id ) const;
	
protected:	
	void DefaultParameters( )
	{
		bias1 = BiasNONE; bias2 = BiasNONE; 
		bf1 = bf2 = 1.0;
		nn1 = nn2 = 2; 
		order = 1;
	}	
    
protected:
	T v1[3], v2[3], v3[3], v4[3];
	MeshBiasType bias1, bias2;
	T bf1, bf2; 
	BASSO_IDTYPE nn1, nn2;
	int order;    
};
	
template<class T>
void Basso_MeshQuadrilateral<T>::EdgeSegs( Basso_Array2D<BASSO_IDTYPE> &conn, int id ) const 
{
    BASSO_IDTYPE cptrn[2];
    
    switch ( id )
    {
        case kRightEdge :  // right edge
            
		    conn.Resize( 2, nn2-1 );
		    cptrn[0]=nn1-1; cptrn[1]=2*nn1-1;	
		    gen_conn_1d( cptrn, 2, nn2-1, conn.Data(), conn.LDA(), nn1 );
            break;
            
        case kTopEdge :  // top edge
            
		    conn.Resize( 2, nn1-1 );
		    cptrn[0]=nn1*nn2-1; cptrn[1]=nn1*nn2-2;	
		    gen_conn_1d( cptrn, 2, nn1-1, conn.Data(), conn.LDA(), -1 );
            break;
            
         case kLeftEdge :  // left edge
            
		    conn.Resize( 2, nn2-1 );
		    cptrn[0]=(nn2-1)*nn1; cptrn[1]=(nn2-2)*nn1;	
		    gen_conn_1d( cptrn, 2, nn2-1, conn.Data(), conn.LDA(), -nn1 );
            break;
            
        case kBottomEdge :  // bottom edge
            
		    conn.Resize( 2, nn1-1 );
		    cptrn[0]=0; cptrn[1]=1;	
		    gen_conn_1d( cptrn, 2, nn1-1, conn.Data(), conn.LDA(), 1 );
            break;
                                 
		default:
		    cout << "Basso_MeshQuadrilateralQuad4<T>::EdgeSegs unknown edge id=" << id << "\n";  
    }
}

template<class T>
void Basso_MeshQuadrilateral<T>::EdgeNodes( set<BASSO_IDTYPE> &nids, int id ) const
{
	Basso_Array2D<BASSO_IDTYPE> conn;
	EdgeSegs(conn,id);
	array2set(conn.Data(),conn.Length(),nids);
}
	
template<class T>
void Basso_MeshQuadrilateral<T>::SetCorners( const T *x1, const T *x2, const T *x3, const T *x4 )
{
	v1[0]=x1[0]; v1[1]=x1[1]; v1[2]=x1[2];
	v2[0]=x2[0]; v2[1]=x2[1]; v2[2]=x2[2];
	v3[0]=x3[0]; v3[1]=x3[1]; v3[2]=x3[2];
	v4[0]=x4[0]; v4[1]=x4[1]; v4[2]=x4[2];	
}

template<class T>
void Basso_MeshQuadrilateral<T>::SetCorners( const T *x1, const T *x7 )
{
	v1[0]=x1[0]; v1[1]=x1[1]; v1[2]=0.0;
	v2[0]=x7[0]; v2[1]=x1[1]; v2[2]=0.0;
	v3[0]=x7[0]; v3[1]=x7[1]; v3[2]=0.0;
	v4[0]=x1[0]; v4[1]=x7[1]; v4[2]=0.0;
}


// ****************************************************************************************** //
/**
\brief Tria3 mesh on a quadrilateral region.
*/
template<class T>
class Basso_MeshQuadrilateralTria3 : public Basso_MeshQuadrilateral<T>
{
public:
                    
	Basso_MeshQuadrilateralTria3() : Basso_MeshQuadrilateral<T>() { this->DefaultParameters( ); }
	
	Basso_MeshQuadrilateralTria3( const T *x1, const T *x2, const T *x3, const T *x4 )
	    {  this->SetCorners(x1,x2,x3,x4); this->DefaultParameters( ); }
		
	Basso_MeshQuadrilateralTria3( const T *x1, const T *x7 ) {  SetCorners(x1,x7); this->DefaultParameters( ); }
	
    virtual BASSO_IDTYPE NumElements() const { return 2*(this->nn1-1)*(this->nn2-1); }
	
	virtual void SetMesh( Basso_Array2D<T> &node, Basso_Array2D<BASSO_IDTYPE> &conn ) const;
    virtual Basso_ElementType MeshElementType() const { return Basso_TRIA3; }
	
};

template<class T>
void Basso_MeshQuadrilateralTria3<T>::SetMesh( Basso_Array2D<T> &node, Basso_Array2D<BASSO_IDTYPE> &conn ) const
{
	node.Resize( 3, this->NumNodes() );
	node_array_2d( this->v1[0], this->v1[1], this->v1[2], this->v2[0], this->v2[1], this->v2[2], 
				this->v3[0], this->v3[1], this->v3[2], this->v4[0], this->v4[1], this->v4[2],
                this->nn1, this->nn2, node.Data(), node.LDA(),
				this->bias1, this->bf1, this->bias2, this->bf2 );
						
	conn.Resize( 3, NumElements() );
    BASSO_IDTYPE *ePtr=conn.Data();
	BASSO_IDTYPE cptrn[3];
	cptrn[0]=0;         cptrn[1]=1;         cptrn[2]=this->nn1;        
	gen_conn_2d( cptrn, 3, this->nn1-1, this->nn2-1, ePtr, conn.LDA() );
	
	cptrn[0]=this->nn1+1;     cptrn[1]=this->nn1;       cptrn[2]=1;
    ePtr = conn.Data() + (this->nn1-1)*(this->nn2-1)*conn.LDA();        
	gen_conn_2d( cptrn, 3, this->nn1-1, this->nn2-1, ePtr, conn.LDA() );
}	

// ****************************************************************************************** //
/**
\brief Quad4 mesh on a quadrilateral region.
*/
template<class T>
class Basso_MeshQuadrilateralQuad4 : public Basso_MeshQuadrilateral<T>
{
public:
	typedef T value_type;

    enum EdgeID { kRightEdge=0, kTopEdge, kLeftEdge, kBottomEdge  };
                    
	Basso_MeshQuadrilateralQuad4() : Basso_MeshQuadrilateral<T>() { DefaultParameters( ); }
	
	Basso_MeshQuadrilateralQuad4( const T *x1, const T *x2, const T *x3, const T *x4 )
	    {  SetCorners(x1,x2,x3,x4); DefaultParameters( ); }
		
	Basso_MeshQuadrilateralQuad4( const T *x1, const T *x7 ) {  SetCorners(x1,x7); DefaultParameters( ); }
	
    virtual BASSO_IDTYPE NumElements() const { return (nn1-1)*(nn2-1); }
	
	virtual void SetMesh( Basso_Array2D<T> &node, Basso_Array2D<BASSO_IDTYPE> &conn ) const;
    virtual Basso_ElementType MeshElementType() const { return Basso_QUAD4; }
	
protected:	
	void DefaultParameters( )
	{
		bias1 = BiasNONE; bias2 = BiasNONE; 
		bf1 = bf2 = 1.0;
		nn1 = nn2 = 2; 
		order = 1;
	}	
    
protected:
	T v1[3], v2[3], v3[3], v4[3];
	MeshBiasType bias1, bias2;
	T bf1, bf2; 
	BASSO_IDTYPE nn1, nn2;
	int order;
};
 
template<class T>
void Basso_MeshQuadrilateralQuad4<T>::SetMesh( Basso_Array2D<T> &node, Basso_Array2D<BASSO_IDTYPE> &conn ) const
{
	node.Resize( 3, this->NumNodes() );
	node_array_2d( v1[0], v1[1], v1[2],   v2[0], v2[1], v2[2], 
				v3[0], v3[1], v3[2],   v4[0], v4[1], v4[2],
                nn1, nn2, node.Data(), node.LDA(),
				bias1, bf1, bias2, bf2 );
						
	conn.Resize( 4, this->NumElements() );
	BASSO_IDTYPE cptrn[4];
	cptrn[0]=0;         cptrn[1]=1;         cptrn[2]=nn1+1;         cptrn[3]=nn1;
	gen_conn_2d( cptrn, 4, nn1-1, nn2-1, conn.Data(), conn.LDA() );
}
	
	
	
	
} // end of namespace

#endif



