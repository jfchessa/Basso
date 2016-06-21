#ifndef _BASSO_FACE_BC_SET_H_
#define _BASSO_FACE_BC_SET_H_

#include <set>
#include <iostream>
#include <iomanip>

#include "Basso_defs.h"
#include "Basso_Array.h"
#include "Basso_Array2D.h"
#include "Basso_iVector.h"
#include "Basso_nVector.h"
#include "Basso_Numeric.h"
#include "Basso_ParentElement.h"

namespace Basso
{
	
class Basso_FaceBCSet
{
public:
	Basso_FaceBCSet( const Basso_Array2D<BASSO_IDTYPE> &fconn, const Basso_ParentElement *eptr,
		 const Basso_nVector &trac, const Basso_iVector &ldofs );
	
	template < class DoFmAp, class RhS >
	void AddForce( const Basso_nMatrix &nodes, const DoFmAp &dofmap, RhS &f ) const;
	
	void Print( std::ostream &out=BASSO_STDOUT ) const;
	
protected:
	Basso_Array2D<BASSO_IDTYPE> fconn_;
	const Basso_ParentElement *eptr_;
	Basso_iVector ldofs_;
	Basso_nVector valv_;

};	

void Basso_FaceBCSet::Print( std::ostream &out ) const
{
	out << "Basso_FaceBC\n";
}

Basso_FaceBCSet::Basso_FaceBCSet( const Basso_Array2D<BASSO_IDTYPE> &fconn, const Basso_ParentElement *eptr,
		 const Basso_nVector &trac, const Basso_iVector &ldofs )
{
	fconn_ = fconn;
	valv_ = trac;
	ldofs_ = ldofs;
	eptr_ = eptr;
}
	

template < class DoFmAp, class RhS >
void Basso_FaceBCSet::AddForce( const Basso_nMatrix &nodes, 
		const DoFmAp &dofmap, RhS &f ) const
{
	int nn=eptr_->NumNodes(), sdim=nodes.NumRows(), edim=eptr_->Dimension(),
		ndof=ldofs_.Length();
	
	Basso_nVector fe(nn*ndof);
	Basso_Array<BASSO_IDTYPE> sctr(nn*ndof);
	Basso_nVector Na(nn);
	Basso_nMatrix ecoord(sdim,nn);
	Basso_QuadratureRule qrule;
	eptr_->Quadrature(1,qrule);
	for ( int e=0; e<fconn_.NumCols(); ++e )
	{
		element_coordinates( nodes, fconn_[e], fconn_.M(), ecoord );
		fe.Zero();
		dofmap.SetScatter(fconn_[e],fconn_.M(),sctr,ldofs_); 
		Basso_QuadratureRule::ConstIterator qitr;
		for ( qitr=qrule.Begin(); qitr!=qrule.End(); ++qitr )
		{
			eptr_->Na( qitr->Pt(), Na );
			Basso_Numeric jac = eptr_->DetJacobian( qitr->Pt(), ecoord, sdim );
			for ( int I=0; I<nn; ++I )
				for ( int i=0; i<ndof; ++i )
					fe[ndof*I+i] += valv_[i]*Na[I]*jac*(qitr->Wt());
		}
		f.Scatter( fe, sctr );
	}
}


std::ostream &operator << ( std::ostream &out, const Basso_FaceBCSet &A )
{
	A.Print( out );
	return out;
}

} // of namespace
#endif
