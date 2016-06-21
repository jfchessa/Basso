/*! \file TBasso_FESolver.h

	\brief Defines the class TBasso_FESolver

This class requires Trilinos Epetra	
	
\author Jack Chessa, jfchessa@utep.edu
\date Thursday, June 8, 2016

*/

#ifndef _TRILINOS_BASSO_FE_SOLVER_H_
#define _TRILINOS_BASSO_FE_SOLVER_H_

// std includes

// mpi includes
#include "mpi.h"

// Trilinos includes
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_LinearProblem.h"
#include "AztecOO_config.h"
#include "AztecOO.h"

// basso includes
#include "Basso_defs.h"
#include "Basso_PointBCSet.h"

// tbasso includes
#include "TBasso_FECrsMatrix.h"
#include "TBasso_FEVector.h"

using namespace std;

namespace Basso
{

/** 
Solves the system K*d=f with essential boundary conditions using AztecOO 
*/

class TBasso_FESolver
{
	public:
	
		TBasso_FESolver( TBasso_FECrsMatrix &K, TBasso_FEVector &d, TBasso_FEVector &f,
			const TBasso_DOFMap &gdof, Basso_PointBCSet &spcs ) 
			:   Kmat(K.EpetraObject()), fvct(f.EpetraObject()), bcs(spcs), dofmap(gdof),
				problem(&K.EpetraObject(),&d.EpetraObject(),&f.EpetraObject()), solver(problem)
			{ 	
				solved_ = false;
				SetDefaultOptions( ); 
			}
		
		int Solve(  );
		
		BASSO_IDTYPE NumEquations() const { return Kmat.NumGlobalRows(); }
		
		void Results( std::ostream &out=BASSO_STDOUT ) const;  
		
		void Print( std::ostream &out=BASSO_STDOUT ) const;
	
	protected:
		
		int EnforceEssentialBCs(  );
		void SetDefaultOptions( );
	
	
	protected:
	
		Basso_PointBCSet &bcs;
		Epetra_FECrsMatrix &Kmat;
		Epetra_FEVector &fvct;
		const TBasso_DOFMap &dofmap;
		
		Epetra_LinearProblem problem;  
		AztecOO solver;  
		
		Basso_Numeric beta;
		long long max_iter;
		Basso_Numeric tol;
		
		bool solved_;
	
};

void TBasso_FESolver::Results( std::ostream &out ) const
{
	if ( !solved_ )
		out << "System not yes solved";
	
	else 
		std::cout << "Solution performed in " << solver.NumIters() 
			<< " iterations out of " << max_iter << " iterations\n"
			<< "Norm of the residual (scaled) is " << solver.ScaledResidual() 
			<< "\nSolved in " << solver.SolveTime() << "\n";

}

void TBasso_FESolver::Print( std::ostream &out ) const 
{
	out << "TBasso_FESolver\n";
	Results(out);
}

void TBasso_FESolver::SetDefaultOptions(  )
{
	solver.SetAztecOption(AZ_solver,AZ_cg);
	solver.SetAztecOption(AZ_precond,AZ_Jacobi);
	
	//solver.SetAztecOption(AZ_precond,AZ_dom_decomp);

	beta = 1.0e6;
	max_iter=2*NumEquations();
	tol=1.0e-6;
}

int TBasso_FESolver::EnforceEssentialBCs( )
{
	
	// apply the essential boundary conditions (penalty method)
	int neq=Kmat.NumGlobalRows();
	Basso_Numeric Kpen = Kmat.GlobalMaxNumEntries()*Kmat.NormInf()*beta, fii;

	Basso_PointBCSet::ConstIterator itr;
	for ( itr=bcs.Begin(); itr!=bcs.End(); ++itr )
	{
		BASSO_IDTYPE iglob = dofmap.GDOF( itr->LNID(), itr->LDOF() );
		Basso_Numeric ival = itr->Value();
		Kmat.ReplaceGlobalValues(1,&iglob,&Kpen);
		fii=ival*Kpen;
		fvct.ReplaceGlobalValues(1,&iglob,&fii);
	}
}
	
int TBasso_FESolver::Solve(  )
{
 	EnforceEssentialBCs();
	BASSO_IDTYPE neq=Kmat.NumGlobalRows();
	
	// now solve using AztecOO 
	int maxSolveAttempts = 100;
	int sstatus = solver.AdaptiveIterate(max_iter,maxSolveAttempts,tol);
	
	solved_ = true;
	
	return sstatus;
} 

std::ostream &operator << ( std::ostream &out, const TBasso_FESolver &A )
{
	A.Print( out );
	return out;
}

} // end of namspace
#endif