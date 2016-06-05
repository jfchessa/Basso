#ifndef _BASSO_LAPACK_WRAPPERS_H_
#define _BASSO_LAPACK_WRAPPERS_H_

extern "C" void sgetrf_( int *, int *, float *, int *, int *, int * );
extern "C" void dgetrf_( int *, int *, double *, int *, int *, int * );

extern "C" void sgetrs_( char *, int *, int *, float *, int *, int *, float *, int *, int * );
extern "C" void dgetrs_( char *, int *, int *, double *, int *, int *, double *, int *, int * );

extern "C" void dgesv_( const int *, const int *, double *, const int *, int *, double *, const int *, int * );
extern "C" void sgesv_( const int *, const int *, float *,  const int *, int *, float *, const int *, int * );

extern "C" void sgbsv_( const int *, const int *, const int *, const int *, float *, const int *, int *, float *, const int *, int * );
extern "C" void dgbsv_( const int *, const int *, const int *, const int *, double *, const int *, int *, double *, const int *, int * );

extern "C" void spbtrs_( const char *, const int *, const int *, const int *, float *, const int *, float *, const int *, int * );
extern "C" void dpbtrs_( const char *, const int *, const int *, const int *, double *, const int *, double *, const int *, int * );

extern "C" void spbsv_( const char *, const int *, const int *, const int *, float *, const int *, float *, const int *, int * );
extern "C" void dpbsv_( const char *, const int *, const int *, const int *, double *, const int *, double *, const int *, int * );

extern "C" void ssyev_( char *, char *, int *, float *, int *, float *, float *, int *, int * );
extern "C" void dsyev_( char *, char *, int *, double *, int *, double *, double *, int *, int * );


//  			     *****************   G B S V    *****************
/**      SUBROUTINE GBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
*
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGBSV computes the solution to a real system of linear equations
*  A * X = B, where A is a band matrix of order N with KL subdiagonals
*  and KU superdiagonals, and X and B are N-by-NRHS matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor A as A = L * U, where L is a product of permutation
*  and unit lower triangular matrices with KL subdiagonals, and U is
*  upper triangular with KL+KU superdiagonals.  The factored form of A
*  is then used to solve the system of equations A * X = B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, and the solution has not been computed.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  M = N = 6, KL = 2, KU = 1:
*
*  On entry:                       On exit:
*
*      *    *    *    +    +    +       *    *    *   u14  u25  u36
*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*
*  Array elements marked * are not used by the routine; elements marked
*  + need not be set on entry, but are required by the routine to store
*  elements of U because of fill-in resulting from the row interchanges.
*
*  =====================================================================
**/
inline int gbsv( int n, int kl, int ku, int nrhs, float *ab, int ldab, int *ipiv, float *b, int ldb )
{
	int info;
	sgbsv_( &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info );
	return info;
}

inline int gbsv( int n, int kl, int ku, int nrhs, double *ab, int ldab, int *ipiv, double *b, int ldb )
{
	int info;
	dgbsv_( &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info );
	return info;
}

//  				*****************   P B T R S    *****************
/*
*	Solves A*x=b for a symetric (Hermitian) banded matrix
*/
inline int pbtrs( const char uplo, const int n, const int kd, const int nrhs, float *ab, const int ldab, float *b, const int ldb )
{
	int info;
	spbtrs_( &uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info );
	return info;
}

inline int pbtrs( const char uplo, const int n, const int kd, const int nrhs, double *ab, const int ldab, double *b, const int ldb )
{
	int info;
	dpbtrs_( &uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info );
	return info;
}

//  				*****************   P B S V    *****************
inline int pbsv( const char uplo, const int n, const int kd, const int nrhs, float *ab, const int ldab, float *b, const int ldb )
{	
	int info;
	spbsv_( &uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info );
	return info;
}

inline int pbsv( const char uplo, const int n, const int kd, const int nrhs, double *ab, const int ldab, double *b, const int ldb )
{	
	int info;
	dpbsv_( &uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info );
	return info;
}

//  				*****************   G E T R F    *****************

inline int getrf( int m, int n, float *A, int lda, int *ipiv )
{
	int info;
	sgetrf_( &m, &n, A, &lda, ipiv, &info );
	return info;
}

inline int getrf( int m, int n, double *A, int lda, int *ipiv )
{
	int info;
	dgetrf_( &m, &n, A, &lda, ipiv, &info );
	return info;
}

//  				*****************   G E T R S    *****************
inline int getrs( char trans, int n, int nrhs, float *A, int lda, int *ipiv, float *b, int ldb  )
{
	int info;
	sgetrs_( &trans, &n, &nrhs, A, &lda, ipiv, b, &ldb, &info );
	return info;
}

inline int getrs( char trans, int n, int nrhs, double *A, int lda, int *ipiv, double *b, int ldb  )
{
	int info;
	dgetrs_( &trans, &n, &nrhs, A, &lda, ipiv, b, &ldb, &info );
	return info;
}

//  				*****************   G E S V    *****************
inline int gesv( const int n, const int nrhs, float *a, const int lda, int *ipiv, float *b, const int ldb )
{
	int info;
	sgesv_( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
	return info;
}

inline int gesv( const int n, const int nrhs, double *a, const int lda, int *ipiv, double *b, const int ldb )
{
	int info;
    dgesv_( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
	return info;
}


//  				*****************   S Y E V    *****************
int syev( char jobz, char uplo, int N, float *A, int lda, float *w, 
               float *work, int lwork )
{
	int info;
	ssyev_( &jobz, &uplo, &N, A, &lda, w, work, &lwork, &info );
	return info;
}

int syev( char jobz, char uplo, int N, double *A, int lda, double *w, 
               double *work, int lwork )
{
	int info;
	dsyev_( &jobz, &uplo, &N, A, &lda, w, work, &lwork, &info );
	return info;
}
#endif

