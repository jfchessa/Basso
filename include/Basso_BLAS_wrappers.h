#ifndef _BASSO_BLAS_WRAPPERS_H_
#define _BASSO_BLAS_WRAPPERS_H_ 

#include <complex>
#define CoMpLeXtYpE std::complex<float>

//----  Level 1

extern "C" void sscal_( const int *, const float *, float *, const int * );

extern "C" void dscal_( const int *, const double *, double *, const int * );

extern "C" void zscal_( const int *, const CoMpLeXtYpE *, CoMpLeXtYpE *, const int * );

//--------

extern "C" float snrm2_( const int *, const float *, const int * );

extern "C" double dnrm2_( const int *, const double *, const int * );

extern "C" CoMpLeXtYpE znrm2_( const int *, const CoMpLeXtYpE *, const int * );

//--------

extern "C" float sdot_( const int *, const float *, const int *, const float *, const int * );

extern "C" double ddot_( const int *, const double *, const int *, const double *, const int * );

extern "C" CoMpLeXtYpE zdot_( const int *, const CoMpLeXtYpE *, const int *, const CoMpLeXtYpE *, const int * );

//--------

extern "C" void saxpy_( const int *, const float *, const float *, 
                        const int *, float *, const int * );
                        
extern "C" void daxpy_( const int *, const double *, const double *, 
                        const int *, double *, const int * );

extern "C" void zaxpy_( const int *, const CoMpLeXtYpE *, const CoMpLeXtYpE *, 
                        const int *, CoMpLeXtYpE *, const int * );
                        
//----  Level 2

extern "C" void sgemv_( const char *, const int *, const int *,  
                        const float *, const float *, const int *, 
                        const float *, const int *, 
                        const float *, float *, const int * );

extern "C" void dgemv_( const char *, const int *, const int *,  
                        const double *, const double *, const int *, 
                        const double *, const int *, 
                        const double *, double *, const int * );

extern "C" void zgemv_( const char *, const int *, const int *,  
                        const CoMpLeXtYpE *, const CoMpLeXtYpE *, const int *, 
                        const CoMpLeXtYpE *, const int *, 
                        const CoMpLeXtYpE *, CoMpLeXtYpE *, const int * );

extern "C" void ssbmv_( const char *, const int *, const int *, const float *, const float *, const int *,
						const float *, const int *, const float *, const float *, const int * );
						
extern "C" void dsbmv_( const char *, const int *, const int *, const double *, const double *, const int *,
						const double *, const int *, const double *, const double *, const int * );
						
extern "C" void zsbmv_( const char *, const int *, const int *, const CoMpLeXtYpE *, const CoMpLeXtYpE *, const int *,
						const CoMpLeXtYpE *, const int *, const CoMpLeXtYpE *, const CoMpLeXtYpE *, const int * );
						
//---- Level 3

extern "C" void sgemm_( const char *, const char *, const int *, const int *, const int *, 
                        const float *, const float *, const int *, 
                        const float *, const int *, 
                        const float *, float *, const int * );

extern "C" void dgemm_( const char *, const char *, const int *, const int *, const int *, 
                        const double *, const double *, const int *, 
                        const double *, const int *, 
                        const double *, double *, const int * );

extern "C" void zgemm_( const char *, const char *, const int *, const int *, const int *, 
        				const CoMpLeXtYpE *, const CoMpLeXtYpE *, const int *, 
        				const CoMpLeXtYpE *, const int *, 
        				const CoMpLeXtYpE *, CoMpLeXtYpE *, const int * );
        				
//  				*****************   S C A L   *****************
/**
      SUBROUTINE SSCAL(N,SA,SX,INCX)
*     .. Scalar Arguments ..
      REAL SA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      REAL SX(*)
*     ..
*
*  Purpose
*  =======
*
*     scales a vector by a constant.
*     uses unrolled loops for increment equal to 1.
*/
inline void scal( const int n, const float a, float *x, const int incx )
{  sscal_( &n, &a, x, &incx ); }

inline void scal( const int n, const double a, double *x, const int incx )
{  dscal_( &n, &a, x, &incx ); }

inline void scal( const int n, const CoMpLeXtYpE a, CoMpLeXtYpE *x, const int incx )
{  zscal_( &n, &a, x, &incx ); }


//  				****************   N O R M  *******************
inline float nrm2( const int n, const float *x, const int incx ) {  return snrm2_( &n, x, &incx ); }

inline double nrm2( const int n, const double *x, const int incx ) {  return dnrm2_( &n, x, &incx ); }

inline CoMpLeXtYpE nrm2( const int n, const CoMpLeXtYpE *x, const int incx ) {  return znrm2_( &n, x, &incx ); }

//  				******************   D O T  *******************
/**
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DDOT forms the dot product of two vectors.
*     uses unrolled loops for increments equal to one.
*/
inline float dot( const int n, const float *x, const int incx, const float *y, const int incy )
{  return sdot_( &n, x, &incx, y, &incy ); }

inline double dot( const int n, const double *x, const int incx, const double *y, const int incy )
{  return ddot_( &n, x, &incx, y, &incy ); }

inline CoMpLeXtYpE dot( const int n, const CoMpLeXtYpE *x, const int incx, const CoMpLeXtYpE *y, const int incy )
{  return zdot_( &n, x, &incx, y, &incy ); }

//  				*****************   A X P Y   *****************
/**
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
*     .. Scalar Arguments ..
      REAL SA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL SX(*),SY(*)
*     ..
*
*  Purpose
*  =======
*
*     SAXPY constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*
*            y <- a*x + y
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (SA.EQ.0.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*
*/
inline void axpy( const int n, const float a, const float *x, const int incx, 
                           float *y, const int incy )
{  saxpy_( &n, &a, x, &incx, y, &incy ); }

inline void axpy( const int n, const double a, const double *x, const int incx, 
                           double *y, const int incy )
{  daxpy_( &n, &a, x, &incx, y, &incy ); }

inline void axpy( const int n, const CoMpLeXtYpE a, const CoMpLeXtYpE *x, const int incx, 
                           CoMpLeXtYpE *y, const int incy )
{  zaxpy_( &n, &a, x, &incx, y, &incy ); }
					
//  				*****************   G E M V   *****************
/*
 GEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*/
inline void gemv( const char transa, const int m, const int n, 
					const float alpha, const float *A, const int lda,
 					const float *x, const int incx, 
					const float beta, float *y, const int incy )
{ sgemv_( &transa, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

inline void gemv( const char transa, const int m, const int n, 
					const double alpha, const double *A, const int lda,
 					const double *x, const int incx, 
					const double beta, double *y, const int incy )
{ dgemv_( &transa, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

inline void gemv( const char transa, const int m, const int n, 
					const CoMpLeXtYpE alpha, const CoMpLeXtYpE *A, const int lda,
 					const CoMpLeXtYpE *x, const int incx, 
					const CoMpLeXtYpE beta, CoMpLeXtYpE *y, const int incy )
{ zgemv_( &transa, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }


//  				*****************   S B M V   *****************
//
// Symmetric banded matrix vector multiplication

inline void sbmv( const char uplo, const int n, const int k, const float alpha, const float *A, const int lda,
						const float *x, const int incx, const float beta, const float *y, const int incy )
{ ssbmv_( &uplo, &n, &k, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }


inline void sbmv( const char uplo, const int n, const int k, const double alpha, const double *A, const int lda,
						const double *x, const int incx, const double beta, const double *y, const int incy )
{ dsbmv_( &uplo, &n, &k, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }


inline void sbmv( const char uplo, const int n, const int k, const CoMpLeXtYpE alpha, const CoMpLeXtYpE *A, const int lda,
						const CoMpLeXtYpE *x, const int incx, const CoMpLeXtYpE beta, const CoMpLeXtYpE *y, const int incy )
{ zsbmv_( &uplo, &n, &k, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }


//  				*****************   G E M M   *****************
/*
    
    C = alpha*op(A)*op(B) + beta*C
    
    op(A) is m by k
    op(B) is k by n
    C is m by n

    op is either 't', 'n', or 'c'
    
    lda, ldb and ldc are the leading dimensions of A, B and C respectively
    
    Only C is changed in this function call.

	C cannot be copied onto A or B (i.e. they must be three separate memory allocations)
    
*/
inline void gemm( const char transa, const char transb, const int m, const int n, const int k, 
					const float alpha, const float *A, const int lda,
 					const float *B, const int ldb, 
					const float beta, float *C, const int ldc )
{ sgemm_( &transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc ); }

inline void gemm( const char transa, const char transb, const int m, const int n, const int k, 
					const double alpha, const double *A, const int lda,
 					const double *B, const int ldb, 
					const double beta, double *C, const int ldc )
{ dgemm_( &transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc ); }

inline void gemm( const char transa, const char transb, const int m, const int n, const int k, 
					const CoMpLeXtYpE alpha, const CoMpLeXtYpE *A, const int lda,
 					const CoMpLeXtYpE *B, const int ldb, 
					const CoMpLeXtYpE beta, CoMpLeXtYpE *C, const int ldc )
{ zgemm_( &transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc ); }




#endif

