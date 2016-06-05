# A make file for building Basso programs
#


#CPPC = mpicxx
CPPC = icc
CPPFLAGS = -O3 

BASSO_DIR = /home1/02736/jfchessa/Basso

# BLAS and LAPACK
BLAS_LAPACK_DIR =  
BLAS_LAPACK_LIB =  -mkl

# Basso
BASSO_INC_DIR = -I$(BASSO_DIR)/include
BASSO_LIB_DIR = 

# trilinos
#TRILINOS_INC_DIR = -I$(TACC_TRILINOS_DIR)/include -I$(TACC_MKL_INC)
#TRILINOS_LIB_DIR = -L$(TACC_TRILINOS_LIB) -L$(TACC_MKL_LIB)
#TRILINOS_LIBS_TRIL = -lepetra -laztecoo 
 
# all together now
INC_DIR = $(BASSO_INC_DIR)  
LIB_DIR = $(BLAS_LAPACK_DIR)
LIBS = $(BLAS_LAPACK_LIB)
 
all:
	make test1 test2 test3
	
test1: test1.cpp
	$(CPPC) $(INC_DIR) $(CPPFLAGS) -o test1 test1.cpp $(LIB_DIR) $(LIBS)
	
test2: test2.cpp
	$(CPPC) $(INC_DIR) $(CPPFLAGS) -o test2 test2.cpp $(LIB_DIR) $(LIBS)
	
test3: test3.cpp
	$(CPPC) $(INC_DIR) $(CPPFLAGS) -o test3 test3.cpp $(LIB_DIR) $(LIBS)
	
test4: test4.cpp
	$(CPPC) $(INC_DIR) $(CPPFLAGS) -o test4 test4.cpp $(LIB_DIR) $(LIBS)
		
clean:
	rm -f test1 test2 test3 test4