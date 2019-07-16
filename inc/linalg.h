
# ifndef LINALG_H
# define LINALG_H

      //                   N      NRHS   A         LDA    IPIV   B         LDB,   info )
   extern "C" void dgesv_( int *, int *, double *, int *, int *, double *, int *, int * );

      //                   N      NRHS   A        LDA    IPIV   B        LDB,   info )
   extern "C" void sgesv_( int *, int *, float *, int *, int *, float *, int *, int * );

# endif
