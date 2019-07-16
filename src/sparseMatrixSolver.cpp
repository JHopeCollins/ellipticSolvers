# include "sparseMatrixSolver.h"
# include "linalg.h"

# include <cmath>
# include <assert.h>

   SparseMatrixSolver::SparseMatrixSolver()
  {
      matrixRank=0;
      maximumSparsity=0;
      jacobiTol=10e10;

      sparsityVector=nullptr;
      connectionVector=nullptr;
      offDiagonalCoefficientVector=nullptr;
      diagonalCoefficientVector=nullptr;
      rhsVector=nullptr;
      currentSolutionVector=nullptr;
      previousSolutionVector=nullptr;
  }

   SparseMatrixSolver::~SparseMatrixSolver()
  {
      delete[] sparsityVector;
      sparsityVector=nullptr;

      delete[] diagonalCoefficientVector;
      diagonalCoefficientVector=nullptr;

      delete[] rhsVector;
      rhsVector=nullptr;

      delete[] currentSolutionVector;
      currentSolutionVector=nullptr;

      delete[] previousSolutionVector;
      previousSolutionVector=nullptr;

      if( connectionVector )
     {
         for( int i=0; i<matrixRank; i++ )
        {
            delete[] connectionVector[i];
            connectionVector[i]=nullptr;
        }
         delete[] connectionVector;
         connectionVector=nullptr;
     }

      if( offDiagonalCoefficientVector )
     {
         for( int i=0; i<matrixRank; i++ )
        {
            delete[] offDiagonalCoefficientVector[i];
            offDiagonalCoefficientVector[i]=nullptr;
        }
         delete[] offDiagonalCoefficientVector;
         offDiagonalCoefficientVector=nullptr;
     }
  }

   void SparseMatrixSolver::mallocVectors()
  {
      assert( matrixRank>0 );
      assert( maximumSparsity>0 );
      int i,j;

      sparsityVector   = new int [ matrixRank ];

      diagonalCoefficientVector = new float [ matrixRank ];

      rhsVector              = new float[ matrixRank ];
      currentSolutionVector  = new float[ matrixRank ];
      previousSolutionVector = new float[ matrixRank ];

      connectionVector = new int*[ matrixRank ];
      for( int i=0; i<matrixRank; i++ )
     {
         connectionVector[i] = new int[ maximumSparsity ];
     }

      offDiagonalCoefficientVector = new float*[ matrixRank ];
      for( int i=0; i<matrixRank; i++ )
     {
         offDiagonalCoefficientVector[i] = new float[ maximumSparsity ];
     }

      for( i=0; i<matrixRank; i++ )
     { 
         sparsityVector[i]=0;
         diagonalCoefficientVector[i]=0;
         rhsVector[i]=0;
         currentSolutionVector[i]=0;
         previousSolutionVector[i]=0;

         for( j=0; j<maximumSparsity; j++ )
        {
            connectionVector[i][j]=0;
            offDiagonalCoefficientVector[i][j]=0;
        }
     }
  }

   void SparseMatrixSolver::constructFullMatrix( float *a )
  {
      float coeff;
      int row,column;
      int k,s,rowSparsity;
      int n=matrixRank;

      for( k=0; k<n*n; k++ ){ a[k]=0.; }

      for( row=0; row<n; row++ )
     {
      // diagonal elements
         column = row;
         coeff  = diagonalCoefficientVector[row];

         k= row + column*n;
         a[k] = coeff;

      // off-diagonal elements
         rowSparsity=sparsityVector[row];
         for( s=0; s<rowSparsity; s++ )
        {
            column = connectionVector[row][s];
            coeff  = offDiagonalCoefficientVector[row][s];

            k= row + column*n;
            a[k]=coeff;
        }
     }
      return;
  }

   void SparseMatrixSolver::jacobiIteration( float &residual )
  {
      int n=matrixRank;

      float next,coeff;
      int row,column,k;

      float res=0;
      for( row=0; row<n; row++ )
     {
         next = rhsVector[row];
         for( k=0; k<sparsityVector[row]; k++ )
        {
            column = connectionVector[row][k];
            coeff  = offDiagonalCoefficientVector[row][k];

            next  -= coeff*previousSolutionVector[column];
        }
         next/=diagonalCoefficientVector[row];

         res += fabs( next - previousSolutionVector[row] );
         currentSolutionVector[row]=next;
     }
      residual=res/n;
  }

   int SparseMatrixSolver::directSolve()
  {
      int n=matrixRank;
      if( !currentSolutionVector ){ currentSolutionVector=new float[n]; }

      float *fullMatrix=new float[ n*n ];
      int   *pivots=new int[n];
      int    info;
      int    nrhs=1;

      constructFullMatrix( fullMatrix );

      for( int i=0; i<n; i++ ){ previousSolutionVector[i]=rhsVector[i]; }

   //         N    NRHS  A           LDA IPIV    B          LDB, info )
      sgesv_( &n, &nrhs, fullMatrix, &n, pivots, previousSolutionVector, &n,  &info );

      for( int i=0; i<n; i++ )
     {
         currentSolutionVector[i]=previousSolutionVector[i];
     }

      delete[] fullMatrix; fullMatrix=nullptr;
      delete[] pivots; pivots=nullptr;
      return info;
  }

   int SparseMatrixSolver::jacobiSolve( int nIterations )
  {
      int n=matrixRank;
      float res=10e10;

      int i=0;
      jacobiIteration( res ); i++;

      while( res>jacobiTol )
     {
         for( int i=0; i<n; i++ )
        {
            previousSolutionVector[i]=currentSolutionVector[i];
        }
         jacobiIteration( res ); i++;
     }
      nIterations=i;
      return 0;
  }
