# ifndef SPARSE_MATRIX_SOLVER_H
# define SPARSE_MATRIX_SOLVER_H

   struct SparseMatrixSolver
  {
      int   matrixRank;
      int   maximumSparsity;
      float jacobiTol;

      int    *sparsityVector;
      int   **connectionVector;

      float **offDiagonalCoefficientVector;
      float  *diagonalCoefficientVector;

      float  *rhsVector;
      float  *currentSolutionVector;
      float  *previousSolutionVector;

      SparseMatrixSolver();
     ~SparseMatrixSolver();

      void    mallocVectors();
      void    constructFullMatrix( float *a );

      void    jacobiIteration( float &residial );

      int     directSolve();
      int     jacobiSolve( int nIterations );
  };

# endif
