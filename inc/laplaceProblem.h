# ifndef LAPLACE_H
# define LAPLACE_H 
# include "sparseMatrixSolver.h"
# include "grid.h"

   class LaplaceProblem : protected SparseMatrixSolver
  {
   private:

      float  leftBoundaryValue;
      float rightBoundaryValue;

   public:

      LaplaceProblem();
     ~LaplaceProblem();

      Grid  grid;

      void scrapeGrid();

      void setDirichletLeft(  float ld );
      void setDirichletRight( float rd );

      float dirichletLeft();
      float dirichletRight();

      int  maxSparsity();
      void constructSparsity();
      int  sparsity( int i );

      void constructConnections();
      int  connection( int i, int c );

      void  constructDiagonalCoefficients();
      float diagonalCoefficient( int i );

      void  constructOffDiagonalCoefficients();
      float offDiagonalCoefficient( int i, int c );

      virtual void constructRHS();
      float rhs( int i );

      void  setJacobiTolerance( float tol );
      float jacobiTolerance();

      using SparseMatrixSolver::directSolve;
      using SparseMatrixSolver::jacobiIteration;
      using SparseMatrixSolver::jacobiSolve;

      float solution( int i );

      void fprintf( char *fileName );
  };

# endif
