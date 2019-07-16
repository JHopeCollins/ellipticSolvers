# ifndef LAPLACE_H
# define LAPLACE_H

# include "grid.h"

   class LaplaceProblem
  {
   private:

      float  leftBoundaryValue;
      float rightBoundaryValue;

      int     maximumSparsity;
      int    *sparsityVector;
      int   **connectionVector;

      float **offDiagonalCoefficientVector;
      float  *diagonalCoefficientVector;

      float *rhsVector;
      float *currentSolutionVector;
      float *previousSolutionVector;

      float jacobiTol;

      void mallocSolutionVectors();
      void constructFullMatrix( float *a );

   public:

      LaplaceProblem();
     ~LaplaceProblem();

      Grid  grid;

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

      void  constructRHS();
      float rhs( int i );

      void  setJacobiTolerance( float tol );
      float jacobiTolerance();

      void  jacobiIteration( float &residial );

      int   directSolve();
      int   jacobiSolve( int nIterations );

      float solution( int i );

      void fprintf( char *fileName );
  };

# endif
