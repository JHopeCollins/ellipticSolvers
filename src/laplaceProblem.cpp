
# include "laplaceProblem.h"
# include "linalg.h"

# include <assert.h>
# include <cstdio>
# include <cmath>

   LaplaceProblem::LaplaceProblem()
  {
      leftBoundaryValue =0.;
      rightBoundaryValue=0.;
      jacobiTol=10e10;
      maximumSparsity=2;
  }

   LaplaceProblem::~LaplaceProblem(){}

   void LaplaceProblem::scrapeGrid()
  {
      matrixRank=grid.nx()-2;
      mallocVectors();
  }

   void LaplaceProblem::setDirichletLeft( float ld )
  {
      leftBoundaryValue=ld;
      return;
  }

   void LaplaceProblem::setDirichletRight( float rd )
  {
      rightBoundaryValue=rd;
      return;
  }

   float LaplaceProblem::dirichletLeft()
  {
      return leftBoundaryValue;
  }

   float LaplaceProblem::dirichletRight()
  {
      return rightBoundaryValue;
  }

   int LaplaceProblem::maxSparsity()
  {
      return maximumSparsity;
  }

   void LaplaceProblem::constructSparsity()
  {
      if( !( grid.isDiscretised() ) ){ grid.discretise(); }

      sparsityVector[0]=1;
      sparsityVector[matrixRank-1]=1;

      for( int i=1; i<matrixRank-1; i++ ){ sparsityVector[i]=2.; }
  }

   void LaplaceProblem::constructConnections()
  {
      if( !( grid.isDiscretised() ) ){ grid.discretise(); }

      for( int i=0; i<matrixRank; i++ )
     {
//       connectionVector[i] = new int[ maximumSparsity ];
         for( int j=0; j<maximumSparsity; j++ )
        {
            connectionVector[i][j]=-1;
        }
     }

      connectionVector[0][0]=1;
      connectionVector[matrixRank-1][0]=matrixRank-2;

   // shifts for boundary nodes
      for( int i=1; i<matrixRank-1; i++ )
     {
         connectionVector[i][0]=i-1;
         connectionVector[i][1]=i+1;
     }
  }

   void LaplaceProblem::constructDiagonalCoefficients()
  {
      if( !( grid.isDiscretised() ) ){ grid.discretise(); }

      for( int i=0; i<matrixRank; i++ ){ diagonalCoefficientVector[i]=-2.; }
  }

   void LaplaceProblem::constructOffDiagonalCoefficients()
  {
      if( !( grid.isDiscretised() ) ){ grid.discretise(); }

      for( int i=0; i<matrixRank; i++ )
     {
         for( int j=0; j<maximumSparsity; j++ )
        {
            offDiagonalCoefficientVector[i][j]=1.;
        }
     }
      offDiagonalCoefficientVector[0][1]=0;
      offDiagonalCoefficientVector[matrixRank-1][1]=0;
  }

   void LaplaceProblem::constructRHS()
  {
      if( !( grid.isDiscretised() ) ){ grid.discretise(); }

      rhsVector[0]= -leftBoundaryValue;
      rhsVector[matrixRank-1]=-rightBoundaryValue;

      for( int i=1; i<matrixRank-1; i++ ){ rhsVector[i]=0.; }
  }

   void LaplaceProblem::setJacobiTolerance( float tol )
  {
      jacobiTol=tol;
      return;
  }

   float LaplaceProblem::jacobiTolerance()
  {
      return jacobiTol;
  }

   int LaplaceProblem::sparsity( int i )
  {
   // shifts for boundary nodes
      assert( (i>0) and (i<grid.nx()-1) );
      return sparsityVector[i-1];
  }

   int LaplaceProblem::connection( int i, int c )
  {
   // shifts for boundary nodes
      assert( (i> 0) and (i<grid.nx()-1) );
      assert( (c>-1) and (c<maximumSparsity) );
      int k=connectionVector[i-1][c];

   // -ve for unconnected
      if( k<0 ){ return k; }
      return k+1;
  }

   float LaplaceProblem::diagonalCoefficient( int i )
  {
   // shifts for boundary nodes
      assert( (i>0) and (i<grid.nx()-1) );
      return diagonalCoefficientVector[i-1];
  }

   float LaplaceProblem::offDiagonalCoefficient( int i, int c )
  {
   // shifts for boundary nodes
      assert( (i> 0) and (i<grid.nx()-1) );
      assert( (c>-1) and (c<maximumSparsity) );
      return offDiagonalCoefficientVector[i-1][c];
  }

   float LaplaceProblem::rhs( int i )
  {
   // shifts for boundary nodes
      assert( (i>0) and (i<grid.nx()-1) );
      return rhsVector[i-1];
  }

   float LaplaceProblem::solution( int i )
  {
      assert( (i>-1) and (i<grid.nx()) );
      if( i==0 )
     {
         return leftBoundaryValue;
     }
      else if( i==grid.nx()-1 )
     {
         return rightBoundaryValue;
     }
      else
     {
      // boundary shift
         return currentSolutionVector[i-1];
     }
  }

   void LaplaceProblem::fprintf( char *fileName )
  {
      FILE *f=fopen( fileName, "w" );
      assert( f );

      for( int i=0; i<grid.nx(); i++ )
     {
       ::fprintf( f, "%f   %f\n",  grid[i], solution(i) );
     }
      fclose( f );
  }

