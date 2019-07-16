
# include "laplaceProblem.h"
# include "linalg.h"

# include <assert.h>
# include <cstdio>
# include <cmath>

   LaplaceProblem::LaplaceProblem()
  {
      leftBoundaryValue =0.;
      rightBoundaryValue=0.;

      maximumSparsity=2;
      offDiagonalCoefficientVector=nullptr;
         diagonalCoefficientVector=nullptr;
            previousSolutionVector=nullptr;
             currentSolutionVector=nullptr;
                  connectionVector=nullptr;
                    sparsityVector=nullptr;
                         rhsVector=nullptr;
  }

   LaplaceProblem::~LaplaceProblem()
  {
      delete[] rhsVector;
      rhsVector=nullptr;

      delete[] sparsityVector;
      sparsityVector=nullptr;

      delete[] diagonalCoefficientVector;
      diagonalCoefficientVector=nullptr;

      if( connectionVector )
     {
         for( int i=0; i<grid.nx()-2; i++ )
        {
            delete[] connectionVector[i];
            connectionVector[i]=nullptr;
        }
     }

      if( offDiagonalCoefficientVector )
     {
         for( int i=0; i<grid.nx()-2; i++ )
        {
            delete[] offDiagonalCoefficientVector[i];
            offDiagonalCoefficientVector[i]=nullptr;
        }
     }

      delete[] connectionVector;
      connectionVector=nullptr;

      delete[] offDiagonalCoefficientVector;
      offDiagonalCoefficientVector=nullptr;

      delete[] currentSolutionVector;
      currentSolutionVector=nullptr;

      delete[] previousSolutionVector;
      previousSolutionVector=nullptr;
  }

   void LaplaceProblem::mallocSolutionVectors()
  {
      int n=grid.nx()-2;
      if( !currentSolutionVector  ){ currentSolutionVector=new float[n]; }
      if( !previousSolutionVector )
     {
         previousSolutionVector = new float[n];
         for( int i=0; i<n; i++ )
        {
            previousSolutionVector[i]=0.;
            currentSolutionVector[i]=0.;
        }
     }
      return;
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
      int n=grid.nx()-2;

      if( !sparsityVector ){ sparsityVector = new int[n]; }

      sparsityVector[  0]=1;
      sparsityVector[n-1]=1;

      for( int i=1; i<n-1; i++ ){ sparsityVector[i]=2.; }
  }

   void LaplaceProblem::constructConnections()
  {
      if( !( grid.isDiscretised() ) ){ grid.discretise(); }
      int n=grid.nx()-2;

      if( !connectionVector ){ connectionVector = new int*[n]; }

      for( int i=0; i<n; i++ )
     {
         connectionVector[i] = new int[ maximumSparsity ];
         for( int j=0; j<maximumSparsity; j++ )
        {
            connectionVector[i][j]=-1;
        }
     }

      connectionVector[  0][0]=1;
      connectionVector[n-1][0]=n-2;

   // shifts for boundary nodes
      for( int i=1; i<n-1; i++ )
     {
         connectionVector[i][0]=i-1;
         connectionVector[i][1]=i+1;
     }
  }

   void LaplaceProblem::constructDiagonalCoefficients()
  {
      if( !( grid.isDiscretised() ) ){ grid.discretise(); }
      int n=grid.nx()-2;

      if( !diagonalCoefficientVector ){ diagonalCoefficientVector = new float[n]; }

      for( int i=0; i<n; i++ ){ diagonalCoefficientVector[i]=-2.; }
  }

   void LaplaceProblem::constructOffDiagonalCoefficients()
  {
      if( !( grid.isDiscretised() ) ){ grid.discretise(); }
      int n=grid.nx()-2;

      if( !offDiagonalCoefficientVector ){ offDiagonalCoefficientVector = new float*[n]; }

      for( int i=0; i<n; i++ )
     {
         offDiagonalCoefficientVector[i] = new float[ maximumSparsity ];
         for( int j=0; j<maximumSparsity; j++ )
        {
            offDiagonalCoefficientVector[i][j]=1.;
        }
     }
      offDiagonalCoefficientVector[  0][1]=0;
      offDiagonalCoefficientVector[n-1][1]=0;
  }

   void LaplaceProblem::constructRHS()
  {
      if( !( grid.isDiscretised() ) ){ grid.discretise(); }

      int n=grid.nx()-2;
      if( !rhsVector ){ rhsVector = new float[n]; }

      rhsVector[  0]= -leftBoundaryValue;
      rhsVector[n-1]=-rightBoundaryValue;

      for( int i=1; i<n-1; i++ ){ rhsVector[i]=0.; }
  }

   void LaplaceProblem::constructFullMatrix( float *a )
  {
      float coeff;
      int row,column;
      int k,s,rowSparsity;
      int n=grid.nx()-2;

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

   int LaplaceProblem::directSolve()
  {
      int n=grid.nx()-2;
      if( !currentSolutionVector ){ currentSolutionVector=new float[n]; }

      float *fullMatrix=new float[ n*n ];
      int   *pivots=new int[n];
      int    info;
      int    nrhs=1;

      constructFullMatrix( fullMatrix );

   //         N    NRHS  A           LDA IPIV    B          LDB, info )
      sgesv_( &n, &nrhs, fullMatrix, &n, pivots, rhsVector, &n,  &info );

      for( int i=0; i<n; i++ )
     {
         currentSolutionVector[i]=rhsVector[i];
     }

      constructRHS();
      delete[] fullMatrix; fullMatrix=nullptr;
      delete[] pivots; pivots=nullptr;
      return info;
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

   void LaplaceProblem::jacobiIteration( float &residual )
  {
      int n=grid.nx()-2;

      float next,coeff;
      int row,column,k;

      mallocSolutionVectors();

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

   int LaplaceProblem::jacobiSolve( int nIterations )
  {
      int n=grid.nx()-2;
      float res=10e10;

      mallocSolutionVectors();

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

