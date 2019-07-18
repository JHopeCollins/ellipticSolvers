
# include "test-SparseMatrixSolver.h"

   void Test_SparseMatrixSolver::setUp(){}

   void Test_SparseMatrixSolver::tearDown(){}

   void Test_SparseMatrixSolver::test_mallocVectors()
  {
      SparseMatrixSolver solver;
      int rank=3;
      int bandwidth=2;

      solver.matrixRank=rank;
      solver.maximumSparsity=bandwidth;

      solver.mallocVectors();

      CPPUNIT_ASSERT( solver.rhsVector                    );
      CPPUNIT_ASSERT( solver.sparsityVector               );
      CPPUNIT_ASSERT( solver.connectionVector             );
      CPPUNIT_ASSERT( solver.currentSolutionVector        );
      CPPUNIT_ASSERT( solver.previousSolutionVector       );
      CPPUNIT_ASSERT( solver.diagonalCoefficientVector    );
      CPPUNIT_ASSERT( solver.offDiagonalCoefficientVector );

      for( int i=0; i<rank; i++ )
     {
         CPPUNIT_ASSERT( solver.connectionVector[i]             );
         CPPUNIT_ASSERT( solver.offDiagonalCoefficientVector[i] );
     }
  }

   void Test_SparseMatrixSolver::test_constructFullMatrix()
  {
      SparseMatrixSolver   solver;
      int   rank=4;
      int   sparse=2;
      float a[rank][rank];
      float b[rank*rank];

      solver.matrixRank=rank;
      solver.maximumSparsity=sparse;
      solver.mallocVectors();

      a[0][0]=2.3; a[0][1]=0.2; a[0][2]=0.0; a[0][3]=0.0;
      a[1][0]=0.8; a[1][1]=3.2; a[1][2]=1.2; a[1][3]=0.0;
      a[2][0]=0.3; a[2][1]=0.0; a[2][2]=1.0; a[2][3]=0.9;
      a[3][0]=0.1; a[3][1]=0.0; a[3][2]=0.0; a[3][3]=0.1;

      solver.diagonalCoefficientVector[0]=a[0][0];
      solver.diagonalCoefficientVector[1]=a[1][1];
      solver.diagonalCoefficientVector[2]=a[2][2];
      solver.diagonalCoefficientVector[3]=a[3][3];

      solver.sparsityVector[0]=1;
      solver.sparsityVector[1]=2;
      solver.sparsityVector[2]=2;
      solver.sparsityVector[3]=1;

      solver.connectionVector[0][0]=1;  solver.offDiagonalCoefficientVector[0][0]=a[0][1];
      solver.connectionVector[1][0]=0;  solver.offDiagonalCoefficientVector[1][0]=a[1][0];
      solver.connectionVector[1][1]=2;  solver.offDiagonalCoefficientVector[1][1]=a[1][2];
      solver.connectionVector[2][0]=0;  solver.offDiagonalCoefficientVector[2][0]=a[2][0];
      solver.connectionVector[2][1]=3;  solver.offDiagonalCoefficientVector[2][1]=a[2][3];
      solver.connectionVector[3][0]=0;  solver.offDiagonalCoefficientVector[3][0]=a[3][0];

      solver.constructFullMatrix( b );

      int i,j,k=0;
      float expected,actual;
      for( i=0; i<rank; i++ )
     {
         for( j=0; j<rank; j++ )
        {
            expected = a[j][i];
            actual = b[k];
            CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, floatError );
            k++;
        }
     }
  }

   void Test_SparseMatrixSolver::test_directSolve()
  {
      SparseMatrixSolver   solver;
      int   rank=4;
      int   sparse=2;
      float a[rank][rank];
      float x[rank];

      solver.matrixRank=rank;
      solver.maximumSparsity=sparse;
      solver.mallocVectors();

      a[0][0]=2.3; a[0][1]=0.2; a[0][2]=0.0; a[0][3]=0.0;
      a[1][0]=0.8; a[1][1]=3.2; a[1][2]=1.2; a[1][3]=0.0;
      a[2][0]=0.3; a[2][1]=0.0; a[2][2]=1.0; a[2][3]=0.9;
      a[3][0]=0.1; a[3][1]=0.0; a[3][2]=0.0; a[3][3]=0.1;

      solver.diagonalCoefficientVector[0]=a[0][0];
      solver.diagonalCoefficientVector[1]=a[1][1];
      solver.diagonalCoefficientVector[2]=a[2][2];
      solver.diagonalCoefficientVector[3]=a[3][3];

      solver.sparsityVector[0]=1;
      solver.sparsityVector[1]=2;
      solver.sparsityVector[2]=2;
      solver.sparsityVector[3]=1;

      solver.connectionVector[0][0]=1;  solver.offDiagonalCoefficientVector[0][0]=a[0][1];
      solver.connectionVector[1][0]=0;  solver.offDiagonalCoefficientVector[1][0]=a[1][0];
      solver.connectionVector[1][1]=2;  solver.offDiagonalCoefficientVector[1][1]=a[1][2];
      solver.connectionVector[2][0]=0;  solver.offDiagonalCoefficientVector[2][0]=a[2][0];
      solver.connectionVector[2][1]=3;  solver.offDiagonalCoefficientVector[2][1]=a[2][3];
      solver.connectionVector[3][0]=0;  solver.offDiagonalCoefficientVector[3][0]=a[3][0];

      solver.rhsVector[0]=4.;
      solver.rhsVector[1]=3.;
      solver.rhsVector[2]=2.;
      solver.rhsVector[3]=1.;

      x[0] =  1315./882.;
      x[1] =  5035./1764.;
      x[2] = -1795./294.;
      x[3] =  7505./882.;
      
      int info = solver.directSolve();
      CPPUNIT_ASSERT_EQUAL( 0, info );

      float expected,actual;
      for( int i=0; i<rank; i++ )
     {
         expected = x[i];
         actual = solver.currentSolutionVector[i];
         CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, floatError );
     }
  }

   void Test_SparseMatrixSolver::test_jacobiIteration()
  {
      SparseMatrixSolver   solver;
      int   rank=4;
      int   sparse=2;
      float a[rank][rank];
      float x0[rank];

      solver.matrixRank=rank;
      solver.maximumSparsity=sparse;
      solver.mallocVectors();

      a[0][0]=2.3; a[0][1]=0.2; a[0][2]=0.0; a[0][3]=0.0;
      a[1][0]=0.8; a[1][1]=3.2; a[1][2]=1.2; a[1][3]=0.0;
      a[2][0]=0.3; a[2][1]=0.0; a[2][2]=1.0; a[2][3]=0.9;
      a[3][0]=0.1; a[3][1]=0.0; a[3][2]=0.0; a[3][3]=0.1;

      solver.diagonalCoefficientVector[0]=a[0][0];
      solver.diagonalCoefficientVector[1]=a[1][1];
      solver.diagonalCoefficientVector[2]=a[2][2];
      solver.diagonalCoefficientVector[3]=a[3][3];

      solver.sparsityVector[0]=1;
      solver.sparsityVector[1]=2;
      solver.sparsityVector[2]=2;
      solver.sparsityVector[3]=1;

      solver.connectionVector[0][0]=1;  solver.offDiagonalCoefficientVector[0][0]=a[0][1];
      solver.connectionVector[1][0]=0;  solver.offDiagonalCoefficientVector[1][0]=a[1][0];
      solver.connectionVector[1][1]=2;  solver.offDiagonalCoefficientVector[1][1]=a[1][2];
      solver.connectionVector[2][0]=0;  solver.offDiagonalCoefficientVector[2][0]=a[2][0];
      solver.connectionVector[2][1]=3;  solver.offDiagonalCoefficientVector[2][1]=a[2][3];
      solver.connectionVector[3][0]=0;  solver.offDiagonalCoefficientVector[3][0]=a[3][0];

      solver.rhsVector[0]=4.;
      solver.rhsVector[1]=3.;
      solver.rhsVector[2]=2.;
      solver.rhsVector[3]=1.;

      float res=0;
      solver.jacobiIteration( res );

      x0[0] = solver.rhsVector[0]/a[0][0];
      x0[1] = solver.rhsVector[1]/a[1][1];
      x0[2] = solver.rhsVector[2]/a[2][2];
      x0[3] = solver.rhsVector[3]/a[3][3];

      float expected,actual;
      for( int i=0; i<rank; i++ )
     {
         expected = x0[i];
         actual = solver.currentSolutionVector[i];
         CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, floatError );
     }
  }

   void Test_SparseMatrixSolver::test_jacobiSolve()
  {
      SparseMatrixSolver   solver;
      int   rank=4;
      int   sparse=2;
      float a[rank][rank];
      float x[rank];

      solver.matrixRank=rank;
      solver.maximumSparsity=sparse;
      solver.mallocVectors();

      a[0][0]=2.3; a[0][1]=0.2; a[0][2]=0.0; a[0][3]=0.0;
      a[1][0]=0.8; a[1][1]=3.2; a[1][2]=1.2; a[1][3]=0.0;
      a[2][0]=0.3; a[2][1]=0.0; a[2][2]=1.0; a[2][3]=0.9;
      a[3][0]=0.1; a[3][1]=0.0; a[3][2]=0.0; a[3][3]=0.1;

      solver.diagonalCoefficientVector[0]=a[0][0];
      solver.diagonalCoefficientVector[1]=a[1][1];
      solver.diagonalCoefficientVector[2]=a[2][2];
      solver.diagonalCoefficientVector[3]=a[3][3];

      solver.sparsityVector[0]=1;
      solver.sparsityVector[1]=2;
      solver.sparsityVector[2]=2;
      solver.sparsityVector[3]=1;

      solver.connectionVector[0][0]=1;  solver.offDiagonalCoefficientVector[0][0]=a[0][1];
      solver.connectionVector[1][0]=0;  solver.offDiagonalCoefficientVector[1][0]=a[1][0];
      solver.connectionVector[1][1]=2;  solver.offDiagonalCoefficientVector[1][1]=a[1][2];
      solver.connectionVector[2][0]=0;  solver.offDiagonalCoefficientVector[2][0]=a[2][0];
      solver.connectionVector[2][1]=3;  solver.offDiagonalCoefficientVector[2][1]=a[2][3];
      solver.connectionVector[3][0]=0;  solver.offDiagonalCoefficientVector[3][0]=a[3][0];

      solver.rhsVector[0]=4.;
      solver.rhsVector[1]=3.;
      solver.rhsVector[2]=2.;
      solver.rhsVector[3]=1.;

      x[0] =  1315./882.;
      x[1] =  5035./1764.;
      x[2] = -1795./294.;
      x[3] =  7505./882.;

      int info, nIterations=0;

      solver.jacobiTol=jacobiError;
      info = solver.jacobiSolve( nIterations );

      CPPUNIT_ASSERT_EQUAL( 0, info );

      float expected,actual;
      for( int i=0; i<rank; i++ )
     {
         expected = x[i];
         actual = solver.currentSolutionVector[i];
         CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, 10*jacobiError );
     }
  }

