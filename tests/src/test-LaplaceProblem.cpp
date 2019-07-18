
# include "test-LaplaceProblem.h"

   void Test_LaplaceProblem::setUp(){}

   void Test_LaplaceProblem::tearDown(){}

   void Test_LaplaceProblem::test_setDirichletLeft()
  {
      LaplaceProblem problem;
      float leftDirichletValue=1.0;

      problem.setDirichletLeft( leftDirichletValue );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( leftDirichletValue, problem.dirichletLeft(), floatError );
  }

   void Test_LaplaceProblem::test_setDirichletRight()
  {
      LaplaceProblem problem;
      float rightDirichletValue=5.0;

      problem.setDirichletRight( rightDirichletValue );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( rightDirichletValue, problem.dirichletRight(), floatError );
  }

   void Test_LaplaceProblem::test_constructRHS_equal_dirichlet_conditions()
  {
      LaplaceProblem problem;
      int   nx=5;
      float  leftBoundary=0.;
      float rightBoundary=1.;

      float  leftDirichlet=1.;
      float rightDirichlet=1.;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      problem.setDirichletLeft(   leftDirichlet );
      problem.setDirichletRight( rightDirichlet );
      problem.constructRHS();

      CPPUNIT_ASSERT_DOUBLES_EQUAL(  -leftDirichlet, problem.rhs(1), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(              0., problem.rhs(2), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( -rightDirichlet, problem.rhs(3), floatError );
  }

   void Test_LaplaceProblem::test_constructRHS_different_dirichlet_conditions()
  {
      LaplaceProblem problem;
      int   nx=5;
      float leftBoundary=0.;
      float rightBoundary=1.;

      float  leftDirichlet=1.;
      float rightDirichlet=4.5;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      problem.setDirichletLeft(   leftDirichlet );
      problem.setDirichletRight( rightDirichlet );
      problem.constructRHS();

      CPPUNIT_ASSERT_DOUBLES_EQUAL(  -leftDirichlet, problem.rhs(1), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(              0., problem.rhs(2), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( -rightDirichlet, problem.rhs(3), floatError );
  }

   void Test_LaplaceProblem::test_sparsity_vector()
  {
      LaplaceProblem problem;
      int    nx=5;
      float  leftBoundary=0.;
      float rightBoundary=1.;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      CPPUNIT_ASSERT_EQUAL( 2, problem.maxSparsity() );

      problem.constructSparsity();

      CPPUNIT_ASSERT_EQUAL( 1, problem.sparsity(1) );
      CPPUNIT_ASSERT_EQUAL( 2, problem.sparsity(2) );
      CPPUNIT_ASSERT_EQUAL( 1, problem.sparsity(3) );

      for( int i=1; i<nx-1; i++ )
     {
         CPPUNIT_ASSERT( problem.sparsity(i)<=problem.maxSparsity() );
     }
  }

   void Test_LaplaceProblem::test_connection_vector()
  {
      LaplaceProblem problem;
      int    nx=6;
      float  leftBoundary=0.;
      float rightBoundary=1.;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      problem.constructConnections();

      CPPUNIT_ASSERT_EQUAL( 2, problem.connection( 1, 0 ) );
      CPPUNIT_ASSERT_EQUAL(-1, problem.connection( 1, 1 ) );

      CPPUNIT_ASSERT_EQUAL( 1, problem.connection( 2, 0 ) );
      CPPUNIT_ASSERT_EQUAL( 3, problem.connection( 2, 1 ) );

      CPPUNIT_ASSERT_EQUAL( 2, problem.connection( 3, 0 ) );
      CPPUNIT_ASSERT_EQUAL( 4, problem.connection( 3, 1 ) );

      CPPUNIT_ASSERT_EQUAL( 3, problem.connection( 4, 0 ) );
      CPPUNIT_ASSERT_EQUAL(-1, problem.connection( 4, 1 ) );
 }

   void Test_LaplaceProblem::test_diagonal_coefficient_vector()
  {
      LaplaceProblem problem;
      int    nx=6;
      float  leftBoundary=0.;
      float rightBoundary=1.;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      problem.constructDiagonalCoefficients();

      for( int i=1; i<nx-1; i++ )
     {
         CPPUNIT_ASSERT_DOUBLES_EQUAL( -2., problem.diagonalCoefficient( i ), floatError );
     }
  }

   void Test_LaplaceProblem::test_offdiagonal_coefficient_vector()
  {
      LaplaceProblem problem;
      int    nx=6;
      float  leftBoundary=0.;
      float rightBoundary=1.;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      problem.constructOffDiagonalCoefficients();

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1, problem.offDiagonalCoefficient( 1, 0 ), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, problem.offDiagonalCoefficient( 1, 1 ), floatError );

      for( int i=2; i<nx-2; i++ )
     {
         CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., problem.offDiagonalCoefficient( i, 0 ), floatError );
         CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., problem.offDiagonalCoefficient( i, 1 ), floatError );
     }

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1, problem.offDiagonalCoefficient( 4, 0 ), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, problem.offDiagonalCoefficient( 4, 1 ), floatError );
  }

   void Test_LaplaceProblem::test_directSolve_gives_constant_for_equal_boundary_values()
  {
      LaplaceProblem problem;
      int    nx=6;
      float  leftBoundary=0.;
      float rightBoundary=1.;

      float  solutionValue=1.;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      problem.setDirichletLeft(  solutionValue );
      problem.setDirichletRight( solutionValue );

      problem.constructSparsity();
      problem.constructConnections();

      problem.constructRHS();
      problem.constructDiagonalCoefficients();
      problem.constructOffDiagonalCoefficients();

      int info=problem.directSolve();
      CPPUNIT_ASSERT_EQUAL( 0, info );

      for( int i=0; i<nx; i++ )
     {
         CPPUNIT_ASSERT_DOUBLES_EQUAL( solutionValue, problem.solution(i), floatError );
     }
  }

   void Test_LaplaceProblem::test_directSolve_gives_linear_for_different_boundary_values()
  {
      LaplaceProblem problem;
      int    nx=6;
      float  leftBoundary=0.;
      float rightBoundary=1.;

      float  leftValue=1.;
      float rightValue=2.;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      problem.setDirichletLeft(   leftValue );
      problem.setDirichletRight( rightValue );

      problem.constructSparsity();
      problem.constructConnections();

      problem.constructRHS();
      problem.constructDiagonalCoefficients();
      problem.constructOffDiagonalCoefficients();

      int info=problem.directSolve();
      CPPUNIT_ASSERT_EQUAL( 0, info );

      for( int i=0; i<nx; i++ )
     {
         float x=problem.grid[i];
         float expected = 1. + x;
         CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, problem.solution(i), floatError );
     }
  }

   void Test_LaplaceProblem::test_jacobiIteration()
  {
      LaplaceProblem problem;
      int    nx=6;
      float  leftBoundary=2.;
      float rightBoundary=1.;

      float  solutionValue=1.;
      float  residual;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      problem.setDirichletLeft(  solutionValue );
      problem.setDirichletRight( solutionValue );

      problem.constructSparsity();
      problem.constructConnections();

      problem.constructRHS();
      problem.constructDiagonalCoefficients();
      problem.constructOffDiagonalCoefficients();

      problem.jacobiIteration( residual );

      CPPUNIT_ASSERT_DOUBLES_EQUAL(     solutionValue, problem.solution(0), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*solutionValue, problem.solution(1), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(                 0, problem.solution(2), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(                 0, problem.solution(3), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5*solutionValue, problem.solution(4), floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(     solutionValue, problem.solution(5), floatError );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( solutionValue/(nx-2), residual, floatError );

  }

   void Test_LaplaceProblem::test_jacobiTolerance()
  {
      LaplaceProblem problem;
      float tol=10e-3;

      problem.setJacobiTolerance( tol );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( tol, problem.jacobiTolerance(), floatError );
  }

   void Test_LaplaceProblem::test_jacobiSolve_gives_constant_for_equal_boundary_values()
  {
      LaplaceProblem problem;
      int    nx=6;
      float  leftBoundary=0.;
      float rightBoundary=1.;

      float  solutionValue=1.;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      problem.setDirichletLeft(  solutionValue );
      problem.setDirichletRight( solutionValue );

      problem.constructSparsity();
      problem.constructConnections();

      problem.constructRHS();
      problem.constructDiagonalCoefficients();
      problem.constructOffDiagonalCoefficients();

      float tol=10e-4;
      problem.setJacobiTolerance( tol );

      int nIterations=0;
      int info=problem.jacobiSolve( nIterations );
      CPPUNIT_ASSERT_EQUAL( 0, info );

      for( int i=0; i<nx; i++ )
     {
         CPPUNIT_ASSERT_DOUBLES_EQUAL( solutionValue, problem.solution(i), 10*tol );
     }
  }

   void Test_LaplaceProblem::test_jacobiSolve_gives_linear_for_different_boundary_values()
  {
      LaplaceProblem problem;
      int    nx=6;
      float  leftBoundary=0.;
      float rightBoundary=1.;

      float  leftValue=1.;
      float rightValue=2.;

      problem.grid.setNx( nx );
      problem.grid.setBounds( leftBoundary, rightBoundary );
      problem.grid.discretise();
      problem.scrapeGrid();

      problem.setDirichletLeft(   leftValue );
      problem.setDirichletRight( rightValue );

      problem.constructSparsity();
      problem.constructConnections();

      problem.constructRHS();
      problem.constructDiagonalCoefficients();
      problem.constructOffDiagonalCoefficients();

      float tol=10e-4;
      problem.setJacobiTolerance( tol );

      int nIterations=0;
      int info=problem.jacobiSolve( nIterations );
      CPPUNIT_ASSERT_EQUAL( 0, info );

      for( int i=0; i<nx; i++ )
     {
         float x=problem.grid[i];
         float expected = 1. + x;
         CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, problem.solution(i), 10*tol );
     }
  }

