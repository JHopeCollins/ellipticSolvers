
# include "test-PoissonProblem.h"

   void Test_PoissonProblem::setUp()
  {
      int    nx=6;
      float  leftBoundary=0.;
      float rightBoundary=1.;

      float  leftValue=1.;
      float rightValue=2.;

      poisson.grid.setNx( nx );
      poisson.grid.setBounds( leftBoundary, rightBoundary );
      poisson.grid.discretise();

      poisson.setDirichletLeft(   leftValue );
      poisson.setDirichletRight( rightValue );
  }

   void Test_PoissonProblem::tearDown(){}

   void Test_PoissonProblem::test_set_number_of_point_sources_to_zero()
  {
      int nSources = 0;
      poisson.setNPointSources( nSources );

      CPPUNIT_ASSERT_EQUAL( nSources, poisson.nPointSources() );
  }

   void Test_PoissonProblem::test_set_number_of_point_sources_to_nonzero()
  {
      int nSources = 4;
      poisson.setNPointSources( nSources );

      CPPUNIT_ASSERT_EQUAL( nSources, poisson.nPointSources() );
  }

   void Test_PoissonProblem::test_set_one_PointSource()
  {
      int   sourceLocation=2;
      float sourceStrength=3.5;

      poisson.setNPointSources( 1 );
      poisson.setPointSource( sourceLocation, sourceStrength );

      float s=poisson.pointSource( sourceLocation );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( sourceStrength, s, floatError );
  }

   void Test_PoissonProblem::test_set_multiple_PointSources()
  {
      int   nsources=2;
      int   sourceLocation0=2;
      float sourceStrength0=3.5;
      int   sourceLocation1=4;
      float sourceStrength1=1.5;

      poisson.setNPointSources( nsources );
      poisson.setPointSource( sourceLocation0, sourceStrength0 );
      poisson.setPointSource( sourceLocation1, sourceStrength1 );

      float s=poisson.pointSource( sourceLocation0 );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( sourceStrength0, s, floatError );

      s=poisson.pointSource( sourceLocation1 );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( sourceStrength1, s, floatError );
  }

   void Test_PoissonProblem::test_constructRHS_no_sources()
  {
      LaplaceProblem laplace;
      float lb,rb;
      poisson.grid.bounds( lb,rb );

      laplace.grid.setNx( poisson.grid.nx() );
      laplace.grid.setBounds( lb,rb );
      laplace.grid.discretise();

      laplace.setDirichletLeft(  poisson.dirichletLeft()  );
      laplace.setDirichletRight( poisson.dirichletRight() );
      laplace.constructRHS();

      poisson.constructRHS();

      for( int i=1; i<poisson.grid.nx()-1; i++ )
     {
         float expected=laplace.rhs(i);
         float actual  =poisson.rhs(i);
         CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, actual, floatError );
     }
  }

   void Test_PoissonProblem::test_constructRHS_one_source()
  {
  }

   void Test_PoissonProblem::test_constructRHS_many_sources()
  {
  }

