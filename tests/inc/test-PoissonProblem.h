# ifndef  TEST_POISSON_H
# define  TEST_POISSON_H

#include "cppunit/TestFixture.h"
#include "cppunit/extensions/HelperMacros.h"

#include "poissonProblem.h"

/*
   Tests for public methods of the PoissonProblem class
*/

   class Test_PoissonProblem : public CppUnit::TestFixture
  {
   private:
      CPPUNIT_TEST_SUITE( Test_PoissonProblem );

         CPPUNIT_TEST( test_set_number_of_point_sources_to_zero );
         CPPUNIT_TEST( test_set_number_of_point_sources_to_nonzero );
         CPPUNIT_TEST( test_set_one_PointSource );
         CPPUNIT_TEST( test_set_multiple_PointSources );
         CPPUNIT_TEST( test_constructRHS_no_sources );
         CPPUNIT_TEST( test_constructRHS_one_source );
         CPPUNIT_TEST( test_constructRHS_many_sources );

      CPPUNIT_TEST_SUITE_END();

      float floatError=10e-6;
      PoissonProblem poisson;

   public:
      void setUp();
      void tearDown();

      void test_set_number_of_point_sources_to_zero();
      void test_set_number_of_point_sources_to_nonzero();
      void test_set_one_PointSource();
      void test_set_multiple_PointSources();
      void test_constructRHS_no_sources();
      void test_constructRHS_one_source();
      void test_constructRHS_many_sources();
  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_PoissonProblem );

# endif
