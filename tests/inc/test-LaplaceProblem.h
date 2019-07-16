# ifndef  TEST_LAPLACE_H
# define  TEST_LAPLACE_H

#include "cppunit/TestFixture.h"
#include "cppunit/extensions/HelperMacros.h"

#include "laplaceProblem.h"

/*
   Tests for public methods of the LaplaceProblem class
*/

   class Test_LaplaceProblem : public CppUnit::TestFixture
  {
   private:
      CPPUNIT_TEST_SUITE( Test_LaplaceProblem );

         CPPUNIT_TEST( test_setDirichletLeft  );
         CPPUNIT_TEST( test_setDirichletRight );

         CPPUNIT_TEST( test_constructRHS_equal_dirichlet_conditions );
         CPPUNIT_TEST( test_constructRHS_different_dirichlet_conditions );

         CPPUNIT_TEST( test_sparsity_vector );
         CPPUNIT_TEST( test_connection_vector );
         CPPUNIT_TEST( test_diagonal_coefficient_vector );
         CPPUNIT_TEST( test_offdiagonal_coefficient_vector );

         CPPUNIT_TEST( test_directSolve_gives_constant_for_equal_boundary_values );
         CPPUNIT_TEST( test_directSolve_gives_linear_for_different_boundary_values );

         CPPUNIT_TEST( test_jacobiTolerance );
         CPPUNIT_TEST( test_jacobiIteration );
         CPPUNIT_TEST( test_jacobiSolve_gives_constant_for_equal_boundary_values );
         CPPUNIT_TEST( test_jacobiSolve_gives_linear_for_different_boundary_values );

      CPPUNIT_TEST_SUITE_END();

      float floatError=10e-6;

   public:
      void setUp();
      void tearDown();

      void test_setDirichletLeft();
      void test_setDirichletRight();

      void test_constructRHS_equal_dirichlet_conditions();
      void test_constructRHS_different_dirichlet_conditions();

      void test_sparsity_vector();
      void test_connection_vector();
      void test_diagonal_coefficient_vector();
      void test_offdiagonal_coefficient_vector();

      void test_directSolve_gives_constant_for_equal_boundary_values();
      void test_directSolve_gives_linear_for_different_boundary_values();

      void test_jacobiTolerance();
      void test_jacobiIteration();
      void test_jacobiSolve_gives_constant_for_equal_boundary_values();
      void test_jacobiSolve_gives_linear_for_different_boundary_values();

  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_LaplaceProblem );

# endif
