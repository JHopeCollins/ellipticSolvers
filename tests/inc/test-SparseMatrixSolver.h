# ifndef  TEST_SparseMatrixSolver_H
# define  TEST_SparseMatrixSolver_H

#include "cppunit/TestFixture.h"
#include "cppunit/extensions/HelperMacros.h"

#include "sparseMatrixSolver.h"

/*
   Tests for public methods of the SparseMatrixSolver class
*/

   class Test_SparseMatrixSolver : public CppUnit::TestFixture
  {
   private:
      CPPUNIT_TEST_SUITE( Test_SparseMatrixSolver );

         CPPUNIT_TEST( test_mallocVectors );
         CPPUNIT_TEST( test_constructFullMatrix );
         CPPUNIT_TEST( test_directSolve );
         CPPUNIT_TEST( test_jacobiIteration );
         CPPUNIT_TEST( test_jacobiSolve );

      CPPUNIT_TEST_SUITE_END();

      float floatError=10e-5;
      float jacobiError=10e-4;

   public:
      void setUp();
      void tearDown();

      void test_mallocVectors();
      void test_constructFullMatrix();
      void test_directSolve();
      void test_jacobiIteration();
      void test_jacobiSolve();
  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_SparseMatrixSolver );

# endif
