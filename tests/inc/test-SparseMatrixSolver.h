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

         CPPUNIT_TEST( test_SparseMatrixSolver );

      CPPUNIT_TEST_SUITE_END();

   public:
      void setUp();
      void tearDown();

      void test_SparseMatrixSolver();
  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_SparseMatrixSolver );

# endif
