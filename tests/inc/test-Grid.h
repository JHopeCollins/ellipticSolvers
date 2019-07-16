# ifndef  TEST_GRID_H
# define  TEST_GRID_H

#include "cppunit/TestFixture.h"
#include "cppunit/extensions/HelperMacros.h"

#include "grid.h"

/*
   Tests for public methods of the Grid class
*/

   class Test_Grid : public CppUnit::TestFixture
  {
   private:
      CPPUNIT_TEST_SUITE( Test_Grid );

         CPPUNIT_TEST( test_set_nx );
         CPPUNIT_TEST( test_set_bounds );
         CPPUNIT_TEST( test_discretise_with_linspace );

      CPPUNIT_TEST_SUITE_END();

      float floatError=10e-6;

   public:
      void setUp();
      void tearDown();

      void test_set_nx();
      void test_set_bounds();
      void test_discretise_with_linspace();
  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_Grid );

# endif
