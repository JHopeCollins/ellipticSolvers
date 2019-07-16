# ifndef  TEST_LINSPACE_H
# define  TEST_LINSPACE_H

#include "cppunit/TestFixture.h"
#include "cppunit/extensions/HelperMacros.h"

#include "linspace.h"

/*
   Tests for linspace function
*/

   class Test_linspace : public CppUnit::TestFixture
  {
   private:
      CPPUNIT_TEST_SUITE( Test_linspace );

         CPPUNIT_TEST( test_n_is_1_returns_lowerLimit );
         CPPUNIT_TEST( test_0to1_spacing );
         CPPUNIT_TEST( test_integer_spacing );
         CPPUNIT_TEST( test_0to0_returns0 );
         CPPUNIT_TEST( test_1to0_reverse_spacing );

      CPPUNIT_TEST_SUITE_END();

      float floatError=10e-6;

   public:
      void setUp();
      void tearDown();

      void test_n_is_1_returns_lowerLimit();
      void test_0to1_spacing();
      void test_integer_spacing();
      void test_0to0_returns0();
      void test_1to0_reverse_spacing();
  };

CPPUNIT_TEST_SUITE_REGISTRATION( Test_linspace );

# endif
