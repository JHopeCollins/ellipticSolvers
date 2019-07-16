
# include "test-linspace.h"

   void Test_linspace::setUp(){}

   void Test_linspace::tearDown(){}

   void Test_linspace::test_n_is_1_returns_lowerLimit()
  {
      int   n=1;
      float l=6;
      float r=1;
      float  dx;
      float *x=new float[n];

      linspace( n, l, r, x, dx );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( l, x[0], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, dx,   floatError );

      delete[] x;

      return;
  }

   void Test_linspace::test_0to1_spacing()
  {
      int   n=5;
      float l=0;
      float r=1;
      float  dx;
      float *x=new float[n];

      linspace( n, l, r, x, dx );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( l,    x[0], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, x[1], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.50, x[2], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.75, x[3], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( r,    x[4], floatError );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.25, dx, floatError );

      delete[] x;

      return;
  }

   void Test_linspace::test_integer_spacing()
  {
      int   n=5;
      float l=0;
      float r=4;
      float  dx;
      float *x=new float[n];

      linspace( n, l, r, x, dx );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( l,  x[0], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., x[1], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 2., x[2], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 3., x[3], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( r,  x[4], floatError );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 1., dx, floatError );

      delete[] x;

      return;
  }

   void Test_linspace::test_0to0_returns0()
  {
      int   n=3;
      float l=0;
      float r=0;
      float  dx;
      float *x=new float[n];

      linspace( n, l, r, x, dx );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0., x[0], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0., x[1], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0., x[2], floatError );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0., dx, floatError );

      delete[] x;

      return;
  }

   void Test_linspace::test_1to0_reverse_spacing()
  {
      int   n= 6;
      float l= 1;
      float r=-1;
      float   dx;
      float *x=new float[n];

      linspace( n, l, r, x, dx );

      CPPUNIT_ASSERT_DOUBLES_EQUAL(  l,   x[0], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(  0.6, x[1], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(  0.2, x[2], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( -0.2, x[3], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( -0.6, x[4], floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(  r,   x[5], floatError );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.4, dx, floatError );

      delete[] x;

      return;
  }
