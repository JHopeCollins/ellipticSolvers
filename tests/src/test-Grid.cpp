
# include "test-Grid.h"
# include "linspace.h"

   void Test_Grid::setUp(){}

   void Test_Grid::tearDown(){}

   void Test_Grid::test_set_nx()
  {
      Grid  grid;
      int   nx=5;

      grid.setNx( nx );

      CPPUNIT_ASSERT_EQUAL( nx, grid.nx() );

      return;
  }

   void Test_Grid::test_set_bounds()
  {
      Grid  grid;
      float lBound=0;
      float rBound=1;
      float l,r;

      grid.setBounds( lBound, rBound );
      grid.bounds( l,r );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( lBound, l, floatError );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( rBound, r, floatError );

      return;
  }

   void Test_Grid::test_discretise_with_linspace()
  {
      Grid  grid;
      int   nx=5;
      float lBound=0;
      float rBound=1;
      float dx;
      float *x=new float[nx];

      grid.setNx( nx );
      grid.setBounds( lBound, rBound );

      CPPUNIT_ASSERT( !(grid.isDiscretised()) );
      grid.discretise();
      CPPUNIT_ASSERT(   grid.isDiscretised()  );

      linspace( nx, lBound, rBound, x, dx );

      for( int i=0; i<nx; i++ )
     {
         CPPUNIT_ASSERT_DOUBLES_EQUAL( x[i], grid[i], floatError );
     }

      CPPUNIT_ASSERT_DOUBLES_EQUAL( dx, grid.dx(), floatError );

      delete[] x;

      return;
  }
