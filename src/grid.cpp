
# include "grid.h"
# include "linspace.h"

# include <assert.h>

   Grid::Grid()
  {
      n=0;
      leftBound =0;
      rightBound=0;
      x=nullptr;
      gridExists=false;
  }

   Grid::~Grid()
  {
      delete[] x;
      x=nullptr;
  }

   void Grid::setNx( int nx )
  {
      n=nx;
      return;
  }

   int Grid::nx()
  {
      return n;
  }

   void Grid::setBounds( float lb, float rb )
  {
      leftBound =lb;
      rightBound=rb;
      return;
  }

   void Grid::bounds( float &l, float &r )
  {
      l= leftBound;
      r=rightBound;
      return;
  }

   void Grid::discretise()
  {
      x=new float[n];
      linspace( n, leftBound, rightBound, x, deltax );
      gridExists=true;
  }

   float Grid::dx()
  {
      return deltax;
  }

   float Grid::operator[]( int i )
  {
      if( isDiscretised() )
     {
         assert( (i>-1) and (i<n) );
         return x[i];
     }
      else
     {
         return 1./0.;
     }
  }

   bool Grid::isDiscretised()
  {
      return gridExists;
  }
