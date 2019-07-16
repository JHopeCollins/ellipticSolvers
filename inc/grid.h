# ifndef GRID_H
# define GRID_H

   class Grid
  {
   private:

      int    n;            // number of elements
      float  leftBound;    // lower boundary coordinate
      float  rightBound;   // upper boundary coordinate
      float  deltax;       // grid spacing
      float *x;            // vector of coordinates
      bool   gridExists;   // has grid been discretised?

   public:

      Grid();
     ~Grid();

      void setNx( int nx );
      int  nx();

      void setBounds( float lb, float rb );
      void bounds(    float &l, float &r );

      void discretise();
      float dx();

      float operator[]( int i );
      bool isDiscretised();
  };

# endif
