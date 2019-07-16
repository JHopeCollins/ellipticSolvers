# include <cmath>

   void linspace( int n, float lowerLimit, float upperLimit, float *x, float &dx )
  {
      if( n==1 or n==0 ){ x[0]=lowerLimit; dx=0.; return; }

      dx= ( upperLimit-lowerLimit )/( n-1 );

      for( int i=0; i<n; i++ )
     {
         x[i] = lowerLimit + i*dx;
     }

      dx=fabs(dx);

      return;
  }
