# ifndef POISSON_H
# define POISSON_H

# include "laplaceProblem.h"

   class PoissonProblem : public LaplaceProblem
  {
   private:

      int    numberPointSources;
      int    pointSourcesSet;
      int   *pointSourceLocations;
      float *pointSourceStrengths;

   public:

      PoissonProblem();
     ~PoissonProblem();

      void setNPointSources( int ns );
      int  nPointSources();

      void  setPointSource( int location, float strength );
      float pointSource( int location );
  };

# endif
