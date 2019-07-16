# include "poissonProblem.h"

# include <assert.h>

   PoissonProblem::PoissonProblem()
  {
      numberPointSources=0;
      pointSourcesSet=0;
      pointSourceLocations=nullptr;
      pointSourceStrengths=nullptr;
  }

   PoissonProblem::~PoissonProblem()
  {
      delete[] pointSourceStrengths;
      delete[] pointSourceLocations;
      pointSourceLocations=nullptr;
      pointSourceStrengths=nullptr;
  }

   void PoissonProblem::setNPointSources( int ns )
  {
      numberPointSources=ns;
      pointSourceLocations = new   int[numberPointSources];
      pointSourceStrengths = new float[numberPointSources];

      for( int i=0; i<numberPointSources; i++ )
     {
         pointSourceLocations[i]=-1;
     }
  }

   int PoissonProblem::nPointSources()
  {
      return numberPointSources;
  }

   void PoissonProblem::setPointSource( int location, float strength )
  {
      assert( pointSourcesSet<numberPointSources );
      pointSourceLocations[ pointSourcesSet ] = location;
      pointSourceStrengths[ pointSourcesSet ] = strength;
      pointSourcesSet++;
  }

   float PoissonProblem::pointSource( int location )
  {
      float strength=0;

      for( int i=0; i<pointSourcesSet; i++ )
     {
         if( pointSourceLocations[i]==location )
        {
            strength=pointSourceStrengths[i];
        }
     }
      return strength;
  }
