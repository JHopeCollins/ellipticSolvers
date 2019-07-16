# include "cppunit/ui/text/TestRunner.h"
# include "cppunit/TestResult.h"

# include "test-PoissonProblem.h"

   int main()
  {
      CppUnit::TextUi::TestRunner   runner;

      runner.addTest( Test_PoissonProblem::suite() );

      bool wasSuccessful = runner.run( "", false );

      return !wasSuccessful;
  }
