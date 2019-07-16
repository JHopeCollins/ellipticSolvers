# include "cppunit/ui/text/TestRunner.h"
# include "cppunit/TestResult.h"

# include "test-LaplaceProblem.h"

   int main()
  {
      CppUnit::TextUi::TestRunner   runner;

      runner.addTest( Test_LaplaceProblem::suite() );

      bool wasSuccessful = runner.run( "", false );

      return !wasSuccessful;
  }
