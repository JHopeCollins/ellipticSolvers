# include "cppunit/ui/text/TestRunner.h"
# include "cppunit/TestResult.h"

# include "test-linspace.h"

   int main()
  {
      CppUnit::TextUi::TestRunner   runner;

      runner.addTest( Test_linspace::suite() );

      bool wasSuccessful = runner.run( "", false );

      return !wasSuccessful;
  }
