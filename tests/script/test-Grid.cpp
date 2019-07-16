# include "cppunit/ui/text/TestRunner.h"
# include "cppunit/TestResult.h"

# include "test-Grid.h"

   int main()
  {
      CppUnit::TextUi::TestRunner   runner;

      runner.addTest( Test_Grid::suite() );

      bool wasSuccessful = runner.run( "", false );

      return !wasSuccessful;
  }
