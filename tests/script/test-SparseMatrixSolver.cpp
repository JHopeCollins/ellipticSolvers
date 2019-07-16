# include "cppunit/ui/text/TestRunner.h"
# include "cppunit/TestResult.h"

# include "test-SparseMatrixSolver.h"

   int main()
  {
      CppUnit::TextUi::TestRunner   runner;

      runner.addTest( Test_SparseMatrixSolver::suite() );

      bool wasSuccessful = runner.run( "", false );

      return !wasSuccessful;
  }
