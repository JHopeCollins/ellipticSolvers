#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/ui/text/TestRunner.h"
#include "cppunit/TestResult.h"

   int main()
  {
      CppUnit::TextUi::TestRunner   runner;

      CppUnit::TestFactoryRegistry  &registry=
         CppUnit::TestFactoryRegistry::getRegistry();

      runner.addTest( registry.makeTest() );

      bool wasSuccessful = runner.run( "", false );

      return !wasSuccessful;
  }
