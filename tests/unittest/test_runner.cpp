/*
 * Execute this program to run the test suite.
 *
 * Usually there's no need to change anything in here.
 */
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include "../../src/arena.hpp"

using namespace CppUnit;

int new_forest_counter = 0;
int delete_forest_counter = 0;
int recombination_counter = 0; // DEBUG
double recomb_opp = 0; // DEBUG
int node_created = 0;
int node_deleted = 0;
/*!
 * Global variable for the memory arena
 */
class Arena* Arena::globalArena;


int main(void) {
	TestResult controller;

	TestResultCollector result;
	controller.addListener(&result);

	BriefTestProgressListener progress;
	controller.addListener(&progress);

	TestRunner runner;
	runner.addTest( TestFactoryRegistry::getRegistry().makeTest() );
	runner.run(controller);

	CompilerOutputter outputter(&result, std::cerr);
	outputter.write();

	return result.wasSuccessful() ? 0 : 1;
}

// EOF
