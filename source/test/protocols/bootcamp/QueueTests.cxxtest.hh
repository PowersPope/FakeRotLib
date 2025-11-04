// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  t/e/s/t///p/r/o/t/o/c/o/l/s///b/o/o/t/c/a/m/p/QueueTests.cxxtest.hh
/// @brief  Building Queue tests during bootcamp
/// @author PowersPope (apowers4@uoregon.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/bootcamp/Queue.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("QueueTests");


class QueueTests : public CxxTest::TestSuite {
	//Define Variables

public:

	protocols::bootcamp::Queue queue;

	void setUp() {
		core_init();

	}

	void tearDown() {

	}



	void test_first() {
		TS_TRACE( "Running my first unit test!" );
		TS_ASSERT( true );
	}

	void test_create_queue() {
		TR << "Testing Create Queue..." << std::endl;
		protocols::bootcamp::Queue queue;

		TR << "Checking No elements are added initially..." << std::endl;
		TS_ASSERT( queue.is_empty() );

// 		TR << "Checking that is_queue returns false here..." << std::endl;
// 		TS_ASSERT( !queue.is_queue() );
	}


	void test_enqueue()	{
		TR <<  "Testing Enqueue functionality..." << std::endl;
		protocols::bootcamp::Queue queue;

		TR << "Check that size has increased..." << std::endl;
		core::Size size_before = queue.size();
		queue.enqueue( "blah" );
		TS_ASSERT( queue.size() == size_before + 1 );

		TR << "Check that a new element is added to the back..." << std::endl;
		queue.enqueue( "new" );

	}

	void test_dequeue() {
		TR << "Testing Dequeue Functionality..." << std::endl;
		protocols::bootcamp::Queue queue;
		queue.enqueue( "A" );
		queue.enqueue( "B" );
// 		core::Size expectedSize = 2;
		TS_ASSERT( queue.size() == 2 );

// 		TS_ASSERT( queue.is_queue() );

		// Now lets remove one of the elements to see if it gets removed correctly
		std::string out;
		out = queue.dequeue();
// 		core::Size removedSize = 1;
		TS_ASSERT( queue.size() == 1 );
		TS_ASSERT( out == "A" );

		// Now Remove the last element and check that the queue is empty
		out = queue.dequeue();
		TS_ASSERT( out == "B" );
		TS_ASSERT( queue.is_empty() );
	}


};
