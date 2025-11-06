// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/boocamp/FoldTreeFromSS.cxxtest.hh
/// @brief
/// @author Andrew Powers (apowers4@uoregon.edu)


// Test headers
#include <cxxtest/TestSuite.h>

#include <devel/init.hh>
#include <test/util/pose_funcs.hh>
#include <test/protocols/init_util.hh>
#include <test/core/init_util.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <test/util/pose_funcs.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMover.fwd.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/moves/Mover.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/deep_copy.hh>
#include <basic/Tracer.hh>

/// Project headers
#include <core/types.hh>

// C++ headers
#include <cmath>

//Auto Headers


static basic::Tracer TR("protocols.bootcamp.BootCampMover.cxxtest");

// --------------- Test Class --------------- //

class BootCampMover : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


	// Shared initialization goes here.
	void setUp() {
		protocols_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------- Function Testing --------------------- //
	

	// --------------- Test Cases --------------- //
	void test_mover() {
		// Grab a MoverOP instantiation and then grab an OP of our mover by using a dynamic pointer cast
		protocols::moves::MoverOP bootcamp( protocols::moves::MoverFactory::get_instance()->newMover( "BootCampMover" ) );
		protocols::bootcamp::BootCampMoverOP bcm_op( utility::pointer::dynamic_pointer_cast< protocols::bootcamp::BootCampMover > ( bootcamp ) );

		TS_ASSERT( bcm_op != nullptr );
	}

};
