// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/boocamp/BootCampMover.cxxtest.hh
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
// #include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <test/util/pose_funcs.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMover.fwd.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/moves/Mover.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/deep_copy.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

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

	// @brief Check to make sure that our class when loaded doesnt make an empty dynamic pointer
	void test_mover() {
		// Grab a MoverOP instantiation and then grab an OP of our mover by using a dynamic pointer cast
		protocols::moves::MoverOP bootcamp( protocols::moves::MoverFactory::get_instance()->newMover( "BootCampMover" ) );
		protocols::bootcamp::BootCampMoverOP bcm_op( utility::pointer::dynamic_pointer_cast< protocols::bootcamp::BootCampMover > ( bootcamp ) );

		// Final test for non nullptr match
		TS_ASSERT( bcm_op != nullptr );
	}

	// @brief This test is specific for making sure that our Mover's setters and getters are working correctly
	void test_get_setters() {
		// Grab a MoverOP instantiation and then grab an OP of our mover by using a dynamic pointer cast
		protocols::moves::MoverOP bootcamp( protocols::moves::MoverFactory::get_instance()->newMover( "BootCampMover" ) );
		protocols::bootcamp::BootCampMoverOP bcm_op( utility::pointer::dynamic_pointer_cast< protocols::bootcamp::BootCampMover > ( bootcamp ) );

		// Set our values to check against
		core::scoring::ScoreFunctionOP score_check( core::scoring::get_score_function( "ref2015" ) );
		core::Size check_iteration( 5 );
		// Set our variables of interest for testing
		bcm_op->set_num_iterations( check_iteration );
		bcm_op->set_scorefxn( score_check );

		// Final equal assert
		TS_ASSERT_EQUALS( bcm_op->get_scorefxn(), score_check );
		TS_ASSERT_EQUALS( bcm_op->get_num_iterations(), check_iteration );
	}

	// @brief Test our XML tag reader for the number of iterations of our MonteCarlo search and set the ScoreFunction
	void test_tag_readers() {
		// init variables
		std::string mover_name( "BootCampMover" );

		// Instantiate our values and datamap
		core::scoring::ScoreFunctionOP sfxn	= core::scoring::get_score_function( "ref2015" );
		basic::datacache::DataMap data;
		data.add( "scorefxns", "testing123", sfxn );

		// Init tag objects to be used in the datamap
		utility::tag::TagOP tag ( utility::pointer::make_shared< utility::tag::Tag >() );
		tag->setName( mover_name );
		tag->setOption< core::Size >( "num_iterations", 4 );
		tag->setOption< std::string >( "scorefxn", "testing123" );

		// Now pass our datamap and tags to our mover
		protocols::moves::MoverOP bootcamp = protocols::moves::MoverFactory::get_instance()->newMover( tag, data );

		// Grab a MoverOP instantiation and then grab an OP of our mover by using a dynamic pointer cast
// 		protocols::moves::MoverOP bootcamp( protocols::moves::MoverFactory::get_instance()->newMover( my_configured_mover ) );
		protocols::bootcamp::BootCampMoverOP bcm_op( utility::pointer::dynamic_pointer_cast< protocols::bootcamp::BootCampMover > ( bootcamp ) );
		TS_ASSERT( bcm_op );

		// Final Assert to check that our scorefunction and passed in variable works
		TS_ASSERT_EQUALS( bcm_op->get_num_iterations(), 4 );
		core::scoring::ScoreFunctionOP stored_scorefxn( bcm_op->get_scorefxn() );
		TS_ASSERT( (*stored_scorefxn).weights() == sfxn->weights() );
	}



};
