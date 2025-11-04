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

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

/// Project headers
#include <core/types.hh>

// C++ headers

//Auto Headers


static basic::Tracer TR("protocols.bootcamp.FoldTreeFromSS.cxxtest");

// --------------- Test Class --------------- //

class FoldTreeFromSS : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------- Function Testing --------------------- //
	utility::vector1< std::pair< core::Size, core::Size > >
	identify_secondary_structure_spans( std::string const & ss_string )
	{
		utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
		core::Size strand_start = -1;
		for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
			if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
				if ( int( strand_start ) == -1 ) {
					strand_start = ii;
				} else if ( ss_string[ii] != ss_string[strand_start] ) {
					ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
					strand_start = ii;
				}
			} else {
				if ( int( strand_start ) != -1 ) {
					ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
					strand_start = -1;
				}
			}
		}
		if ( int( strand_start ) != -1 ) {
			// last residue was part of a ss-eleemnt                                                                                                                                
			ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
		}
		for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
			std::cout << "SS Element " << ii << " from residue "
				<< ss_boundaries[ ii ].first << " to "
				<< ss_boundaries[ ii ].second << std::endl;
		}
		return ss_boundaries;
	}

	// --------------- Test Cases --------------- //
	void test_case1() {
		TR << "Running First Test Case For SS Position Extraction (Should be 7 elements)..." << std::endl;
		utility::vector1< std::pair< core::Size, core::Size > > vector_of_ss_pairs, control_pairs;
		vector_of_ss_pairs = identify_secondary_structure_spans("   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH");

		// Generate our control pairs vector to check against
		control_pairs.push_back({4, 8});
		control_pairs.push_back({12, 19});
		control_pairs.push_back({22, 26});
		control_pairs.push_back({36, 41});
		control_pairs.push_back({45, 55});
		control_pairs.push_back({58, 62});
		control_pairs.push_back({65, 68});

		// Iterate through pairs and check to see that they mathc what we expect
		for (long unsigned int i; i<=vector_of_ss_pairs.size(); i++) {
			TS_ASSERT( vector_of_ss_pairs[i] == control_pairs[i] );
		}

	}

	void test_case2() {
		TR << "Running Second Test Case For SS Position Extraction (Should be 7 elements)..." << std::endl;
		utility::vector1< std::pair< core::Size, core::Size > > vector_of_ss_pairs, control_pairs;
		vector_of_ss_pairs = identify_secondary_structure_spans("HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH ");

		// Generate our control pairs vector to check against
		control_pairs.push_back({1, 7});
		control_pairs.push_back({11, 22});
		control_pairs.push_back({29, 40});
		control_pairs.push_back({41, 50});
		control_pairs.push_back({51, 57});
		control_pairs.push_back({59, 62});
		control_pairs.push_back({63, 65});

		// Iterate through pairs and check to see that they mathc what we expect
		for (long unsigned int i; i<=vector_of_ss_pairs.size(); i++) {
			TS_ASSERT( vector_of_ss_pairs[i] == control_pairs[i] );
		}

	}

	void test_case3() {
		TR << "Running Third Test Case For SS Position Extraction (Should be 9 elements)..." << std::endl;
		utility::vector1< std::pair< core::Size, core::Size > > vector_of_ss_pairs, control_pairs;
		vector_of_ss_pairs = identify_secondary_structure_spans("EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE");

		// Generate our control pairs vector to check against
		control_pairs.push_back({1, 9});
		control_pairs.push_back({11, 18});
		control_pairs.push_back({20, 28});
		control_pairs.push_back({30, 30});
		control_pairs.push_back({32, 36});
		control_pairs.push_back({38, 38});
		control_pairs.push_back({40, 40});
		control_pairs.push_back({42, 42});
		control_pairs.push_back({44, 51});

		// Iterate through pairs and check to see that they mathc what we expect
		for (long unsigned int i; i<=vector_of_ss_pairs.size(); i++) {
			TS_ASSERT( vector_of_ss_pairs[i] == control_pairs[i] );
		}

	}
};
