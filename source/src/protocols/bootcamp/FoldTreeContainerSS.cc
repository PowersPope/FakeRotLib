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

// Main Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

/// Project headers
#include <core/types.hh>

// C++ headers
#include <iostream>
#include <cmath>

//Auto Headers


namespace protocols {
namespace bootcamp {

class FoldTreeFromSS {
	public:
		// @brief Construct our class
// 		FoldTreeFromSS( std::string const & ssstring );
		FoldTreeFromSS( core::Size const & stringSize ) : 
		loop_for_residue_( stringSize, 0 ) {};
		// @brief Deconstructor
		~FoldTreeFromSS(){};

		//@brief return the foldtree that was computed on start up
		core::kinematics::FoldTree const & fold_tree() const { return ft_; }

		// @brief grab the Loop constructed for passed in resi index
		protocols::loops::Loop const & loop( core::Size index ) const { return loop_vector_[ index ]; }

		// @brief return the correct index for the loop vector given a specific 
		// sequence position passed in
		core::Size loop_for_residue( core::Size seqpose ) const {
			return loop_for_residue_[ seqpose ];
		}

		// @brief add and an element with the correct index for our loop
		// ss struct information that is available
		void add_loop_to_vector( protocols::loops::Loop const & loop, core::Size const & index ) {
			if ( loop_vector_.size() + 1 == index ) { loop_vector_.push_back( loop ); }
		}

		// @brief add element to loop_for_residue_
		void add_reference_to_loop_for_residue( core::Size const & ref, core::Size index ) { loop_for_residue_[ index ] = ref; }

		// @brief add a fold tree
		void add_fold_tree( core::kinematics::FoldTree ft ) { ft_ = ft; }

	private:
		core::kinematics::FoldTree ft_; // FoldTree output we produce and can use on pose
		utility::vector1< protocols::loops::Loop > loop_vector_; // Vector1 that holds our Loop objects and can be accessed by indexing
		utility::vector1< core::Size > loop_for_residue_; // The indices needed to match up with our loop_vector_, so you grab the right Loop mover.
};

} // bootcamp
} // protocols

