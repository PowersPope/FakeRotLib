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
#include <protocols/loops/Loop.hh>
#include <core/scoring/dssp/Dssp.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

/// Project headers
#include <core/types.hh>

// C++ headers
#include <cmath>
#include <iostream>

//Auto Headers

namespace protocols {
namespace bootcamp {

class FoldTreeFromSS {
	public:
// 		FoldTreeFromSS( std::string const & ssstring );
		FoldTreeFromSS( core::Size const & stringSize );
		~FoldTreeFromSS();

		core::kinematics::FoldTree const & fold_tree() const;

		protocols::loops::Loop const & loop( core::Size index ) const;

		core::Size loop_for_residue( core::Size seqpose ) const;

		void add_loop_to_vector( 
				protocols::loops::Loop const & loop, 
				core::Size const & index 
				);

		void add_reference_to_loop_for_residue( 
				core::Size const & ref, 
				core::Sixe index
				);

		void add_fold_tree( core::kinematics::FoldTree );

	private:
		core::kinematics::FoldTree ft_;
		utility::vector1< protocols::loops::Loop > loop_vector_;
		utility::vector1< core::Size > loop_for_residue_;
};

} // bootcamp
} // protocols
