// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/boocamp/FoldTreeContainerSS.cxxtest.hh
/// @brief
/// @author Andrew Powers (apowers4@uoregon.edu)

#ifndef INCLUDED_protocols_bootcamp_fold_tree_from_ss_hh
#define INCLUDED_protocols_bootcamp_fold_tree_from_ss_hh

// Main Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/bootcamp/FoldTreeContainerSS.hh>

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

// @brief Calculate where the secondary structures are within our passed in string
utility::vector1< std::pair< core::Size, core::Size > >
identify_secondary_structure_spans( std::string const & ss_string );

// @brief This takes in a pose and returns a FoldTree for the given pose.
FoldTreeFromSS fold_tree_from_ss( core::pose::Pose inpose );

// @brief take in a dssp based string and return a FoldTree that can be passed to a pose.
FoldTreeFromSS fold_tree_from_dssp_string( std::string const & in_dssp );

// @brief determine the middle residue of our range
core::Size determine_middle_residue( core::Size start, core::Size end );

} // bootcamp
} // protocols


#endif
