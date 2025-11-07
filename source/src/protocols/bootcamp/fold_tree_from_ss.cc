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

// @brief Calculate where the secondary structures are within our passed in string
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

// @brief This takes in a pose and returns a FoldTree for the given pose.
// This returns three different values
// 1. A New Fold Tree (built by secondary structure elements)
// 2. Additionally a pair/tuple that contains a Loop as the first element
// 3. As a second element in the pair/tuple is the bool of the vector1.
std::pair< core::kinematics::FoldTree, std::pair< utility::vector1< protocols::loops::Loop >, utility::vector1< core::Size > > >
fold_tree_from_ss( core::pose::Pose inpose ) {
  // Init our return variable
  core::kinematics::FoldTree out_fold_tree;
  utility::vector1< protocols::loops::Loop > loops;
	utility::vector1< core::Size > disorder_sections; 

  // Initialize a dssp, and extract out our secondary structure sequence
  core::scoring::dssp::Dssp dssp( inpose );
  std::string dssp_string = dssp.get_dssp_secstruct();

  // Pass our dssp string to fold_tree_from_dssp_string to get our new FoldTree
  out_fold_tree = fold_tree_from_dssp_string( dssp_string );

  return out_fold_tree;
}

// @brief take in a dssp based string and return a FoldTree that can be passed to a pose.
core::kinematics::FoldTree
fold_tree_from_dssp_string( std::string const & in_dssp ) {
  core::kinematics::FoldTree ft;
  utility::vector1< protocols::loops::Loop > loops;
	utility::vector1< core::Size > disorder_sections(in_dssp.size(), 0); // A vector of amino acid size and filled with zeros

// 		std::cout << "Length of our DSSP String Input: " << in_dssp.size() << std::endl;
  // Pass our sequence to our original function to extract out a vector1< pair< Size, Size >>
  utility::vector1< std::pair< core::Size, core::Size >> vector_dssp_pairs = identify_secondary_structure_spans( in_dssp );

  // Iter through our ranges and grab the middle 
  utility::vector1< core::Size > middle_vector;
  for ( unsigned long int k=1; k<=vector_dssp_pairs.size(); k++ ) {
    middle_vector.push_back( determine_middle_residue( vector_dssp_pairs[k].first, vector_dssp_pairs[k].second ) );
  }

  // Store variables to go through our loop
  core::Size first_middle, previous_end, current_middle, current_start, current_end, jump;
  first_middle = middle_vector[1];
  previous_end = vector_dssp_pairs[1].second;
  jump = 1;

  // Form the initial edge, as our initial edge will always start at 1
  ft.add_edge( first_middle, 1, core::kinematics::Edge::PEPTIDE );
  ft.add_edge( first_middle, vector_dssp_pairs[1].second, core::kinematics::Edge::PEPTIDE );

  // Generate our FoldTree by iterating through our two vector1s
  core::Size skipped = 0;
  for ( unsigned long int j=2; j<=vector_dssp_pairs.size(); j++ ) {
// 			std::cout << "Values Out: " << vector_dssp_pairs[j].first << " " << vector_dssp_pairs[j].second << std::endl;
    // Assign our current values first
    current_middle = middle_vector[j];
    current_start = vector_dssp_pairs[j].first;
    current_end = vector_dssp_pairs[j].second;

    // Skip if the edge is only 1 residue long and add an additional residue to start for the next loop
    if ( current_start == current_end ) { skipped++; continue; }
    if ( skipped ) { current_start--; skipped = 0; }
    // If the the previous_end and current_start dont follow then the previous region that was skipped was diordered and we 
    // should assign that region as an edge and a jump as well.
    if ( previous_end + 1 != current_start ) {
      // First we add the jump from the previous foldtree to this new center
      core::Size disordered_middle = determine_middle_residue( previous_end + 1, current_start - 1);
      ft.add_edge( first_middle, disordered_middle, jump );
      ft.add_edge( disordered_middle, previous_end + 1, core::kinematics::Edge::PEPTIDE );
      ft.add_edge( disordered_middle, current_start - 1, core::kinematics::Edge::PEPTIDE );
      jump++;
    }

    // First we add the jump from the previous foldtree to this new center
    ft.add_edge( first_middle, current_middle, jump );
    ft.add_edge( current_middle, current_start, core::kinematics::Edge::PEPTIDE );
    if ( j == vector_dssp_pairs.size() ) {
      ft.add_edge( current_middle, in_dssp.size(), core::kinematics::Edge::PEPTIDE );
    } else {
      ft.add_edge( current_middle, current_end, core::kinematics::Edge::PEPTIDE );
    }

// 			previous_start = vector_dssp_pairs[j].first;
    previous_end = vector_dssp_pairs[j].second;
    jump++;
  }

// 		std::cout << "Fold Tree output: " << ft.to_string() << std::endl;
  return ft;
}

// @brief determine the middle residue of our range
core::Size determine_middle_residue( core::Size start, core::Size end ) {
  // Init our vector to hold this information
  utility::vector1< core::Size > resi_vector;
  int middle;

  // First make our vector to get the size and to index later on
  for ( core::Size j=start; j<=end; j++ ) {
    resi_vector.push_back(j);
  }
    
  // Determine if our size is even or not
  if ( resi_vector.size() % 2 == 0 ) {
    middle = std::floor( resi_vector.size() / 2 );
  } else {
    middle = std::floor( resi_vector.size() / 2 ) + 1;
  }
  core::Size anchor_residue = resi_vector[middle];
  return anchor_residue;
}

} // bootcamp
} // protocols

