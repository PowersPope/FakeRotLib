// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <iostream>
#include <src/basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <src/core/pose/Pose.fwd.hh>
#include <src/core/pose/Pose.hh>
#include <src/core/import_pose/import_pose.hh>
#include <utility/pointer/owning_ptr.hh>
#include <src/core/scoring/ScoreFunctionFactory.hh>
#include <src/core/scoring/ScoreFunction.hh>
#include <src/numeric/random/random.hh>
#include <src/numeric/random/uniform.hh>
#include <src/protocols/moves/MonteCarlo.hh>
#include <src/core/conformation/Residue.hh>

int main( int argc, char ** argv) {
	// Init Rosetta as a whole
	devel::init( argc, argv );
	// Init our scorefunction
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function("ref2015"); 
	// Set our random number generator
	numeric::random::rg().set_seed("mt19937", 192); 


	// Check if we get a pdb
	utility::vector1< std::string > filenames = basic::options::option[basic::options::OptionKeys::in::file::s ].value();

	if ( filenames.size() > 0 ) {
		std::cout << "You entered: " << filenames[1] << " as the PDB file to be read" << std::endl;
	} else {
		std::cout << "You didn't provide a PDB file with the -in::file::s option" << std::endl;
		return 1;
	}

	// Load our file into a pose object
	core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );

	// Setup our MonteCarlo object
	protocols::moves::MonteCarlo  mc( *mypose, *scorefxn, 100 );

	// Score our pose
// 	core::Real score = scorefxn->score( * mypose );

	// Output our score
// 	std::cout << "File Score: " << score << std::endl;

	
	// Monte Carlo change residues
	for ( int i=0; i<10; i++ ) {
		core::Size size_pose = mypose->size();
		core::Size randres = ( numeric::random::rg().uniform() * size_pose + 1 );
		if ( mypose->residue( randres ).is_protein() ) {
			core::Real pert1 = numeric::random::rg().gaussian();
			core::Real pert2 = numeric::random::rg().gaussian();
			core::Real orig_phi = mypose->phi( randres );
			core::Real orig_psi = mypose->psi( randres );
			mypose->set_phi( randres, orig_phi + pert1 );
			mypose->set_psi( randres, orig_psi + pert2 );
			mc.boltzmann( * mypose );
			std::cout << "MC MOVEMENT!!! Round: " << i << std::endl;
		}
	}

// 	std::cout << "Random Res: " << randres << " Pert1: " << orig_phi + pert1 << std::endl;
// 	std::cout << "Random Res: " << randres << " Pert2: " << orig_psi + pert2 << std::endl;

	return 0;
}

