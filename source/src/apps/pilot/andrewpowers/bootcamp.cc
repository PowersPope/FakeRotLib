// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <basic/Tracer.hh>

#include <iostream>

#include <src/basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <protocols/jd2/JobDistributor.hh>
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
#include <src/core/pack/pack_rotamers.hh>
#include <src/core/pack/task/PackerTask.hh>
#include <src/core/pack/task/TaskFactory.hh>
#include <src/core/kinematics/MoveMap.hh>
#include <src/core/kinematics/FoldTree.hh>
#include <src/core/optimization/MinimizerOptions.hh>
#include <src/core/optimization/AtomTreeMinimizer.hh>
#include <core/pose/variant_util.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMover.fwd.hh>


static basic::Tracer TR( "apps.pilot.andrewpowers.bootcamp" );

int main( int argc, char ** argv) {
	// Init Rosetta as a whole
	devel::init( argc, argv );
	// Init our scorefunction
// 	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function("ref2015"); 
// 	scorefxn->set_weight( core::scoring::ScoreType::linear_chainbreak, 1.0 );
// 	// Set our random number generator
// 	numeric::random::rg().set_seed("mt19937", 192); 
// 	// Setup our Movemap
// 	core::kinematics::MoveMap mm;
// 	mm.set_bb( true );
// 	mm.set_chi( true );
// 	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
// 	core::optimization::AtomTreeMinimizer atm;
// 
// 
// 	// Check if we get a pdb
// 	utility::vector1< std::string > filenames = basic::options::option[basic::options::OptionKeys::in::file::s ].value();
// 
// 	if ( filenames.size() > 0 ) {
// 		TR << "You entered: " << filenames[1] << " as the PDB file to be read" << std::endl;
// 	} else {
// 		TR << "You didn't provide a PDB file with the -in::file::s option" << std::endl;
// 		return 1;
// 	}
// 
// 	// Load our file into a pose object
// 	core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
// 
// 	// Generate our dssp sequence
// 	core::kinematics::FoldTree ft = protocols::bootcamp::fold_tree_from_ss( *mypose );
// 	mypose->fold_tree( ft );
// 	correctly_add_cutpoint_variants( *mypose );
// 
// 	// Setup our MonteCarlo object
// 	protocols::moves::MonteCarlo  mc( *mypose, *scorefxn, 0.8 ); // 100 temp before
// 
// 	// Copy the pose and init values
// 	core::pose::Pose copy_pose;
// 	bool bool_out;
// 	core::Real accepted_moves = 0.0;
// 	core::Real score_overtime = 0.0;
// 	int steps = 4;
// 
// 	// Monte Carlo change residues 
// 	for ( int i=0; i<steps; i++ ) {
// 		core::Size size_pose = mypose->size();
// 		core::Size randres = ( numeric::random::rg().uniform() * size_pose + 1 );
// 		if ( mypose->residue( randres ).is_protein() ) {
// 			// Determine our perterbation amounts
// 			core::Real pert1 = numeric::random::rg().gaussian();
// 			core::Real pert2 = numeric::random::rg().gaussian();
// 			core::Real orig_phi = mypose->phi( randres );
// 			core::Real orig_psi = mypose->psi( randres );
// 			mypose->set_phi( randres, orig_phi + pert1 );
// 			mypose->set_psi( randres, orig_psi + pert2 );
// 
// 			// Pack our pose
// 			core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
// 			repack_task->restrict_to_repacking();
// 			core::pack::pack_rotamers( *mypose, *scorefxn, repack_task );
// 
// 			// Copy and Minimize the pose
// 			copy_pose = *mypose;
// 			atm.run( copy_pose, mm, *scorefxn, min_opts );
// 			*mypose = copy_pose;
// 			bool_out = mc.boltzmann( *mypose );
// 			accepted_moves += bool_out;
// 			score_overtime += mc.last_score();
// 			TR << "MonteCarlo Round: " << i << " Score of Move: " << mc.last_score() << " Move Accepted: " << bool_out << std::endl;
// 		}
// 	}
// 
// 	TR << "Final Score output: " << mc.last_accepted_score() << std::endl;
// 	TR << "Proportion of accepted moves: " << accepted_moves / steps << std::endl;
// 	TR << "Average Score of protein: " << score_overtime / steps << std::endl;
	protocols::bootcamp::BootCampMoverOP bootcamp( new protocols::bootcamp::BootCampMover );
	protocols::jd2::JobDistributor::get_instance()->go(bootcamp);

	return 0;
}

