// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/bootcamp/BootCampMover.cc
/// @brief Perform iterative rounds of packing, minimization, and MonteCarlo based search to find a minimized structure for a protein
/// @author PowersPope (apowers4@uoregon.edu)

// Unit headers
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMoverCreator.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pose/variant_util.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/bootcamp/fold_tree_from_ss.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "protocols.bootcamp.BootCampMover" );

namespace protocols {
namespace bootcamp {

	/////////////////////
	/// Constructors  ///
	/////////////////////

/// @brief Default constructor
BootCampMover::BootCampMover():
	protocols::moves::Mover( BootCampMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BootCampMover::~BootCampMover(){}

////////////////////////////////////////////////////////////////////////////////
	/// Mover Methods ///
	/////////////////////

/// @brief Take in a pose and apply iterative rounds of random Phi/Psi 
/// perturbations, which are followed by packing, minimization, and then MonteCarlo
/// boltzmann movement and pass/fail progression.
void
BootCampMover::apply( core::pose::Pose& mypose ){
	// Init our scorefunction and reweight our score
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function("ref2015"); 
	scorefxn->set_weight( core::scoring::ScoreType::linear_chainbreak, 1.0 );
	// Set our random number generator
	numeric::random::rg().set_seed("mt19937", 192); 
	// Setup our Movemap and minimization method
	core::kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_chi( true );
	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
	core::optimization::AtomTreeMinimizer atm;

	// Generate our dssp sequence
	core::kinematics::FoldTree ft = protocols::bootcamp::fold_tree_from_ss( mypose );
	mypose.fold_tree( ft );
	correctly_add_cutpoint_variants( mypose );

	// Setup our MonteCarlo object
	protocols::moves::MonteCarlo  mc( mypose, *scorefxn, 0.8 ); // 100 temp before

	// Copy the pose and init values
	core::pose::Pose copy_pose;
	bool bool_out;
	core::Real accepted_moves = 0.0;
	core::Real score_overtime = 0.0;
	int steps = 4;

	// Monte Carlo change residues 
	for ( int i=0; i<steps; i++ ) {
		core::Size size_pose = mypose.size();
		core::Size randres = ( numeric::random::rg().uniform() * size_pose + 1 );
		if ( mypose.residue( randres ).is_protein() ) {
			// Determine our perterbation amounts
			core::Real pert1 = numeric::random::rg().gaussian();
			core::Real pert2 = numeric::random::rg().gaussian();
			core::Real orig_phi = mypose.phi( randres );
			core::Real orig_psi = mypose.psi( randres );
			mypose.set_phi( randres, orig_phi + pert1 );
			mypose.set_psi( randres, orig_psi + pert2 );

			// Pack our pose
			core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( mypose );
			repack_task->restrict_to_repacking();
			core::pack::pack_rotamers( mypose, *scorefxn, repack_task );

			// Copy and Minimize the pose
			copy_pose = mypose;
			atm.run( copy_pose, mm, *scorefxn, min_opts );
			mypose = copy_pose;
			bool_out = mc.boltzmann( mypose );
			accepted_moves += bool_out;
			score_overtime += mc.last_score();
			TR << "MonteCarlo Round: " << i << " Score of Move: " << mc.last_score() << " Move Accepted: " << bool_out << std::endl;
		}
	}

	TR << "Final Score output: " << mc.last_accepted_score() << std::endl;
	TR << "Proportion of accepted moves: " << accepted_moves / steps << std::endl;
	TR << "Average Score of protein: " << score_overtime / steps << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
BootCampMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
BootCampMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void BootCampMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Perform iterative rounds of packing, minimization, and MonteCarlo based search to find a minimized structure for a protein", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::fresh_instance() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::clone() const
{
	return utility::pointer::make_shared< BootCampMover >( *this );
}

std::string BootCampMover::get_name() const {
	return mover_name();
}

std::string BootCampMover::mover_name() {
	return "BootCampMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
BootCampMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

std::string
BootCampMoverCreator::keyname() const
{
	return BootCampMover::mover_name();
}

void BootCampMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BootCampMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns PowersPope as its author.
void
BootCampMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"BootCampMover", basic::citation_manager::CitedModuleType::Mover,
		"PowersPope",
		"University of Oregon & The Flatiron Institute",
		"apowers4@uoregon.edu",
		"Wrote the BootCampMover."
		)
	);
}


////////////////////////////////////////////////////////////////////////////////
	/// private methods ///
	///////////////////////


std::ostream &
operator<<( std::ostream & os, BootCampMover const & mover )
{
	mover.show(os);
	return os;
}


} //bootcamp
} //protocols
