/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin;

    // !!! PHYSIMESS CODE BLOCK START !!! //
    // load cells from your CSV file (if enabled)
    load_cells_from_pugixml();

    // new fibre related parameters and bools
    bool isFibreFromFile = false;
    bool fibreanisotropy = parameters.bools("anisotropic_fibres");
    double fibre_length = parameters.doubles("fibre_length");
    double length_normdist_sd = parameters.doubles("length_normdist_sd");
    double fibre_radius = parameters.doubles("fibre_radius");
    double fibre_angle = parameters.doubles("fibre_angle");
    double angle_normdist_sd = parameters.doubles("angle_normdist_sd");

    for( int i=0; i < (*all_cells).size(); i++ ){

        // initialise the following parameters for all cells regardless of type
        (*all_cells)[i]->parameters.mCellVelocityMaximum = parameters.doubles("cell_velocity_max");
        (*all_cells)[i]->parameters.mVelocityAdhesion = parameters.doubles("vel_adhesion");
        (*all_cells)[i]->parameters.mVelocityContact = parameters.doubles("vel_contact");


        (*all_cells)[i]->parameters.fibreDegradationRate = parameters.doubles("fibre_deg_rate");
        (*all_cells)[i]->parameters.stuck_threshold = parameters.doubles("fibre_stuck");
        (*all_cells)[i]->parameters.fibre_degradation = parameters.bools("fibre_degradation");

        (*all_cells)[i]->parameters.fibre_pushing = parameters.bools("fibre_pushing");
        (*all_cells)[i]->parameters.fibre_rotation = parameters.bools("fibre_rotation");
        (*all_cells)[i]->parameters.mFibreStickiness = parameters.doubles("fibre_sticky");

        (*all_cells)[i]->state.crosslink_point.resize(3,0.0);

        const auto agentname = std::string((*all_cells)[i]->type_name);
        const auto ecm = std::string("ecm");
        const auto matrix = std::string("matrix");
        const auto fiber = std::string("fiber");
        const auto fibre = std::string("fibre");
        const auto rod = std::string("rod");

        if (agentname.find(ecm) != std::string::npos ||
            agentname.find(matrix) != std::string::npos ||
            agentname.find(fiber) != std::string::npos ||
            agentname.find(fibre) != std::string::npos ||
            agentname.find(rod) != std::string::npos) {
            /* fibre positions are given by csv
               assign fibre orientation and test whether out of bounds */
            isFibreFromFile = true;

            (*all_cells)[i]->parameters.mLength = NormalRandom(fibre_length, length_normdist_sd) / 2.0;
            (*all_cells)[i]->parameters.mRadius = fibre_radius;

            //assign fibre orientation as a random vector from points on unit sphere/circle
            (*all_cells)[i]->assign_orientation();
            if (default_microenvironment_options.simulate_2D) {
                if (fibreanisotropy){
                    double theta = NormalRandom(fibre_angle,angle_normdist_sd);
                    (*all_cells)[i]->state.orientation[0] = cos(theta);
                    (*all_cells)[i]->state.orientation[1] = sin(theta);
                }
                else{
                    (*all_cells)[i]->state.orientation = UniformOnUnitCircle();
                }
                (*all_cells)[i]->state.orientation[2] = 0.0;
            }
            else {
                (*all_cells)[i]->state.orientation = UniformOnUnitSphere();
            }

            //###########################################//
            //   this bit a hack for PacMan and maze	 //
            //###########################################//
            if ((*all_cells)[i]->type_name == "fibre_vertical") {
                (*all_cells)[i]->state.orientation[0] = 0.0;
                (*all_cells)[i]->state.orientation[1] = 1.0;
                (*all_cells)[i]->state.orientation[2] = 0.0;
            }
            if ((*all_cells)[i]->type_name == "fibre_horizontal") {
                (*all_cells)[i]->state.orientation[0] = 1.0;
                (*all_cells)[i]->state.orientation[1] = 0.0;
                (*all_cells)[i]->state.orientation[2] = 0.0;
            }
            //###########################################//

            // start and end points of a fibre are calculated from fibre center
            double xs = (*all_cells)[i]->position[0] -
                        (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[0];
            double xe = (*all_cells)[i]->position[0] +
                        (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[0];
            double ys = (*all_cells)[i]->position[1] -
                        (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[1];
            double ye = (*all_cells)[i]->position[1] +
                        (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[1];
            double zs = 0.0;
            double ze = 0.0;
            if (default_microenvironment_options.simulate_2D) {
                /*std::cout << " fibre endpoints in 2D are " << xs << " " << ys <<
                               " and " << xe << " " << ye << std::endl; */
            }
            else if (!default_microenvironment_options.simulate_2D) {
                zs = (*all_cells)[i]->position[2] -
                     (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[2];
                ze = (*all_cells)[i]->position[2] +
                     (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[2];
                /*std::cout << " fibre endpoints in 3D are " << xs << " " << ys << " " << zs <<
                               " and " << xe << " " << ye << " " << ze << std::endl; */
            }

            /* check whether a fibre end point leaves the domain and if so initialise fibre again
                     assume user placed the centre of fibre within the domain so reinitialise orientation,
                     break after 10 failures
                     It needs re-writing at some stage to handle the 3D case properly */

            if (fibreanisotropy) {
                if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                    ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax) {
                    (*all_cells)[i]->parameters.fail_count = 10;
                }
            }
            else{
                if (default_microenvironment_options.simulate_2D) {
                    while ((*all_cells)[i]->parameters.fail_count < 10) {
                        if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                            ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax) {
                            (*all_cells)[i]->parameters.fail_count++;
                            (*all_cells)[i]->state.orientation = UniformOnUnitCircle();
                            xs = (*all_cells)[i]->position[0] -
                                 (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[0];
                            xe = (*all_cells)[i]->position[0] +
                                 (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[0];
                            ys = (*all_cells)[i]->position[1] -
                                 (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[1];
                            ye = (*all_cells)[i]->position[1] +
                                 (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[1];
                        }
                        else {
                            break;
                        }
                    }
                }

                if (!default_microenvironment_options.simulate_2D) {
                    while ((*all_cells)[i]->parameters.fail_count < 10) {
                        if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                            ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax ||
                            zs < Zmin || ze > Zmax || ze < Xmin || zs > Xmax) {
                            (*all_cells)[i]->parameters.fail_count++;
                            (*all_cells)[i]->state.orientation = UniformOnUnitSphere();
                            xs = (*all_cells)[i]->position[0] -
                                 (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[0];
                            xe = (*all_cells)[i]->position[0] +
                                 (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[0];
                            ys = (*all_cells)[i]->position[1] -
                                 (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[1];
                            ye = (*all_cells)[i]->position[1] +
                                 (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[1];
                            zs = (*all_cells)[i]->position[2] -
                                 (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[2];
                            ze = (*all_cells)[i]->position[2] +
                                 (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[2];
                        }
                        else {
                            break;
                        }
                    }
                }
            }

            // relabel so that the rest of the code works (HACK)
            (*all_cells)[i]->type_name = "fibre";

        }
        else
        {
            // type is a normal cell
        }
    }

    /* agents have not been added from the file but do want them
       create some of each agent type */

    if(!isFibreFromFile){
        Cell* pC;
        std::vector<double> position = {0, 0, 0};

        for( int k=0; k < cell_definitions_by_index.size() ; k++ ) {

            Cell_Definition *pCD = cell_definitions_by_index[k];
            std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
            
            const auto agentname = std::string(pCD->name);
            const auto ecm = std::string("ecm");
            const auto matrix = std::string("matrix");
            const auto fiber = std::string("fiber");
            const auto fibre = std::string("fibre");
            const auto rod = std::string("rod");

            if (agentname.find(ecm) == std::string::npos &&
                agentname.find(matrix) == std::string::npos &&
                agentname.find(fiber) == std::string::npos &&
                agentname.find(fibre) == std::string::npos &&
                agentname.find(rod) == std::string::npos){
                for (int n = 0; n < parameters.ints("number_of_cells"); n++) {

                    position[0] = Xmin + UniformRandom() * Xrange;
                    position[1] = Ymin + UniformRandom() * Yrange;
                    position[2] = Zmin + UniformRandom() * Zrange;

                    pC = create_cell(*pCD);

                    pC->parameters.mCellVelocityMaximum = parameters.doubles("cell_velocity_max");
                    pC->parameters.mVelocityAdhesion = parameters.doubles("vel_adhesion");
                    pC->parameters.mVelocityContact = parameters.doubles("vel_contact");

                    pC->parameters.fibreDegradationRate = parameters.doubles("fibre_deg_rate");
                    pC->parameters.stuck_threshold = parameters.doubles("fibre_stuck");
                    pC->parameters.fibre_degradation = parameters.bools("fibre_degradation");

                    pC->parameters.fibre_pushing = parameters.bools("fibre_pushing");
                    pC->parameters.fibre_rotation = parameters.bools("fibre_rotation");
                    pC->parameters.mFibreStickiness = parameters.doubles("fibre_sticky");

                    pC->state.crosslink_point.resize(3, 0.0);

                    pC->assign_position(position);
                }
            }

            if (agentname.find(ecm) != std::string::npos ||
                agentname.find(matrix) != std::string::npos ||
                agentname.find(fiber) != std::string::npos ||
                agentname.find(fibre) != std::string::npos ||
                agentname.find(rod) != std::string::npos){
                for ( int nf = 0 ; nf < parameters.ints("number_of_fibres") ; nf++ ) {

                    position[0] = Xmin + UniformRandom() * Xrange;
                    position[1] = Ymin + UniformRandom() * Yrange;
                    position[2] = Zmin + UniformRandom() * Zrange;

                    pC = create_cell(*pCD);

                    pC->parameters.mCellVelocityMaximum = parameters.doubles("cell_velocity_max");
                    pC->parameters.mVelocityAdhesion = parameters.doubles("vel_adhesion");
                    pC->parameters.mVelocityContact = parameters.doubles("vel_contact");

                    pC->parameters.fibreDegradationRate = parameters.doubles("fibre_deg_rate");
                    pC->parameters.stuck_threshold = parameters.doubles("fibre_stuck");
                    pC->parameters.fibre_degradation = parameters.bools("fibre_degradation");

                    pC->parameters.fibre_pushing = parameters.bools("fibre_pushing");
                    pC->parameters.fibre_rotation = parameters.bools("fibre_rotation");
                    pC->parameters.mFibreStickiness = parameters.doubles("fibre_sticky");

                    pC->state.crosslink_point.resize(3,0.0);

                    pC->parameters.mLength = NormalRandom(fibre_length, length_normdist_sd) / 2.0;
                    pC->parameters.mRadius = fibre_radius;

                    //assign fibre orientation as a random vector from points on unit sphere/circle.
                    pC->assign_orientation();
                    if( default_microenvironment_options.simulate_2D) {
                        if (fibreanisotropy) {
                            double theta = NormalRandom(fibre_angle, angle_normdist_sd);
                            pC->state.orientation[0] = cos(theta);
                            pC->state.orientation[1] = sin(theta);
                        }
                        else {
                            pC->state.orientation = UniformOnUnitCircle();
                        }
                        pC->state.orientation[2] = 0.0;
                    }
                    else{
                        pC->state.orientation = UniformOnUnitSphere();
                    }

                    // start and end points of a fibre are calculated from fibre center
                    double xs = position[0] - pC->parameters.mLength*pC->state.orientation[0];
                    double xe = position[0] + pC->parameters.mLength*pC->state.orientation[0];
                    double ys = position[1] - pC->parameters.mLength*pC->state.orientation[1];
                    double ye = position[1] + pC->parameters.mLength*pC->state.orientation[1];
                    double zs = 0.0;
                    double ze = 0.0;
                    if( !default_microenvironment_options.simulate_2D) {
                        zs = position[2] - pC->parameters.mLength * pC->state.orientation[2];
                        ze = position[2] + pC->parameters.mLength * pC->state.orientation[2];
                    }

                    /* check whether a fibre end point leaves the domain and if so initialise fibre again
                     assume user placed the centre of fibre within the domain so reinitialise orientation,
                     break after 10 failures
                     It needs re-writing at some stage to handle the 3D case properly */

                    if (fibreanisotropy) {
                        if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                            ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax) {
                            pC->parameters.fail_count = 10;
                        }
                    }
                    else {
                        if (default_microenvironment_options.simulate_2D) {
                            while (pC->parameters.fail_count < 10) {
                                if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                                    ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax) {
                                    pC->parameters.fail_count++;
                                    pC->state.orientation = UniformOnUnitCircle();
                                    xs = position[0] -
                                         pC->parameters.mLength * pC->state.orientation[0];
                                    xe = position[0] +
                                         pC->parameters.mLength * pC->state.orientation[0];
                                    ys = position[1] -
                                         pC->parameters.mLength * pC->state.orientation[1];
                                    ye = position[1] +
                                         pC->parameters.mLength * pC->state.orientation[1];
                                }
                                else {
                                    break;
                                }
                            }
                        }

                        if (!default_microenvironment_options.simulate_2D) {
                            while (pC->parameters.fail_count < 10) {
                                if (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax ||
                                    ys < Ymin || ye > Ymax || ye < Ymin || ys > Ymax ||
                                    zs < Zmin || ze > Zmax || ze < Zmin || zs > Zmax) {
                                    pC->parameters.fail_count++;
                                    pC->state.orientation = UniformOnUnitSphere();
                                    xs = position[0] -
                                         pC->parameters.mLength * pC->state.orientation[0];
                                    xe = position[0] +
                                         pC->parameters.mLength * pC->state.orientation[0];
                                    ys = position[1] -
                                         pC->parameters.mLength * pC->state.orientation[1];
                                    ye = position[1] +
                                         pC->parameters.mLength * pC->state.orientation[1];
                                    zs = position[2] -
                                         pC->parameters.mLength * pC->state.orientation[2];
                                    ze = position[2] +
                                         pC->parameters.mLength * pC->state.orientation[2];
                                }
                                else {
                                    break;
                                }
                            }
                        }
                    }

                    pC->assign_position(position);

                    // relabel so that the rest of the code works (HACK)
                    pC->type_name = "fibre";

                }
            }

        }

    }

    int number_of_agents = (*all_cells).size();
    for( int i=0; i < number_of_agents; i++ ){
        if ((*all_cells)[i]->parameters.fail_count >= 10) {
            std::cout << "I failed to place " << (*all_cells)[i]->type_name << " " <<
                      (*all_cells)[i]->ID << " in the domain - I am deleting agent " << std::endl;
            delete_cell((*all_cells)[i]);
        }
    }
    // !!! PHYSIMESS CODE BLOCK END !!! //

    std::cout << std::endl;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 