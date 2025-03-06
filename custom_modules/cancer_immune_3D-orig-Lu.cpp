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

#include "./cancer_immune_3D.h"

Cell_Definition* pM1; 
Cell_Definition* pM2;

void update_cancer_phenotype(Cell* pCell, Phenotype& phenotype, double dt)
{
	// if cell is dead, don't bother with future phenotype changes
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}

    // Find the proliferation factor concentration in the cell's environment
    static int proliferation_factor_index = microenvironment.find_density_index("proliferation factor");
	// interalized
	// double prolif_factor = pCell->phenotype.molecular.internalized_total_substrates[proliferation_factor_index];
	double prolif_factor = pCell->nearest_density_vector()[proliferation_factor_index];

    // std::cout << "Vector elements: ";
    // for (double value : pCell->phenotype.molecular.internalized_total_substrates) {
    //     std::cout << value << " ";
    // }
    // std::cout << std::endl;

	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	double base_prolif_rate = parameters.doubles("base_prolif_rate"); 
	
	// check if prolif saturation has been reached
	static int prolif_saturation_i = pCell->custom_data.find_variable_index( "prolif_saturation" ); 
	double prolif_saturation = pCell->custom_data[prolif_saturation_i];

	// max out
	if (prolif_factor > prolif_saturation) {
		prolif_factor = prolif_saturation;
	}
	// if (prolif_factor != 0) {
	// 		std::cout << "prolif_factor: " << prolif_factor << std::endl; 
	// }

	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);

	// Update the cell cycle transition rate based on the proliferation factor
	phenotype.cycle.data.transition_rate( cycle_start_index, cycle_end_index ) = base_prolif_rate * (1 + prolif_factor); 
	// std::cout << "transition_rate: " << phenotype.cycle.data.transition_rate( cycle_start_index, cycle_end_index ) << std::endl; 

	return;
}

void create_macrophage_cell_type( void )
{
	pM1 = find_cell_definition( "M1" ); 
	pM2 = find_cell_definition( "M2" ); 
	
	static int oxygen_ID = microenvironment.find_density_index( "oxygen" );  
	
	// reduce o2 uptake 
	pM1->phenotype.secretion.uptake_rates[oxygen_ID] *= 
		parameters.doubles("immune_o2_relative_uptake");  
	
	pM1->phenotype.mechanics.cell_cell_adhesion_strength *= 
		parameters.doubles("immune_relative_adhesion"); 
	pM1->phenotype.mechanics.cell_cell_repulsion_strength *= 
		parameters.doubles("immune_relative_repulsion"); 

	pM2->phenotype.secretion.uptake_rates[oxygen_ID] *= 
		parameters.doubles("immune_o2_relative_uptake");  
	
	pM2->phenotype.mechanics.cell_cell_adhesion_strength *= 
		parameters.doubles("immune_relative_adhesion"); 
	pM2->phenotype.mechanics.cell_cell_repulsion_strength *= 
		parameters.doubles("immune_relative_repulsion"); 
		
	// figure out mechanics parameters 	
	pM1->phenotype.mechanics.relative_maximum_attachment_distance 
		= pM1->custom_data["max_attachment_distance"] / pM1->phenotype.geometry.radius ; 
		
	pM1->phenotype.mechanics.attachment_elastic_constant 
		= pM1->custom_data["elastic_coefficient"]; 		
	
	pM1->phenotype.mechanics.relative_detachment_distance 
		= pM1->custom_data["max_attachment_distance" ] / pM1->phenotype.geometry.radius ; 		

	pM2->phenotype.mechanics.relative_maximum_attachment_distance 
		= pM2->custom_data["max_attachment_distance"] / pM2->phenotype.geometry.radius ; 
		
	pM2->phenotype.mechanics.attachment_elastic_constant 
		= pM2->custom_data["elastic_coefficient"]; 		
	
	pM2->phenotype.mechanics.relative_detachment_distance 
		= pM2->custom_data["max_attachment_distance" ] / pM2->phenotype.geometry.radius ; 	
	
	// set functions 

	pM1->functions.update_phenotype = m1_update_phenotype; 
	pM1->functions.custom_cell_rule = m1_cell_rule; 
	pM1->functions.update_migration_bias = m1_motility;
	pM1->functions.contact_function = adhesion_contact_function; 

	pM2->functions.update_phenotype = m2_update_phenotype; 
	pM2->functions.custom_cell_rule = m2_cell_rule;
	pM2->functions.update_migration_bias = m2_motility;
	pM2->functions.contact_function = adhesion_contact_function; 
	
	return; 
}

void m1_motility( Cell* pCell, Phenotype& phenotype, double dt)
{
    static int tnfa_index = microenvironment.find_density_index("tnfa");

    // Get local TNF-α concentration
    double tnfa_concentration = pCell->nearest_density_vector()[tnfa_index];

    // // Define the TNF-α threshold for stopping motility
	// static int tnfa_motility_saturation_i = pCell->custom_data.find_variable_index( "tnfa_motility_saturation" );
	// double tnfa_motility_saturation = pCell->custom_data[tnfa_motility_saturation_i];

    // // If TNF-α exceeds the threshold, stop moving
    // if (tnfa_concentration >= tnfa_motility_saturation)
    // {
    //     phenotype.motility.is_motile = true;
    //     phenotype.motility.migration_bias = 0;
		
	// 	// // Debugging output
	// 	std::cout << "TNF-a: " << tnfa_concentration << std::endl;
    // }
    // else if (pCell->state.attached_cells.size() == 0) // If not docked and TNF-a is below threshold, move normally
    if (pCell->state.attached_cells.size() == 0) // If not docked and TNF-a is below threshold, move normally
    {
        phenotype.motility.is_motile = true; 
        phenotype.motility.migration_bias = 0.5;
        // Follow the TNF-a gradient
        phenotype.motility.migration_bias_direction = pCell->nearest_gradient(tnfa_index);
        normalize(&(phenotype.motility.migration_bias_direction));
    }
    else
    {
        // If docked to another cell, remain immobile
        phenotype.motility.is_motile = false;
    }

    return;
}

void m1_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int attach_lifetime_i = pCell->custom_data.find_variable_index( "attachment_lifetime" ); 
	
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}
	
	// if I'm docked
	if( pCell->state.number_of_attached_cells() > 0 )
	{
		// attempt to kill my attached cell
		
		bool detach_me = false; 
		
        if (pCell->state.attached_cells[0]->phenotype.death.dead == false)
		{
			trigger_apoptosis( pCell, pCell->state.attached_cells[0] ); 
			detach_me = true; 
		}
		
		// decide whether to detach 
		if( UniformRandom() < dt / ( pCell->custom_data[attach_lifetime_i] + 1e-15 ) )
		{ detach_me = true; }
		
		// if I dettach, resume motile behavior 
		
		if( detach_me )
		{
			detach_cells( pCell, pCell->state.attached_cells[0] ); 
			phenotype.motility.is_motile = true; 
		}
		return; 
	}
	
	// I'm not docked, look for cells nearby and try to docked
	// if this returns non-NULL, we're now attached to a cell 
	if( check_neighbors_for_attachment( pCell , dt) )
	{
		// set motility off 
		phenotype.motility.is_motile = false; 
		return; 
	}
	phenotype.motility.is_motile = true; 
	
	return; 
}

void m1_update_phenotype(Cell* pCell, Phenotype& phenotype, double dt) {
	static int polarization_lifetime_i = pCell->custom_data.find_variable_index( "polarization_lifetime" ); 
	static int prob_polarize_i = pCell->custom_data.find_variable_index( "prob_polarize" ); 

	static int il10_threshold_i = pCell->custom_data.find_variable_index( "il10_threshold" ); 
	double il10_threshold = pCell->custom_data[il10_threshold_i];

	// static int il10_index = microenvironment.find_density_index( "il10" );
    // double il10 = pCell->nearest_density_vector()[il10_index];
	static int il10_i = microenvironment.find_density_index( "il10" );
	double il10 = pCell->phenotype.molecular.internalized_total_substrates[il10_i];

	// polarize if applicable after buffer time
	if( UniformRandom() < dt / ( pCell->custom_data[polarization_lifetime_i] + 1e-15 ) ) {
		if( il10 > il10_threshold && UniformRandom() < pCell->custom_data[prob_polarize_i]) {
			std::cout << "il10: " << il10 << std::endl;
			while (!pCell->state.attached_cells.empty())
			{
				Cell* attached_cell = pCell->state.attached_cells.back();
				detach_cells(pCell, attached_cell); // Detach the cell
			}
			// Convert the M1 macrophage into an M2 macrophage
			pCell->convert_to_cell_definition(*find_cell_definition("M2"));
			static int tnfa_i = microenvironment.find_density_index( "tnfa" );
			pCell->phenotype.molecular.internalized_total_substrates[tnfa_i] = 0;
			pCell->phenotype.molecular.internalized_total_substrates[il10_i] = 0;
		}
	}	
	return;
}

void m2_motility( Cell* pCell, Phenotype& phenotype, double dt)
{	

    static int il10_index = microenvironment.find_density_index("il10");
    static int oxygen_index = microenvironment.find_density_index("oxygen");

	// get density
	double il10_concentration = pCell->nearest_density_vector()[il10_index];
    double oxygen_concentration = pCell->nearest_density_vector()[oxygen_index];
	// std::cout << "oxygen_concentration: " << oxygen_concentration << std::endl;

	static int hypoxic_motility_saturation_i = pCell->custom_data.find_variable_index( "hypoxic_motility_saturation" ); 
	double hypoxic_motility_saturation = pCell->custom_data[hypoxic_motility_saturation_i];
	static int il10_motility_saturation_i = pCell->custom_data.find_variable_index( "il10_motility_saturation" ); 
	double il10_motility_saturation = pCell->custom_data[il10_motility_saturation_i];
	
	// stop moving if ive reached enough il10 or enough hypoxia
	if (il10_concentration >= il10_motility_saturation_i && oxygen_concentration <= hypoxic_motility_saturation)
    {
        phenotype.motility.is_motile = false;
        return;
    }
	else {
        phenotype.motility.is_motile = true;
	}

	// Get gradients
	std::vector<double> il10_gradient = pCell->nearest_gradient(il10_index);
	std::vector<double> oxygen_gradient = pCell->nearest_gradient(oxygen_index);

	// Reverse oxygen gradient since cells move towards **low oxygen** (hypoxia)
	for (int i = 0; i < 3; i++)
	{
		oxygen_gradient[i] *= -1.0; // Move toward lower oxygen levels (hypoxia)
	}

	// Define weighting factors
	static int hypoxia_il10_motility_bias_i = pCell->custom_data.find_variable_index( "hypoxia_il10_motility_bias" ); 
	double hypoxia_il10_motility_bias = pCell->custom_data[hypoxia_il10_motility_bias_i];

	// Compute weighted migration direction
	for (int i = 0; i < 3; i++)
	{
		phenotype.motility.migration_bias_direction[i] = 
			hypoxia_il10_motility_bias * oxygen_gradient[i] + (1 - hypoxia_il10_motility_bias) * il10_gradient[i];
	}

	// Normalize direction vector
	normalize(&(phenotype.motility.migration_bias_direction));

    return;
}

void m2_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	// std::cout << "pM2: " << pM2->functions.update_migration_bias << std::endl;
	// std::cout << "pCell: " << pCell->functions.update_migration_bias << std::endl;
	m2_motility(pCell, phenotype, dt);
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}
}

void m2_update_phenotype( Cell* pCell, Phenotype& phenotype, double dt)
{
	
	static int polarization_lifetime_i = pCell->custom_data.find_variable_index( "polarization_lifetime" ); 

	static int inhib_threshold_i = pCell->custom_data.find_variable_index( "il10_inhibit_polarization_threshold" ); 
	double il10_inhibit_polarization_threshold = pCell->custom_data[inhib_threshold_i];
	
	// get internal il10
	static int il10_i = microenvironment.find_density_index( "il10" );
	double internal_il10 = pCell->phenotype.molecular.internalized_total_substrates[il10_i];
	// std::cout << "internal_il10: " << internal_il10 << std::endl;

	// dont polarize if too much internal il10
	if (internal_il10 > il10_inhibit_polarization_threshold) {
		return;
	}

	// check for tnfa for polarization to m1
	static int tnfa_threshold_i = pCell->custom_data.find_variable_index( "tnfa_threshold" ); 
	double tnfa_threshold = pCell->custom_data[tnfa_threshold_i];
	
	static int tnfa_i = microenvironment.find_density_index( "tnfa" );
	double tnfa = pCell->nearest_density_vector()[tnfa_i];

	// polarize if applicable after buffer time
	if( UniformRandom() < dt / ( pCell->custom_data[polarization_lifetime_i] + 1e-15 ) ) {
    // If the stimulatory factor exceeds a threshold, differentiate into the new cell type
		if( tnfa > tnfa_threshold ) {
			while (!pCell->state.attached_cells.empty())
			{
				Cell* attached_cell = pCell->state.attached_cells.back();
				detach_cells(pCell, attached_cell); // Detach the cell
			}
			// Convert to M1 macrophage
			pCell->convert_to_cell_definition(*find_cell_definition("M1"));
			static int il10_i = microenvironment.find_density_index( "il10" );
			pCell->phenotype.molecular.internalized_total_substrates[tnfa_i] = 0;
			pCell->phenotype.molecular.internalized_total_substrates[il10_i] = 0;
			return;
		}
	}	
}

bool trigger_apoptosis( Cell* pAttacker, Cell* pTarget )
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );	
	
	// if the Target cell is already dead, don't bother!
	if( pTarget->phenotype.death.dead == true )
	{ return false; }

	pTarget->start_death( apoptosis_model_index );
	return true; 
}

Cell* check_neighbors_for_attachment( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container(); 
    Cell_Definition* pCancerCell = find_cell_definition("cancer cell");

	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to kill yourself 
		if( nearby[i] != pAttacker && nearby[i]->type == pCancerCell->type)
		{
			if( attempt_attachment( pAttacker, nearby[i] , dt ) )
			{ return nearby[i]; }
		}
		i++; 
	}
	
	return NULL; 
}

bool attempt_attachment( Cell* pAttacker, Cell* pTarget , double dt )
{
	static int attach_rate_i = pAttacker->custom_data.find_variable_index( "attachment_rate" ); 

	double max_attachment_distance = 
		pAttacker->custom_data["max_attachment_distance"];   
	double min_attachment_distance = 
		pAttacker->custom_data["min_attachment_distance"];   
	double attachment_difference = max_attachment_distance - min_attachment_distance; 
	
	if(pTarget->phenotype.death.dead == false )
	{
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm( displacement ); 
		if( distance_scale > max_attachment_distance )
		{ return false; } 

		
		distance_scale *= -1.0; 
		distance_scale += max_attachment_distance; 
		distance_scale /= attachment_difference; 
		if( distance_scale > 1.0 )
		{ distance_scale = 1.0; } 
		
		if( UniformRandom() < pAttacker->custom_data[attach_rate_i] * dt * distance_scale )
		{
//			std::cout << "\t attach!" << " " << pTarget->custom_data[oncoprotein_i] << std::endl; 
			attach_cells( pAttacker, pTarget ); 
		}
		
		return true; 
	}
	
	return false; 
}

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	// housekeeping 
	
	initialize_default_cell_definition();
	
	
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	cell_defaults.parameters.o2_reference = 38.0; 

	static int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0 
	// static int immuno_ID = microenvironment.find_density_index( "immunostimulatory factor" ); // 1
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/

	initialize_cell_definitions_from_pugixml(); 
	
	// // change the max cell-cell adhesion distance 
	// cell_defaults.phenotype.mechanics.relative_maximum_attachment_distance = 
	// 	cell_defaults.custom_data["max_attachment_distance"] / cell_defaults.phenotype.geometry.radius;
		
	// cell_defaults.phenotype.mechanics.relative_detachment_distance 
	// 	= cell_defaults.custom_data["max_attachment_distance"] / cell_defaults.phenotype.geometry.radius ; 
		
	// cell_defaults.phenotype.mechanics.attachment_elastic_constant 
	// 	= cell_defaults.custom_data[ "elastic_coefficient" ];	
		
	cell_defaults.functions.update_phenotype = update_cancer_phenotype;
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = adhesion_contact_function; 
	cell_defaults.functions.update_migration_bias = NULL; 

	// create the immune cell type 
	// create_immune_cell_type(); 

	// create the monocyte and macrophage cell type 
	// create_monocyte_cell_type();
	create_macrophage_cell_type();

	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 


	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	
	// if( default_microenvironment_options.simulate_2D == true )
	// {
	// 	std::cout << "Warning: overriding 2D setting to return to 3D" << std::endl; 
	// 	default_microenvironment_options.simulate_2D = false; 
	// }
	
	initialize_microenvironment(); 	

	return; 
}	

void setup_tissue( int init_cell_count, int initial_num_immune_cells )
{
	// place a cluster of tumor cells at the center 
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	
	// double tumor_radius = 
	// 	parameters.doubles("tumor_radius"); // 250.0; 

	// int init_cell_count = 
	// 	parameters.ints("init_cell_count"); // 250.0; 
	
	Cell* pCell = NULL; 
	
	// std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius, tumor_radius); 
	// std::vector<std::vector<double>> positions = create_cell_circle_positions(cell_radius, tumor_radius); 
	std::vector<std::vector<double>> positions = create_cell_circle_positions_by_count(init_cell_count); 
	
	std::cout << "creating " << positions.size() << " closely-packed tumor cells ... " << std::endl; 
		
	for( int i=0; i < positions.size(); i++ )
	{
		pCell = create_cell(); // tumor cell 
		pCell->assign_position( positions[i] );
	}

	double tumor_radius = -9e9; // 250.0; 
	double temp_radius = 0.0; 
	
	// for the loop, deal with the (faster) norm squared 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		temp_radius = norm_squared( (*all_cells)[i]->position ); 
		if( temp_radius > tumor_radius )
		{ tumor_radius = temp_radius; }
	}
	// now square root to get to radius 
	tumor_radius = sqrt( tumor_radius ); 

    double proportion_m1_macrophages = 0;
    double proportion_m2_macrophages = 0;
	
	// for( int i=0; i < initial_num_immune_cells; i++ )
	// {
	// 	double theta = UniformRandom() * 6.283185307179586476925286766559; 
	// 	double phi = acos( 2.0*UniformRandom() - 1.0 );  
	// 	// double radius = NormalRandom( mean_radius, std_radius );
	// 	double r = tumor_radius * sqrt(UniformRandom()); // Random distance within the circle/sphere
	// 	double random_value = UniformRandom();

	// 	// Create m1
		// if (random_value < proportion_m1_macrophages) {
		// 	pCell = create_cell(*pM1);
		// 	// pCell->assign_position(tumor_radius * cos(theta) * sin(phi), tumor_radius * sin(theta) * sin(phi), tumor_radius * cos(phi)); // Place the macrophage at the calculated position
		// 	pCell->assign_position(tumor_radius * cos(theta) * sin(phi), tumor_radius * sin(theta) * sin(phi), 0); // Place the macrophage at the calculated position
		// }
		// // Create m2
		// else {
		// 	pCell = create_cell(*pM2);
		// 	// pCell->assign_position(tumor_radius * cos(theta) * sin(phi), tumor_radius * sin(theta) * sin(phi), tumor_radius * cos(phi)); // Place the M1 macrophage at the calculated position
		// 	pCell->assign_position(tumor_radius * cos(theta) * sin(phi), tumor_radius * sin(theta) * sin(phi), 0); // Place the M1 macrophage at the calculated position
		// }
	// }
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

	for( int n = 0 ; n < initial_num_immune_cells ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		position[0] = Xmin + UniformRandom()*Xrange; 
		position[1] = Ymin + UniformRandom()*Yrange; 
		position[2] = Zmin + UniformRandom()*Zrange; 
		
		if (UniformRandom() < proportion_m1_macrophages) {
			pCell = create_cell(*pM1);
			pCell->assign_position(position); // Place the macrophage at the calculated position
		}
		// Create m2
		else {
			pCell = create_cell(*pM2);
			pCell->assign_position(position); // Place the M1 macrophage at the calculated position
		}
	}
	
	return; 
}

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);
	
	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;
				
				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
			
		}
	}
	return cells;	
}

#include <vector>
#include <cmath>

std::vector<std::vector<double>> create_cell_circle_positions(double cell_radius, double circle_radius)
{
    std::vector<std::vector<double>> cells;
    int xc = 0, yc = 0;
    
    double x_spacing = cell_radius * sqrt(3);
    double y_spacing = cell_radius * 2;

    std::vector<double> tempPoint(2, 0.0);

    for (double x = -circle_radius; x < circle_radius; x += x_spacing, xc++)
    {
        for (double y = -circle_radius; y < circle_radius; y += y_spacing, yc++)
        {
            tempPoint[0] = x + (yc % 2) * 0.5 * cell_radius;  // Offset x if row index is odd
            tempPoint[1] = y + (xc % 2) * cell_radius;        // Offset y if column index is odd
			tempPoint[2] = 0;

            // Check if point is inside the circle
            if (std::sqrt(tempPoint[0] * tempPoint[0] + tempPoint[1] * tempPoint[1]) < circle_radius)
            {
                cells.push_back(tempPoint);
            }
        }
    }

    return cells;
}

std::vector<std::vector<double>> create_cell_circle_positions_by_count(int cell_count) {
	Cell* pC;
    std::vector<std::vector<double>> cells;

	double theta = 0.0; 
	std::vector<double> position(2, 0.0);
	
	double ds = 15; 
	
	int i = 0; 
	double r = 10.0; 
	double two_pi = 6.283185307179586476925286766559;
	while( i < cell_count )
	{
		
		double arclength = r * 6.283185307179586476925286766559; 
		int num_cells = (int) floor( arclength / ds ); 
		double d_theta = 6.283185307179586476925286766559 /num_cells; 
		
		theta = 0.0;
		while( theta < 6.283185307179586476925286766559 && i < cell_count )
		{
			position[0] = r*cos(theta); 
			position[1] = r*sin(theta); 
			position[2] = 0; 
            cells.push_back(position);

			theta += d_theta; 
			i++; 
		}
		
		r += ds; 
	}
	return cells;
}

std::vector<std::string> cancer_immune_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4, "black" ); 

	// m1 are purple
	if( pCell->type == 1 )
	{ 
		output[0] = "rgb(255, 0, 255)";
		output[1] = "rgb(255, 0, 255)";
		output[2] = "rgb(128, 0, 128)";
		return output;
	} 

	// m2 are green
	if( pCell->type == 2 )
	{ 
		output[0] = "rgb(0, 255, 0)";
		output[1] = "rgb(0, 255, 0)";
		output[2] = "rgb(0, 128, 0)";
		return output;
	} 

	// if I'm under attack, color me 
	// if( pCell->state.attached_cells.size() > 0 )
	// {
	// 	output[0] = "darkcyan"; // orangered // "purple"; // 128,0,128
	// 	output[1] = "black"; // "magenta"; 
	// 	output[2] = "cyan"; // "magenta"; //255,0,255
	// 	return output; 
	// }
	
	// live cells are green, but shaded by prolif factor value 
	if( pCell->phenotype.death.dead == false )
	{
		static int prolif_index = microenvironment.find_density_index( "proliferation factor" );
		// double prolif_factor = pCell->phenotype.molecular.internalized_total_substrates[prolif_index];
		double prolif_factor = pCell->nearest_density_vector()[prolif_index];

		// check if prolif saturation has been reached
		// static int prolif_saturation_i = pCell->custom_data.find_variable_index( "prolif_saturation" ); 
		// double prolif_saturation = pCell->custom_data[prolif_saturation_i];
		// if (prolif_factor > prolif_saturation) {
		// 	prolif_factor = prolif_saturation;
		// }

		// if (prolif_factor != 0) {
		// 	std::cout << "prolif factor: " << prolif_factor <<  std::endl; 
		// }
		
		int proliferation = (int) round(prolif_factor * 255.0 ); 
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", proliferation, proliferation, 255-proliferation );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
		output[2].assign( szTempString );
		
		return output; 
	}

	// if not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	
	return output; 
}

void adhesion_contact_function( Cell* pActingOn, Phenotype& pao, Cell* pAttachedTo, Phenotype& pat , double dt )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	
	static double max_elastic_displacement = pao.geometry.radius * pao.mechanics.relative_detachment_distance; 
	static double max_displacement_squared = max_elastic_displacement*max_elastic_displacement; 
	
	// detach cells if too far apart 
	
	if( norm_squared( displacement ) > max_displacement_squared )
	{
		detach_cells( pActingOn , pAttachedTo );
		return; 
	}
	
	axpy( &(pActingOn->velocity) , pao.mechanics.attachment_elastic_constant , displacement ); 
	
	return; 
}

