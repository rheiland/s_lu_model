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

Cell_Definition* pMacro;

void update_cancer_phenotype(Cell* pCell, Phenotype& phenotype, double dt)
{
	// supported cycle models:
		// advanced_Ki67_cycle_model= 0;
		// basic_Ki67_cycle_model=1
		// live_cells_cycle_model = 5; 
    
    //rwh
    static double Xmin = microenvironment.mesh.bounding_box[0]; 
	static double Ymin = microenvironment.mesh.bounding_box[1]; 

	static double Xmax = microenvironment.mesh.bounding_box[3]; 
	static double Ymax = microenvironment.mesh.bounding_box[4]; 

	if ((pCell->position[0] <= Xmin) || (pCell->position[0] >= Xmax)) return;
	if ((pCell->position[1] <= Ymin) || (pCell->position[1] >= Ymax)) return;



	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}
	
	// set up shortcuts to find the Q and K(1) phases (assuming Ki67 basic or advanced model)
	static bool indices_initiated = false; 
	static int start_phase_index; // Q_phase_index; 
	static int end_phase_index; // K_phase_index;
	static int necrosis_index; 
	
	static int oxygen_substrate_index = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	
	if( indices_initiated == false )
	{
		// Ki67 models
		
		if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || 
			phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_negative );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			
			if( phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
			{
				end_phase_index = 
					phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_positive );
				indices_initiated = true; 
			}
			if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model )
			{
				end_phase_index = 
					phenotype.cycle.model().find_phase_index( PhysiCell_constants::Ki67_positive_premitotic );
				indices_initiated = true; 
			}
		}
		
		// live model 
			
		if( phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
			indices_initiated = true; 
		}
		
		// cytometry models 
		
		if( phenotype.cycle.model().code == PhysiCell_constants::flow_cytometry_cycle_model || 
			phenotype.cycle.model().code == PhysiCell_constants::flow_cytometry_separated_cycle_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::G0G1_phase );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::S_phase );
			indices_initiated = true; 
		}	

		if( phenotype.cycle.model().code == PhysiCell_constants::cycling_quiescent_model )
		{
			start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::quiescent );
			necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
			end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::cycling );
			indices_initiated = true; 
		}
		
	}
	
	// don't continue if we never "figured out" the current cycle model. 
	if( indices_initiated == false )
	{
		return; 
	}

	// sample the microenvironment to get the pO2 value 
	double pO2 = (pCell->nearest_density_vector())[oxygen_substrate_index]; // PhysiCell_constants::oxygen_index]; 
	int n = pCell->phenotype.cycle.current_phase_index(); 
	
	// this multiplier is for linear interpolation of the oxygen value 
	double multiplier = 1.0;
	if( pO2 < pCell->parameters.o2_proliferation_saturation )
	{
		multiplier = ( pO2 - pCell->parameters.o2_proliferation_threshold ) 
			/ ( pCell->parameters.o2_proliferation_saturation - pCell->parameters.o2_proliferation_threshold );
	}
	if( pO2 < pCell->parameters.o2_proliferation_threshold )
	{ 
		multiplier = 0.0; 
	}

    // Find the proliferation factor concentration in the cell's environment
    static int proliferation_factor_index = microenvironment.find_density_index("proliferation factor");
	double prolif_factor = pCell->nearest_density_vector()[proliferation_factor_index];

	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	double base_prolif_rate = pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(cycle_start_index,cycle_end_index);
	
	// check if prolif saturation has been reached
	double prolif_saturation = pCell->custom_data["prolif_saturation"];

	// max out
	if (prolif_factor > prolif_saturation) {
		prolif_factor = prolif_saturation;
	}

	// now, update the appropriate cycle transition rate 
	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = (multiplier) * (1 + prolif_factor) *
		pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index);
	
	// Update necrosis rate
	multiplier = 0.0;
	if( pO2 < pCell->parameters.o2_necrosis_threshold )
	{
		multiplier = ( pCell->parameters.o2_necrosis_threshold - pO2 ) 
			/ ( pCell->parameters.o2_necrosis_threshold - pCell->parameters.o2_necrosis_max );
	}
	if( pO2 < pCell->parameters.o2_necrosis_max )
	{ 
		multiplier = 1.0; 
	}	
	
	// now, update the necrosis rate 
	pCell->phenotype.death.rates[necrosis_index] = multiplier * pCell->parameters.max_necrosis_rate; 
	
	// check for deterministic necrosis 
	if( pCell->parameters.necrosis_type == PhysiCell_constants::deterministic_necrosis && multiplier > 1e-16 )
	{ pCell->phenotype.death.rates[necrosis_index] = 9e99; } 

	// adjust adhesion based on mmp
	static int mmp_i = microenvironment.find_density_index("mmp");
	double mmp = pCell->nearest_density_vector()[mmp_i];
	double mmp_saturation = pCell->custom_data["mmp_saturation"];

	// max out
	if (mmp > mmp_saturation) {
		mmp = mmp_saturation;
	}
	phenotype.mechanics.cell_cell_adhesion_strength = (1 - mmp) * 
		pCell->parameters.pReference_live_phenotype->mechanics.cell_cell_adhesion_strength;

	// update dirichlet nodes
	microenvironment.remove_dirichlet_node(microenvironment.nearest_voxel_index( pCell->position));
	return; 
}

void create_macrophage_cell_type( void )
{
	pMacro = find_cell_definition( "macrophage" ); 
	
	static int oxygen_ID = microenvironment.find_density_index( "oxygen" );  
	
	// reduce o2 uptake 
	pMacro->phenotype.secretion.uptake_rates[oxygen_ID] *= 
		parameters.doubles("immune_o2_relative_uptake");  

	pMacro->phenotype.mechanics.relative_maximum_attachment_distance 
		= pMacro->custom_data["max_attachment_distance"] / pMacro->phenotype.geometry.radius ; 
		
	pMacro->phenotype.mechanics.attachment_elastic_constant 
		= pMacro->custom_data["elastic_coefficient"]; 		
	
	pMacro->phenotype.mechanics.relative_detachment_distance 
		= pMacro->custom_data["max_attachment_distance" ] / pMacro->phenotype.geometry.radius ; 	
	
	// set functions 
	pMacro->functions.update_phenotype = macro_update_phenotype; 
	pMacro->functions.custom_cell_rule = macro_cell_rule;
	pMacro->functions.update_migration_bias = NULL; // macro_motility;
	pMacro->functions.contact_function = adhesion_contact_function; 
	
	return; 
}

bool trigger_apoptosis( Cell* pAttacker, Cell* pTarget )
{	
	// std::cout << "trigger_apoptosis" << std::endl;
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

void macro_motility( Cell* pCell, Phenotype& phenotype, double dt)
{
	// // get densities
    // static int il10_index = microenvironment.find_density_index("il10");
    // static int tnfa_index = microenvironment.find_density_index("tnfa");
    // static int oxygen_index = microenvironment.find_density_index("oxygen");
	// double il10_concentration = pCell->nearest_density_vector()[il10_index];
	// double tnfa_concentration = pCell->nearest_density_vector()[tnfa_index];
    // double oxygen_concentration = pCell->nearest_density_vector()[oxygen_index];

	// // get polarizations
	// double m1_polarization = pCell->custom_data["m1_polarization"]; 
	// double m2_polarization = pCell->custom_data["m2_polarization"]; 
	// double tam_polarization = 1.0 - (m1_polarization + m2_polarization);

	// // m1 behaviors
	// double rand = UniformRandom(); 
	// if (rand < m1_polarization) {
	// 	// towards tnfa
	// }

	// // m2 behaviors
	// double rand = UniformRandom(); 
	// if (rand < m2_polarization) {
	// 	// towards il10 and hypoxia
	// }

	// // tam behaviors
	// double rand = UniformRandom(); 
	// if (rand < tam_polarization) {
	// 	// towards il10 and hypoxia
	// }

    return;
}

void macro_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}
	
	// get polarizations
	double m1_polarization = pCell->custom_data["m1_polarization"]; 
	double m2_polarization = pCell->custom_data["m2_polarization"]; 
	double tam_polarization = 1.0 - (m1_polarization + m2_polarization);

    static int il10_index = microenvironment.find_density_index("il10");
	double il10_concentration = pCell->nearest_density_vector()[il10_index];

	// m1 behaviors
	if (UniformRandom() < m1_polarization) {
		// induce apoptosis, inhibited by microenv il10
		double il10_inhibition_saturation = pCell->custom_data["il10_inhibition_saturation"]; 
		// std::cout << "il10: " << il10_concentration << std::endl;
		double inhib = 1 - std::min(1.0, il10_concentration / il10_inhibition_saturation);
		// std::cout << "inhib: " << inhib << std::endl;
		// I'm not docked, look for cells nearby and try to docked
		if( UniformRandom() < pCell->custom_data["attachment_prob"] * inhib ) {
			// if this returns non-NULL, we're now attached to a cell 
			if( check_neighbors_for_attachment( pCell , dt) ) {
				// set motility off 
				phenotype.motility.is_motile = false; 
			} else {
				phenotype.motility.is_motile = true; 
			}
		}
		// if I'm docked
		if( pCell->state.number_of_attached_cells() > 0 ) {
			bool detach_me = false; 
			// attempt to kill my attached cell
       		if (pCell->state.attached_cells[0]->phenotype.death.dead == false) {
				trigger_apoptosis( pCell, pCell->state.attached_cells[0] ); 
				detach_me = true; 
			}
			// if I dettach, resume motile behavior 
			if( detach_me ) {
				detach_cells( pCell, pCell->state.attached_cells[0] ); 
				phenotype.motility.is_motile = true; 
			}
		}
		// secrete tnfa
    	static int tnfa_index = microenvironment.find_density_index("tnfa");
		// double tnfa = pCell->nearest_density_vector()[tnfa_index];
		// std::cout << "tnfa: " << tnfa << std::endl;
		phenotype.secretion.secretion_rates[tnfa_index] = m1_polarization * 
			pCell->parameters.pReference_live_phenotype->secretion.secretion_rates[tnfa_index];
	}

	// m2 behaviors
	if (UniformRandom() < m2_polarization) {
	// promote proliferation
    	static int proliferation_factor_i = microenvironment.find_density_index("proliferation factor");
		phenotype.secretion.secretion_rates[proliferation_factor_i] = m2_polarization * 
			pCell->parameters.pReference_live_phenotype->secretion.secretion_rates[proliferation_factor_i];
	// secrete il10
	    static int il10_i = microenvironment.find_density_index("il10");
		phenotype.secretion.secretion_rates[il10_i] = m2_polarization * 
			pCell->parameters.pReference_live_phenotype->secretion.secretion_rates[il10_i];
	}

	// tam behaviors
	if (UniformRandom() < tam_polarization) {
		double tam_relative_il10_secretion = pCell->custom_data["tam_relative_il10_secretion"]; 
	// promote proliferation
		static int proliferation_factor_i = microenvironment.find_density_index("proliferation factor");
		phenotype.secretion.secretion_rates[proliferation_factor_i] = tam_polarization * tam_relative_il10_secretion * 
			pCell->parameters.pReference_live_phenotype->secretion.secretion_rates[proliferation_factor_i];
	// secrete mmp (weaken adhesion)
		static int mmp_i = microenvironment.find_density_index("mmp");
		phenotype.secretion.secretion_rates[mmp_i] = tam_polarization * 
			pCell->parameters.pReference_live_phenotype->secretion.secretion_rates[mmp_i];
	// secrete more il10 (stronger than m2)
		phenotype.secretion.secretion_rates[il10_index] = tam_polarization * tam_relative_il10_secretion * 
			pCell->parameters.pReference_live_phenotype->secretion.secretion_rates[il10_index];
	}
}

void macro_update_phenotype( Cell* pCell, Phenotype& phenotype, double dt)
{
	// update dirichlet
	microenvironment.remove_dirichlet_node(microenvironment.nearest_voxel_index( pCell->position));

	// polarize if applicable after buffer time 
	if( UniformRandom() < dt / ( pCell->custom_data[ "polarization_lifetime"] + 1e-15 ) ) {
		// get densities
		static int il10_index = microenvironment.find_density_index("il10");
		static int tnfa_index = microenvironment.find_density_index("tnfa");
		double il10 = pCell->nearest_density_vector()[il10_index];
		double tnfa = pCell->nearest_density_vector()[tnfa_index];
		double prop_tnfa = tnfa / (il10 + tnfa);
		double prop_il10 = il10 / (il10 + tnfa);
		std::vector<double> delta_polarization = {prop_tnfa, 0, prop_il10};

		// get current polarization
		double m1_polarization = pCell->custom_data["m1_polarization"]; 
		double m2_polarization = pCell->custom_data["m2_polarization"]; 
		double tam_polarization = 1.0 - (m1_polarization + m2_polarization);
		std::vector<double> curr_polarization = {m1_polarization, m2_polarization, tam_polarization};

		// change polarization based on delta
		double max_polarization_step = pCell->custom_data[ "max_polarization_step" ]; 
		for (int i = 0; i < 3; i++) {
			curr_polarization[i] += max_polarization_step * (delta_polarization[i] - curr_polarization[i]);
		}
		// normalize
		double total_polarization = curr_polarization[0] + curr_polarization[1] + curr_polarization[2];
		for (int i = 0; i < 3; i++) {
			curr_polarization[i] /= total_polarization;
		}
		// update data
		pCell->custom_data["m1_polarization"] = curr_polarization[0];
		pCell->custom_data["m2_polarization"] = curr_polarization[1];
	}
	return;
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

	static int oxygen_ID = microenvironment.find_density_index( "oxygen" );
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/

	initialize_cell_definitions_from_pugixml(); 
		
	cell_defaults.functions.update_phenotype = update_cancer_phenotype;
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = adhesion_contact_function; 
	cell_defaults.functions.update_migration_bias = NULL; 

	// create the immune cell type 
	create_macrophage_cell_type();

	build_cell_definitions_maps(); 

	setup_signal_behavior_dictionaries(); 	

	setup_cell_rules(); 

	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	initialize_microenvironment(); 	

	std::vector< double > position = {0.0,0.0,0.0};
	for( unsigned int n=0; n < microenvironment.number_of_voxels() ; n++ ) {
		microenvironment.add_dirichlet_node(n,default_microenvironment_options.Dirichlet_condition_vector);
	}

	return;
}	

void setup_tissue( int init_cell_count, double macro_prop )
{	
	// place a cluster of tumor cells at the center 
	Cell* pCell = NULL; 
	std::vector<std::vector<double>> positions;
	int initial_num_immune_cells = (int) round(init_cell_count * macro_prop);
	
	positions = create_cell_circle_positions_by_count(init_cell_count + initial_num_immune_cells); 
	shuffle(positions);

    // double proportion_m2_macrophages = 0.5;
	
	std::cout << "creating " << positions.size() << " closely-packed cells ... " << std::endl; 
		
	for( int i=0; i < positions.size(); i++ )
	{	
		// create immune cell
		if (i < initial_num_immune_cells) {
			pCell = create_cell(*pMacro);
			double m1, m2, tam;
    		// randomPolarization(m1, m2, tam);
			m1 = UniformRandom(); 
			m2 = 1 - m1; 
			pCell->custom_data["m1_polarization"] = m1; 
			pCell->custom_data["m2_polarization"] = m2; 
        } else {
            // Create a tumor cell
            pCell = create_cell();
        }
		pCell->assign_position( positions[i] );
	}
	microenvironment.remove_dirichlet_node(microenvironment.nearest_voxel_index( pCell->position));
	
	return; 
}

void randomPolarization(double &p0, double &p1, double &p2) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::gamma_distribution<double> gammaDist(1.0, 1.0);

    // Generate three gamma-distributed values
    double x0 = gammaDist(gen);
    double x1 = gammaDist(gen);
    double x2 = gammaDist(gen);
    
    // Normalize to sum to 1
    double sum = x0 + x1 + x2;
    p0 = x0 / sum;
    p1 = x1 / sum;
    p2 = x2 / sum;
}

void shuffle(std::vector<std::vector<double>>& positions) {
    srand(time(0));  // Seed the random number generator

    for (size_t i = positions.size() - 1; i > 0; i--)
    {
        // Pick a random index between 0 and i
        size_t j = std::rand() % (i + 1);

        // Swap elements
        std::swap(positions[i], positions[j]);
    }
}

void setup_tissue_from_csv( void )
{
    // load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();
	return; 
}

std::vector<std::vector<double>> create_cell_sphere_positions(double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;

	double cell_radius = cell_defaults.phenotype.geometry.radius; 
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

std::vector<std::vector<double>> create_cell_circle_positions(double circle_radius)
{
    std::vector<std::vector<double>> cells;
    int xc = 0, yc = 0;
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
    
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

std::vector<std::vector<double>> create_cell_circle_positions_by_count( int cell_count ) {
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

std::vector<std::vector<double>> create_cell_sphere_positions_by_count(int cell_count) {
    std::vector<std::vector<double>> cells;
    
    double ds = 15;  // Approximate cell spacing
    double r = 10.0; // Initial radius
    double pi = 3.14159265358979323846;
    double two_pi = 2.0 * pi;

    int i = 0;
    while (i < cell_count) {
        double phi = 0.0;  // Latitude angle (0 to π)
        double d_phi = pi / sqrt(cell_count);  // Step size for latitude

        while (phi < pi && i < cell_count) {
            double theta = 0.0;  // Longitude angle (0 to 2π)
            int num_cells = (int)floor(two_pi * sin(phi) * r / ds);  // Cells per latitude ring
            double d_theta = two_pi / num_cells;

            while (theta < two_pi && i < cell_count) {
                double x = r * sin(phi) * cos(theta);
                double y = r * sin(phi) * sin(theta);
                double z = r * cos(phi);

                cells.push_back({x, y, z});

                theta += d_theta;
                i++;
            }

            phi += d_phi;
        }

        r += ds;  // Expand outward
    }

    return cells;
}

std::vector<std::string> cancer_immune_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4, "black" ); 
	// macros are purple
	if( pCell->type == 1 )
	{ 
		// TODO: color based on polarization and tam expression
		double m1 = pCell->custom_data["m1_polarization"];
    	double m2 = pCell->custom_data["m2_polarization"];
    	double tam = 1.0 - (m1 + m2);

		// Determine the highest polarization value and color based on it
        if(m1 >= m2 && m1 >= tam) {
            output[0] = "rgb(255, 0, 0)";  // Red for M1
        } else if(m2 >= m1 && m2 >= tam) {
            output[0] = "rgb(0, 255, 0)";  // Green for M2
        } else {
            output[0] = "rgb(255, 0, 255)";  // Magenta for TAM
        }
        output[1] = output[0];
		char szTempString [128];
		// sprintf( szTempString , "rgb(%u,%u,%u)", r, g, b );
		// output[0].assign( szTempString );
		// output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
		output[2].assign( szTempString );
		return output;
	} 
	
	// live cells are green, but shaded yellow by prolif factor value 
	if( pCell->phenotype.death.dead == false )
	{
		static int prolif_index = microenvironment.find_density_index( "proliferation factor" );
		double prolif_factor = pCell->nearest_density_vector()[prolif_index];
		double prolif_saturation = pCell->custom_data["prolif_saturation"];
		
		int proliferation = (int) round(std::min(1.0, prolif_factor/prolif_saturation) * 255.0 ); 

		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", proliferation, proliferation, 255 - proliferation );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/2.0) , (int)round(output[0][1]/2.0) , (int)round(output[0][2]/2.0) );
		output[2].assign( szTempString );
		
		return output; 
	}

	// if not, dead colors 	
	// if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	// {
	// 	output[0] = "rgb(255,0,0)";
	// 	output[2] = "rgb(125,0,0)";
	// }
	
	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic)
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

