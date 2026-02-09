
#ifndef _CONSTANTS
#define _CONSTANTS

// Geometric constants.
#define fiftyseven (180.0/M_PI)
#define fiftyseventh (M_PI/180)
#define triangular (M_PI/1.5)
#define tetrahedral acos(-1.0/3)
#define square (M_PI/2)
#define hexagonal (M_PI/3)
#define circle M_PI*2
#define packed_sphere_ligancy 12

// Physical and chemical constants.
#define _kcal_per_kJ 0.239006
// #define coplanar_threshold 0.5
#define coplanar_threshold 2.5
#define Avogadro 6.02214076e+23
#define kB 1.380649e-23
#define kB_kJmol (kB/1000 * Avogadro)
#define temperature 310.2
#define R 0.00831446261815324
#define kJ_mol_per_summed_bond_strain (91.0 / 4.332)
#define pH 6.0
#define hydrophilicity_cutoff 0.28
#define charge_attraction 60.0
#define amide_zwitterionic_amount 0.25
#define water_molecule_size 2.8

// Warning - increasing these constants significantly above the maximal 35.0, 60.0 values
// will cause docking fails in the unit tests.
#define polar_repulsion 35.0
#define charge_repulsion 60.0

// Programmatic constants.
#define _def_atc 100
#define _MAX_NUM_FORCES 65536
#define any_element -5141
#define ATOM_NOT_OF_AMINO_ACID 0x907aa
#define BOND_DEF_NOT_FOUND 0xbadb09d
#define memsanity 0x10000000
#define NOT_ATOM_RECORD 0xb1207e19
#define PROT_MAX_RGN 40
#define SPHREACH_MAX 256
#define VALENCE_EXCEEDS_GEOMETRY 22

#if defined(__linux__) || defined(__sun) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__APPLE__)
    #define CMD_CHECK_INSTALLED_3P_SMILES_PARSER "which obabel"
    #define CMD_CALL_3P_SMILES_PARSER "obabel -:'%s' --gen3D -osdf 2> /dev/null"
#elif _WIN32
    #define CMD_CHECK_INSTALLED_3P_SMILES_PARSER "dir \"C:\\Program Files\\OpenBabel*\""
    #define CMD_CALL_3P_SMILES_PARSER "obabel.exe -:'%s' --gen3D -osdf"
#else
    #error It appears your operating system is not supported yet - would you be willing to add it and submit a pull request?
#endif

// Periodic table constants.
#define ALKMETAL 1
#define ALKEARTH 2
#define TRIEL 3
#define TETREL 4
#define PNICTOGEN 5
#define CHALCOGEN 6
#define HALOGEN 7
#define NOBLEGAS 8
#define HEAVYMETAL 30
#define LANTHANIDE 60
#define ACTINIDE 90
#define _ATOM_Z_LIMIT 200

// Atomic and molecular spatial constants.
#define _can_clash_angle (120.0 * fiftyseventh)
#define _SANOM_BOND_ANGLE_WEIGHT 10.0
#define _SANOM_BOND_RAD_WEIGHT 30.0
#define _SANOM_BOND_RADIUS_WEIGHT 50.0
#define _SANOM_CLASHES_WEIGHT 10.0
#define _shield_angle (130.0 * fiftyseventh)
#define _shield_angle_pi (100.0 * fiftyseventh)
#define vdw_surface_density 15
#define water_surface_H 0.59

// Atomic similarity constants
#define similarity_charged 1.0
#define similarity_polar 1.0
#define similarity_antipolar 0.75
#define similarity_nonpolar 0.333
#define similarity_pi 0.13

// Molecular generation and conformation constants.
#define _default_mutation_rate 0.1
#define _default_population_size 100
#define echo_non_coplanar_pi_atoms 0
#define _evolution_atom_displacement 4.0
#define _evolution_default_generations 2000
#define stretch_out_molecules_by_interatomic_distance 1

// Protein structure constants.
#define chain_length_per_aa 3.9
#define helix_hbond_cutoff 2.8
#define hx_tight_cutoff 3.2
#define peptide_bond_length 1.32
#define unconnected_residue_mindist 4.82

// Torsion angles for various helices and for beta strands.
#define ALPHA_PHI fiftyseventh*-57.8
#define ALPHA_PSI fiftyseventh*-47.0
#define BETA_PHI fiftyseventh*-140
#define BETA_PSI fiftyseventh*130
#define _310_PHI fiftyseventh*-49
#define _310_PSI fiftyseventh*-26
#define PI_PHI fiftyseventh*-55
#define PI_PSI fiftyseventh*-70
#define POLYPRO2_PHI fiftyseventh*-75
#define POLYPRO2_PSI fiftyseventh*145
#define POLYPRO2_OMEGA fiftyseventh*124
#define POLYPRO1_PHI fiftyseventh*-75
#define POLYPRO1_PSI fiftyseventh*160
#define POLYPRO1_OMEGA fiftyseventh*113

// Docking basics.
#define _def_ang_momentum (fiftyseventh*5)
#define _def_bnd_momentum (fiftyseventh*15)
#define _def_lin_momentum 0.1
#define _DEFAULT_INTERA_R_CUTOFF 6
#define _INTER_TYPES_LIMIT 10
#define clash_limit_per_aa 7.0
#define clash_limit_per_atom 5.0
#define double_hydrogen_clash_allowance_multiplier 1.5
#define generate_new_sidechain_clashes 0
#define global_clash_allowance 0.4
#define ignore_double_hydrogen_clashes 0
#define ignore_nonpolar_hydrogen_clashes 0
#define iteration_additional_clash_coefficient 10.0
#define iteration_internal_clash_coefficient 0.25
#define Lennard_Jones_epsilon 1.0
#define Lennard_Jones_epsilon_x4 Lennard_Jones_epsilon*4
#define linear_motion_probability 0.1
#define lmpull 0.1
#define lmpush 0.3
#define lmsteps 5
#define metropolis_criterion 1
#define speed_limit 1.5

// Docking features.
#define allow_abhor_vacuum 1
#define allow_axial_tumble 1
#define allow_bond_rots 1
#define _allow_conditional_basicity 1
#define _allow_conditional_basicity_with_acid_ligand 1
#define allow_iter_cb 1
#define allow_linear_motion 1
#define _ALLOW_PROTONATE_PNICTOGENS 0
#define _allow_Schiff_base_formation 1
#define allow_stay_close_flexions 0
#define allow_tethered_rotations 0
#define auto_pK_protonation 0
#define conj_charge_as_polarity 0
#define _dock_result_in_iter 1
#define flexion_selection 1
#define include_eclipses 0
#define include_residue_eclipses 0
#define limit_flexions_by_mclash 0
#define make_thiolates_in_scoring 0
#define mclashables_as_residue_nearbys 1
#define no_zero_flexions 1
#define nodes_no_ligand_360_flex 1
#define nodes_no_ligand_360_tumble 1
#define pocketcen_from_reach_atoms 0
#define pocketcen_is_loneliest 1
#define point_hydro_side_chains_inward 0
#define prevent_ligand_360_on_activate 1
#define recapture_ejected_ligand 0
#define reuse_best_pose 1
#define stay_close_conform_iters 1
#define stays_rotation 1
#define stays_rotation_verbose 0
#define summed_missed_connections 1
#define _teleport_dissatisfied_waters 0
#define use_exclusions 1
#define warn_orphan_atoms 0
// Auto hydroxy makes geraniol fail in OR1A1. So does pre-rotate side chains.
// Auto hydroxy was supposed to always point polar hydrogens towards nearby H-bond acceptors.
// Prerot side chains was supposed to minimize clashes between the ligand and side chains after
// choosing an initial candidate conformer from tumble spheres.
// Obviously, both are DISMAL FAILURES just like everything else this project's main author has ever done.
#define allow_auto_hydroxy 0
#define prerot_sidechains_from_ligand 0

// Flexion settings.
#define _ALLOW_FLEX_RINGS 0
#define _allow_fullrot 1
#define _fullrot_every 5
#define _fullrot_stepdeg 3
#define _fullrot_steprad (fiftyseventh*_fullrot_stepdeg)
#define attempt_to_connect_hydrogen_bonds_to_ligand 1
#define default_pre_ligand_flex_radius 10
#define default_pre_ligand_multimol_radius 15
#define flexion_maxangle hexagonal
#define flexion_probability_multiplier 1
#define flexion_sub_iterations_ligand 15
#define flexion_sub_iterations_sidechain 1
#define fullrot_flex_first_subiter_only 0
#define fullrot_flex_residues_only 0
#define fullrot_flex_unfavorable_energy_only 0
#define fullrot_forbid_residues 1
#define sidechain_flexion_frequency 0.333
#define stochastic_flexion_of_clashing_residues 1

// Dynamic modification constants.
#define _water_satisfaction_threshold -5
#define _water_teleport_tries 200
#define cond_bas_hbond_distance_threshold 2.9
#define cond_bas_hbond_energy_threshold 8.0
#define homology_long_axis_rotations 0
#define homology_phi_psi_rotations 0
#define homology_region_optimization 0
#define hydrogenate_add_missing_heavy_atoms 1

// Ligand repositioning constants.
#define multimol_stays_enforcement 0.5
#define multimol_stays_allow_revert_worsening 0
#define secondary_stays_effect 0.5
#define stays_tolerance_factor .15
#define reuse_pose_probability 0.333
#define best_pose_reset_frequency 0.25
// Drift pulls the ligand towards the loneliest point if it encounters clashes.
// Turning it off can cause the ligand to be ejected from the protein.
#define allow_drift 0
#define initial_drift 0.3
#define drift_decay_rate 0.01
#define drift_energy_threshold 1000

// Tumble spheres constants.
#define tumble_spheres_include_vdW 0

// Best-Binding constants.
#define bb_avoid_eclipsing_contacts 1
#define bb_clash_avoidance_threshold 5e8
#define bb_disqualification_energy 1000
#define bb_eclipsing_divisor 10
#define bb_eclipsing_limit 2
#define bb_enable_residue_disqualifications 0
#define bb_find_empty_space 1
#define bb_group_distance_cutoff 4.3
#define bb_lonely_radius 2
#define bb_pocket_res_extra_spacing 2.5
#define bb_pocket_res_spacing_allowance 1
#define bb_pullaway_allowance 1.0
#define bb_realign_amount 0.3
#define bb_realign_iters 1
#define bb_realign_only_hydro 1
#define bb_secondary_must_be_farthest_from_primary 1
#define bb_secondary_must_be_at_least_threshold_distance_from_primary 0
#define bb_secondary_threshold_distance_times_farthest 0.8
#define bb_stochastic 0.333
#define bb_stochastic_A 1.5
#define bb_stochastic_ligflex 0
#define best_binding_stochastic 0.3
#define enable_bb_scooch 1
#define enforce_no_bb_pullaway 0

// Cavity constants.
#define default_cavity_stuffing 0.03
#define min_cavmatch_ctainmt 0.6
#define cavity_intersect_threshold 1e2
#define min_cvty_ctnmt 0.25

// Soft docking constants.
#define initial_soft_contact_elasticity 0.1
#define move_ligand_with_soft_motion 1
#define progressively_increase_elasticity 1
#define repack_soft_clashes false
#define soft_repack_iterations 20
#define soft_contact_elasticity_progress 0.01
#define soft_iter_min 10
#define soft_push_multiplier 0.01

// Vestibule feature (non-functional)
#define enable_vestibules 0
#define MAX_VESTIBULE 10

// Optimization constants.
#define contact_energy_allowance_for_optimization 0.3
#define optimize_energies_in_contact_anomaly_check 0
#define contact_energy_optimization_iterations 5
#define optimize_internal_contacts_post_iterations 1

// Performance constants.
#define iter_lostreturns_threshold 0.05
#define max_iters_without_ligand_change 5

// Scoring constants.
#define compute_lsrb 0
#define lsrb_vdw_multiplier 1.0
#define priority_weight_group 25
#define ts_priority_coefficient 10
#define ignore_invalid_partial 1
#define include_mcr_in_dock_energy 0
#define compute_missed_connections 1
#define metal_anisotropy 1
#define mandatory_coordination_threshold 5
// Extra weight given to sidechain-ligand binding strengths during conformer search.
#define dock_ligand_bias 1.0
#define limit_interactions_by_hydrophilicity 0
#define polar_sat_influence_for_dock 30
#define polar_sat_influence_for_bb 30
#define polar_sat_influence_for_scoring 0
#define per_atom_occlusions 0

// Path-based docking.

// Amount to reduce momenta for path nodes beyond zero. Since the point of path based
// docking is to keep as closely as possible the same ligand pose and move it through
// the protein, we want to minimize the ligand's conformational changes from node to
// node. Linear momenta are not affected by this number. Angular momenta are multiplied
// by this number, and bond rotation momenta by the square of this number.
#define internode_momentum_mult 0.25
#define internode_momentum_only_on_activation 1
#define recenter_ligand_each_node 0


// Obsolete stuff.
#define compute_vdw_repulsion 0
#define compute_clashdirs 0
#define redo_tumble_spheres_every_node 0
#define accommodate_ligand_in_post 0


// Debugging stuff.

// Can be up to 3:
#define _bb_max_grp 3

// For auditing binding energies between individual atoms:
#define _peratom_audit 0
#define _peratom_audit_nans 0

// A short-term feature that can be used in case the bond reciprocity problem recurs. See issue #423.
#define bond_reciprocity_fix 0

// Should normally be false or zero:
#define _dbg_259 0
#define _dbg_415 0
#define _dbg_51e2_ionic 0
#define _dbg_A100 0
#define _dbg_allow_excessive_aa_clashes 0
#define _dbg_anemia 0
#define _dbg_appears 0
#define _dbg_asunder_atoms 0
#define _dbg_atom_mov_to_clash 0
#define _dbg_atom_pointing 0
#define _dbg_bb_clash_avoidance 0
#define _dbg_bb_contact_lonely 0
#define _dbg_bb_pairs 0
#define _dbg_bb_pullaway 0
#define _dbg_bb_realign 0
#define _dbg_bb_rots 0
#define _dbg_bb_scoring 0
#define _dbg_bconstr 0
#define _dbg_bridges 0
#define _dbg_can_rotate 0
#define _dbg_chirality_detection 0
#define _dbg_cond_basic 0
#define _dbg_cond_basic_acd_lig 0
#define _dbg_conj_chg 0
#define _dbg_conjugation 0
#define _dbg_contact_anomaly 0
#define _dbg_cs_pairing 0
#define _dbg_cvty_pose_filter 0
#define _dbg_eclipsing_contacts 0
#define _dbg_eclipses 0
#define _dbg_fitness_plummet 0
#define _dbg_flexion_selection 0
#define _dbg_groupsalign 0
#define _dbg_groupsel 0
#define _DBG_H2O_TELEPORT 0
#define _DBG_HISFLIP 0
#define _dbg_homology 0
#define _dbg_Huckel 0
#define _dbg_hxrax 0
#define _dbg_hydrogenate 0
#define _dbg_icreshape 0
#define _dbg_identify_rings 0
#define _dbg_imidazole_check 0
#define _dbg_improvements_only_residue 158
#define _dbg_improvements_only_rule 0
#define _dbg_intera_applicable 0
#define _dbg_interatomic_forces 0
#define _dbg_internal_energy 0
#define _dbg_linear_motion 0
#define _DBG_LONELINESS 0
#define _dbg_lsrb 0
#define _dbg_mand_conn 0
#define _DBG_MAX_CLASHES 0
#define _dbg_moieties 0
#define _DBG_MOLBB 0
#define _dbg_mol_flexion 0
#define _dbg_mol_frames 0
#define _dbg_molstruct_evolutions 0
#define _dbg_molstruct_evolution_bond_lengths 0
#define _dbg_monaxial 0
#define _dbg_multiflex 0
#define _dbg_null_flexions 0
#define _dbg_octant_occlusion 0
#define _dbg_optimize_hydrogens 0
#define _dbg_optimize_hydrogens_resno 260
#define _dbg_path_search 0
#define _dbg_peptide_bond_formation 0
#define _dbg_pocket_DeltaG_solv 0
#define _dbg_point_avg 0
#define _dbg_polar_calc 0
#define _dbg_polsat 0
#define _dbg_rap_resolve_special_atom 0
#define _dbg_repack 0
#define _DBG_RESBMULT 0
#define _dbg_residue_poses 0
#define _dbg_rock_pic 0
#define _dbg_rshpm_apply 0
#define _dbg_soft 0
#define _dbg_softpivot 0
#define _dbg_stays_assignment 0
#define _dbg_stays_enforce 0
#define _DBG_SPACEDOUT 0
#define _DBG_STEPBYSTEP 0
#define _DBG_TOOLARGE_DIFFNUMS 0
#define _DBG_TUMBLE_SPHERES 0
#define _dbg_unreciprocated_bonds 0
#define _dbg_worst_energy 0
#define _debug_active_bond_rot 0
#define debug_stop_after_tumble_sphere 0
#define _DORESPHRES 0
#define _dummy_atoms_for_debug 0
#define output_tumble_debug_docs 0
#define _show_cond_bas_hbond_energy 0
#define _show_final_group_pairs 0

#endif











