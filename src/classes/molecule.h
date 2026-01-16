
#include "intera.h"

#ifndef _MOLECULE
#define _MOLECULE

#include <vector>
#include <memory>

#define _solve_nonpol 2.1 / 84.5 / _kcal_per_kJ
#define _solve_np2pol -6.9 / 19.9 / _kcal_per_kJ

struct SMILES_Parenthetical
{
    Atom* startsfrom=0;
    char* smilesstr=0;
};

struct HistidineFlip
{
    Atom* H;
    Atom* C;
    Atom* N1;
    Atom* N2;
};

enum MovabilityType
{
    MOV_CAN_RECEN   =   0x8000,           // Molecule can move through space.
    MOV_CAN_AXIAL   =    0x800,           // Whole molecule can rotate in space.
    MOV_MC_AXIAL    =    0x400,           // Monte Carlo whole molecule rotations allowed.
    MOV_MUST_FLEX   =     0x80,           // Molecule's rotatable bonds are guaranteed free to rotate.
    MOV_MC_FLEX     =     0x40,           // Monte Carlo flexion is allowed.
    MOV_CAN_FLEX    =     0x20,           // Molecule's rotatable bonds may rotate if selected for flexion.
    MOV_ALL			=   0xfff0,
    MOV_NORECEN		=   0x0ff0,
    MOV_NOAXIAL     =   0xf0f0,
    MOV_FORCEFLEX   =     0xf0,
    MOV_FLEXONLY	=     0x70,
    MOV_PINNED      =     0x04,
    MOV_FLXDESEL    =     0x02,
    MOV_FORBIDDEN   =     0x0f,
    MOV_NONE		=     0x00,
    MOV_BKGRND      = 0x020000,
    MOV_CAN_CLASH   = 0x100000
};

enum MoleculeType
{
    MOLTYP_UNKNOWN,
    MOLTYP_LIGAND,
    MOLTYP_WATER,
    MOLTYP_AMINOACID
};

class Pose
{
public:
    Pose();
    Pose(Molecule* from_mol);
    Pose(const Pose& p);
    Pose& operator=(const Pose& p);
    Pose(Pose&& p) noexcept;
    Pose& operator=(Pose&& p) noexcept;
    ~Pose();
    void copy_state(Molecule* from_mol);
    void restore_state(Molecule* to_mol);
    void restore_state_relative(Molecule* to_mol, const char* rel_atom_name = "");
    float total_atom_motions();
    void reset();
    bool has_data();
    void deallocate();

protected:
    int sz = 0;
    Point* saved_atom_locs = nullptr;
    int* saved_atom_Z = nullptr;
    char** saved_atom_name = nullptr;
    int nsaved_atom = 0;
    Molecule* saved_from = nullptr;
};


class Molecule
{
    friend class Pose;

public:
    Molecule();
    Molecule(const char* name);
    Molecule(const char* name, Atom** collection);
    Molecule(const Molecule& copyfrom);
    Molecule(Molecule** monomers);
    virtual ~Molecule();
    void set_name(const char* new_name)
    {
        name = new char[strlen(new_name)+2];
        strcpy(name, new_name);
    }
    void add_monomer(Molecule* toadd);
    void make_multimer(int n);
    // TODO: linkages between monomers.
    Molecule* get_monomer(int i);

    // Load and save functions.
    int from_sdf(const char* sdf_dat);		// returns number of atoms loaded.
    bool save_sdf(FILE* outf);
    bool save_sdf(FILE* outf, Molecule** included_ligands);
    void save_pdb(FILE* outf, int atomno_offset=0, bool endpdb = true);
    int from_pdb(FILE* inf, bool het_only = false);  // returns number of atoms loaded.
    void identify_acidbase();				// called within every load.
    bool from_smiles(char const * smilesstr, bool use_parser = true);
    void clear_cache();

    // Getters.
    const char* get_name() const
    {
        return name;
    }
    int get_atom_count() const
    {
        return atcount;
    }
    int get_heavy_atom_count() const;
    int get_bond_count(bool unidirectional) const;
    float get_molecular_wt();
    Atom* get_nearest_atom(Point loc) const;
    Atom* get_nearest_atom(Point loc, intera_type capable_of, int sgn_polarity = any_element) const;
    Atom* get_nearest_atom_to_line(Point A, Point B) const;
    Atom* get_farthest_atom(Point loc) const;
    Point get_bounding_box() const;				// Return the +x+y+z vertex of a bounding box, including vdW radii, if center={0,0,0}.
    float get_volume();
    float get_surface_area(bool polaronly = false, bool overwrite_atom_areas = true);
    float get_exposed_surface_area(Molecule** neighbors, bool polaronly = false, bool overwrite_atom_areas = true);
    float get_charge() const;
    int is_residue();
    bool is_thiol();
    bool is_water();
    float pi_stackability(bool include_backbone = false);

    // Spatial functions.
    Point get_barycenter(bool bond_weighted = false) const;
    virtual void move(Vector move_amt, bool override_residue = false);
    virtual void move(Point move_amt, bool override_residue = false);
    virtual void recenter(Point new_location);
    const Point* obtain_vdW_surface(float density);
    Atom** get_vdW_vertex_atoms() { return vdw_vertex_atom; }
    int get_vdW_vertex_count() { return vdw_vertex_count; }
    int get_atom_vdW_vertex_count(Atom* a);
    void wipe_vdw_surface();
    void rotate(Vector* Vector, float theta, bool bond_weighted = false);
    void rotate(LocatedVector vec, float theta);
    void rotate(Rotation rot);
    void rotate(Rotation rot, Point origin);
    bool shielded(Atom* a, Atom* b) const;
    float correct_structure(int iters = 500);
    float close_loop(Atom** path, float closing_bond_cardinality);
    float sum_interatomic_distances();
    bool is_chiral();
    void mirror();
    float total_eclipses();
    void crumple(float theta);					// Randomly rotate all rotatable bonds by +/- the specified angle.
    float distance_to(Molecule* other_mol);
    std::vector<Atom*> longest_dimension();
    float get_atom_bond_length_anomaly(Atom* atom, Atom* ignore = nullptr);
    float refine_structure(int generations = _evolution_default_generations, float mutation_rate = _default_mutation_rate, int pop_size = _default_population_size);
    int atoms_inside_sphere(Sphere container, bool* byindex, float radius_multiplier = 1);     // If byindex is not null, sets byindex[n] to true for atoms inside the sphere, but does not set to false.
    float surface_occlusion(Molecule** ligands);
    float surface_occlusion(Molecule* ligand);
    float octant_occlusion(Molecule** ligands, bool ignore_polar = false);
    float octant_occlusion(Molecule* ligand, bool ignore_polar = false);

    // Atom functions.
    Atom* add_atom(const char* elemsym, const char* aname, Atom* bond_to, const float bcard);
    Atom* add_atom(const char* elemsym, const char* aname, const Point* location, Atom* bond_to, const float bcard, const int charge = 0);
    void add_existing_atom(Atom* to_add);
    char** get_atom_names() const;
    Atom* get_atom(const char* aname) const;
    Atom* get_atom(const int a_idx) const
    {
        return atoms[a_idx];
    }
    int count_atoms_by_element(const char* esym);
    Point get_atom_location(const char* aname);
    Point get_atom_location(int idx);
    int atom_idx_from_ptr(Atom* a);
    void delete_atom(Atom* a);
    void delete_all_atoms();
    int get_hydrogen_count();
    virtual void hydrogenate(bool steric_only = false);
    virtual void dehydrogenate();
    void clear_atom_binding_energies();
    int has_hbond_donors();
    int has_hbond_acceptors();                    // N+ is not an h-bond acceptor.
    int has_pi_atoms(bool include_backbone = false);
    bool protonate();
    bool deprotonate();
    void propagate_stays();

    // Bond functions.
    Bond** get_rotatable_bonds(bool include_can_flip = true);
    Bond** get_all_bonds(bool unidirectional);
    void clear_all_bond_caches();					// Call this any time you add or remove an atom.
    bool rotate_bond(const Bond* rot8b, const float angle);
    void do_histidine_flip(HistidineFlip* hf);
    void identify_conjugations();
    bool check_Greek_continuity();

    // Ring functions.
    int identify_rings();
    int get_num_rings() const;
    int add_ring(Atom** atoms);
    bool ring_is_coplanar(int ringid);
    bool ring_is_aromatic(int ringid) const;
    Point get_ring_center(int ringid);
    Vector get_ring_normal(int ringid);
    Atom** get_ring_atoms(int ringid);
    int get_ring_num_atoms(int ringid);
    void identify_cages();

    // Interaction functions.
    float get_internal_clashes(bool subtract_baseline = false);
    void minimize_internal_clashes();
    float get_base_clashes() { return base_internal_clashes; }
    float get_intermol_clashes(Molecule* ligand);
    float get_intermol_clashes(Molecule** ligands);
    static float total_intermol_clashes(Molecule** ligands);
    float get_worst_clash();
    Interaction get_intermol_binding(Molecule* ligand, bool subtract_clashes = true, bool priority_boost = false);
    Interaction get_intermol_binding(Molecule** ligands, bool subtract_clashes = true, bool priority_boost = false);
    float get_intermol_potential(Molecule* ligand, bool disregard_distance = false);
    float get_intermol_potential(Molecule** ligands, bool disregard_distance = false);
    float hydrophilicity();
    float solvent_free_energy(float epsilon = 80.0, float kappa = 0.316, bool compute_surface_area = true);
    float solvent_bound_energy(Molecule** neighbors);
    float get_intermol_polar_sat(Molecule* ligand);
    float get_intermol_contact_area(Molecule* ligand, bool hydrophobic_only = false);
    void mutual_closest_atoms(Molecule* mol2, Atom** atom1, Atom** atom2);
    void mutual_closest_hbond_pair(Molecule* mol2, Atom** atom1, Atom** atom2);
    float get_total_mclashes();
    Interaction optimize_intermol_contact(Molecule* ligand);

    #if compute_vdw_repulsion
    float get_vdW_repulsion(Molecule* ligand);
    #endif

    float bindability_by_type(intera_type type, bool include_backbone = false);

    static Interaction total_intermol_binding(Molecule** ligands);

    static void conform_molecules(Molecule** molecules, int iterations = 50,
        void (*callback)(int, Molecule**) = nullptr,
        void (*progress)(float) = nullptr,
        int min_iter = 0
        );
    
    static void conform_molecules(Molecule** molecules, Molecule** background, int iterations = 50,
        void (*callback)(int, Molecule**) = nullptr,
        void (*progress)(float) = nullptr
        );

    void quick_conform(Molecule** background, int iterations = 25);

    void conform_atom_to_location(int atom_idx, Point target, int iterations = 50, float optimal_distance = 0);
    void conform_atom_to_location(const char* atom_name, Point target, int iterations = 20, float optimal_distance = 0);
    void conform_atom_to_location(Atom* a, Atom* target, int iterations = 50, float optimal_distance = 0);
    Vector motion_to_optimal_contact(Molecule* ligand);

    // Returns the sum of all possible atom-molecule interactions if all distances and anisotropies were somehow optimal.
    float get_atom_mol_bind_potential(Atom* a);
    float find_mutual_max_bind_potential(Molecule* other);
    bool check_stays();
    bool check_stays_dry();
    void enforce_stays(float amount=1, void (*stepscb)(std::string mesg) = nullptr);

    float get_springy_bond_satisfaction();

    void reset_conformer_momenta();
    Atom** get_most_bindable(int max_num = 3);						// Return the atoms with the greatest potential intermol binding.
    Atom** get_most_bindable(int max_num, Atom* for_atom);

    void allocate_mandatory_connections(int mcmax);
    void add_mandatory_connection(Molecule* addmol);
    void remove_mandatory_connection(Molecule* rmvmol);
    void zero_mandatory_connection_cache();
    void delete_mandatory_connections();

    // Debug stuff.
    #if debug_break_on_move
    void set_atoms_break_on_move(bool break_on_move)
    {
        if (atoms)
        {
            int i;
            for (i=0; atoms[i]; i++) atoms[i]->break_on_move = break_on_move;
        }
    }
    #endif

    bool echo_iters = false;
    MovabilityType movability = MOV_ALL;
    float lastbind = 0;
    float lastbind_history[10];
    float lastshielded = 0;
    float lastmc = 0;
    HistidineFlip** hisflips = nullptr;
    Bond* springy_bonds = nullptr;
    int springy_bondct = 0;
    bool been_flexed = false;
    bool priority = false;
    Atom *coordmtl = nullptr;
    Molecule** mclashables = nullptr;
    Atom *clash1 = nullptr, *clash2 = nullptr;
    Atom *best_intera = nullptr, *best_other_intera = nullptr;
    Molecule* best_interactor = nullptr;
    int eclipse_hash = 0;
    Atom *stay_close_mine = nullptr, *stay_close_other = nullptr;
    Atom *stay_close2_mine = nullptr, *stay_close2_other = nullptr;
    Molecule *stay_close_water = nullptr, *stay_close_mol = nullptr, *stay_close2_mol = nullptr;
    float stay_close_tolerance = 0, stay_close_optimal = 2, stay_close2_optimal = 2;
    bool is_ic_res = false;

protected:

    Atom** atoms = 0;
    int atcount = 0;
    char* name = 0;
    char* smiles = 0;
    Molecule** monomers = nullptr;
    int nmonomers = 0;
    float clash_worst = 0;
    Point* vdw_surface = nullptr;
    Atom** vdw_vertex_atom = nullptr;
    int vdw_vertex_count = 0;
    Atom*** paths = nullptr;
    Ring** rings = nullptr;
    Bond** rotatable_bonds = nullptr;
    bool immobile = false;
    bool doing_bkbend = false;
    float base_internal_clashes = 0;					// Baseline computed internal clashes due to unavoidably close atoms.
    float base_eclipses = 0;
    std::string sdfgen_aboutline = "";
    Molecule** mandatory_connection = nullptr;
    float* last_mc_binding = nullptr;
    Atom** most_bindable = nullptr;
    Pose* iterbegan = nullptr;
    int iters_without_change = 0;

    // For intermol conformer optimization:
    float lmx=0,lmy=0,lmz=0;			// Linear momentum xyz.
    float amx=0,amy=0,amz=0;			// Angular momentum xyz.

    bool from_smiles(char const * smilesstr, Atom* ipreva);
    int smlen = 0;
    SMILES_Parenthetical* paren;
    int spnum = 0;
    MoleculeType mol_typ = MOLTYP_UNKNOWN;

    void find_paths();
    int path_contains_atom(int path_idx, Atom* a);
    int path_get_length(int path_idx);
    Atom* path_get_terminal_atom(int path_idx);
    void copy_path(int old_idx, int new_idx);
    bool path_is_subset_of(int short_path, int long_path);
    void echo_path(int idx);
    int aidx(Atom* a);
    void reallocate();
    void atoms_from_multimers();
    float fsb_lsb_anomaly(Atom* first, Atom* last, float lcard, float bond_length);
    void make_coplanar_ring(Atom** ring_members, int ringid);
    void recenter_ring(int ringid, Point new_ring_cen);
    void rotate_ring(int ringid, Rotation rot);
    bool in_same_ring(Atom* a, Atom* b);
    float get_atom_error(int atom_idx, LocatedVector* best_lv, bool hemispherical = true);
    Interaction intermol_bind_for_multimol_dock(Molecule* othermol, bool allow_clash);
    static Interaction cfmol_multibind(Molecule* mol, Molecule** nearby_mols);
    bool faces_any_ligand(Molecule** ligands);

    public:
    const int& num_monomers = nmonomers;
};

float g_total_mclash(void* mol);
bool mclash_delta(void* mol, float previous_mclashes);

extern float conformer_momenta_multiplier;
extern float conformer_tumble_multiplier;
extern bool allow_ligand_360_tumble;
extern bool allow_ligand_360_flex;
extern float cavity_stuffing;
extern float clash_fleeing;
extern float _momentum_rad_ceiling;
extern Molecule *worst_clash_1, *worst_clash_2;
extern float worst_mol_clash;
extern Molecule global_water;
extern FILE* audit;
extern bool cfmols_have_metals;

#if _dbg_improvements_only_rule
extern Molecule** check_mols;
extern Molecule* check_ligand;
extern bool excuse_deterioration;
#endif

#endif

