
#include "aminoacid.h"
#include "progress.h"

#ifndef _PROTEIN
#define _PROTEIN

#include <string>
#include <climits>

#define resplc_rgnstart -1
#define resplc_rgnend -99

enum region_source
{
    rgn_none,
    rgn_pdb,
    rgn_manual
};

class Region
{
    public:
    int start=0;
    int end=0;
    std::string name="";

    Point start_CA_location(Protein* p);
    Point end_CA_location(Protein* p);
};

class BallesterosWeinstein
{
    public:
    int helix_no = 0;
    int member_no = 0;

    BallesterosWeinstein() { ; }
    BallesterosWeinstein(int b, int w) { helix_no = b; member_no = w; }
    BallesterosWeinstein(const char* fromstr) { this->from_string(fromstr); }

    void from_string(const char* inpstr);
    std::string to_string();
};

class ResiduePlaceholder
{
    public:
    int node = 0;
    int resno = 0;
    int hxno = 0;
    int bwpos = 0;
    std::string bw;
    std::string allowed_aas;

    void set(const char* str);
    void resolve_resno(Protein* prot);

    protected:
    Protein* m_prot = nullptr;
};

class ResidueAtomPlaceholder : public ResiduePlaceholder
{
    protected:
    std::string aname;

    public:
    void set(const char* str);
    Point loc();
    Atom* atom();
    bool resolve_special_atom(Protein* p, Point rel);
    std::string get_aname();
    std::string get_orig_aname() { return aname; }
};

struct MCoord
{
    int Z = 29;
    int charge = 2;
    Atom* mtl = nullptr;
    Point mtl_original_location;
    ResiduePlaceholder coordres[16];
    int ncoordres = 0;
};

struct AARenumber
{
    AminoAcid* aa = nullptr;
    AminoAcid* replace_with = nullptr;
    int new_resno = 0;
};

class Protein
{
public:
    // Constructors.
    Protein();
    Protein(const char* name);
    ~Protein();

    // Build functions.
    bool add_residue(const int resno, const char aaletter);
    bool add_sequence(const char* sequence);
    bool add_residue(const char* pdbdata);
    void set_clashables(int resno = -1, bool recursed = false);
    void delete_residue(int resno);
    void delete_sidechain(int resno);
    void delete_residues(int startres, int endres);
    void delete_sidechains(int startres, int endres);
    MCoord* coordinate_metal(MCoord* mtlcoords, int count);
    void set_region(std::string name, int start, int end);
    void set_bw50(int helixno, int resno);
    void renumber_residues(int startres, int endres, int new_startres);
    bool disulfide_bond(int resno1, int resno2);

    // Serialization.
    void set_name_from_pdb_name(const char* pdb_name);
    int load_pdb(FILE* infile, int resno_offset = 0, char chain = 'A');				// Returns number of residues loaded.
    void save_pdb(FILE* outfile, Molecule* ligand = nullptr);
    void end_pdb(FILE* outfile);
    void revert_to_pdb();
    void copy_mcoords(Protein* copyfrom);

    // Getters.
    int get_seq_length();
    int get_start_resno();
    int get_end_resno();
    std::string get_sequence();
    Molecule* metals_as_molecule();
    int get_metals_count();
    bool has_bwnos();
    AminoAcid* get_residue(int resno);
    AminoAcid* get_residue(BallesterosWeinstein bw);
    AminoAcid* get_residue_bw(int helixno, int bwno);
    AminoAcid* get_residue_bw(const char* bwno);
    AminoAcid** clashable_residues(int resno) { return res_can_clash[resno]; }
    BallesterosWeinstein get_bw_from_resno(int resno);
    Region get_region(std::string name);
    const Region* get_regions() { return regions; }
    int get_region_end(std::string name);
    int get_region_end(int hxno);
    int get_region_start(std::string name);
    int get_region_start(int hxno);
    bool aa_ptr_in_range( AminoAcid* aaptr );
    Atom* get_atom(int resno, const char* aname);
    std::string get_name()
    {
        return std::string(name);
    }
    Point get_atom_location(int resno, const char* aname);
    Atom* get_nearest_atom(Point pt, int start_resno = 1, int end_resno = RAND_MAX);
    void add_remark(const char* remark);
    void add_remark(std::string new_remark);
    std::vector<std::string> get_remarks(std::string search_for = "");
    int get_bw50(int helixno);
    int search_sequence(const int start_resno, const int end_resno, const char* search_for, const int threshold = -1, int* similarity = nullptr);

    char get_pdb_chain() const { return pdbchain; }
    char set_pdb_chain(char chain);

    // Metrics functions.
    float get_internal_clashes(int start_resno = 0, int end_resno = 0, bool repack = false, int repack_iters = 10);
    float get_rel_int_clashes();
    Interaction get_internal_binding();
    float get_intermol_clashes(Molecule* ligand);
    float get_intermol_clashes(Molecule* ligand, int startres, int endres);
    Interaction get_intermol_binding(Molecule* ligand);
    AminoAcid** get_residues_can_clash(int resno);
    std::vector<AminoAcid*> get_residues_can_clash(int start_resno, int end_resno);
    int get_residues_can_clash_ligand
    (	AminoAcid** reaches_spheroid,
        Molecule* ligand,
        const Point nodecen,
        const Point search_size,
        const int* addl_resno = nullptr,
        bool ignore_priority = false,
        const int* excl_resno = nullptr
    );

    int fetch_residues_near(Point pt, float max_distance, AminoAcid** results, bool facing=true);
    std::vector<AminoAcid*> get_contact_residues(Protein* other_prot, float contact_distance = 2.5);
    Molecule** all_residues_as_molecules();
    Molecule** all_residues_as_molecules_except(Molecule** mm);
    Point get_region_center(int startres, int endres);
    Vector get_region_axis(int startres, int endres);
    float get_helix_orientation(int startres, int endres);
    float aa_angle_about_helical_axis(int resno);
    Point find_loneliest_point(Point search_center, Point search_size = Point(_INTERA_R_CUTOFF*2, _INTERA_R_CUTOFF*2, _INTERA_R_CUTOFF*2), float* pminr = nullptr);
    float total_mclashes();
    Point estimate_pocket_size(AminoAcid** ba);
    Point estimate_pocket_size(std::vector<AminoAcid*> ba);
    Interaction binding_to_nearby_residues(int resno);
    void minimize_residue_clashes(int resno);
    float region_can_move(int startres, int endres, Vector direction, bool repack = false, int ignore_startres = 0, int ignore_endres = 0);
    float region_can_rotate(int startres, int endres, LocatedVector axis, bool repack = false, float extra_clash_allowance = 0, int ignore_startres = 0, int ignore_endres = 0);     // Searches positive theta.
    void region_optimal_positioning(int startres, int endres, Vector* output_transformation, Rotation* output_rotation, Protein** other_strands = nullptr);
    void set_conditional_basicities();
    float A100();
    Atom* region_pivot_atom(Region region, Atom** other_atom = nullptr);
    Point get_region_bounds(int startres, int endres);
    float optimize_hydrogens(int start_resno = 1, int end_resno = INT_MAX, int* force_resnos = nullptr);

    // Motion functions
    void upright();
    void move_piece(int start_res, int end_res, Point new_center);
    void move_piece(int start_res, int end_res, Vector move_amt);
    LocRotation rotate_piece(int start_res, int end_res, int align_res, Point align_target, int pivot_res = 0);		// If no pivot res, rotate about the center.
    LocRotation rotate_piece(int start_res, int end_res, Rotation rot, int pivot_res);
    LocRotation rotate_piece(int start_res, int end_res, Point origin, Vector axis, float theta);
    LocRotation rotate_piece(int start_res, int end_res, LocRotation lr);

    void rotate_backbone(int residue_no, bb_rot_dir direction, float angle);
    void conform_backbone(int startres, int endres, Atom* a, Point target, int iters = 50);
    void rotate_backbone_partial(int startres, int endres, bb_rot_dir direction, float angle);
    void conform_backbone(int startres, int endres, int iters = 50, bool backbone_atoms_only = false);
    void conform_backbone(int startres, int endres, Atom* a1, Point target1, Atom* a2, Point target2, Atom* a3, Point target3, int iters = 50);
    void conform_backbone(int startres, int endres, Atom* a1, Point target1, Atom* a2, Point target2, int iters = 50, bool backbone_atoms_only = false);
    void conform_backbone(int startres, int endres,
                          Atom* a1, Point target1,
                          Atom* a2, Point target2,
                          Atom* a3, Point target3,
                          int iters, bool backbone_atoms_only
                         );

    void backconnect(int startres, int endres);
    void find_residue_initial_bindings();
    void undo();

    // Secondary structure
    void make_helix(int startres, int endres, float phi, float psi);
    void make_helix(int startres, int endres, int stopat, float phi, float psi);
    float orient_helix
    (	int startres, int endres,						// Boundaries of helix.
        int stopat,										// Last residue to move with helix.
        float angle,									// 0 = horizontal; positive = ascending (+Y) with increasing resno.
        int iterations
    );

    // Homology
    void homology_conform(Protein* target_structure, Protein* reference_structure);
    void bridge(int resno1, int resno2);
    int replace_side_chains_from_other_protein(Protein* other);

    // Pre-placement
    float tumble_ligand_inside_pocket(Molecule* ligand, Point pocketcen,
        float clashmult = 1, Progressbar* pgb = nullptr, Point* ligcen = nullptr);

    int mcoord_resnos[32];

    Vector last_uprighted_xform;
    LocRotation last_uprighted_A, last_uprighted_B;
    Vector last_int_clash_dir;

    AminoAcid *stop1, *stop2;
    Atom *stop1a, *stop2a;
    int last_saved_atom_number = 0;
    Point pocketcen;

protected:
    Atom** ca = nullptr;
    int arrlimit = 0;
    std::string name;
    char* sequence = nullptr;
    AminoAcid** residues = nullptr;
    AminoAcid*** res_can_clash = nullptr;
    float* res_reach = nullptr;
    Atom* metals[16];
    int metcount = 0;
    Star aaptrmin, aaptrmax;
    float initial_int_clashes = 0;
    Region regions[PROT_MAX_RGN];
    region_source regions_from = rgn_none;
    char** remarks = nullptr;
    int remarksz = 0;
    MCoord m_mcoords[16];
    int nm_mcoords = 0;
    int Ballesteros_Weinstein[79];
    std::vector<AABridge> aabridges;
    std::vector<Bond*> connections;
    std::vector<Pose> origpdb_residues;
    char pdbchain = ' ';
    Pose** undo_poses = nullptr;
    bool mass_undoable = false;

    int* get_residues_in_reach(int resno);
    float get_coord_anomaly(Atom* metal, AminoAcid* coord_res);
    void allocate_undo_poses();
    void save_undo_state();
};

std::ostream& operator<<(std::ostream& os, const BallesterosWeinstein& bw);

extern float *g_rgnxform_r, *g_rgnxform_theta, *g_rgnxform_y, *g_rgnrot_alpha, *g_rgnrot_w, *g_rgnrot_u;

#endif

