#include "protein.h"
#include "space.h"

#ifndef _CAVITY
#define _CAVITY

class CPartial : public SPartial
{
    public:
    // int resno = 0;
    std::string resnos_as_string(Protein* p);
    int resnos_as_array(Protein* p, int* output);
    int from_cvty_line(char* lndata);               // Returns the cavity number from the first column.
    void write_cvty_line(char* outdata, int cno, Protein* protein);
};

class Cavity : public Space
{
    public:
    static int scan_in_protein(Protein* p, Cavity* results, int results_max, Progressbar* pgb = nullptr);
    float molecule_inside_pocket(Molecule* m, bool match_attributes = false);
    float cavity_filling(Molecule* m);
    float containment_violations(Molecule* m, float stop_if_more_than = -1);
    float find_best_containment(Molecule* m, bool match_binding_types = false);
    float match_ligand(Molecule* ligand, Atom** match_atom = nullptr, CPartial** match_partial = nullptr, Protein* prot = nullptr);
    int resnos(Protein* p, AminoAcid** result);
    std::string resnos_as_string(Protein* p);
    int resnos_as_array(Protein* p, int* output);
    Protein* prot = nullptr;
    float cavity_intersection(Cavity* other);
    void unify(Cavity* cavfrom);
    int estimate_multiplicity(Molecule* ligand);
    CPartial* get_partial_by_idx(int idx) { return &partials[idx]; }

    protected:
    CPartial* partials = nullptr;
};

extern float cav_xmax, cav_xmin, cav_ymax, cav_ymin, cav_zmax, cav_zmin, cav_xyrlim, cav_xzrlim, cav_yzrlim;
extern int cav_resmin, cav_resmax;

#endif

