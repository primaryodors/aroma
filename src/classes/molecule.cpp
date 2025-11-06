
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <ctime>
#include <regex>
#include <stdint.h>
#include "molecule.h"
#include "aminoacid.h"

#define _DBGCLSLOOP 0

using namespace std;

float conformer_momenta_multiplier = 1;
float conformer_tumble_multiplier = 1;

float cavity_stuffing = default_cavity_stuffing;
float clash_fleeing = lmpush;

bool allow_ligand_360_tumble = true;
bool allow_ligand_360_flex = true;

bool cfmols_have_metals = false;

Molecule *worst_clash_1 = nullptr, *worst_clash_2 = nullptr;
float worst_mol_clash = 0;
FILE* audit = nullptr;

Molecule global_water("H2O");

#if _dbg_improvements_only_rule
Molecule** check_mols = nullptr;
Molecule* check_ligand = nullptr;
bool excuse_deterioration = false;
#endif

float _momentum_rad_ceiling = fiftyseventh * 30;

Molecule::Molecule(char const* lname)
{
    name = new char[strlen(lname)+1];
    strcpy(name, lname);

    atoms = nullptr;
    smiles = nullptr;
    rings = nullptr;
    atcount = 0;
    reset_conformer_momenta();
    rotatable_bonds = nullptr;

    int j;
    for (j=0; j<10; j++) lastbind_history[j] = 0;
}

Molecule::Molecule()
{
    atoms = nullptr;
    smiles = nullptr;
    rings = nullptr;
    atcount = 0;
    reset_conformer_momenta();
    rotatable_bonds = 0;
    paren = nullptr;

    int j;
    for (j=0; j<10; j++) lastbind_history[j] = 0;
}

Molecule::~Molecule()
{
    if (atoms)
    {
        delete[] atoms;
    }
    if (smiles) delete[] smiles;
    if (rings) delete[] rings;
    if (rotatable_bonds) delete[] rotatable_bonds;
    if (vdw_surface) delete[] vdw_surface;
    if (vdw_vertex_atom) delete[] vdw_vertex_atom;
}

int length(Atom** array)
{
    int numAtoms;
    for (numAtoms=0; array[numAtoms]; numAtoms++);	// Get count.
    return numAtoms;
}

bool hasAtoms(Atom** array)
{
    if (array == nullptr)
        return false;
    return array[0] != nullptr;
}

bool noAtoms(Atom** array)
{
    return !hasAtoms(array);
}

Molecule::Molecule(char const* lname, Atom** collection)
{
    if (!collection)
    {
        cout << "Temporary molecule creation attempted from nullptr atom pointer array." << endl;
        throw 0xbadca22;
    }

    int numAtoms = length(collection);
    atoms = new Atom*[numAtoms + 1];
    for (int j=0; j < numAtoms; j++)
        atoms[j] = collection[j];
    atoms[numAtoms] = 0;
    atcount = numAtoms;

    name = new char[strlen(lname)+1];
    strcpy(name, lname);

    smiles = nullptr;
    rings = nullptr;
    reset_conformer_momenta();
    identify_conjugations();
    rotatable_bonds = 0;
}

Molecule::Molecule(const Molecule &o)
{
    if (!o.atcount) return;

    echo_iters = o.echo_iters;
    movability = o.movability;
    priority = o.priority;
    is_ic_res = o.is_ic_res;
    if (o.smiles)
    {
        smiles = new char[strlen(o.smiles)+4];
        strcpy(smiles, o.smiles);
    }
    immobile = o.immobile;
    base_internal_clashes = o.base_internal_clashes;
    base_eclipses = o.base_eclipses;
    sdfgen_aboutline = o.sdfgen_aboutline;

    atoms = new Atom*[o.atcount+4];
    int i;
    for (i=0; i<o.atcount; i++)
    {
        atoms[i] = new Atom(*o.atoms[i]);
    }
    atoms[i] = nullptr;
    atcount = i;

    for (i=0; i<o.atcount; i++)
    {
        int n = o.atoms[i]->get_geometry();
        int j;
        for (j=0; j<n; j++)
        {
            Bond* b = o.atoms[i]->get_bond_by_idx(j);
            if (!b || !b->atom2) continue;
            Atom* a = get_atom(b->atom2->name);
            if (!a) continue;
            if (atoms[i]->is_bonded_to(a)) continue;
            atoms[i]->bond_to(a, b->cardinality);
        }
    }

    if (o.hisflips)
    {
        int n;
        for (n=0; o.hisflips[n]; n++);
        hisflips = new HistidineFlip*[n+4];
        for (i=0; o.hisflips[i]; i++)
        {
            hisflips[i] = new HistidineFlip();
            Atom* a = o.hisflips[i]->C;
            Atom* b = get_atom(a->name);
            hisflips[i]->C = b;
            a = o.hisflips[i]->H;
            b = get_atom(a->name);
            hisflips[i]->H = b;
            a = o.hisflips[i]->N1;
            b = get_atom(a->name);
            hisflips[i]->N1 = b;
            a = o.hisflips[i]->N2;
            b = get_atom(a->name);
            hisflips[i]->N2 = b;
        }
    }

    identify_acidbase();
    identify_rings();
    identify_cages();
    identify_conjugations();
}

Molecule::Molecule(Molecule **m)
{
    if (!m) return;
    int i, n;
    for (n=0; m[n]; n++);
    monomers = new Molecule*[n+4];
    for (i=0; i<n; i++) monomers[i] = m[i];
    monomers[i] = nullptr;
    nmonomers = n;
}

Pose::Pose()
{
    reset();
}

Pose::Pose(Molecule* m)
{
    reset();
    copy_state(m);
}

Pose::Pose(const Pose &p)
{
    sz = p.sz;
    saved_atom_locs = new Point[sz+4];
    saved_atom_Z = new int[sz+4];
    saved_atom_name = new char*[sz+4];
    int i;
    for (i=0; i<sz; i++)
    {
        saved_atom_locs[i] = p.saved_atom_locs[i];
        saved_atom_Z[i] = p.saved_atom_Z[i];
        saved_atom_name[i] = new char[16];
        strcpy(saved_atom_name[i], p.saved_atom_name[i]);
    }
    nsaved_atom = p.nsaved_atom;
    saved_from = p.saved_from;
}

Pose::Pose(Pose &&p) noexcept
{
    sz = p.sz;
    saved_atom_locs = new Point[sz+4];
    saved_atom_Z = new int[sz+4];
    saved_atom_name = new char*[sz+4];
    int i;
    for (i=0; i<sz; i++)
    {
        saved_atom_locs[i] = p.saved_atom_locs[i];
        saved_atom_Z[i] = p.saved_atom_Z[i];
        saved_atom_name[i] = new char[16];
        strcpy(saved_atom_name[i], p.saved_atom_name[i]);
    }
    nsaved_atom = p.nsaved_atom;
    saved_from = p.saved_from;
}

Pose &Pose::operator=(const Pose &p)
{
    if (this != &p)
    { 
        sz = p.sz;
        saved_atom_locs = new Point[sz+4];
        saved_atom_Z = new int[sz+4];
        saved_atom_name = new char*[sz+4];
        int i;
        for (i=0; i<sz; i++)
        {
            saved_atom_locs[i] = p.saved_atom_locs[i];
            saved_atom_Z[i] = p.saved_atom_Z[i];
            saved_atom_name[i] = new char[16];
            strcpy(saved_atom_name[i], p.saved_atom_name[i]);
        }
        nsaved_atom = p.nsaved_atom;
        saved_from = p.saved_from;
    }
    return *this;
}

Pose &Pose::operator=(Pose &&p) noexcept
{
    if (this != &p)
    { 
        sz = p.sz;
        saved_atom_locs = new Point[sz+4];
        saved_atom_Z = new int[sz+4];
        saved_atom_name = new char*[sz+4];
        int i;
        for (i=0; i<sz; i++)
        {
            saved_atom_locs[i] = p.saved_atom_locs[i];
            saved_atom_Z[i] = p.saved_atom_Z[i];
            saved_atom_name[i] = new char[16];
            strcpy(saved_atom_name[i], p.saved_atom_name[i]);
        }
        nsaved_atom = p.nsaved_atom;
        saved_from = p.saved_from;
    }
    return *this;
}

Pose::~Pose()
{
    deallocate();
}

void Pose::deallocate()
{
    if (saved_atom_locs) delete[] saved_atom_locs;
    saved_atom_locs = nullptr;
    if (saved_atom_Z) delete[] saved_atom_Z;
    saved_atom_Z = nullptr;
    if (saved_atom_name)
    {
        int i;
        for (i=0; i<sz; i++)
        {
            if (saved_atom_name[i]) delete[] saved_atom_name[i];
        }
        delete[] saved_atom_name;
    }
    saved_atom_name = nullptr;
}

void Pose::reset()
{
    sz = 0;
    nsaved_atom = 0;
    saved_from = nullptr;
}

bool Pose::has_data()
{
    return nsaved_atom;
}

void Pose::copy_state(Molecule* m)
{
    if (sz != m->get_atom_count()) deallocate();
    sz = m->get_atom_count();
    int i;
    if (saved_from != m) deallocate();
    if (!nsaved_atom || saved_from != m)
    {
        saved_from = m;
        if (!m || !m->atoms) return;

        saved_atom_locs = new Point[sz+4];
        saved_atom_Z = new int[sz+4];
        saved_atom_name = new char*[sz+4];
        for (i=0; i<sz; i++)
        {
            saved_atom_name[i] = new char[16];
            saved_atom_name[i][0] = 0;
        }
    }

    for (i=0; m->atoms[i] && i<sz; i++)
    {
        saved_atom_locs[i] = m->atoms[i]->loc;
        saved_atom_locs[i].weight = m->atoms[i]->get_atomic_weight();       // Important for restore_state_relative.
        saved_atom_Z[i] = m->atoms[i]->Z;
        strcpy(saved_atom_name[i], m->atoms[i]->name);
    }
    sz = i;
}

void Pose::restore_state(Molecule* m)
{
    if (!m || !m->atoms || !sz) return;
    int i, n;
    if (m != saved_from)
    {
        n = nsaved_atom;
        for (i=0; i<n; i++)
        {
            if (i == n-1 && !m->atoms[i] && (saved_atom_Z[i] == 1)) break;
            if (!m->atoms[i] && !saved_atom_Z[i]) break;
            if (/*n != sz ||*/ !m->atoms[i] || (saved_atom_Z[i] != m->atoms[i]->Z) || strcmp(saved_atom_name[i], m->atoms[i]->name) )
            {
                cout << "Attempt to restore pose to incompatible molecule (from " << saved_from->name << " to " << m->name << ")." << endl;
                if (m->is_residue())
                {
                    Star s;
                    s.pmol = m;
                    if (s.paa->conditionally_basic()) return;
                }
                throw -4;
            }
        }

        saved_from = m;
    }

    for (i=0; i<sz && m->atoms[i]; i++)
    {
        m->atoms[i]->move(saved_atom_locs[i]);
    }
}

void Pose::restore_state_relative(Molecule *m, const char *ran)
{
    if (!m || !m->atoms || !sz) return;
    int i, n;    
    if (m != saved_from)
    {
        n = nsaved_atom;
        for (i=0; i<n; i++)
        {
            if (i == n-1 && !m->atoms[i] && (saved_atom_Z[i] == 1)) break;
            if (!m->atoms[i] && !saved_atom_Z[i]) break;
            if (/*n != sz ||*/ !m->atoms[i] || (saved_atom_Z[i] != m->atoms[i]->Z) || strcmp(saved_atom_name[i], m->atoms[i]->name) )
            {
                cout << "Attempt to restore pose to incompatible molecule (from " << saved_from->name << " to " << m->name << ")." << endl;
                if (m->is_residue())
                {
                    Star s;
                    s.pmol = m;
                    if (s.paa->conditionally_basic()) return;
                }
                throw -4;
            }
        }

        saved_from = m;
    }

    Point relcen, savcen;
    Atom* acen = nullptr;
    int censaidx = -1;
    if (ran)
    {
        acen = m->get_atom(ran);
        for (i=0; i<sz && m->atoms[i]; i++)
        {
            if (!strcmp(m->atoms[i]->name, ran)) censaidx = i;
            if (censaidx >= 0) break;
        }
    }

    if (acen) relcen = acen->loc;
    else relcen = m->get_barycenter();

    if (censaidx >= 0) savcen = saved_atom_locs[censaidx];
    else savcen = average_of_points(saved_atom_locs, sz);           // Points are weighted so average will be barycenter.

    for (i=0; i<sz && m->atoms[i]; i++)
    {
        Point pt = saved_atom_locs[i];
        pt = pt.subtract(savcen);
        pt = pt.add(relcen);
        m->atoms[i]->move(pt);
    }
}

float Pose::total_atom_motions()
{
    if (!saved_from || !saved_from->atoms || !sz) return 0;
    int i;
    float result = 0;

    for (i=0; i<sz && saved_from->atoms[i]; i++)
    {
        float r = saved_from->atoms[i]->loc.get_3d_distance(saved_atom_locs[i]);
        result += r;
    }

    return result;
}

void Molecule::delete_atom(Atom* a)
{
    if (!a) return;
    paths = nullptr;

    int i, j;

    // Delete from monomers too.
    if (nmonomers && monomers)
    {
        for (i=0; i<nmonomers; i++)
        {
            if (monomers[i]->get_atom(a->name) == a)
            {
                monomers[i]->delete_atom(a);
            }
        }
    }

    if (!atoms[0]->residue)
    {
        i = 0;
    }

    if (hasAtoms(atoms))
    {
        for (i=0; atoms[i]; i++)
        {
            if (atoms[i] == a)
            {
                a->unbond_all();
                for (j=i+1; atoms[j]; j++) atoms[j-1] = atoms[j];
                atoms[j-1] = nullptr;
                rotatable_bonds = nullptr;
                for (atcount=0; atoms[atcount]; atcount++);     // Get count.
                return;
            }
        }
    }

    cout << "Attempt to delete atom " << a->name << " from a molecule it is not part of." << endl << flush;
    throw 0xbada70b;
}

void Molecule::delete_all_atoms()
{
    if (!atoms) return;
    if (paths) delete[] paths;
    paths = nullptr;

    int i;
    for (i=0; atoms[i]; i++)
    {
        delete atoms[i];
    }
    
    delete[] atoms;
    atoms = nullptr;
    atcount = 0;
}

void Molecule::reset_conformer_momenta()
{
    srand(time(nullptr));

    lmx = _def_lin_momentum * sgn(0.5-(rand()&1));
    lmy = _def_lin_momentum * sgn(0.5-(rand()&1));
    lmz = _def_lin_momentum * sgn(0.5-(rand()&1));
    amx = _def_ang_momentum * conformer_momenta_multiplier * conformer_tumble_multiplier * sgn(0.5-(rand()&1));
    amy = _def_ang_momentum * conformer_momenta_multiplier * conformer_tumble_multiplier * sgn(0.5-(rand()&1));
    amz = _def_ang_momentum * conformer_momenta_multiplier * conformer_tumble_multiplier * sgn(0.5-(rand()&1));

    Bond** b = get_rotatable_bonds();
    int i;

    if (b)
    {
        for (i=0; b[i]; i++)
        {
            b[i]->angular_momentum = _def_bnd_momentum * conformer_momenta_multiplier * sgn(0.5-(rand()&1));
        }
    }
}

void Molecule::reallocate()
{
    if (!(atcount % _def_atc))
    {
        int oac = atcount;
        int ac1 = oac + _def_atc;
        Atom** latoms = new Atom*[ac1+10];

        // if (atcount) cout << "Molecule " << (name?name:"(no name)") << " has " << atcount << " atoms." << endl;

        int i;
        if (atoms && oac)
        {
            for (i=0; i<oac; i++) latoms[i] = atoms[i];
        }
        for (i=oac; i<ac1; i++) latoms[i] = nullptr;
        // delete[] atoms;
        atoms = latoms;
    }
    rotatable_bonds = nullptr;
}

void Molecule::atoms_from_multimers()
{
    if (!nmonomers) return;
    if (atoms) delete[] atoms;
    atoms = nullptr;
    atcount = 0;
    int i, j, l, n=1;
    for (i=0; i<nmonomers; i++) for (j=0; monomers[i]->atoms[j]; j++)
    {
        add_existing_atom(monomers[i]->atoms[j]);
        char* aname = monomers[i]->atoms[j]->name;
        l = strlen(aname);
        if (aname[l-1] >= '0' && aname[l-1] <= '9')
        {
            std::string str = aname;
            std::regex patt("[0-9]");
            str = std::regex_replace(str, patt, "");
            str = str + std::to_string(n++);
            strcpy(monomers[i]->atoms[j]->name, str.c_str());
        }
    }
}

void Molecule::add_existing_atom(Atom* a)
{
    if (!atoms) atcount = 0;
    else for (atcount=0; atoms[atcount]; atcount++);     // Get count.

    reallocate();
    atoms[atcount++] = a;
    atoms[atcount] = nullptr;

    a->mol = reinterpret_cast<void*>(this);
    preflex_cb = &g_total_mclash;
    postflex_cb = &mclash_delta;
    if (atcount > 1)
    {
        strcpy(a->aa3let, atoms[0]->aa3let);
        a->residue = atoms[0]->residue;
        a->aaletter = atoms[0]->aaletter;
    }

    clear_all_bond_caches();
}

Atom* Molecule::add_atom(char const* elemsym, char const* aname, const Point* location, Atom* bond_to, const float bcard, const int charge)
{
    paths = nullptr;

    Atom* a = new Atom(elemsym, location, charge);
    a->name = new char[strlen(aname)+1];
    a->residue = 0;
    strcpy(a->name, aname);
    if (bond_to && bcard) a->bond_to(bond_to, bcard);
    a->mol = reinterpret_cast<void*>(this);
    preflex_cb = &g_total_mclash;
    postflex_cb = &mclash_delta;

    reallocate();
    atoms[atcount++] = a;
    atoms[atcount] = nullptr;

    return a;
}

Atom* Molecule::add_atom(char const* elemsym, char const* aname, Atom* bondto, const float bcard)
{
    if (!bondto || !bcard)
    {
        Point pt;
        return add_atom(elemsym, aname, &pt, bondto, bcard);
    }

    paths = nullptr;
    for (atcount=0; atoms[atcount]; atcount++);

    if (bondto->is_pi()) bondto->aromatize();
    Vector v = bondto->get_next_free_geometry(bcard);
    Atom* a = new Atom(elemsym);
    a->name = new char[strlen(aname)+1];
    a->residue = 0;
    strcpy(a->name, aname);
    a->mol = reinterpret_cast<void*>(this);
    preflex_cb = &g_total_mclash;
    postflex_cb = &mclash_delta;

    try
    {
        v.r = InteratomicForce::covalent_bond_radius(a, bondto, bcard);
        Point pt(&v);
        Point loc = bondto->loc;
        pt = pt.add(&loc);
        a->move(&pt);

        reallocate();
        atoms[atcount++] = a;
        atoms[atcount] = nullptr;
        a->bond_to(bondto, bcard);

        clear_all_bond_caches();

        if ((atcount & 1) && bondto->get_bonded_atoms_count() == 2)
        {
            Bond* b = bondto->get_bond_between(a);
            if (b && b->can_rotate)
            {
                b->rotate(M_PI);
            }
        }

        return a;
    }
    catch (int err)
    {
        if (err == BOND_DEF_NOT_FOUND)
            cout << "Covalent bond parameters not found for " << elemsym
                 << " to " << bondto->get_elem_sym()
                 << " cardinality " << bcard
                 << ". Please check bindings.dat." << endl;
        throw 0xbad6160;
    }
}

void Molecule::clear_all_bond_caches()
{
    if (!atoms) return;
    int i, j;
    for (i=0; atoms[i]; i++)
    {
        Bond* b[16];
        atoms[i]->fetch_bonds(b);
        if (b[0])
        {
            for (j=0; b[j]; j++) b[j]->clear_moves_with_cache();
        }
        atoms[i]->used = 0;
    }
}

int Molecule::is_residue()
{
    if (noAtoms(atoms)) return 0;
    int i;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->Z < 9 && atoms[i]->get_family() >= TETREL && atoms[i]->residue) return atoms[i]->residue;     // Rule out metals.
    }
    return 0;
}

int Molecule::get_hydrogen_count()
{
    if (noAtoms(atoms)) return 0;
    int i, retval=0;

    for (i=0; atoms[i]; i++)
        if (atoms[i]->Z == 1)
            retval++;
    
    return retval;
}

int Molecule::count_atoms_by_element(const char* esym)
{
    if (noAtoms(atoms)) return 0;
    int findZ = Atom::Z_from_esym(esym);
    int i, retval=0;

    for (i=0; atoms[i]; i++)
        if (atoms[i]->Z == findZ)
            retval++;
    
    return retval;
}

void Molecule::hydrogenate(bool steric_only)
{
    if (noAtoms(atoms)) return;
    int i, j;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_metal()) continue;					// No hydrides.
        if (atoms[i]->dnh) continue;
        if (atoms[i]->is_pi() && atoms[i]->get_bonded_atoms_count() > 2) continue;

        float valence = atoms[i]->get_valence();
        if (valence > 4) valence = 8 - valence;

        #if _dbg_hydrogenate
        cout << atoms[i]->name << " has valence " << valence << endl;
        #endif

        float bcardsum = 0;

        Bond* aib[16];
        atoms[i]->fetch_bonds(aib);
        int db = 0;
        if (aib[0])
        {
            for (j=0; aib[j]; j++)
            {
                if (aib[j]->atom2) bcardsum += aib[j]->cardinality;
                if (aib[j]->cardinality > 1) db++;
            }
        }
        if (!db && steric_only) continue;
        if (atoms[i]->aaletter == 'W' && !strcmp(atoms[i]->name, "NE1")) bcardsum--;

        if (bcardsum && atoms[i]->Z == 1) continue;
        bcardsum = ceil(bcardsum);

        #if _dbg_hydrogenate
        cout << " minus existing bonds " << bcardsum;
        #endif

        bcardsum -= atoms[i]->get_charge();
        if (atoms[i]->is_backbone && !strcmp(atoms[i]->name, "N")) bcardsum = 1;

        #if _dbg_hydrogenate
        cout << " given charge makes " << bcardsum << endl;
        #endif

        int fam = atoms[i]->get_family();
        if (fam == PNICTOGEN || fam == CHALCOGEN)
        {
            Atom* C = atoms[i]->is_bonded_to_pi(TETREL, true);
            if (!C) C = atoms[i]->is_bonded_to_pi(PNICTOGEN, true);
            if (C) atoms[i]->aromatize();
        }

        int h_to_add = atoms[i]->is_conjugated_to_charge() > 0 ? ceil(valence - bcardsum) : floor(valence - bcardsum);

        // These two lines are bad kludges and the code ought to be improved to not require them.
        if (atoms[i]->aaletter == 'H' && !strcmp(atoms[i]->name, "NE2")) h_to_add++;
        if (atoms[i]->aaletter == 'R' && !strcmp(atoms[i]->name, "NE")) h_to_add = 1;

        for (j=0; j<h_to_add; j++)
        {
            if (atoms[i]->is_pi() && atoms[i]->get_bonded_atoms_count() > 2) break;

            char hname[15];
            sprintf(hname, "H%d", atcount+1);
            Atom* H = add_atom("H", hname, atoms[i], 1);
            #if _dbg_hydrogenate
            cout << "Adding " << hname << " to " << atoms[i]->name << " whose valence is " << valence << " and has " << bcardsum << " bonds already." << endl;
            #endif

            if (atoms[i]->get_geometry() == 3)
            {
                Bond* aib = atoms[i]->get_bond_by_idx(1);
                int k=0;
                Bond* b0 = atoms[i]->get_bond_by_idx(k++);
                if (!b0->atom2 || b0->atom2 == H) b0 = atoms[i]->get_bond_by_idx(k++);
                Bond* b1 = atoms[i]->get_bond_by_idx(k++);
                if (!b1->atom2 || b1->atom2 == H) b1 = atoms[i]->get_bond_by_idx(k++);

                if (b0->atom2 && b1->atom2 && (abs((__int64_t)(b1) - (__int64_t)b1->atom2) < memsanity))
                {
                    Point source = atoms[i]->loc;
                    Point movable = b1->atom2->loc;
                    Vector axis = compute_normal(source, b0->atom2->loc, b1->atom2->loc);
                    Point plus  = rotate3D(movable, source, axis,  triangular);
                    Point minus = rotate3D(movable, source, axis, -triangular);

                    float rp = plus.get_3d_distance(b0->atom2->loc);
                    float rm = minus.get_3d_distance(b0->atom2->loc);

                    Point pt = ((rp > rm) ? plus : minus).subtract(source);
                    pt.scale(InteratomicForce::covalent_bond_radius(atoms[i], H, 1));
                    H->move(pt.add(source));
                }
            }
        }
    }

    for (atcount=0; atoms[atcount]; atcount++);     // Update count.

    clear_all_bond_caches();
}

bool Molecule::protonate()
{
    if (!atoms) return false;
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->is_amide()) continue;
        int fam = atoms[i]->get_family();
        if (fam != PNICTOGEN && fam != CHALCOGEN) continue;
        int g = atoms[i]->get_geometry();
        int nb = atoms[i]->get_bonded_atoms_count();
        if (nb >= g) continue;

        char aname[6];
        if (atoms[i]->get_Greek())
        {
            strcpy(aname, atoms[i]->name);
            if (aname[0] < 'A') aname[1] = 'H';
            else aname[0] = 'H';
        }
        else
        {
            sprintf(aname, "H%i", get_atom_count());
        }

        Atom* H = add_atom("H", aname, atoms[i], 1);
        H->residue = atoms[i]->residue;
        strcpy(H->aa3let, atoms[i]->aa3let);
        H->aaletter = atoms[i]->aaletter;
        H->conjugation = atoms[i]->conjugation;
        H->pdbchain = atoms[i]->pdbchain;
        H->region = atoms[i]->region;
        atoms[i]->increment_charge(1);
        return true;
    }
    return false;
}

bool Molecule::deprotonate()
{
    if (!atoms) return false;
    int i;
    for (i=0; atoms[i]; i++);
    i--;
    for (; i>=0; i--)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->is_amide()) continue;
        int Z = atoms[i]->Z;
        if (Z > 1) continue;

        Bond* b = atoms[i]->get_bond_by_idx(0);
        if (!b) continue;
        Atom* heavy = b->atom2;
        if (!heavy) continue;

        int fam = heavy->get_family();
        if (fam != PNICTOGEN && fam != CHALCOGEN) continue;

        delete_atom(atoms[i]);
        heavy->increment_charge(-1);
        return true;
    }
    return false;
}

void Molecule::propagate_stays()
{
    if (!nmonomers) return;
    int i, j;
    for (i=0; i<nmonomers; i++)
    {
        for (j=0; monomers[i]->atoms[j]; j++)
        {
            if (monomers[i]->atoms[j] == stay_close_mine)
            {
                monomers[i]->stay_close_mol = stay_close_mol;
                monomers[i]->stay_close_mine = stay_close_mine;
                monomers[i]->stay_close_other = stay_close_other;
                monomers[i]->stay_close_optimal = stay_close_optimal;
                monomers[i]->stay_close_tolerance = stay_close_tolerance;
                monomers[i]->stay_close2_mol = stay_close2_mol;
                monomers[i]->stay_close2_mine = stay_close2_mine;
                monomers[i]->stay_close2_other = stay_close2_other;
                monomers[i]->stay_close2_optimal = stay_close2_optimal;
            }
        }
    }
}

void Molecule::dehydrogenate()
{
    if (!atoms) return;

    Atom* tmp[atcount];
    int i, j=0;
    for (i=0; i<atcount; i++)
    {
        if (!atoms[i]) continue;
        if (atoms[i]->Z < 2)
        {
            atoms[i]->get_heavy_atom()->unbond(atoms[i]);
            continue;
        }
        tmp[j++] = atoms[i];
    }

    for (i=0; i<j; i++) atoms[i] = tmp[i];
    atoms[j] = nullptr;
}

char** Molecule::get_atom_names() const
{
    int i;
    char** retval = new char*[atcount+1];

    for (i=0; atoms[i]; i++) retval[i] = atoms[i]->name;
    retval[atcount] = 0;

    return retval;
}

Atom* Molecule::get_atom(char const* aname) const
{
    if (noAtoms(atoms)) return 0;

    int i;
    for (i=0; atoms[i]; i++)
        if (!strcmp(atoms[i]->name, aname)) return atoms[i];

    return 0;
}

int Molecule::atom_idx_from_ptr(Atom* a)
{
    if (!atoms) return -1;
    int i;
    for (i=0; atoms[i]; i++) if (atoms[i] == a) return i;

    return -1;		// If not found.
}

Point Molecule::get_atom_location(char const * aname)
{
    if (noAtoms(atoms))
    {
        Point pt;
        return pt;
    }
    Atom* a = get_atom(aname);
    if (!a) return get_barycenter();
    return a->loc;
}

Point Molecule::get_atom_location(int i)
{
    if (noAtoms(atoms))
    {
        Point pt;
        return pt;
    }
    Atom* a = get_atom(i);
    if (!a) return get_barycenter();
    return a->loc;
}

Atom* Molecule::get_nearest_atom(Point loc) const
{
    if (noAtoms(atoms)) return 0;

    int i, j;
    float best, r;
    for (i=0; atoms[i]; i++)
    {
        r = loc.get_3d_distance(atoms[i]->loc);
        if (!i || r < best)
        {
            j=i;
            best=r;
        }
    }

    return atoms[j];
}

Atom* Molecule::get_nearest_atom(Point loc, intera_type capable_of, int sp) const
{
    if (noAtoms(atoms)) return 0;

    int i, j=-1;
    float best=0, r;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;             // TODO: What if I'm a glycine? Or I'm a MAILV but caller wants an hbond?
        if (!InteratomicForce::atom_is_capable_of(atoms[i], capable_of)) continue;
        if (sp != any_element)
        {
            float apol = atoms[i]->is_polar();
            if (fabs(apol) < hydrophilicity_cutoff) apol = 0;
            if (sgn(sp) != sgn(apol)) continue;
        }
        r = loc.get_3d_distance(atoms[i]->loc);
        if (!i || !best || r < best)
        {
            j=i;
            best=r;
        }
    }

    if (j<0) return nullptr;
    return atoms[j];
}

Atom *Molecule::get_nearest_atom_to_line(Point A, Point B) const
{
    if (noAtoms(atoms)) return 0;

    int i, j;
    float best, r;
    for (i=0; atoms[i]; i++)
    {
        r = atoms[i]->loc.get_distance_to_line(A, B);
        if (!i || r < best)
        {
            j=i;
            best=r;
        }
    }

    return atoms[j];
}

std::vector<Atom*> Molecule::longest_dimension()
{
    std::vector<Atom*> retval;
    if (!atoms) return retval;
    int i, j;
    float rmax = 0;

    for (i=0; atoms[i]; i++)
    {
        for (j=i+1; atoms[j]; j++)
        {
            float r = atoms[i]->distance_to(atoms[j]);
            if (r > rmax)
            {
                retval.clear();
                retval.push_back(atoms[i]);
                retval.push_back(atoms[j]);
                rmax = r;
            }
        }
    }

    return retval;
}

float Molecule::total_eclipses()
{
    if (!atoms) return 0;
    float result = 0;
    int i, j, k, l, m, n;
    if (!rotatable_bonds) get_rotatable_bonds();
    if (!rotatable_bonds) return 0;
    while (!eclipse_hash) eclipse_hash = rand();

    for (i=0; rotatable_bonds[i]; i++)
    {
        if (rotatable_bonds[i]->eclipse_hash == eclipse_hash)
        {
            result += rotatable_bonds[i]->eclipse_partial;
            // cout << " " << name << flush;
            continue;
        }
        else
        {
            rotatable_bonds[i]->eclipse_partial = 0;
        }

        Bond *abt[16], *bbt[16];
        rotatable_bonds[i]->atom1->fetch_bonds(abt);
        rotatable_bonds[i]->atom2->fetch_bonds(bbt);
        m = rotatable_bonds[i]->atom1->get_geometry();
        n = rotatable_bonds[i]->atom2->get_geometry();

        Point b2loc;
        float avdw;
        Vector* bgeov;
        Point normal;
        Point antinormal;

        if (n == 3)
        {
            b2loc = rotatable_bonds[i]->atom2->loc;
            avdw = rotatable_bonds[i]->atom2->vdW_radius;
            bgeov = rotatable_bonds[i]->atom2->get_geometry_aligned_to_bonds();
            normal = compute_normal(bgeov[0], bgeov[1], bgeov[2]);
            normal.scale(avdw*1.3);
            antinormal = Point(0,0,0).subtract(normal);
        }

        for (j=0; j<m; j++)
        {
            if (!abt[j]) continue;
            if (!abt[j]->atom1) continue;
            if (!abt[j]->atom2) continue;
            if (abt[j]->atom2 == rotatable_bonds[i]->atom2) continue;
            Point abtjloc = abt[j]->atom2->loc;

            if (n != 3) for (l=0; l<n; l++)
            {
                if (!bbt[l]) continue;
                if (!bbt[l]->atom1) continue;
                if (!bbt[l]->atom2) continue;
                if (bbt[l]->atom2 == rotatable_bonds[i]->atom1) continue;

                // This doesn't fit exactly the values from https://www.sas.upenn.edu/~kimg/mcephome/chem502/ethbutconform/ethbutmm2.html
                // but for carbon-carbon energies, it's close enough.
                float sigma = abt[j]->atom2->vdW_radius + bbt[l]->atom2->vdW_radius;
                float r = abt[j]->atom2->distance_to(bbt[l]->atom2);
                float sigma_r = sigma / r;

                float f = Lennard_Jones_epsilon_x4 * pow(sigma_r, 12) / 2;
                result += f;
                rotatable_bonds[i]->eclipse_partial += f;
            }

            if (n == 3 && rotatable_bonds[i]->atom2->is_pi())
            {
                // Pi orbitals.
                float sigma = abt[j]->atom2->vdW_radius + avdw*1.5;
                float r = abtjloc.get_3d_distance(b2loc.add(normal));
                float sigma_r = sigma / r;

                float f = Lennard_Jones_epsilon_x4 * pow(sigma_r, 12) / 2;
                result += f;
                rotatable_bonds[i]->eclipse_partial += f;

                r = abtjloc.get_3d_distance(b2loc.add(antinormal));
                sigma_r = sigma / r;

                f = Lennard_Jones_epsilon_x4 * pow(sigma_r, 12) / 2;
                result += f;
                rotatable_bonds[i]->eclipse_partial += f;
            }
        }

        rotatable_bonds[i]->eclipse_hash = eclipse_hash;
    }

    if (base_eclipses > result) base_eclipses = result;
    return result - base_eclipses;
}

int Molecule::atoms_inside_sphere(Sphere s, bool* bi, float rm)
{
    if (!atoms) return 0;
    int numfound = 0;
    int i;
    for (i=0; atoms[i]; i++)
    {
        float r = s.center.get_3d_distance(atoms[i]->loc);
        if (r <= s.radius*rm + atoms[i]->vdW_radius)
        {
            numfound++;
            if (bi) bi[i] = true;
        }
    }

    return numfound;
}

float Molecule::surface_occlusion(Molecule **ligands)
{
    int i, j, l;
    float total_occlusions = 0;
    float surf = get_surface_area();

    for (i=0; ligands[i]; i++)
    {
        if (ligands[i] == this) continue;
        Atom* a;
        for (j=0; a = ligands[i]->atoms[j]; j++)
        {
            if (!j && !a->molsurf_area) ligands[i]->get_surface_area();
            for (l=0; atoms[l]; l++)
            {
                float r = a->distance_to(atoms[l]);
                if (r > water_molecule_size) continue;

                r /= (a->vdW_radius / atoms[l]->vdW_radius)/2;
                float partial = a->molsurf_area * atoms[l]->molsurf_area / (r*r) / 2;

                total_occlusions += partial;
            }
        }
    }

    return total_occlusions / surf;
}

float Molecule::surface_occlusion(Molecule *ligand)
{
    Molecule* tmp[2];
    tmp[0] = ligand;
    tmp[1] = nullptr;
    return surface_occlusion(tmp);
}

float Molecule::octant_occlusion(Molecule *ligand)
{
    Molecule* tmp[2];
    tmp[0] = ligand;
    tmp[1] = nullptr;
    return octant_occlusion(tmp);
}

float Molecule::octant_occlusion(Molecule **ligands)
{
    if (!atoms) return 0;
    int h, i, j, l;

    #if per_atom_occlusions
    float worst_occlusion = 1835;
    float total_occlusion = 0;
    int occl_samples = 0;
    for (h=0; atoms[h]; h++)
    {
        if (atoms[h]->is_backbone) continue;
        Sphere octant_atoms[8];
        for (i=0; i<8; i++)
        {
            octant_atoms[i].center = Point(0,0,0);
            octant_atoms[i].radius = 0;
        }
        int ag = atoms[h]->get_geometry();
        for (j=0; j<ag; j++)
        {
            Bond* b = atoms[h]->get_bond_by_idx(j);
            if (!b) continue;
            if (!b->atom2) continue;

            Vector rel = b->atom2->loc.subtract(atoms[h]->loc);
            int octi = rel.octant_idx();
            float oim = octant_atoms[octi].center.magnitude();

            float r = rel.r;
            if (!oim || r < oim)
            {
                octant_atoms[octi].center = rel;
                octant_atoms[octi].radius = atoms[h]->vdW_radius + b->atom2->vdW_radius;
            }
        }

        for (j=0; ligands[j]; j++)
        {
            Atom* b = ligands[j]->get_nearest_atom(atoms[h]->loc);
            if (!b) continue;

            Vector rel = b->loc.subtract(atoms[h]->loc);
            Atom* a = get_nearest_atom(b->loc);
            if (atoms[h]->get_heavy_atom()->is_bonded_to(a->get_heavy_atom()) 
                || atoms[h]->get_heavy_atom()->shares_bonded_with(a->get_heavy_atom())) 
                rel = b->loc.subtract(a->loc);
            int octi = rel.octant_idx();
            float oim = octant_atoms[octi].center.magnitude();

            float r = rel.r;
            if (!oim || r < oim)
            {
                octant_atoms[octi].center = rel;
                octant_atoms[octi].radius = atoms[h]->vdW_radius + b->vdW_radius;
            }
        }

        #if _dbg_octant_occlusion
        cout << "Atom " << atoms[h]->name;
        #endif
        float local_occlusions = 0;
        Point zero(0,0,0);
        float local_worst = 1;
        for (i=0; i<4; i++)
        {
            int j = 7-i;
            #if _dbg_octant_occlusion
            if (!i) cout << endl;
            cout << "Octant " << i << " vs. " << j << ":";
            #endif
            if (octant_atoms[i].radius && octant_atoms[j].radius)
            {
                float theta = find_3d_angle(octant_atoms[i].center, octant_atoms[j].center, zero);
                #if _dbg_octant_occlusion
                cout << " theta " << (theta*fiftyseven);
                #endif
                float partial = pow(cos(fmin(theta-M_PI, square)), 0.333);
                #if _dbg_octant_occlusion
                cout << " partial " << partial;
                #endif
                float r = fmax(octant_atoms[i].center.magnitude()-octant_atoms[i].radius, 0) / octant_atoms[i].radius + 1;
                partial /= r;
                #if _dbg_octant_occlusion
                cout << " / " << r;
                #endif
                r = fmax(octant_atoms[j].center.magnitude()-octant_atoms[j].radius, 0) / octant_atoms[i].radius + 1;
                partial /= r;
                #if _dbg_octant_occlusion
                cout << " / " << r << " = " << partial << endl;
                #endif
                local_occlusions += partial/3;          // it only takes 3 out of 4.
                if (partial < local_worst) local_worst = partial;
            }
            #if _dbg_octant_occlusion
            else cout << " (vacant)" << endl;
            #endif
        }
        local_occlusions = fmin(1, local_occlusions-local_worst/3);

        #if _dbg_octant_occlusion
        cout << "Local: " << local_occlusions << endl << endl;
        #endif

        if (local_occlusions < worst_occlusion) worst_occlusion = local_occlusions;
        total_occlusion += local_occlusions;
        occl_samples++;
    }

    // return worst_occlusion;
    return occl_samples ? (total_occlusion / occl_samples) : 0;
    #else
    float total_occlusions = 0;
    Sphere octant_atoms[8];

    for (i=0; i<8; i++) octant_atoms[i].radius = 0;

    #if _dbg_octant_occlusion
    cout << endl;
    #endif
    for (i=0; ligands[i]; i++)
    {
        if (ligands[i] == this) continue;
        Atom* a;
        for (j=0; a = ligands[i]->atoms[j]; j++)
        {
            if (!j && !a->molsurf_area) ligands[i]->get_surface_area();
            Atom* b = get_nearest_atom(a->loc);
            Vector rel = a->loc.subtract(b->loc);
            int octi = rel.octant_idx();
            float oim = octant_atoms[octi].center.magnitude();

            float r = rel.r;
            if (!oim || r < oim)
            {
                octant_atoms[octi].center = rel;
                octant_atoms[octi].radius = a->vdW_radius + b->vdW_radius;
                #if _dbg_octant_occlusion
                cout << "Atom " << ligands[i]->name << ":" << a->name << " for octant " << octi << endl;
                #endif
            }
        }
    }

    Point zero(0,0,0);
    for (i=0; i<4; i++)
    {
        int j = 7-i;
        #if _dbg_octant_occlusion
        if (!i) cout << endl;
        cout << "Octant " << i << " vs. " << j << ":";
        #endif
        if (octant_atoms[i].radius && octant_atoms[j].radius)
        {
            float theta = find_3d_angle(octant_atoms[i].center, octant_atoms[j].center, zero);
            #if _dbg_octant_occlusion
            cout << " theta " << (theta*fiftyseven);
            #endif
            float partial = pow(cos(fmin(theta-M_PI, square)), 0.333);
            #if _dbg_octant_occlusion
            cout << " partial " << partial;
            #endif
            float r = fmax(octant_atoms[i].center.magnitude()-octant_atoms[i].radius, 0) + 1;
            partial /= r;
            #if _dbg_octant_occlusion
            cout << " / " << r;
            #endif
            r = fmax(octant_atoms[j].center.magnitude()-octant_atoms[j].radius, 0) + 1;
            partial /= r;
            #if _dbg_octant_occlusion
            cout << " / " << r << " = " << partial << endl;
            #endif
            total_occlusions += partial/4;
        }
    }
    #if _dbg_octant_occlusion
    cout << endl << endl;
    #endif

    return total_occlusions;
    #endif
}

float Molecule::bindability_by_type(intera_type t, bool ib)
{
    if (!atoms) return 0;

    float result = 0;
    int i, heavy = 0;
    float c;
    for (i=0; atoms[i]; i++)
    {
        if (!ib && atoms[i]->is_backbone) continue;
        if (atoms[i]->Z > 1) heavy++;
        switch (t)
        {
            case covalent:
            // TODO
            ;
            break;

            case ionic:
            c = atoms[i]->get_charge();
            if (c) result += c;
            break;

            case hbond:
             if (atoms[i]->Z > 1)
             {
                 result += fabs(atoms[i]->is_polar());
                 result += fabs(atoms[i]->get_charge());
             }
            break;

            case pi:
            result += fabs(atoms[i]->is_pi());
            break;

            case polarpi:
            result += fabs(atoms[i]->is_polar()) + fabs(atoms[i]->is_pi());
            break;

            case mcoord:
            if (atoms[i]->is_metal())
            {
                result += atoms[i]->get_charge();
            }
            else
            {
                int fam = atoms[i]->get_family();
                int Z = atoms[i]->Z;

                if (fam == PNICTOGEN && Z <= 15)
                    result -= (atoms[i]->is_pi() ? 0.5 : 1);

                if (fam == CHALCOGEN && Z <= 35)
                    result -= (atoms[i]->is_pi() ? 0.25 : 1);
            }

            break;

            case vdW:
            default:
            result += 1;
            break;
        }
    }

    if (heavy) result /= heavy;
    if (t == hbond)
    {
        if (result < 0.3) result = 0;
    }

    return result;
}

int Molecule::has_hbond_donors()
{
    if (!atoms) return 0;

    int i, result=0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->get_heavy_atom()->get_family() == TETREL) continue;
        if (atoms[i]->is_polar() >= hydrophilicity_cutoff)
        {
            result++;
            // if (get_charge() < 0.5) cout << name << ":" << atoms[i]->name << " is an hbond donor." << endl;
        }
        if (atoms[i]->get_family() == HALOGEN && atoms[i]->is_bonded_to(TETREL)) result++;
    }

    return result;
}

int Molecule::has_hbond_acceptors()
{
    if (!atoms) return 0;

    int i, result=0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->get_family() == TETREL) continue;
        if (atoms[i]->get_family() == HALOGEN && atoms[i]->is_bonded_to(TETREL)) continue;
        if (atoms[i]->get_bonded_atoms_count() > 3) continue;
        if (atoms[i]->is_polar() <= -hydrophilicity_cutoff) result++;
    }

    return result;
}

int Molecule::has_pi_atoms(bool ib)
{
    if (!atoms) return 0;

    int i, result=0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->Z < 2) continue;
        if (!ib && atoms[i]->is_backbone) continue;
        if (atoms[i]->is_pi())
        {
            result++;
        }
    }

    return result;
}

void Molecule::add_monomer(Molecule *toadd)
{
    if (!nmonomers) make_multimer(1);
    int i;
    Molecule** temp = monomers;
    monomers = new Molecule*[nmonomers+4];
    if (temp) for (i=0; i<nmonomers; i++) monomers[i] = temp[i];
    monomers[i++] = toadd;
    monomers[i] = nullptr;
    nmonomers++;
    if (temp) delete[] temp;

    for (i=0; toadd->atoms[i]; i++) add_existing_atom(toadd->atoms[i]);
}

void Molecule::make_multimer(int n)
{
    int i, j, l;
    Molecule* source = this;
    if (nmonomers) source = monomers[0];

    Molecule** old = monomers;
    monomers = new Molecule*[n+4];
    for (i=0; i<n; i++)
    {
        std::string lname = name ? name : "ligand";
        lname += std::to_string(i+1);
        monomers[i] = new Molecule(lname.c_str());

        for (j=0; source->atoms[j]; j++)
        {
            Bond* b = source->atoms[j]->get_bond_by_idx(0);
            char esym[5];
            strcpy(esym, Atom::esym_from_Z(source->atoms[j]->Z));
            monomers[i]->add_atom(esym, source->atoms[j]->name, &source->atoms[j]->loc,
                nullptr, 0, source->atoms[j]->get_orig_charge());
        }

        for (j=0; source->atoms[j]; j++)
        {
            Atom* a1 = monomers[i]->get_atom(source->atoms[j]->name);
            for (l=0; l<source->atoms[j]->get_valence(); l++)
            {
                Bond* b = source->atoms[j]->get_bond_by_idx(l);
                if (!b) continue;
                if (!b->atom2) continue;
                if (b->cardinality < 1) continue;
                Atom* a2 = monomers[i]->get_atom(b->atom2->name);
                if (!a1->is_bonded_to(a2)) a1->bond_to(a2, b->cardinality);
            }
        }

        srand(i);
        Vector mov(frand(2, _INTERA_R_CUTOFF), frand(0, M_PI*2), frand(0, M_PI*2));
        monomers[i]->movability = MOV_ALL;
        monomers[i]->move(mov);
        monomers[i]->rotate(mov, frand(0, M_PI*2));
        monomers[i]->reset_conformer_momenta();
        monomers[i]->identify_conjugations();
    }
    monomers[i] = nullptr;
    nmonomers = n;
    if (old) delete[] old;
    atoms_from_multimers();
}

Molecule *Molecule::get_monomer(int i)
{
    if (i<0 || i>=nmonomers) return this;
    return monomers[i];
}

int Molecule::from_sdf(char const *sdf_dat)
{
    if (!sdf_dat) return 0;
    char const* lines[8192];
    int i,j=0,lncount;

    immobile = false;

    // cout << sdf_dat << endl;

    lines[j++] = sdf_dat;
    for (i=0; sdf_dat[i]; i++)
    {
        if (sdf_dat[i] == '\n')
        {
            lines[j++] = (sdf_dat+i+1);
        }
        if (j > 8190) break;
    }
    lines[lncount = j] = 0;

    int na, nb;
    int added=0;
    char** words;

    for (j=3; j<lncount; j++)
    {
        char line[1024];
        strncpy(line, lines[j], 1023);
        words = chop_spaced_words(line);

        if (!words || !words[0] || !words[1]) break;
        if (!strcmp(words[1], "END")) break;

        if (!strcmp(words[1], "CHG"))
        {
            for (i=3; words[i] && words[i+1]; i+=2)
            {
                int aidx = atoi(words[i]);
                atoms[aidx-1]->increment_charge(atof(words[i+1]));
            }
        }

        if (j == 3)
        {
            na = atoi(words[0]);
            nb = atoi(words[1]);

            atoms = new Atom*[na+4];
        }
        else if (added < na)
        {
            Point* loc = new Point(atof(words[0]), atof(words[1]), atof(words[2]));
            if (words[3][0] >= 'a' && words[3][0] <= 'z') words[3][0] -= 0x20;
            Atom* a = new Atom(words[3], loc);
            delete loc;
            a->name = new char[16];
            sprintf(a->name, "%s%d", words[3], added+1);
            a->residue = 0;
            atoms[atcount++] = a;
            atoms[atcount] = nullptr;
            added++;
        }
        else
        {
            int a1i = atoi(words[0]);
            int a2i = atoi(words[1]);

            if (!a1i || !a2i) break;
            atoms[a1i-1]->bond_to(atoms[a2i-1], atof(words[2]));
        }

        if (words) delete[] words;
    }
    atoms[atcount] = 0;
    if (words) delete[] words;

    identify_conjugations();
    identify_rings();
    identify_cages();
    identify_acidbase();
    return added;
}

void Molecule::identify_conjugations()
{
    if (!atoms) return;
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_pi() && !atoms[i]->conjugation)
        {
            Conjugation* conj = new Conjugation(atoms[i]);          // Don't have to save this pointer because the atoms will save it.
        }
    }
}

bool Molecule::check_Greek_continuity()
{
    if (!atoms) return true;
    int i, j, k, l, m;
    for (i=0; atoms[i]; i++)
    {
        if (!atoms[i]->check_Greek_continuity()) return false;

        continue;

        if (!atoms[i]->residue) continue;
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->Z < 2) continue;

        int n = atoms[i]->get_geometry();
        if (!n) continue;
        Bond* bb[n+8];
        for (j=0; j<n+8; j++) bb[j] = nullptr;
        atoms[i]->fetch_bonds(bb);

        int ig = greek_from_aname(atoms[i]->name);
        if (ig < 0) continue;
        for (j=i+1; atoms[j]; j++)
        {
            if (atoms[j]->Z < 5) continue;
            int jg = greek_from_aname(atoms[j]->name);
            if (jg > ig)
            {
                bool found = false;
                for (k=0; bb[k]; k++)
                {
                    Atom* mwb[256];
                    for (l=0; l<256; l++) mwb[l] = nullptr;
                    bb[k]->fetch_moves_with_atom2(mwb);
                    for (l=0; mwb[l]; l++)
                    {
                        if (mwb[l] == atoms[j])
                        {
                            found = true;
                            goto _exit_mwbsearch;
                        }
                    }
                }
                if (!found) return false;
            }
            _exit_mwbsearch:
            ;
        }
    }

    return true;
}

int Molecule::from_pdb(FILE* is, bool het_only)
{
    /*
    ATOM     55  SG  CYS     4       6.721  -8.103   4.542  1.00001.00           S
    */
    char buffer[1024];
    int added=0;

    while (!feof(is))
    {
        char* got = fgets(buffer, 1003, is);
        if (!got) break;

        int charge = 0, offset = (buffer[21] != ' ' && buffer[22] == ' ') ? 1 : 0;
        char chgstr[3];
        chgstr[0] = 0;
        if (buffer[78] && buffer[78] > ' ')
        {
            chgstr[0] = (buffer[78] > ' ') ? buffer[78] : 0;
            chgstr[1] = (buffer[79] > ' ') ? buffer[79] : 0;
            chgstr[2] = 0;

            if 		(!strcmp(chgstr, "+")) charge = 1;
            else if	(!strcmp(chgstr, "++")) charge = 2;
            else if	(!strcmp(chgstr, "-")) charge = -1;
            else if	(!strcmp(chgstr, "--")) charge = -2;
            else if (atoi(chgstr)) charge = atoi(chgstr);
        }

        if (buffer[0] == 'H' && buffer[1] == 'E' && buffer[2] == 'T' && buffer[3] == 'A' && buffer[4] == 'T' && buffer[5] == 'M'
            && buffer[6] >= '0' && buffer[6] <= '9')
        {
            char buftmp[1024];
            strcpy(buftmp, buffer+6);
            buffer[6] = ' ';
            strcpy(buffer+7, buftmp);
        }
        char** words = chop_spaced_words(buffer);

        if (words)
        {
            if (
                  (!strcmp(words[0], "ATOM") && !het_only)
                  ||
                  !strcmp(words[0], "HETATM")
               )
            {
                try
                {
                    char esym[7];
                    if (words[2][0] >= '0' && words[2][0] <= '9')
                        strcpy(esym, &words[2][1]);
                    else
                        strcpy(esym, words[2]);

                    int i;
                    for (i=1; i<6; i++)
                    {
                        if (!esym[i+1]) esym[i] = 0;
                        if (esym[i+1] >= '0' && esym[i+1] <= '9') esym[i]=0;
                        if (i>1) esym[i] = 0;
                        if (!esym[i]) break;
                    }
                    esym[1] &= 0x5f;

                    // cout << buffer[21] << buffer[22] << " " << offset << endl;
                    Point aloc(atof(words[5+offset]), atof(words[6+offset]),atof(words[7+offset]));

                    Atom* a = add_atom(esym, words[2], &aloc, 0, 0, charge);
                    added++;

                    // a->residue = atoi(words[4]);

                    for (i=0; atoms[i]; i++)
                    {
                        if (atoms[i] == a) continue;
                        float r = aloc.get_3d_distance(atoms[i]->loc);
                        if (r < 5)
                        {
                            if (r < 1.05* InteratomicForce::covalent_bond_radius(a, atoms[i], 1))
                                a->bond_to(atoms[i], 1);
                        }
                    }

                }
                catch (int ex)
                {
                    ;
                }
            }
        }
        buffer[0] = 0;

        delete words;
    }

    identify_conjugations();
    return added;
}

int Molecule::get_bond_count(bool unidirectional) const
{
    int i, j, bc=0;

    if (noAtoms(atoms)) return 0;
    for (i=0; i<atcount && atoms[i]; i++)
    {
        Bond* b[16];
        atoms[i]->fetch_bonds(b);
        if (!b[0]) continue;

        for (j=0; b[j]; j++)
        {
            if (b[j]->atom2 > atoms[i] || !unidirectional) bc++;
        }
    }

    return bc;
}

int Molecule::aidx(Atom* a)
{
    if (!a) return -1;
    int i;
    if (noAtoms(atoms)) return -1;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i] == a) return i;
    }
    return -1;
}

bool Molecule::save_sdf(FILE* os)
{
    return save_sdf(os, 0);
}

bool Molecule::save_sdf(FILE* os, Molecule** lig)
{
    if (!os) return false;
    fprintf(os, "%s\n", name);

    time_t now = time(0);
    tm *gmtm = gmtime(&now);

    // If we used obabel or another third party app, give due credit.
    if (sdfgen_aboutline.length()) fprintf(os, "%s\n", sdfgen_aboutline.c_str());
    else
    {
        fprintf(os, "  AROMA-%02d%02d%02d%02d%02d%02d3D\n", gmtm->tm_year % 100, gmtm->tm_mon+1, gmtm->tm_mday,
                gmtm->tm_hour, gmtm->tm_min, gmtm->tm_sec
               );

        fprintf(os, "https://github.com/primaryodors/aroma\n");
    }

    int ac, bc, chargeds=0;
    ac = get_atom_count();
    bc = get_bond_count(true);

    int i, j, k, l;
    Atom* latoms[65536];
    Bond* lbonds[65536];

    if (hasAtoms(atoms))
        for (j=0; atoms[j]; j++)
            latoms[j] = atoms[j];

    Bond** b = get_all_bonds(true);
    if (b)
        for (l=0; b[l]; l++)
            lbonds[l] = b[l];
    if (b) delete[] b;

    if (lig)
    {
        for (i=0; lig[i] && lig[i]->atoms; i++)
        {
            if (lig[i] == this) continue;
            ac += lig[i]->get_atom_count();
            bc += lig[i]->get_bond_count(true);

            for (k=0; lig[i]->atoms[k]; k++)
                latoms[j++] = lig[i]->atoms[k];

            b = lig[i]->get_all_bonds(true);

            for (k=0; b[k]; k++)
                lbonds[l++] = b[k];

            if (b) delete[] b;
        }
    }

    fprintf(os, " %d %d  0     0  0  0  0  0  0999 V2000\n", ac, bc );

    for (i=0; i<ac; i++)
    {
        Point p = latoms[i]->loc;
        char const* esym = latoms[i]->get_elem_sym();
        if (!esym) continue;

        if (latoms[i]->get_charge()) chargeds++;

        // This avoids some weird "-0.0000" nonsense that messes up the alignment and corrupts the SDF.
        if (!p.x) p.x=0;
        if (!p.y) p.y=0;
        if (!p.z) p.z=0;

        /*fprintf(os, "   %s%5.4f   %s%5.4f   %s%5.4f %s%s  0  0  0  0  0  0  0  0  0  0  0  0\n",
        			(p.x<0)?"":" ",p.x, (p.y<0)?"":" ",p.y, (p.z<0)?"":" ",p.z, esym, esym[1]?"":" "
        	   );*/
        char buffer[256];
        string str;
        sprintf(buffer, "%5.4f", p.x);
        str = buffer;
        str = str_pad(str, 10, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        sprintf(buffer, "%5.4f", p.y);
        str = buffer;
        str = str_pad(str, 10, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        sprintf(buffer, "%5.4f", p.z);
        str = buffer;
        str = str_pad(str, 10, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        fprintf(os, " %s%s  0  0  0  0  0  0  0  0  0  0  0  0\n", esym, esym[1]?"":" ");
    }


    for (i=0; i<bc; i++)
    {
        int laidx=0, lbidx=0;

        for (j=0; j<ac; j++)
        {
            if (latoms[j] == lbonds[i]->atom1) laidx = j+1;
            if (latoms[j] == lbonds[i]->atom2) lbidx = j+1;
        }

        // fprintf(os, " %d %d  %d  0  0  0  0\n", laidx, lbidx, (int)lbonds[i]->cardinality);
        char buffer[256];
        string str;
        sprintf(buffer, "%d", laidx);
        str = buffer;
        str = str_pad(str, 3, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        sprintf(buffer, "%d", lbidx);
        str = buffer;
        str = str_pad(str, 3, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        sprintf(buffer, " %3.2f", lbonds[i]->cardinality);
        str = buffer;
        str = str_pad(str, 3, " ", STR_PAD_LEFT);
        fprintf(os, "%s", str.c_str());
        fprintf(os, "  0  0  0  0\n");
    }

    if (chargeds)
    {
        int thisline = min(chargeds, 8);
        fprintf(os, "M  CHG  %d  ", thisline);		// TODO: Multiline if chargeds>8.
        k = 0;
        for (i=0; i<ac; i++)
        {
            float chg = latoms[i]->get_charge();
            if (!chg) continue;
            int ichg = (chg<0) ? floor(chg) : ceil(chg);
            fprintf(os, "%d  %d  ", i+1, ichg);
            k++;
            if (k > 7)
            {
                fprintf(os, "\n");
                chargeds -= k;
                if (chargeds <= 0) break;
                thisline = min(chargeds, 8);
                fprintf(os, "M  CHG  %d  ", thisline);
            }
        }
        fprintf(os, "\n");
    }

    fprintf(os, "M  END\n\n");

    return true;
}

void Molecule::save_pdb(FILE* os, int atomno_offset, bool endpdb)
{
    int i;

    for (i=0; atoms[i]; i++)
    {
        atoms[i]->save_pdb_line(os, i+1+atomno_offset);
    }

    if (endpdb) fprintf(os, "\nTER\nEND\n\n");
    else fprintf(os, "\n\n");
}

int Molecule::add_ring(Atom** atoms)
{
    int i, j, l, n, ringcount;

    l = 0;
    Atom* min_atom = nullptr;
    for (i=0; atoms[i]; i++)
    {
        if (!min_atom || atoms[i] < min_atom)
        {
            min_atom = atoms[i];
            l = i;
        }
    }
    n = i;

    i = l-1; if (i<0) i += n;
    j = l+1; if (j>=n) j -= n;
    bool reversed = (atoms[j] > atoms[i]);

    Atom* atoms_ordered[n+4];
    for (i=0; i<n; i++)
    {
        if (reversed) j = l - i;
        else j = i + l;
        if (j < 0) j += n;
        if (j >= n) j -= n;
        atoms_ordered[i] = atoms[j];
    }
    atoms_ordered[n] = nullptr;


    int m;
    for (m=0; atoms_ordered[m]; m++);

    if (rings)
    {
        bool already_exists = false;
        for (i=0; rings[i]; i++)
        {
            if (rings[i]->get_overlap_count(atoms_ordered) == m) already_exists = true; 
        }
        ringcount = i;
        if (already_exists)
        {
            return ringcount;
        }
    }
    else
    {
        ringcount=0;
    }

    Ring** ringstmp = new Ring*[ringcount+4];

    if (rings)
    {
        for (i=0; i<ringcount; i++) ringstmp[i] = rings[i];
        delete[] rings;
    }

    Ring* r = new Ring(atoms_ordered);
    ringstmp[ringcount++] = r;
    ringstmp[ringcount] = nullptr;
    rings = ringstmp;

    return ringcount-1;
}

int Molecule::identify_rings()
{
    find_paths();
    // return 0;

    Atom** ringstmp[256];
    int ringcount;
    Atom *a;
    int chainlen[256];
    bool is_ring[256];
    int found_rings=0, chains=0, cnvchain, active, i, j, k, l, m, n, p;
    Atom *cnva, *cnvb;
    Atom *ra, *rb;

    if (!rings) return 0;

    for (ringcount = 0; rings[ringcount]; ringcount++)
    {
        #if _dbg_identify_rings
        cout << "Ring number " << ringcount << *rings[ringcount] << ": coplanar? " << rings[ringcount]->is_coplanar() << ", conjugated? " << rings[ringcount]->is_conjugated() << endl;
        #endif

        if (rings[ringcount]->is_coplanar() && rings[ringcount]->is_conjugated())
        {
            rings[ringcount]->aromatize();
            #if _dbg_identify_rings
            cout << "Aromatized." << endl;
            #endif
        }
    }
    return ringcount;
}

int Molecule::path_contains_atom(int path_idx, Atom* a)
{
    if (!paths) return 0;
    if (!paths[path_idx]) return 0;

    int i;
    for (i=0; paths[path_idx][i]; i++)
        if (paths[path_idx][i] == a) return i+1;

    return 0;
}

int Molecule::path_get_length(int path_idx)
{
    if (!paths) return 0;
    if (!paths[path_idx]) return 0;

    int i;
    for (i=0; paths[path_idx][i]; i++);

    return i;
}

Atom* Molecule::path_get_terminal_atom(int path_idx)
{
    if (!paths) return nullptr;
    if (!paths[path_idx]) return nullptr;

    int n = path_get_length(path_idx);
    if (!n) return nullptr;
    return paths[path_idx][n-1];
}

void Molecule::copy_path(int old_idx, int new_idx)
{
    if (!paths) return;
    if (!paths[old_idx]) return;
    if (abs((__int64_t)paths[old_idx][0] - (__int64_t)paths[old_idx]) >= memsanity) return;

    if (!paths[new_idx]) paths[new_idx] = new Atom*[get_atom_count()];
    int i;
    for (i=0; paths[old_idx][i]; i++)
        paths[new_idx][i] = paths[old_idx][i];
    
    paths[new_idx][i] = nullptr;
}

bool Molecule::path_is_subset_of(int short_path, int long_path)
{
    if (!paths) return false;
    if (!paths[long_path]) return false;
    if (!paths[short_path]) return true;

    int i;
    for (i=0; paths[short_path][i] && paths[long_path][i]; i++)
    {
        if (paths[long_path][i] != paths[short_path][i]) return false;
    }

    return true;
}

void Molecule::echo_path(int i)
{
    int j;
    cout << i << ":";
    for (j=0; paths[i][j]; j++)
        cout << " " << paths[i][j]->name;
    
    cout << endl;
}

void Molecule::find_paths()
{
    if (!atoms || !atoms[0]) return;

    if (paths)
    {
        #if _dbg_path_search
        cout << "Paths already set; skipping." << endl;
        #endif
        return;
    }
    #if _dbg_path_search
    cout << "Searching for paths..." << endl;
    #endif

    int h, i, j, k, l, m, n = 0, p, q, limit;
    Bond* b[16];
    for (i=0; atoms[i]; i++)
    {
        n += atoms[i]->get_bonded_heavy_atoms_count();
    }

    // paths = new Atom**[n];
    limit = n*n;
    Atom** lpaths[limit];
    paths = lpaths;
    for (i=0; i<limit; i++) paths[i] = nullptr;

    atcount = get_atom_count();

    Atom* a = atoms[0];
    a->fetch_bonds(b);
    if (!b[0]) return;
    n=0;
    for (i=0; b[i]; i++)
    {
        if (!b[i]->atom2) continue;
        if (b[i]->atom2->Z < 2) continue;
        if (b[i]->atom2->residue && b[i]->atom2->residue != a->residue) continue;
        paths[n] = new Atom*[atcount];
        for (q=0; q<atcount; q++) paths[n][q] = nullptr;
        paths[n][0] = a;
        paths[n][1] = b[i]->atom2;
        paths[n][2] = nullptr;
        n++;
    }

    int num_added, iter;
    for (iter=0; iter<1000; iter++)
    {
        num_added = 0;
        p = n;

        for (i=0; i<p; i++)
        {
            m = path_get_length(i);
            if (!m) continue;
            a = path_get_terminal_atom(i);
            if (!a) continue;
            if (!a->name) continue;
            if (abs((__int64_t)(atoms[0]) - (__int64_t)a) >= memsanity) continue;
            a->fetch_bonds(b);
            if (!b[0]) continue;

            k=0;
            for (j=0; b[j]; j++)
            {
                if (abs((__int64_t)(a) - (__int64_t)b[j]) > memsanity) break;
                if (abs((__int64_t)(b[j]) - (__int64_t)b[j]->atom2) > memsanity) break;
                if (!b[j]->atom2) continue;
                if (b[j]->atom2->Z < 2) continue;
                if (b[j]->atom2->get_bonded_heavy_atoms_count() < 2) continue;
                if (b[j]->atom2->residue && b[j]->atom2->residue != a->residue) continue;

                #if _dbg_path_search
                cout << "Trying " << b[j]->atom2->name << "... ";
                #endif

                l = path_contains_atom(i, b[j]->atom2);
                if (l > 0)
                {
                    if ((m-l) > 1)
                    {
                        Atom* ring_atoms[m];
                        for (h=l; h<m; h++)
                        {
                            ring_atoms[h-l] = paths[i][h];
                        }
                        ring_atoms[h-l] = b[j]->atom2;
                        h++;
                        ring_atoms[h-l] = nullptr;

                        add_ring(ring_atoms);
                        #if _dbg_path_search
                        cout << "Created ring from ";
                        Atom::dump_array(ring_atoms);
                        #endif
                    }
                }
                else
                {
                    if (paths[n]) delete[] paths[n];
                    paths[n] = new Atom*[atcount];
                    for (q=0; q<atcount; q++) paths[n][q] = nullptr;
                    copy_path(i, n);
                    paths[n][m] = b[j]->atom2;
                    paths[n][m+1] = nullptr;

                    #if _dbg_path_search
                    cout << "Created ";
                    echo_path(n);
                    #endif

                    n++;
                    if (n >= limit) goto _exit_paths;
                    k++;
                    num_added++;
                }
            }

            #if _dbg_path_search
            cout << endl;
            #endif
        }

        for (j=n-2; j>=0; j--)
        {
            if (path_is_subset_of(j, n-1))
            {
                n--;
                int plj = path_get_length(j);
                int pln = path_get_length(n);
                if (plj == pln) num_added--;

                #if _dbg_path_search
                cout << plj << "/" << pln << " ";
                cout << "Replacing ";
                echo_path(j);
                cout << "...with ";
                echo_path(n);
                cout << endl;
                #endif

                copy_path(n, j);
                delete[] paths[n];
                paths[n] = nullptr;
            }
        }

        if (!num_added) break;
    }

    _exit_paths:
    #if _dbg_path_search
    cout << "Paths:" << endl;
    for (i=0; i<limit && paths[i]; i++) echo_path(i);
    #else
    ;
    #endif
}

void Molecule::identify_acidbase()
{
    if (noAtoms(atoms)) return;

    // For every atom in the molecule:
    int i, j, k, l;
    Bond* b[16];

    for (i=0; atoms[i]; i++)
    {
        // If it is a carbon, pi-bonded to a chalcogen, not bonded to a pnictogen,
        // Or if it is a heavier tetrel, a pnictogen, a chalcogen, or a halogen,
        // And it is single-bonded to at least one chalcogen,
        // And the single-bonded chalcogen is either bonded to a hydrogen or carries a negative charge:
        // Then all chalcogens bonded to the carbon are acidic.
        int sbOH = 0;
        int fama = atoms[i]->get_family();
        bool carbon = false;

        if (fama == TETREL || fama == PNICTOGEN || fama == CHALCOGEN || fama == HALOGEN)
        {
            carbon = (atoms[i]->Z == 6);
            //cout << "Atom " << atoms[i]->name << " is of family " << fama << endl;
            atoms[i]->fetch_bonds(b);
            if (!b[0]) goto _not_acidic;

            int nb2 = atoms[i]->get_bonded_atoms_count();
            if ((fama == TETREL || fama == PNICTOGEN) && !atoms[i]->is_pi() && nb2 < 4)
            {
                if (nb2 >= 2)
                {
                    l=0;
                    Point planarity_check[4];
                    planarity_check[l++] = atoms[i]->loc;
                    for (j=0; j<nb2; j++)
                    {
                        if (!b[j]) break;
                        if (!b[j]->atom2) break;
                        if (!b[j]->atom2->Z) break;
                        planarity_check[l++] = b[j]->atom2->loc;
                    }

                    bool lplanar = false;
                    if (l > 3)
                    {
                        float coplanarity = are_points_planar(planarity_check[0], planarity_check[1], planarity_check[2], planarity_check[3]);
                        if (coplanarity < coplanar_threshold) lplanar = true;
                    }
                    else if (l == 3)
                    {
                        float theta = find_3d_angle(planarity_check[1], planarity_check[2], planarity_check[0]);
                        if (theta > square && theta < M_PI && fabs(theta-triangular) < fabs(theta - tetrahedral)) lplanar = true;
                        #if _dbg_internal_energy
                        if (!atoms[i]->residue) cout << atoms[i]->name << " " << (theta*fiftyseven) << (lplanar ? " pi" : "") << endl;
                        #endif
                    }

                    if (lplanar)
                    {
                        atoms[i]->aromatize();
                        int chalcogens = 0;
                        for (j=0; j<l; j++)
                        {
                            if (!b[j]) break;
                            if (!b[j]->atom2) break;
                            if (!b[j]->atom2->Z) break;
                            int bfam = b[j]->atom2->get_family();
                            if (bfam == PNICTOGEN || bfam == CHALCOGEN)
                            {
                                b[j]->atom2->aromatize();
                                b[j]->cardinality = 1.5;

                                if (bfam == CHALCOGEN && b[j]->atom2->get_bonded_heavy_atoms_count() < 2)
                                {
                                    chalcogens++;
                                    if (chalcogens > 1 && !b[j]->atom2->get_charge())
                                    {
                                        b[j]->atom2->increment_charge(-1);
                                        chalcogens = -65536;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (carbon)
            {
                if (!atoms[i]->is_pi() || !atoms[i]->is_bonded_to_pi(CHALCOGEN, true))
                {
                    goto _not_acidic;
                }
                for (j=0; b[j]; j++)
                {
                    if (!b[j]->atom2) continue;
                    if (b[j]->cardinality == 2)
                    {
                        int fam = b[j]->atom2->get_family();
                        if (fam != CHALCOGEN)
                        {
                            goto _not_acidic;
                        }
                    }
                }
            }
            for (j=0; b[j]; j++)
            {
                if (!b[j]->atom2) continue;
                int fam = b[j]->atom2->get_family();
                if (carbon && fam == PNICTOGEN)
                {
                    goto _not_acidic;
                }
                //cout << "Fam: " << fam << endl;
                if (fam == CHALCOGEN && b[j]->cardinality < 2)
                {
                    if (b[j]->atom2->get_charge() < 0)
                    {
                        sbOH++;
                        break;
                    }
                    else
                    {
                        Bond* b1[16];
                        b[j]->atom2->fetch_bonds(b1);
                        if (!b1[0])
                        {
                            goto _not_acidic;
                        }
                        for (k=0; b1[k]; k++)
                        {
                            if (!b1[k]->atom2) continue;
                            if (b1[k]->atom2->Z == 1)
                            {
                                sbOH++;
                                //cout << "OH" << endl;
                                break;
                            }
                        }
                    }
                }
            }
        }
        if (sbOH)
        {
            atoms[i]->fetch_bonds(b);
            for (j=0; b[j]; j++)
            {
                if (!b[j]->atom2) continue;
                int fam = b[j]->atom2->get_family();
                if (fam == CHALCOGEN)
                {
                    b[j]->atom2->set_acidbase(-1);
                    //cout << "Atom " << b[j]->atom2->name << " is acidic." << endl;
                }
            }
        }
    _not_acidic:
        //if (b) delete[] b;

        // If it is a pnictogen,
        // And it is not bonded to a chalcogen or a halogen,
        // And it is not part of an amide,
        // Or if it has a positive charge,
        // TODO: Or if any hydrogen bound to it has a positive charge,
        // Then it is basic.
        if (fama == PNICTOGEN)
        {
            Atom* c;
            c = atoms[i]->is_bonded_to("C");
            if (c)
            {
                Atom* bto = c->is_bonded_to("O");
                if (bto)
                {
                    if (amide_zwitterionic_amount)
                    {
                        // Amides are weakly zwitterionic.
                        float arity = c->is_bonded_to(bto);
                        if (arity >= 1.5)
                        {
                            atoms[i]->increment_charge(amide_zwitterionic_amount);
                            bto->increment_charge(-amide_zwitterionic_amount);
                        }
                    }

                    goto _not_basic;
                }
            }
            if (atoms[i]->get_charge() <= 0)
            {
                if (atoms[i]->is_bonded_to(CHALCOGEN)) goto _not_basic;
                if (atoms[i]->is_bonded_to(HALOGEN)) goto _not_basic;
            }

            atoms[i]->set_acidbase(1);

        }
    _not_basic:
        ;
    }
}

Bond** Molecule::get_rotatable_bonds(bool icf)
{
    if (noAtoms(atoms)) return 0;
    if (rotatable_bonds) return rotatable_bonds;
    if (mol_typ == MOLTYP_AMINOACID)
    {
        // TODO: There has to be a better way.
        Star s;
        s.pmol = this;
        if (!rotatable_bonds) rotatable_bonds = s.paa->get_rotatable_bonds();
        else if (rotatable_bonds[0] && rotatable_bonds[1] && rotatable_bonds[0]->atom1 == rotatable_bonds[1]->atom1) rotatable_bonds = s.paa->get_rotatable_bonds();
        else if (rotatable_bonds[0] && abs(rotatable_bonds[0]->atom1 - rotatable_bonds[0]->atom2) >= 524288) rotatable_bonds = s.paa->get_rotatable_bonds();
        return rotatable_bonds;
    }
    // cout << name << " Molecule::get_rotatable_bonds()" << endl << flush;

    Bond* btemp[65536];

    int i,j, bonds=0;
    if (!immobile)
        for (i=0; atoms[i]; i++)
        {
            Bond* lb[32];
            atoms[i]->fetch_bonds(lb);
            int g = atoms[i]->get_geometry();
            for (j=0; j<g && lb[j]; j++)
            {
                if (!lb[j]->atom1 || !lb[j]->atom2) continue;

                if (lb[j]->cardinality > 1 && (!lb[j]->atom1->is_pi() || !lb[j]->atom2->is_pi()) && lb[j]->cardinality < 3)
                    lb[j]->cardinality = 1;

                bool pia = lb[j]->atom1->is_pi(),
                     pib = lb[j]->atom2->is_pi();

                int fa = lb[j]->atom1->get_family(),
                    fb = lb[j]->atom2->get_family();

                if (lb[j]->atom1->in_same_ring_as(lb[j]->atom2))
                {
                    #if _ALLOW_FLEX_RINGS
                    lb[j]->can_rotate = false;
                    lb[j]->compute_flip_capability();
                    #else
                    lb[j]->can_rotate = lb[j]->can_flip = false;
                    continue;
                    #endif
                }

                // Generally, a single bond from a pi atom to an amino group cannot rotate.
                if (pia && pib)
                {
                    lb[j]->can_rotate = false;
                    lb[j]->compute_flip_capability();
                }

                // If atoms a and b are pi, and a-b cannot rotate, then a-b can flip.
                if (!lb[j]->can_rotate
                    && lb[j]->atom2->is_bonded_to(CHALCOGEN)
                    && fa != CHALCOGEN
                    && !(lb[j]->atom1->in_same_ring_as(lb[j]->atom2))
                    )
                {
                    lb[j]->compute_flip_capability();
                }

                if (lb[j]->atom2
                        &&
                        lb[j]->atom1 < lb[j]->atom2
                        &&
                        (lb[j]->can_rotate || (icf && lb[j]->can_flip))
                   )
                {
                    btemp[bonds++] = lb[j];
                    btemp[bonds] = 0;
                }
            }
        }
    else
        for (i=0; atoms[i]; i++)
        {
            Bond* lb[16];
            atoms[i]->fetch_bonds(lb);
            int g = atoms[i]->get_geometry();
            for (j=0; j<g; j++)
            {
                // Generally, a single bond between pi atoms cannot rotate.
                // Same if pi atom bonded to a pnictogen or chalcogen without a single bond to other atoms.
                if (lb[j]->atom1 && lb[j]->atom2
                    &&  (
                            (
                                lb[j]->atom1->is_pi() &&
                                (   lb[j]->atom2->is_pi()
                                    ||
                                    (   
                                        !lb[j]->atom2->is_bonded_to_pi(TETREL, false)
                                        &&
                                        (
                                            lb[j]->atom2->get_family() == PNICTOGEN
                                            ||
                                            lb[j]->atom2->get_family() == CHALCOGEN
                                        )
                                    )
                                )
                            )
                            ||
                            (
                                lb[j]->atom2->is_pi() &&
                                (   lb[j]->atom1->is_pi()
                                    ||
                                    (   
                                        !lb[j]->atom1->is_bonded_to_pi(TETREL, false)
                                        &&
                                        (
                                            lb[j]->atom1->get_family() == PNICTOGEN
                                            ||
                                            lb[j]->atom1->get_family() == CHALCOGEN
                                        )
                                    )
                                )
                            )
                        )
                   )
                    lb[j]->can_rotate = false;

                if ((lb[j]->can_rotate || (icf && lb[j]->can_flip))
                        &&
                        lb[j]->atom1 && lb[j]->atom2
                        &&
                        (!lb[j]->atom1->is_backbone || !strcmp(lb[j]->atom1->name, "CA"))
                        &&
                        !lb[j]->atom2->is_backbone
                        &&
                        greek_from_aname(lb[j]->atom1->name) == (greek_from_aname(lb[j]->atom2->name)-1)
                        &&
                        lb[j]->atom2->Z > 1
                   )
                {
                    btemp[bonds++] = lb[j];
                    btemp[bonds] = 0;
                }
            }
        }

        rotatable_bonds = new Bond*[bonds+1];
        for (i=0; i<=bonds; i++) rotatable_bonds[i] = btemp[i];
        rotatable_bonds[bonds] = 0;

        return rotatable_bonds;
}

void Molecule::crumple(float theta)
{
    Bond** b = get_rotatable_bonds();
    if (!b) return;

    float int_clsh = get_internal_clashes();

    int i;
    for (i=0; b[i]; i++)
    {
        if (b[i]->can_rotate)
        {
            float ltheta = frand(-theta, theta)*pow(frand(0,1),2);
            b[i]->rotate(ltheta);
            if (get_internal_clashes() > clash_limit_per_aa*4) b[i]->rotate(-ltheta);
        }
        else if (b[i]->can_flip)
        {
            float sint = sin(theta);
            if (frand(0,2) <= sint)
            {
                b[i]->rotate(b[i]->flip_angle);
                if (get_internal_clashes() > clash_limit_per_aa*4) b[i]->rotate(b[i]->flip_angle);
            }
        }
    }
}

void Molecule::clear_cache()
{
    rotatable_bonds = nullptr;
}

// TODO: There has to be a better way.
Bond** AminoAcid::get_rotatable_bonds()
{
    if (rotatable_bonds) return rotatable_bonds;

    // cout << name << " AminoAcid::get_rotatable_bonds()" << endl << flush;
    // Return ONLY side chain bonds, from lower to higher Greek. E.g. CA-CB but NOT CB-CA.
    // Exclude CA-N and CA-C as these will be managed by the Protein class.
    if (noAtoms(atoms)) return 0;

    // TODO: Something is overwriting the cached rotatable_bonds, causing segfaults.
    // So the cache is unusable for amino acids until the problem gets fixed.
    // if (rotatable_bonds) return rotatable_bonds;
    if (aadef && aadef->proline_like)
    {
        // cout << "Proline-like! No rotbonds!" << endl;
        return nullptr;
    }
    Bond* btemp[65536];

    int i,j, bonds=0;
    for (i=0; i<65536; i++) btemp[i] = nullptr;
    if (aadef && aadef->aabonds)
    {
        for (i=0; aadef->aabonds[i] && aadef->aabonds[i]->Za && aadef->aabonds[i]->Zb; i++)
        {
            // cout << (name ? name : "(no name)") << "." << *(aadef->aabonds[i]) << endl;
            if (    (
                        aadef->aabonds[i]->cardinality == 1
                        &&
                        (aadef->aabonds[i]->can_rotate || aadef->aabonds[i]->can_flip)
                    )
                    ||
                    (
                        (aadef->aabonds[i]->cardinality < 2 && aadef->aabonds[i]->Za == 6 && aadef->aabonds[i]->Zb == 8)
                        ||
                        (aadef->aabonds[i]->cardinality < 2 && aadef->aabonds[i]->Za == 8 && aadef->aabonds[i]->Zb == 6)
                    )
               )
            {
                Atom* la = get_atom(aadef->aabonds[i]->aname);
                if (la
                        && (!la->is_backbone || !strcmp(la->name, "CA"))
                   )
                {
                    Bond* lb = la->get_bond_between(aadef->aabonds[i]->bname);
                    if (!lb)
                    {
                        // TODO: Add the missing bond if possible.
                        cout << "Warning: No bond between " << la->residue << ":" << la->name
                             << " and " << aadef->aabonds[i]->bname
                             << endl << flush;

                        Bond* lbb[32];
                        la->fetch_bonds(lbb);
                        if (lbb[0])
                        {
                            cout << la->name << " is bonded to:";
                            int o, ag = la->get_geometry();
                            for (o=0; o<ag; o++)
                                if (lbb[o]
                                    && (abs((__int64_t)(this) - (__int64_t)lbb[o]) < memsanity)
                                    && lbb[o]->atom2
                                    && (abs((__int64_t)(this) - (__int64_t)lbb[o]->atom2) < memsanity)
                                    )
                                    cout << " " << lbb[o]->atom2->name;
                            cout << "." << endl;
                        }

                        Atom* lba = get_atom(aadef->aabonds[i]->bname);
                        if (lba)
                        {
                            lba->fetch_bonds(lbb);
                            if (lbb[0])
                            {
                                cout << lba->name << " is bonded to:";
                                int o, ag = lba->get_geometry();
                                for (o=0; o<ag; o++) if (lbb[o]->atom2) cout << " " << lbb[o]->atom2->name;
                                cout << "." << endl;
                            }
                        }
                        else cout << aadef->aabonds[i]->bname << " not found." << endl;
                    }
                    else
                    {
                        // cout << (name ? name : "(no name)") << ":" << *(lb) << endl;
                        // Generally, a single bond between pi atoms cannot rotate.
                        if (lb->atom1->is_pi() && lb->atom2 && lb->atom2->is_pi())
                        {
                            lb->can_rotate = false;
                            lb->compute_flip_capability();
                            lb->flip_angle = M_PI;
                        }

                        lb->can_rotate = aadef->aabonds[i]->can_rotate;

                        if ((!la->is_backbone || !strcmp(la->name, "CA"))
                                &&
                                la->Z > 1
                                &&
                                (	greek_from_aname(la->name) == (greek_from_aname(lb->atom2->name)+1)
                                    ||
                                    greek_from_aname(la->name) == (greek_from_aname(lb->atom2->name)-1)
                                )
                           )
                        {
                            // cout << "Included." << endl;

                            if (greek_from_aname(la->name) < greek_from_aname(lb->atom2->name))
                                btemp[bonds] = la->get_bond_between(lb->atom2);
                            else
                                btemp[bonds] = lb->atom2->get_bond_between(la);

                            btemp[++bonds] = 0;

                            // cout << (name ? name : "(no name)") << ":" << *(btemp[bonds-1]) << endl;
                        }
                    }
                }
            }
        }

        goto _found_aadef;
    }

    // cout << name << " ";
    for (i=0; atoms[i]; i++)
    {
        Bond* lb[16];
        atoms[i]->fetch_bonds(lb);
        int g = atoms[i]->get_geometry();
        for (j=0; j<g; j++)
        {
            if (!lb[j]) break;
            if (lb[j]->can_rotate
                    &&
                    lb[j]->atom1 && lb[j]->atom2
                    &&
                    (!lb[j]->atom1->is_backbone || !strcmp(lb[j]->atom1->name, "CA"))
                    &&
                    !lb[j]->atom2->is_backbone
                    &&
                    greek_from_aname(lb[j]->atom1->name) == (greek_from_aname(lb[j]->atom2->name)-1)
                    &&
                    lb[j]->atom2->Z > 1
               )
            {
                // cout << *lb[j] << " ";
                btemp[bonds++] = lb[j];
                btemp[bonds] = 0;
            }
            else lb[j]->can_rotate = false;
        }
    }
    // cout << endl;

_found_aadef:
    Bond** retval = new Bond*[bonds+8];
    for (i=0; i<=bonds; i++) retval[i] = btemp[i];
    retval[i] = nullptr;
    rotatable_bonds = retval;

    return retval;
}

float Molecule::hydrophilicity()
{
    int i, count = 0;
    float total = 0;
    for (i=0; atoms[i]; i++)
    {
        int Z = atoms[i]->Z;
        if (Z==1) continue;

        total += atoms[i]->hydrophilicity_rule();
        count++;
    }
    return count ? (total / count) : 0;
}

float Molecule::solvent_free_energy(float epsilon, float kappa, bool csa)
{
    if (!atoms) return 0;

    // Sources:
    // https://doi.org/10.1002/prot.20033
    // https://doi.org/10.1126/science.2011744
    // https://en.wikipedia.org/wiki/Implicit_solvation#Generalized_Born_model
    float DeltaGsurf = get_surface_area(false, csa) * _solve_nonpol;
    double DeltaGGB = 0;

    #if 1
    double polsurf = get_surface_area(true, csa);
    DeltaGGB = polsurf * (epsilon/80) * _solve_np2pol;
    #else
    // Hoping someone with a solid math background can find out why the following is not giving accurate numbers at all.
    // Nobody is going to help the uneducated disgraced former code monkey.
    int i, j;
    for (i=0; atoms[i]; i++)
    {
        float Ri = atoms[i]->vdW_radius;
        for (j=0; atoms[j]; j++)
        {
            if (j==i) continue;
            if (!atoms[i]->is_bonded_to(atoms[j])) continue;

            float rij = atoms[i]->distance_to(atoms[j]);
            float Rj = atoms[j]->vdW_radius;
            float aij = sqrt(Ri*Rj);
            float D = pow(rij/(2*aij), 2);
            double fGB = sqrt(rij*rij + aij*aij * exp(-D));

            double qi = fabs(atoms[i]->is_polar());
            double qj = fabs(atoms[j]->is_polar());
            double lDeltaGGB = -0.5 * (qi*qj / fGB) * (1.0 - exp(-kappa*fGB)/epsilon) * 48.4;
            #if 1
            if (fabs(lDeltaGGB) > 0.1) cout << atoms[i]->name << " (" << qi << ") "
                << "~"
                << atoms[j]->name << " (" << qj << ") "
                << ": " << (lDeltaGGB * _kcal_per_kJ) << endl;
            #endif
            DeltaGGB += lDeltaGGB;
        }
    }
    #endif

    // cout << "# DeltaGsurf = " << (DeltaGsurf * _kcal_per_kJ) << endl;
    // cout << "# DeltaGGB = " << (DeltaGGB * _kcal_per_kJ) << endl;
    return DeltaGGB + DeltaGsurf;
}

float Molecule::solvent_bound_energy(Molecule **neighbors)
{
    float DeltaGsurf = get_exposed_surface_area(neighbors, false, true) * _solve_nonpol;
    double DeltaGGB = get_exposed_surface_area(neighbors, true, true) * _solve_np2pol;
    return 0.0f;
}

Bond** Molecule::get_all_bonds(bool unidirectional)
{
    if (noAtoms(atoms)) return 0;
    Bond* btemp[65536];

    int i,j, bonds=0;
    for (i=0; atoms[i]; i++)
    {
        Bond* lb[16];
        atoms[i]->fetch_bonds(lb);
        int g = atoms[i]->get_geometry();
        for (j=0; j<g; j++)
        {
            if (!lb[j]) break;
            if (lb[j]->atom1 < lb[j]->atom2
                    ||
                    !unidirectional
               )
            {
                btemp[bonds++] = lb[j];
                btemp[bonds] = 0;
            }
        }
    }

    Bond** retval = new Bond*[bonds+1];
    for (i=0; i<=bonds; i++) retval[i] = btemp[i];

    return retval;
}


float Molecule::get_internal_clashes(bool sb)
{
    int i, j;
    float r;
    Atom* a, *b;
    float clash = 0;

    if (noAtoms(atoms)) return 0;

    for (i=0; atoms[i]; i++)
    {
        Point pta = atoms[i]->loc;
        float avdW = atoms[i]->vdW_radius;
        for (j=i+1; atoms[j]; j++)
        {
            if (atoms[i]->residue && atoms[i]->residue == atoms[j]->residue && !strcmp(atoms[i]->name, atoms[j]->name)) continue;
            if (atoms[i]->is_bonded_to(atoms[j]) /* || atoms[j]->is_bonded_to(atoms[i]) */)
            {
                Bond* ab = atoms[i]->get_bond_between(atoms[j]);
                if (atoms[i] < atoms[j] && atoms[i]->is_pi() && atoms[j]->is_pi() && !atoms[i]->in_same_ring_as(atoms[j]))
                {
                    Atom* c = atoms[i]->get_heaviest_bonded_atom_that_isnt(atoms[j]);
                    Atom* d = atoms[j]->get_heaviest_bonded_atom_that_isnt(atoms[i]);

                    if (c && d)
                    {
                        Vector axis = atoms[j]->loc.subtract(atoms[i]->loc);
                        float theta = find_angle_along_vector(c->loc, d->loc, atoms[i]->loc, axis);
                        float cpartial = 13.5 - 13.5 * cos(theta*2);
                        #if _dbg_internal_energy
                        if (!is_residue())
                            cout << "Conjugated " << atoms[i]->name << "-" << atoms[j]->name
                                << " " << ab->cardinality << " bond theta = " << (theta*fiftyseven) << "deg."
                                << " adding " << cpartial << " kJ/mol."
                                << endl;
                        #endif
                        clash += cpartial;
                    }
                }
                continue;
            }

            Point ptb = atoms[j]->loc;
            float bvdW = atoms[j]->vdW_radius;

            r = pta.get_3d_distance(&ptb);
            if (r >= 0.9 && atoms[i]->shares_bonded_with(atoms[j])) continue;

            if (!r) r += 10e-15;
            if (r < avdW + bvdW)
            {
                float lclash = fmax(InteratomicForce::Lennard_Jones(atoms[i], atoms[j]), 0); // sphere_intersection(avdW, bvdW, r);
                clash += lclash;

                if (false && lclash > 3)
                {
                    cout << atoms[i]->name << " clashes with " << atoms[j]->name << " by " << lclash << " cu. A. resulting in " << clash << endl;
                    int g = atoms[i]->get_geometry();
                    cout << "Geometry: " << g << endl;
                    Bond* b[16];
                    atoms[i]->fetch_bonds(b);
                    int k;
                    for (k=0; k<g; k++)
                        cout << atoms[i]->name << " is bonded to " << hex << b[k]->atom2 << dec << " "
                             << (b[k]->atom2 ? b[k]->atom2->name : "") << "." << endl;
                }
            }
        }
    }

    return clash - (sb ? base_internal_clashes : 0);
}

#if compute_vdw_repulsion
float Molecule::get_vdW_repulsion(Molecule* ligand)
{
    if (!ligand) return 0;
    if (ligand == this) return 0;
    if (!atoms || !ligand->atoms) return 0;

    int i, j;
    float retval = 0;

    for (i=0; atoms[i]; i++)
    {
        float achg = atoms[i]->get_charge();
        bool api = atoms[i]->is_pi();

        for (j=0; ligand->atoms[j]; j++)
        {
            float bchg = ligand->atoms[j]->get_charge();
            bool bpi = ligand->atoms[j]->is_pi();

            if (!achg || !bchg)
            {
                // TODO: Hard coded values get from bindings.dat instead.
                float rlim = 4, kJmol = 0.4;
                if (api && bpi)
                {
                    rlim = 3.87;
                    kJmol = 2;
                }
                float halflim = rlim/2;
                float asphere = 4.0/3 * M_PI * halflim * halflim * halflim;

                float r = atoms[i]->distance_to(ligand->atoms[j]);
                if (r < rlim)
                {
                    retval += fabs(sphere_intersection(halflim, halflim, r) * kJmol / asphere);
                }
            }
        }
    }

    return retval;
}
#endif

float Molecule::get_intermol_clashes(Molecule* ligand)
{
    Molecule * ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    return get_intermol_clashes(ligands);
}

float Molecule::get_intermol_clashes(Molecule** ligands)
{
    int i, j, l;
    float r;
    Atom* a, *b;
    float clash = 0;

    clash1 = clash2 = nullptr;
    float worst = 0;

    if (noAtoms(atoms)) return 0;
    if (!ligands) return 0;
    if (!ligands[0]) return 0;

    for (i=0; atoms[i]; i++)
    {
        Atom* a = atoms[i];
        Point pta = a->loc;
        float avdW = a->vdW_radius;
        for (l=0; ligands[l]; l++)
        {
            if (!ligands[l]->atoms) continue;
            if (ligands[l] == this) continue;
            for (j=0; ligands[l]->atoms[j]; j++)
            {
                Atom* b = ligands[l]->atoms[j];
                // if (atoms[i]->is_bonded_to(ligands[l]->atoms[j])) continue;
                // if (atoms[i]->shares_bonded_with(ligands[l]->atoms[j])) continue;

                if (a->is_backbone && b->is_backbone && abs(a->residue - b->residue) < 2) continue;

                float f = fmax(InteratomicForce::Lennard_Jones(a, b), 0);

                if (a->residue != b->residue)
                {
                    if (f > worst)
                    {
                        worst = f;
                        clash1 = a;
                        clash2 = b;
                    }

                    if (f > worst_mol_clash)
                    {
                        worst_mol_clash = f;
                        worst_clash_1 = this;
                        worst_clash_2 = ligands[l];
                    }
                }

                clash += f;
                continue;
            }
        }
    }

    return clash;
}

float Molecule::total_intermol_clashes(Molecule** ligands)
{
    if (!ligands) return 0;
    int i;
    float clash = 0;
    for (i=0; ligands[i]; i++)
    {
        clash += ligands[i]->get_intermol_clashes(ligands);
    }
    return clash;
}

float Molecule::get_worst_clash()
{
    if (nmonomers && monomers)
    {
        clash_worst = 0;
        int i;
        for (i=0; i<nmonomers; i++) clash_worst += monomers[i]->clash_worst;
    }
    return clash_worst;
}

float Molecule::get_total_mclashes()
{
    if (!mclashables) return 0;
    return total_intermol_clashes(mclashables);
}

Interaction Molecule::optimize_intermol_contact(Molecule *ligand)
{
    Interaction result;

    Atom *a, *b;
    mutual_closest_atoms(ligand, &a, &b);
    if (!a || !b) return result;
    a = a->get_heavy_atom();
    b = b->get_heavy_atom();

    Bond *b1, *b2;
    b1 = a->get_bond_by_idx(0);
    if (b1) b1 = b1->get_reversed();
    b2 = b->get_bond_by_idx(0);
    if (b2) b2 = b2->get_reversed();

    bool dorot1 = !a->is_pi() && a->get_bonded_heavy_atoms_count() <= 1;
    bool dorot2 = !b->is_pi() && b->get_bonded_heavy_atoms_count() <= 1;
    float step1 = b1->can_rotate ? hexagonal/2 : (b1->can_flip ? b1->flip_angle : M_PI*2);
    float step2 = b2->can_rotate ? hexagonal/2 : (b2->can_flip ? b2->flip_angle : M_PI*2);
    float theta1, theta2, opt1=0, opt2=0;


    Interaction e;
    for (theta1=0; theta1<M_PI*2; theta1+=step1)
    {
        for (theta2=0; theta2<M_PI*2; theta2+=step2)
        {
            e = get_intermol_binding(ligand);
            if (e.improved(result))
            {
                result = e;
                opt1 = theta1;
                opt2 = theta2;
            }
            if (!step2) break;
            if (dorot2) b2->rotate(step2);
            else break;
        }
        if (!step1) break;
        if (dorot1) b1->rotate(step1);
        else break;
    }

    if (dorot1 && opt1 && (b1->can_rotate || b1->can_flip)) b1->rotate(opt1);
    if (dorot2 && opt2 && (b2->can_rotate || b2->can_flip)) b2->rotate(opt2);

    if (b1->can_rotate || b2->can_rotate)
    {
        step1 /= 2;
        step2 /= 2;
        while (step1 > 0.1*fiftyseven && step2 > 0.1*fiftyseven)
        {
            if (b1->can_rotate)
            {
                b1->rotate(-step1);
                e = get_intermol_binding(ligand);
                if (e.improved(result)) result = e;
                else
                {
                    b1->rotate(step1*2);
                    e = get_intermol_binding(ligand);
                    if (e.improved(result)) result = e;
                    else b1->rotate(-step1);
                }
            }
            if (b2->can_rotate)
            {
                b2->rotate(-step2);
                e = get_intermol_binding(ligand);
                if (e.improved(result)) result = e;
                else
                {
                    b2->rotate(step2*2);
                    e = get_intermol_binding(ligand);
                    if (e.improved(result)) result = e;
                    else b2->rotate(-step2);
                }
            }

            step1 *= 0.666;
            step2 *= 0.666;
        }
    }

    return result;
}

void Molecule::mutual_closest_atoms(Molecule* mol, Atom** a1, Atom** a2)
{
    if (!a1 || !a2) return;

    *a1 = *a2 = nullptr;

    int i, j, m, n;
    Atom *a, *b;
    float rbest = Avogadro;

    m = get_atom_count();
    n = mol->get_atom_count();
    for (i=0; i<m; i++)
    {
        a = get_atom(i);
        if (!a) break;
        int aZ = a->Z;
        for (j=0; j<n; j++)
        {
            b = mol->get_atom(j);
            if (!b) break;
            int bZ = b->Z;
            if (aZ == 1 && bZ == 1) continue;

            float r = a->distance_to(b);
            if (r < rbest)
            {
                rbest = r;
                *a1 = a;
                *a2 = b;
            }
        }
    }
}

void Molecule::move(Vector move_amt, bool override_residue)
{
    if (noAtoms(atoms)) return;
    if (immobile)
    {
        cout << "Warning: Attempt to move \"immobile\" molecule " << name << endl;
        return;
    }
    int i;

    #if _dbg_improvements_only_rule
    Interaction before;
    if (check_ligand && check_mols && !excuse_deterioration && is_residue() == _dbg_improvements_only_residue) before = cfmol_multibind(check_ligand, check_mols);
    #endif

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->residue && !override_residue) return;
        Point loc = atoms[i]->loc;
        loc = loc.add(&move_amt);
        atoms[i]->move(&loc);
    }

    #if _dbg_improvements_only_rule
    if (check_ligand && check_mols && !excuse_deterioration && is_residue() == _dbg_improvements_only_residue) 
    {
        Interaction after = cfmol_multibind(check_ligand, check_mols);
        if (!after.improved(before)) throw 0xb0661;
    }
    excuse_deterioration = false;
    #endif
}

void Molecule::move(Point move_amt, bool override_residue)
{
    if (noAtoms(atoms)) return;
    if (immobile)
    {
        cout << "Warning: Attempt to move \"immobile\" molecule " << name << endl;
        return;
    }
    int i;
    int vvc = vdw_vertex_count;

    #if _dbg_improvements_only_rule
    Interaction before;
    if (check_ligand && check_mols && !excuse_deterioration && is_residue() == _dbg_improvements_only_residue) before = cfmol_multibind(check_ligand, check_mols);
    #endif

    for (i=0; atoms[i]; i++)
    {
        // cout << atoms[i]->name << " ";
        if (atoms[i]->residue && !override_residue) return;
        Point loc = atoms[i]->loc;
        loc = loc.add(&move_amt);
        atoms[i]->move(&loc);
    }

    #if _dbg_improvements_only_rule
    if (check_ligand && check_mols && !excuse_deterioration && is_residue() == _dbg_improvements_only_residue) 
    {
        Interaction after = cfmol_multibind(check_ligand, check_mols);
        if (!after.improved(before)) throw 0xb0661;
    }
    excuse_deterioration = false;
    #endif

    vdw_vertex_count = vvc;
}

Point Molecule::get_barycenter(bool bond_weighted) const
{
    if (noAtoms(atoms))
    {
        Point pt;
        return pt;
    }

    Point locs[atcount];
    int i;

    for (i=0; atoms[i]; i++)
    {
        locs[i] = atoms[i]->loc;
        locs[i].weight = atoms[i]->get_atomic_weight();
        #if allow_tethered_rotations
        if (bond_weighted)
        {
            locs[i].weight -= fmin(0, atoms[i]->last_bind_energy);
        }
        #endif
    }

    return average_of_points(locs, atcount);
}

float Molecule::get_charge() const
{
    int i;
    float charge=0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->Z == 1) continue;
        charge += atoms[i]->get_charge();
    }
    return charge;
}

void Molecule::recenter(Point nl)
{
    if (movability <= MOV_NORECEN) return;
    Point loc = get_barycenter();
    Point rel = nl.subtract(&loc);
    Vector v(&rel);
    move(v);

    if (vdw_surface && vdw_vertex_count)
    {
        int i;
        for (i=0; i<vdw_vertex_count; i++)
        {
            Point loc = vdw_surface[i];
            Point rel = nl.subtract(&loc);
            vdw_surface[i] = rel;
        }
    }
}

void Molecule::rotate(Vector* v, float theta, bool bond_weighted)
{
    if (noAtoms(atoms)) return;
    // cout << name << " Molecule::rotate()" << endl;

    if (movability <= MOV_FLEXONLY) return;
    if (movability <= MOV_NORECEN) bond_weighted = false;
    Point cen = get_barycenter(bond_weighted);

    LocatedVector lv = *v;
    lv.origin = cen;

    rotate(lv, theta);
}

void Molecule::rotate(LocatedVector lv, float theta)
{
    if (noAtoms(atoms)) return;

    if (movability <= MOV_FLEXONLY) return;
    if (movability <= MOV_NORECEN) lv.origin = get_barycenter();

    #if _dbg_improvements_only_rule
    Interaction before;
    if (check_ligand && check_mols && !excuse_deterioration && is_residue() == _dbg_improvements_only_residue) before = cfmol_multibind(check_ligand, check_mols);
    #endif

    int i;
    if (fabs(theta) > hexagonal)
    {
        i = 0;
    }
    for (i=0; atoms[i]; i++)
    {
        // if (atoms[i]->residue) return;
        Point loc = atoms[i]->loc;
        Point nl  = rotate3D(&loc, &lv.origin, &lv, theta);
        atoms[i]->move(&nl);
    }

    if (vdw_surface && vdw_vertex_count)
    {
        for (i=0; i<vdw_vertex_count; i++)
        {
            Point loc = vdw_surface[i];
            Point nl  = rotate3D(&loc, &lv.origin, &lv, theta);
            vdw_surface[i] = nl;
        }
    }

    #if _dbg_improvements_only_rule
    if (check_ligand && check_mols && !excuse_deterioration && is_residue() == _dbg_improvements_only_residue) 
    {
        Interaction after = cfmol_multibind(check_ligand, check_mols);
        if (!after.improved(before)) throw 0xb0661;
    }
    excuse_deterioration = false;
    #endif
}

void Molecule::rotate(Rotation rot)
{
    LocatedVector lv = rot.v;
    lv.origin = get_barycenter(true);
    rotate(lv, rot.a);
}

void Molecule::rotate(Rotation rot, Point origin)
{
    LocatedVector lv = rot.v;
    lv.origin = origin;
    rotate(lv, rot.a);
}

bool Molecule::shielded(Atom* a, Atom* b) const
{
    int i;
    float r = a->distance_to(b);
    float r6 = r*1.26, r125 = 1.25*r;
    float va = a->vdW_radius, vb = b->vdW_radius;
    if (r < 2) return false;

    a->shielding_angle = b->shielding_angle = 0;

    Point aloc = a->loc, bloc = b->loc;
    for (i=0; atoms[i]; i++)
    {
        Atom* ai = atoms[i];
        if (!ai) break;
        if (ai == a || ai == b) continue;
        float rai = ai->distance_to(a);
        if (rai > r6) continue;
        float rbi = ai->distance_to(b);
        if (rbi > r6) continue;
        if ((rai+rbi) > r125) continue;
        float vs = ai->vdW_radius;
        if (rai < va+vs) return true;
        if (rbi < vb+vs) return true;
        Point sloc = ai->loc;
        float f3da = find_3d_angle(&aloc, &bloc, &sloc);
        if (f3da > a->shielding_angle) a->shielding_angle = b->shielding_angle = f3da;
        if (f3da > _shield_angle)
        {
            if (last_iter && (a->residue == 114 || b->residue == 114) && ((a->residue + b->residue) == 114))
            {
                /*cout << ai->name << " shields "
                	 << a->residue << ":" << a->name << "..."
                	 << b->residue << ":" << b->name
                	 << " angle " << (f3da*fiftyseven)
                	 << endl;*/
                return true;
            }
        }
    }

    return false;
}

float Molecule::pi_stackability(bool ib)
{
    if (noAtoms(atoms)) return 0;
    int i, j=0;
    float result = 0;

    for (i=0; atoms[i]; i++)
    {
        if (!ib && atoms[i]->is_backbone) continue;
        if (atoms[i]->is_pi()) result += 1;
        j++;
    }

    if (j) result /= j;
    return result;
}

float Molecule::get_atom_mol_bind_potential(Atom* a)
{
    if (noAtoms(atoms)) return 0;
    int i, j, n;

    float hydro = 0;
    j = 0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        if (atoms[i]->Z == 1) continue;

        hydro += fabs(atoms[i]->is_polar());
        j++;
    }

    if (j) hydro /= j;

    float retval=0;
    n = 0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;

        InteratomicForce* ifs[32];
        InteratomicForce::fetch_applicable(a, atoms[i], ifs);
        if (!ifs[0]) continue;

        for (j=0; ifs[j]; j++)
        {
            float partial;

            if (hydro > hydrophilicity_cutoff && ifs[j]->get_type() == vdW) continue;

            if (ifs[j]->get_type() == ionic)
            {
                if (sgn(a->get_charge()) != -sgn(atoms[i]->get_charge())) continue;
                partial = 60;
            }
            else
            {
                partial = ifs[j]->get_kJmol();
            }

            if (ifs[j]->get_type() == hbond)
            {
                if (fabs(a->is_polar()) < hydrophilicity_cutoff || fabs(atoms[i]->is_polar()) < hydrophilicity_cutoff) continue;
                if (sgn(a->is_polar()) > 0 && !this->has_hbond_acceptors()) continue;
                if (sgn(a->is_polar()) < 0 && !this->has_hbond_donors()) continue;
                partial *= fmin(fabs(a->is_polar()), fabs(atoms[i]->is_polar()));
            }
            
            if (ifs[j]->get_type() == pi && fabs(atoms[i]->is_polar()) >= hydrophilicity_cutoff) partial /= 3;
            if (ifs[j]->get_type() == pi && fabs(this->get_charge()) >= hydrophilicity_cutoff) partial /= 3;

            if (ifs[j]->get_type() == polarpi) partial /= 6;            // Config is for benzene rings.

            if (ifs[j]->get_type() == mcoord)
            {
                partial *= (1.0 + 1.0 * cos((a->get_electronegativity() + atoms[i]->get_electronegativity()) / 2 - 2.25));
            }

            if (ifs[j]->get_type() == vdW)
            {
                if (fabs(a->is_polar()) >= hydrophilicity_cutoff || fabs(atoms[i]->is_polar()) >= hydrophilicity_cutoff) continue;
            }

            retval += partial;
            n++;
        }
    }

    if (n) retval /= n;
    return retval;
}

float Molecule::find_mutual_max_bind_potential(Molecule* other)
{
    int i, j, m = get_atom_count(), n = other->get_atom_count();
    if (!m || !n) return 0;

    float best_potential = 0, nextbest_potential = 0;
    for (i=0; i<m; i++)
    {
        if (atoms[i]->is_backbone) continue;
        int heavy1 = atoms[i]->get_bonded_heavy_atoms_count();
        for (j=0; j<n; j++)
        {
            if (other->atoms[j]->is_backbone) continue;
            int heavy2 = other->atoms[j]->get_bonded_heavy_atoms_count();
            float b = InteratomicForce::potential_binding(atoms[i], other->atoms[j]) * frand(0.8, 1.3);
            b /= (heavy1+heavy2);
            // float r = atoms[i]->distance_to(other->atoms[j]);
            // b /= sqrt(r);
            // cout << atoms[i]->name << " ~ " << other->atoms[j]->name << " = " << b;
            if (b > best_potential)
            {
                best_potential = b;
                stay_close2_mol = stay_close_mol;
                stay_close2_mine = stay_close_mine;
                stay_close2_other = stay_close_other;
                stay_close2_optimal = stay_close_optimal;
                stay_close_mol = other;
                stay_close_mine = atoms[i];
                stay_close_other = other->atoms[j];
                stay_close_optimal = InteratomicForce::optimal_distance(stay_close_mine, stay_close_other);
                stay_close_tolerance = stay_close_optimal*stays_tolerance_factor;

                // cout << " *";
                #if _dbg_cs_pairing || _dbg_stays_assignment
                cout << "Stay-close atoms ";
                if (stay_close_mine->residue)
                    cout << stay_close_mine->aa3let << stay_close_mine->residue << ":";
                else cout << "LIG:";
                cout << stay_close_mine->name << " (polarity " << stay_close_mine->is_polar() << ") and ";
                if (stay_close_other->residue)
                    cout << stay_close_other->aa3let << stay_close_other->residue << ":";
                else cout << "LIG:";
                cout << stay_close_other->name << " (polarity " << stay_close_other->is_polar() << ")"
                    << " have potential binding energy of " << -b << " kJ/mol at " << stay_close_optimal << "Å."
                    << endl;
                #endif
            }
            else if (b > nextbest_potential)
            {
                stay_close2_mol = other;
                stay_close2_mine = atoms[i];
                stay_close2_other = other->atoms[j];
                stay_close2_optimal = InteratomicForce::optimal_distance(stay_close2_mine, stay_close2_other);
            }
            // cout << endl;
        }
    }

    #if _dbg_cs_pairing
    cout << endl;
    #endif

    return best_potential;
}

bool Molecule::check_stays()
{
    if (!stay_close_mine || !stay_close_other) return true;
    if (!stay_close_water) return check_stays_dry();

    Atom* nm = stay_close_water->get_nearest_atom(stay_close_mine->loc, hbond, -sgn(stay_close_mine->is_polar()));
    Atom* no = stay_close_water->get_nearest_atom(stay_close_other->loc, hbond, -sgn(stay_close_other->is_polar()));
    if (!nm || !no)
    {
        stay_close_water = nullptr;
        return check_stays();             // RECURSION!
    }

    float r1 = stay_close_mine->distance_to(nm);
    float r2 = stay_close_other->distance_to(no);
    float optimal1 = InteratomicForce::optimal_distance(stay_close_mine, nm);
    float optimal2 = InteratomicForce::optimal_distance(stay_close_other, no);
    return (r1 <= optimal1+stay_close_tolerance && r2 <= optimal2+stay_close_tolerance);
}

bool Molecule::check_stays_dry()
{
    float r = stay_close_other->distance_to(stay_close_mine);
    bool result = (r <= stay_close_optimal+stay_close_tolerance);
    if (!result) return result;
    r = stay_close2_other->distance_to(stay_close2_mine);
    return (r <= stay_close2_optimal+stay_close_tolerance);
}

void Molecule::enforce_stays(float amt, void (*stepscb)(std::string mesg))
{
    if (!stay_close_mine || !stay_close_other) return;
    if (is_residue())
    {
        // TODO: flex side chain
        return;
    }

    if (stay_close_other->vanished)
    {
        stay_close_mine = stay_close_other = nullptr;
        return;
    }
    if (stay_close2_other && stay_close2_other->vanished)
    {
        stay_close2_mine = stay_close2_other = nullptr;
    }

    Rotation rot;
    LocatedVector lv;
    MovabilityType wasmov;

    if (stay_close_mol && !stay_close_other->is_metal())
    {
        intera_type btyp = vdW;
        float polar = stay_close_mine->is_polar();
        if (stay_close_mine->get_charge() && stay_close_other->get_charge()) btyp = ionic;
        else if (polar && stay_close_other->is_polar()) btyp = hbond;
        else if (stay_close_mine->is_pi() && stay_close_other->is_pi()) btyp = pi;

        Atom* sco = stay_close_mol->get_nearest_atom(get_barycenter(), btyp, fabs(polar) >= hydrophilicity_cutoff ? -sgn(polar) : 0);
        if (sco) stay_close_other = sco;
    }

    if (stay_close_water)
    {
        Atom* nm = stay_close_water->get_nearest_atom(stay_close_mine->loc, hbond, -sgn(stay_close_mine->is_polar()));
        Atom* no = stay_close_water->get_nearest_atom(stay_close_other->loc, hbond, -sgn(stay_close_other->is_polar()));
        if (!nm || !no)
        {
            stay_close_water = nullptr;
            enforce_stays(amt);             // RECURSION!
            if (stepscb) stepscb("forcing dry contact");
            return;
        }

        Vector movamt1 = nm->loc.subtract(stay_close_mine->loc);
        Vector movamt2 = stay_close_other->loc.subtract(no->loc);
        float optimal1 = InteratomicForce::optimal_distance(stay_close_mine, nm);
        float optimal2 = InteratomicForce::optimal_distance(stay_close_other, no);  

        movamt1.r -= (optimal1+stay_close_tolerance);
        movamt2.r -= (optimal2+stay_close_tolerance);

        if (movamt2.r > 0)
        {
            rot = align_points_3d(no->loc, stay_close_other->loc, stay_close_water->get_barycenter());
            lv = rot.v;
            lv.origin = stay_close_water->get_barycenter();
            stay_close_water->rotate(lv, rot.a*amt);
            if (stepscb) stepscb("rotated water");
        }

        movamt1 = nm->loc.subtract(stay_close_mine->loc);
        movamt2 = stay_close_other->loc.subtract(no->loc);
        movamt1.r -= (optimal1+stay_close_tolerance);
        movamt2.r -= (optimal2+stay_close_tolerance);

        if (movamt2.r > 0)
        {
            movamt2.r *= amt;
            wasmov = stay_close_water->movability;
            stay_close_water->movability = MOV_ALL;
            stay_close_water->move(movamt2);
            stay_close_water->movability = wasmov;
            if (stepscb) stepscb("moved water");
        }

        movamt1 = nm->loc.subtract(stay_close_mine->loc);
        movamt1.r -= (optimal1+stay_close_tolerance);

        if (movamt1.r > 0)
        {
            rot = align_points_3d(nm->loc, stay_close_mine->loc, no->loc);
            lv = rot.v;
            lv.origin = no->loc;
            stay_close_water->rotate(lv, rot.a*amt);
            if (stepscb) stepscb("rotated water to ligand");
        }

        movamt1 = nm->loc.subtract(stay_close_mine->loc);
        movamt1.r -= (optimal1+stay_close_tolerance);

        wasmov = movability;
        movability = MOV_ALL;

        if (movamt1.r > 0)
        {
            rot = align_points_3d(stay_close_mine->loc, nm->loc, get_barycenter());
            lv = rot.v;
            lv.origin = get_barycenter();
            rotate(lv, rot.a*amt);
            if (stepscb) stepscb("rotated ligand");
        }

        movamt1 = nm->loc.subtract(stay_close_mine->loc);
        movamt1.r -= (optimal1+stay_close_tolerance);

        if (movamt1.r > 0)
        {
            movamt1.r *= amt;
            move(movamt1);
            if (stepscb) stepscb("moved ligand");
        }

        movability = wasmov;
        return;
    }

    wasmov = movability;
    if (!is_residue()) movability = MOV_ALL;

    if (stay_close_mine->Z < 2)
    {
        Atom* heavy = stay_close_mine->get_heaviest_bonded_atom_that_isnt(nullptr);
        Bond* b = heavy->get_bond_by_idx(0);
        if (b) b = b->get_reversed();
        if (b && b->can_rotate && b->atom1)
        {
            LocatedVector lv = b->get_axis();
            lv.origin = b->atom1->loc;
            float theta = find_angle_along_vector(stay_close_mine->loc, stay_close_other->loc,
                lv.origin, lv);
            float was = stay_close_mine->distance_to(stay_close_other);
            b->rotate(theta);
            float isnow = stay_close_mine->distance_to(stay_close_other);
            if (isnow > was) b->rotate(theta*-2);
        }
    }

    #if allow_stay_close_flexions
    This causes side chain clashes.
    if (stay_close_mol)
    {
        if (stay_close_mol->movability & MOV_CAN_FLEX)
            stay_close_mol->conform_atom_to_location(stay_close_other->name, stay_close_mine->loc,
                stay_close_conform_iters, stay_close_optimal+stay_close_tolerance);
    }
    #endif

    Vector movamt = stay_close_other->loc.subtract(stay_close_mine->loc);
    movamt.r -= (stay_close_optimal+stay_close_tolerance);
    if (movamt.r < 0) movamt.r = 0;

    Interaction e;
    if (stay_close_mol)
    {
        e = get_intermol_binding(stay_close_mol);
        if (e.clash >= clash_limit_per_aa)
        {
            if (!clash1 || !clash2) return;
            Vector rplamt = clash1->loc.subtract(clash2->loc);
            rplamt.r = clash1->vdW_radius + clash2->vdW_radius - rplamt.r;          // Sum of vdW radii minus already distance.
            movamt = movamt.add(rplamt);
        }
    }

    if (!movamt.r) return;

    rot = align_points_3d(stay_close_mine->loc, stay_close_other->loc, get_barycenter());
    rot.a = fmin(fabs(rot.a), fiftyseventh*22.5)*sgn(rot.a);
    lv = rot.v;
    lv.origin = get_barycenter();
    #if _dbg_improvements_only_rule
    excuse_deterioration = true;
    #endif
    rotate(lv, fmin(fabs(rot.a), hexagonal)*sgn(rot.a)*amt);

    movamt = stay_close_other->loc.subtract(stay_close_mine->loc);
    movamt.r -= (stay_close_optimal+stay_close_tolerance);
    if (movamt.r > 0)
    {
        movamt.r *= amt;

        if (fabs(movamt.r) > speed_limit) movamt.r = sgn(movamt.r) * speed_limit;

        #if _dbg_improvements_only_rule
        excuse_deterioration = true;
        #endif
        move(movamt);
    }

    if (stay_close2_mine && stay_close2_other)
    {
        Vector v = stay_close2_other->loc.subtract(stay_close2_mine->loc);
        v.r -= (stay_close2_optimal+stay_close_tolerance);
        rot = align_points_3d(stay_close2_mine->loc, stay_close2_other->loc, stay_close_mine->loc);
        rot.a *= amt;
        lv = rot.v;
        lv.origin = stay_close2_mine->loc;
        rotate(lv, rot.a);
    }

    movability = wasmov;

    #if _dbg_stays_enforce
    movamt = stay_close_other->loc.subtract(stay_close_mine->loc);
    cout << stay_close_mine->name << " - " << stay_close_other->residue << ":" << stay_close_other->name << " = " << movamt << endl << endl;
    #endif
}

Interaction Molecule::get_intermol_binding(Molecule* ligand, bool subtract_clashes, bool priority_boost)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    return get_intermol_binding(ligands, subtract_clashes, priority_boost);
}

void Molecule::clear_atom_binding_energies()
{
    if (!atoms) return;
    int i;
    for (i=0; atoms[i]; i++)
        atoms[i]->last_bind_energy = 0;
}

float Molecule::get_intermol_potential(Molecule* ligand, bool pure)
{
    Molecule* ligands[4];
    ligands[0] = ligand;
    ligands[1] = nullptr;
    return get_intermol_potential(ligands, pure);
}

float Molecule::get_intermol_potential(Molecule** ligands, bool pure)
{
    if (!atoms) return 0;
    if (!ligands) return 0;
    if (!ligands[0]) return 0;
    int i, j, l, n;
    float kJmol = 0;

    for (i=0; atoms[i]; i++)
    {
        if (!atoms[i]) continue;
        Point aloc = atoms[i]->loc;
        for (l=0; ligands[l]; l++)
        {
            if (ligands[l] == this) continue;
            for (j=0; j<ligands[l]->atcount; j++)
            {
                if (!ligands[l]->atoms[j]) continue;
                float r = ligands[l]->atoms[j]->loc.get_3d_distance(&aloc);
                float f = 1.0 / r;		// Regular invert rather than inv square so that actual bonding will take over at short range.
                InteratomicForce* iff[32];
                InteratomicForce::fetch_applicable(atoms[i], ligands[l]->atoms[j], iff);

                if (iff) for (n=0; iff[n]; n++)
                {
                    if (iff[n]->get_type() == vdW) continue;

                    if (pure || r < iff[n]->get_distance())
                        kJmol += iff[n]->get_kJmol();
                    else
                        kJmol += iff[n]->get_kJmol()*f;
                }
            }
        }
    }

    return kJmol;
}

Interaction Molecule::get_intermol_binding(Molecule** ligands, bool subtract_clashes, bool priority_boost)
{
    if (!ligands) return 0;
    if (!ligands[0]) return 0;
    int i, j, l;
    if (nmonomers)
    {
        bool self_among_ligands = false;
        for (i=0; ligands[i]; i++)
        {
            if (ligands[i] == this) self_among_ligands = true;
            if (self_among_ligands) break;
        }

        if (self_among_ligands)
        {
            Interaction result = 0;
            for (i=0; i<nmonomers; i++)
            {
                if (monomers[i]->nmonomers) throw 0xbadc0de;
                for (j=i; j<nmonomers; j++)
                {
                    // cout << i << "," << j << endl;
                    result += monomers[i]->get_intermol_binding(monomers[j], subtract_clashes, priority_boost);
                }
                result += monomers[i]->get_intermol_binding(ligands, subtract_clashes, priority_boost);
            }
            return result;      // RECURSION.
        }
    }

    Interaction kJmol;
    if (!atoms) return 0;
    if (subtract_clashes) kJmol.clash += get_internal_clashes();
    if (!check_stays()) kJmol.stays_met = false;

    lastmc = 0;
    lastshielded = 0;
    clash1 = clash2 = nullptr;
    float best_atom_energy = 0;

    #if _dbg_internal_energy
    cout << (name ? name : "") << " base internal clashes: " << base_internal_clashes << "; final internal clashes " << kJmol.summed() << endl;
    #endif

    for (i=0; atoms[i]; i++)
    {
        atoms[i]->last_bind_energy = 0;
        atoms[i]->strongest_bind_energy = 0;
        atoms[i]->strongest_bind_atom = nullptr;
    }

    clash_worst = 0;
    for (i=0; atoms[i]; i++)
    {
        Point aloc = atoms[i]->loc;
        for (l=0; ligands[l]; l++)
        {
            bool skip = false;
            if (ligands[l]->nmonomers)
            {
                for (j=0; j<ligands[l]->nmonomers; j++)
                {
                    if (ligands[l]->monomers[j] == this) skip = true;
                }
            }
            if (skip) continue;

            #if _dbg_51e2_ionic
            if (!is_residue() && atoms[i]->get_family() == CHALCOGEN && ligands[l]->is_residue() == 262)
            {
                j = 0;
            }
            #endif

            if (ligands[l] == this) continue;
            for (j=0; j<ligands[l]->atcount; j++)
            {
                // TODO: Fix this in the hydrogenate function, but for now we'll fix it here and hope for the best. 🤞🏼
                if (!ligands[l]->atoms[j])
                {
                    ligands[l]->atcount = j;
                    break;
                }
                if (atoms[i]->is_backbone && ligands[l]->atoms[j]->is_backbone
                    &&
                    (	(	atoms[i]->residue == ligands[l]->atoms[j]->residue - 1
                            &&
                            !strcmp(atoms[i]->name, "C")
                            &&
                            !strcmp(ligands[l]->atoms[j]->name, "N")
                        )
                        ||
                        (	atoms[i]->residue == ligands[l]->atoms[j]->residue + 1
                            &&
                            !strcmp(atoms[i]->name, "N")
                            &&
                            !strcmp(ligands[l]->atoms[j]->name, "C")
                        )
                    )) continue;			// kludge to prevent adjacent residue false clashes.
                float r = ligands[l]->atoms[j]->loc.get_3d_distance(&aloc);
                if (r < _INTERA_R_CUTOFF)
                {
                    if (	!shielded(atoms[i], ligands[l]->atoms[j])
                            &&
                            !ligands[l]->shielded(atoms[i], ligands[l]->atoms[j])
                       )
                    {
                        missed_connection.r = -Avogadro;
                        mc_bpotential = 0;

                        // cout << ligands[l]->atoms[j]->loc.subtract(atoms[i]->loc) << ": ";
                        Interaction abind;
                        if (atoms[i]->residue && atoms[i]->get_family() == CHALCOGEN && atoms[i]->Z > 10
                            && ligands[l]->atoms[j]->residue && ligands[l]->atoms[j]->get_family() == CHALCOGEN && ligands[l]->atoms[j]->Z > 10
                            && atoms[i]->residue != ligands[l]->atoms[j]->residue
                            && atoms[i]->is_bonded_to(ligands[l]->atoms[j]))
                        {
                            abind.attractive = 251;     // Cystine kludge.
                            abind.clash = 0;
                        }
                        else abind = InteratomicForce::total_binding(atoms[i], ligands[l]->atoms[j]);

                        if (((atoms[i]->residue && !ligands[l]->atoms[j]->residue) || (!atoms[i]->residue && ligands[l]->atoms[j]->residue))
                            && (coordmtl || ligands[l]->coordmtl)
                            )
                        {
                            total_binding_by_type[ionic - covalent] -= abind.repulsive;
                            abind.repulsive = 0;        // Mcoord thiolated ligand/thiolated cysteine kludge.
                        }

                        float asum = abind.summed();
                        #if _dbg_internal_energy
                        if (ligands[l] == this)
                        {
                            cout << "Energy between " << atoms[i]->name << "..." << ligands[l]->atoms[j]->name
                                << " = " << abind << " kJ/mol." << endl;
                        }
                        #endif

                        if (compute_interall)
                        {
                            int h;
                            Atom* heavy1 = atoms[i]->get_heavy_atom();
                            Atom* heavy2 = ligands[l]->atoms[j]->get_heavy_atom();

                            if (!heavy1->residue && !heavy2->residue
                                && heavy1->get_family() == TETREL && heavy2->get_family() == TETREL
                                && asum < 0)
                            {
                                cout << atoms[i]->name << " (" << heavy1->name << ") ~ "
                                    << ligands[l]->atoms[j]->name << " (" << heavy2->name << "): " << asum << endl;
                            }

                            for (h=0; h<ninterall; h++)
                            {
                                if ((heavy1 == interall_a1[h] && heavy2 == interall_a2[h])
                                    || (heavy1 == interall_a2[h] && heavy2 == interall_a1[h]))
                                {
                                    interall[h] += asum;

                                    if (!heavy1->residue && !heavy2->residue
                                        && heavy1->get_family() == TETREL && heavy2->get_family() == TETREL
                                        && interall[h] < -10)
                                    {
                                        throw 0xbadc0de;
                                    }

                                    goto _updated_interall;
                                }
                            }

                            interall_a1[ninterall] = heavy1;
                            interall_a2[ninterall] = heavy2;
                            interall[ninterall++] = asum;

                            if (ninterall >= MAX_INTERALL-1)
                            {
                                ninterall = MAX_INTERALL-1;
                                compute_interall = false;
                            }

                            _updated_interall:
                            ;
                        }

                        if ((abind.attractive || abind.repulsive || abind.clash) && !isnan(asum) && !isinf(asum))
                        {
                            if (priority_boost && asum < 0 && minimum_searching_aniso && ligands[l]->priority) abind.attractive *= 1.5;
                            kJmol += abind;

                            if (asum < best_atom_energy)
                            {
                                best_atom_energy = asum;
                                best_intera = atoms[i];
                                best_interactor = ligands[l];
                                best_other_intera = ligands[l]->atoms[j];
                            }

                            atoms[i]->last_bind_energy += asum;
                            if (abind.attractive > atoms[i]->strongest_bind_energy)
                            {
                                atoms[i]->strongest_bind_energy = abind.attractive;
                                atoms[i]->strongest_bind_atom = ligands[l]->atoms[j];
                            }

                            if (asum > 0 && abind.clash > clash_worst)
                            {
                                clash_worst = abind.clash;
                                if (atoms[i]->residue != ligands[l]->atoms[j]->residue)
                                {
                                    clash1 = atoms[i];
                                    clash2 = ligands[l]->atoms[j];
                                }
                            }

                            if (asum > 0 && ligands[l]->is_residue() && movability >= MOV_ALL)
                            {
                                Point ptd = aloc.subtract(ligands[l]->atoms[j]->loc);
                                ptd.multiply(fmin(fabs(asum) / 1000, 1));
                                lmx += lmpush * sgn(ptd.x);
                                lmy += lmpush * sgn(ptd.y);
                                lmz += lmpush * sgn(ptd.z);
                            }
                        }

                        if (missed_connection.r > 0)
                        {
                            lastmc += mc_bpotential;
                            Point mc = missed_connection;
                            // cout << mc << endl;
                            if (ligands[l]->priority) mc_bpotential *= 20;
                            float lc = ligands[l]->atoms[j]->get_charge();
                            if (lc && sgn(lc) == -sgn(atoms[i]->get_charge())) mc_bpotential *= 3.333;
                            float mcrr = missed_connection.r * missed_connection.r;
                            lmx += lmpull * mc.x * mc_bpotential / mcrr;
                            lmy += lmpull * mc.y * mc_bpotential / mcrr;
                            lmz += lmpull * mc.z * mc_bpotential / mcrr;
                        }
                    }
                    else lastshielded -= InteratomicForce::total_binding(atoms[i], ligands[l]->atoms[j]).summed();
                }
            }
        }
    }

    return kJmol;
}

float Molecule::get_intermol_contact_area(Molecule* ligand, bool hpho)
{
    if (!ligand) return 0;
    if (!atoms) return 0;
    int i, j;
    float result = 0;

    for (i=0; atoms[i]; i++)
    {
        Point aloc = atoms[i]->loc;
        float avdw = atoms[i]->vdW_radius;

        if (hpho && atoms[i]->is_polar()) continue;

        for (j=0; j<ligand->atcount; j++)
        {
            if (!ligand->atoms[j])
            {
                ligand->atcount = j;
                break;
            }

            if (hpho && ligand->atoms[j]->is_polar()) continue;

            float r = ligand->atoms[j]->loc.get_3d_distance(&aloc);
            if (!r) continue;
            float bvdw = ligand->atoms[j]->vdW_radius;
            if (r > avdw+bvdw) continue;

            float f = sphere_inter_area(avdw, bvdw, r);
            if (!isnan(f)) result += f;
        }
    }

    return result;
}

float Molecule::distance_to(Molecule* om)
{
    if (!atoms || !om || !om->atoms) return nanf("Bad molecule.");

    int i, j;
    float minr = Avogadro;

    for (i=0; atoms[i]; i++)
    {
        for (j=0; om->atoms[j]; j++)
        {
            float r = atoms[i]->distance_to(om->atoms[j]);

            if (r < minr) minr = r;
        }
    }

    return minr;
}

float Molecule::get_intermol_polar_sat(Molecule* ligand)
{
    if (!ligand) return 0;
    if (!atoms) return 0;
    int i, j, l, n;
    float result = 0;

    for (i=0; atoms[i]; i++)
    {
        Point aloc = atoms[i]->loc;
        float aapol = fabs(atoms[i]->is_polar());
        for (j=0; j<ligand->atcount; j++)
        {
            float r = ligand->atoms[j]->loc.get_3d_distance(&aloc);
            if (r > 5) continue;
            float abpol = fabs(ligand->atoms[j]->is_polar());
            float f = (aapol >= 0.5 ? 2 : 1) + (abpol >= 0.5 ? 2 : 1);
            if (!f) continue;
            f = 1.0 / pow(fmax(1,r/f),2);

            if (aapol < 0.2)
            {
                if (abpol < 0.2) result += f;
                else if (abpol >= 0.5) result -= f;
            }
            else if (aapol >= 0.5)
            {
                if (abpol < 0.2) result -= f;
                else if (abpol >= 0.5) result += f;
            }
        }
    }

    return result;
}

void Molecule::minimize_internal_clashes()
{
    // cout << (name ? name : "(no name)");
    if (noAtoms(atoms)) return;
    base_internal_clashes = base_eclipses = 0;

    int i, j, iter;
    float clash = get_internal_clashes() + total_eclipses();

    if (!clash) return;		// Already zero, nothing to decrease to.

    Bond** b = get_rotatable_bonds();
    if (!b || !b[0])
    {
        base_internal_clashes = get_internal_clashes();
        base_eclipses = total_eclipses();
        return;		// No bonds to rotate.
    }

    int numrb = 0;
    for (i=0; b[i]; i++) numrb = i;

    float angle[numrb];
    for (i=0; i<numrb; i++) angle[i] = M_PI;

    for (i=0; i<numrb; i++)
    {
        float step = hexagonal / 20;
        float theta = 0;
        for (j=0; step*j < M_PI*2; j++)
        {
            b[i]->rotate(step);
            float clash1 = get_internal_clashes() + total_eclipses();
            // cout << (step*j*fiftyseven) << "deg " << clash1 << endl;

            if (clash1 < clash)
            {
                clash = clash1;
                theta = step*j;
            }
        }
        b[i]->rotate(theta);
        angle[i] = step/2;
    }

    clash = get_internal_clashes() + total_eclipses();
    for (iter=0; iter<100; iter++)
    {
        for (i=0; i<numrb; i++)
        {
            b[i]->rotate(angle[i]);
            float clash1 = get_internal_clashes() + total_eclipses();

            if (clash1 <= clash)
            {
                clash = clash1;
                angle[i] *= 0.99;
            }
            else
            {
                b[i]->rotate(-angle[i]);		// Put it back.
                angle[i] *= -0.5;
            }
        }
    }

    base_internal_clashes = get_internal_clashes();
    base_eclipses = total_eclipses();
    // cout << name << " base internal clashes: " << base_internal_clashes << endl;
}

void Molecule::do_histidine_flip(HistidineFlip* hf)
{
    if (hf->N1->get_charge() >= 0.2 || hf->N2->get_charge() >= 0.2) return;

    Point ptC  = hf->C->loc;
    Point ptN1 = hf->N1->loc;
    Point ptN2 = hf->N2->loc;
    Point ptH  = hf->H->loc;

    Point arr[2] = {ptN1, ptN2};
    Point Navg = average_of_points(arr, 2);

    Point newloc = rotate3D(ptH, Navg, ptC.subtract(Navg), M_PI);
    hf->H->move(newloc);

    Atom* was_bonded = hf->H->get_bond_by_idx(0)->atom2;
    hf->H->unbond(was_bonded);
    float rN1 = newloc.get_3d_distance(ptN1), rN2 = newloc.get_3d_distance(ptN2);
    Atom* new_bonded = (rN1 > rN2) ? hf->N2 : hf->N1;
    hf->H->bond_to(new_bonded, 1);
    strcpy(hf->H->name+1, new_bonded->name+1);

    #if _DBG_HISFLIP
    cout << hf->H->name << " moved from " << ptH << " to " << newloc << endl;
    #endif
}

float Molecule::get_springy_bond_satisfaction()
{
    if (!springy_bonds) return 0;

    float retval = 0;
    int i;
    for (i=0; i<springy_bondct; i++)
    {
        if (!springy_bonds[i].atom1 || !springy_bonds[i].atom2 || !springy_bonds[i].optimal_radius) continue;
        float r = springy_bonds[i].atom1->loc.get_3d_distance(springy_bonds[i].atom2->loc);
        r /= springy_bonds[i].optimal_radius;
        if (r < 1) retval += r;
        else retval += 1.0/(r*r);
    }

    return retval;
}

void Molecule::allocate_mandatory_connections(int mcmax)
{
    delete_mandatory_connections();
    mandatory_connection = new Molecule*[mcmax+4];
    mandatory_connection[0] = nullptr;
    last_mc_binding = new float[mcmax+4];
    int i;
    for (i=0; i<mcmax; i++) last_mc_binding[i] = 0;
}

void Molecule::add_mandatory_connection(Molecule* addmol)
{
    int i;
    for (i=0; mandatory_connection[i]; i++)                 // get count.
    {
        if (mandatory_connection[i] == addmol)
        {
            #if _dbg_mand_conn
            cout << "Already have mandatory connection " << addmol->name;
            if (last_mc_binding) cout << " last binding " << last_mc_binding[i];
            cout << endl;
            #endif
            return;      // already have it.
        }
    }
    mandatory_connection[i] = addmol;
    mandatory_connection[i+1] = nullptr;
    if (last_mc_binding) last_mc_binding[i] = 0;
    #if _dbg_mand_conn
    cout << "Add mandatory connection " << addmol->name << endl;
    #endif
}

void Molecule::zero_mandatory_connection_cache()
{
    if (!last_mc_binding) return;
    if (!mandatory_connection) return;
    int i;
    for (i=0; mandatory_connection[i]; i++) last_mc_binding[i] = 0;
}

void Molecule::remove_mandatory_connection(Molecule* rmvmol)
{
    if (!last_mc_binding) return;
    if (!mandatory_connection) return;
    int i, j;
    for (i=0; mandatory_connection[i]; i++)
    {
        if (mandatory_connection[i] == rmvmol)
        {
            for (j=i; mandatory_connection[j]; j++)
            {
                mandatory_connection[j] = mandatory_connection[j+1];
                last_mc_binding[j] = last_mc_binding[j+1];
            }
            #if _dbg_mand_conn
            cout << "Remove mandatory connection " << rmvmol->name << endl;
            #endif
            return;
        }
    }
}

void Molecule::delete_mandatory_connections()
{
    if (last_mc_binding) delete last_mc_binding;
    last_mc_binding = nullptr;
    if (mandatory_connection) delete mandatory_connection;
    mandatory_connection = nullptr;
    #if _dbg_mand_conn
    cout << "No more mandatory connections." << endl;
    #endif
}

Interaction Molecule::intermol_bind_for_multimol_dock(Molecule* om, bool is_ac)
{
    float lbias = 1.0 + (sgn(is_residue()) == sgn(om->is_residue()) ? 0 : dock_ligand_bias);
    Interaction rawbind = get_intermol_binding(om, !is_ac, true);
    Interaction lbind = rawbind * lbias;
    lbind.attractive += get_intermol_contact_area(om, true) * cavity_stuffing;
    lbind.clash *= iteration_additional_clash_coefficient;
    lbind.clash += (get_internal_clashes() + total_eclipses()) * iteration_internal_clash_coefficient;

    int i;

    #if _dbg_415
    if (lbind.clash > clash_limit_per_aa)
    {
        int resno = is_residue(), rno1 = om->is_residue();
        if ((resno == 104 && rno1 == 158) || (resno == 158 && rno1 == 104))
        {
            cout << "\033[A\033[K" << name << " ! " << om->name << " " << lbind.clash << endl << endl << flush;
        }
    }
    #endif

    if (mandatory_connection && rawbind.summed() < 0)                   // Only enforce pullaway if mols are not clashing.
    {
        if (!last_mc_binding)
        {
            for (i=0; mandatory_connection[i]; i++);            // Get count.
            last_mc_binding = new float[i+4];
            for (i=0; mandatory_connection[i]; i++) last_mc_binding[i] = -Avogadro;
        }
        for (i=0; mandatory_connection[i]; i++)
        {
            if (mandatory_connection[i] == om)
            {
                if (lbind.summed() > last_mc_binding[i] && lbind.summed() > mandatory_coordination_threshold)
                {
                    return -1e9;
                }
                else last_mc_binding[i] = lbind.summed();
            }
        }
    }

    #if limit_interactions_by_hydrophilicity
    else
    {
        float apol = fabs(this->hydrophilicity()) + fabs(this->get_charge());
        float bpol = fabs(om->hydrophilicity()) + fabs(om->get_charge());

        float factor = 1.0 / fmax(1.0, fabs(apol-bpol));
        lbind.attractive *= factor;
    }
    #endif

    return lbind;
}

Interaction Molecule::cfmol_multibind(Molecule* a, Molecule** nearby)
{
    if (a->is_residue() && ((AminoAcid*)a)->conditionally_basic()) ((AminoAcid*)a)->set_conditional_basicity(nearby);

    Interaction result = -a->total_eclipses();
    if (a->is_residue()) result += reinterpret_cast<AminoAcid*>(a)->initial_eclipses;

    int i, j;
    if (a->mclashables)
    {
        for (j=0; a->mclashables[j]; j++)
        {
            Interaction f = a->intermol_bind_for_multimol_dock(a->mclashables[j], false);
            result += f;
        }
    }
    else if (nearby) for (j=0; nearby[j]; j++)
    {
        Interaction f = a->intermol_bind_for_multimol_dock(nearby[j], false);
        result += f;
    }

    if (cfmols_have_metals && nearby) for (i=0; nearby[i]; i++)
    {
        float ichg;
        if (nearby[i]->coordmtl && (ichg = nearby[i]->get_charge()))
        {
            for (j=0; nearby[j]; j++)
            {
                if (j==i) continue;
                float jchg = nearby[j]->get_charge();
                if (!jchg) continue;
                if (sgn(ichg) == sgn(jchg)) result.clash += 60.0 / pow(fmax(1, nearby[i]->distance_to(nearby[j]) / 2.0), 2);
            }
        }
    }

    return result;
}

void Molecule::conform_molecules(Molecule** mm, Molecule** bkg, int iters, void (*cb)(int, Molecule**),
    void (*progress)(float))
{
    int m, n;

    if (!mm) m=0;
    else for (m=0; mm[m]; m++);         // Get count.

    if (!bkg) n=0;
    else for (n=0; bkg[n]; n++);        // Get count.

    Molecule* all[m+n+8];
    int i, j, l=0;

    for (i=0; i<m; i++)
    {
        bool duplicate = false;
        for (j=0; j<l; j++)
        {
            if (all[j] == mm[i]) duplicate = true;
        }
        if (duplicate) continue;

        all[l++] = mm[i];
        mm[i]->movability = static_cast<MovabilityType>(static_cast<int>(mm[i]->movability & !MOV_BKGRND));
    }

    for (i=0; i<n; i++)
    {
        bool duplicate = false;
        for (j=0; j<l; j++)
        {
            if (all[j] == bkg[i]) duplicate = true;
        }
        if (duplicate) continue;

        all[l++] = bkg[i];
        bkg[i]->movability = static_cast<MovabilityType>(static_cast<int>(bkg[i]->movability | MOV_BKGRND));
    }

    all[l] = nullptr;

    conform_molecules(all, iters, cb);

    for (i=0; i<n; i++)
    {
        bkg[i]->movability = static_cast<MovabilityType>(static_cast<int>(bkg[i]->movability & !MOV_BKGRND));
    }
}

void Molecule::quick_conform(Molecule **bkg, int iters)
{
    Molecule* mols[2];
    mols[0] = this;
    mols[1] = nullptr;
    conform_molecules(mols, bkg, iters);
}

void Molecule::conform_atom_to_location(const char* an, Point t, int iters, float od)
{
    if (!atoms) return;
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (!strcmp(atoms[i]->name, an))
        {
            conform_atom_to_location(i, t, iters, od);
            return;
        }
    }
}

void Molecule::conform_atom_to_location(Atom *a, Atom *target, int iters, float od)
{
    if (!(movability & MOV_CAN_FLEX)) return;
    if (!a) return;

    int iter, j, l, circdiv = 144;
    Bond** b = get_rotatable_bonds();
    if (!b) return;

    float oc = get_internal_clashes();

    Pose best(this);
    float bestr = Avogadro;
    for (iter = 0; iter < iters; iter++)
    {
        float r;
        for (j=0; b[j]; j++)
        {
            if (!b[j]->can_rotate)
            {
                if (b[j]->can_flip) circdiv = 2;
                else continue;
            }
            float step = M_PI*2.0/circdiv;
            for (l=0; l<=circdiv; l++)
            {
                b[j]->rotate(step, false, true);
                r = a->distance_to(target);
                if (od) r = fabs(r-od);

                #if _dbg_atom_pointing
                cout << " " << bestr;
                #endif

                float c = get_internal_clashes();
                if (r < bestr && c < oc+2.5*clash_limit_per_aa)
                {
                    bestr = r;
                    best.copy_state(this);
                }
            }
            best.restore_state(this);

            #if _dbg_atom_pointing
            cout << endl << name << "." << iter << ": " << circdiv << "|" << bestr << endl;
            #endif
        }
    }
    best.restore_state(this);
}

#define _dbg_atom_pointing 0

void Molecule::conform_atom_to_location(int i, Point t, int iters, float od)
{
    if (!(movability & MOV_CAN_FLEX)) return;

    int iter, j, l, circdiv = 144;
    Bond** b = get_rotatable_bonds();
    if (!b) return;
    Atom* a = atoms[i];
    if (!a) return;

    float oc = get_internal_clashes();

    Pose best(this);
    float bestr = Avogadro;
    for (iter = 0; iter < iters; iter++)
    {
        float r;
        for (j=0; b[j]; j++)
        {
            if (!b[j]->can_rotate)
            {
                if (b[j]->can_flip) circdiv = 2;
                else continue;
            }
            float step = M_PI*2.0/circdiv;
            for (l=0; l<=circdiv; l++)
            {
                b[j]->rotate(step, false, true);
                r = a->loc.get_3d_distance(t);
                if (od) r = fabs(r-od);

                #if _dbg_atom_pointing
                cout << " " << bestr;
                #endif

                float c = get_internal_clashes();
                if (r < bestr && c < oc+2.5*clash_limit_per_aa)
                {
                    bestr = r;
                    best.copy_state(this);
                }
            }
            best.restore_state(this);

            #if _dbg_atom_pointing
            cout << endl << name << "." << iter << ": " << circdiv << "|" << bestr << endl;
            #endif
        }
    }
    best.restore_state(this);
}

Interaction Molecule::total_intermol_binding(Molecule** l)
{
    int i;
    Interaction f;

    for (i=0; l[i]; i++)
    {
        f += l[i]->get_intermol_binding(l);
    }

    return f;
}

void Molecule::conform_molecules(Molecule** mm, int iters, void (*cb)(int, Molecule**), void (*progress)(float), int mi)
{
    if (!mm) return;
    int i, imer, j, l, n, iter;
    
    int mmmono = 0;
    for (i=0; mm[i]; i++) if (mm[i]->nmonomers) mmmono += mm[i]->nmonomers;
    if (mmmono)
    {
        Molecule* lmm[i+mmmono*2+4];
        j=0;
        for (i=0; mm[i]; i++)
        {
            if (mm[i]->nmonomers)
            {
                for (l=0; l<mm[i]->nmonomers; l++)
                {
                    if (mm[i]->monomers[l]) lmm[j++] = mm[i]->monomers[l];
                }
            }
            else lmm[j++] = mm[i];
        }
        lmm[j] = nullptr;
        return conform_molecules(lmm, iters, cb, progress, mi);             // RECURSION!
    }

    end_iterations = false;
    minimum_searching_aniso = 0.5;

    #if _dbg_improvements_only_rule
    check_mols = mm;
    #endif

    for (n=0; mm[n]; n++);
    int nmm = n;
    Pose absolute_best[nmm+4];
    Interaction abe = mm[0]->get_intermol_binding(mm);          // Absolute Best Energy.
    Interaction abtest;
    int abc;                                                    // Absolute Best Counter.
    char triedchange[1024];
    strcpy(triedchange, "");

    for (abc=0; mm[abc]; abc++) absolute_best[abc].copy_state(mm[abc]);

    #define test_and_update_absolute_best_poses abtest = mm[0]->get_intermol_binding(mm); \
                if (abtest.improved(abe)) \
                { \
                    if (audit) fprintf(audit, "Iter %d accepted %s %s from %g (%g/%g) to %g (%g/%g).\n", iter, a->get_name(), triedchange, abe.summed(), abe.attractive, abe.clash, abtest.summed(), abtest.attractive, abtest.clash); \
                    for (abc=0; mm[abc]; abc++) absolute_best[abc].copy_state(mm[abc]); \
                    abe = abtest; \
                }
    
    for (i=0; mm[i]; i++) if (mm[i]->coordmtl) { cfmols_have_metals = true; break;}

    for (iter=0; iter<iters; iter++)
    {
        if (frand(0,1) < best_pose_reset_frequency) for (abc=0; mm[abc]; abc++) absolute_best[abc].restore_state(mm[abc]);
        for (i=0; mm[i]; i++) mm[i]->lastbind = 0;

        Molecule* nearby[2048];
        memset(nearby, 0, sizeof(Molecule*) * 2048);
        bool do_full_rotation = _allow_fullrot && ((iter % _fullrot_every) == 0);

        for (i=0; mm[i]; i++)
        {
            Molecule* a = mm[i];
            if (a->movability & MOV_PINNED) continue;
            bool flipped_rings = false;
            int ares = a->is_residue();

            if (a->movability & MOV_BKGRND) continue;

            #if _dbg_asunder_atoms
            if (!a->check_Greek_continuity()) throw 0xbadc0de;
            #endif

            Point aloc = a->get_barycenter();

            Interaction benerg = 0;
            if (!ares)
            {
                l = 0;
                for (j=0; mm[j]; j++)
                {
                    if (j==i) continue;
                    Molecule* b = mm[j];

                    #if mclashables_as_residue_nearbys
                    if (ares && b->is_residue()) continue;
                    #endif

                    Point bloc = b->get_barycenter();

                    Atom* na = a->get_nearest_atom(bloc);
                    if (!na) continue;
                    float r = na->distance_to(b->get_nearest_atom(aloc));
                    if (r > _INTERA_R_CUTOFF+2.5) continue;
                    nearby[l++] = b;
                }
                #if mclashables_as_residue_nearbys
                if (ares)
                {
                    for (j=0; a->mclashables[j]; j++) nearby[l++] = a->mclashables[j];
                }
                #endif
                nearby[l] = 0;
            }
            benerg = cfmol_multibind(a, nearby);

            #if _dbg_fitness_plummet
            if (!i) cout << "# mol " << a->name << " iter " << iter << ": initial " << -benerg << " ";
            #endif

            #if _dbg_asunder_atoms
            if (!a->check_Greek_continuity()) throw 0xbadc0de;
            #endif

            if (!a->iterbegan) a->iterbegan = new Pose(a);
            a->iterbegan->copy_state(a);
            if (!iter) a->iters_without_change = 0;

            Interaction tryenerg;
            Pose pib;
            // pib.copy_state(a);

            /**** Linear Motion ****/
            #if allow_linear_motion
            if ((a->movability & MOV_CAN_RECEN) && !(a->movability & MOV_FORBIDDEN))
            {
                Point motion(a->lmx, a->lmy, a->lmz);
                if (motion.magnitude() > speed_limit) motion.scale(speed_limit);

                benerg = cfmol_multibind(a, nearby);
                if (motion.magnitude() > 0.01*speed_limit)
                {
                    pib.copy_state(a);
                    motion.scale(motion.magnitude()/lmsteps);
                    for (l=0; l<lmsteps; l++)
                    {
                        #if _dbg_improvements_only_rule
                        excuse_deterioration = true;
                        #endif

                        if (audit) sprintf(triedchange, "lmpush/lmpull linear motion [%f,%f,%f]", motion.x, motion.y, motion.z);

                        a->move(motion);
                        tryenerg = cfmol_multibind(a, nearby);

                        if (tryenerg.improved(benerg))
                        {
                            benerg = tryenerg;
                            pib.copy_state(a);
                            test_and_update_absolute_best_poses
                        }

                        a->enforce_stays(1.0/lmsteps);
                        tryenerg = cfmol_multibind(a, nearby);

                        if (tryenerg.improved(benerg))
                        {
                            benerg = tryenerg;
                            pib.copy_state(a);
                            test_and_update_absolute_best_poses
                        }
                    }
                    pib.restore_state(a);
                    benerg = cfmol_multibind(a, nearby);
                }
                pib.copy_state(a);

                #if _dbg_linear_motion
                if (!ares) cout << iter << "! ";
                #endif

                a->lmx = a->lmy = a->lmz = 0;

                int xyz;
                for (xyz=0; xyz<3; xyz++)
                {
                    motion.scale(0);

                    switch(xyz)
                    {
                        case 0: motion.x = frand(-speed_limit,speed_limit); break;
                        case 1: motion.y = frand(-speed_limit,speed_limit); break;
                        case 2: motion.z = frand(-speed_limit,speed_limit); break;
                        default:
                        ;
                    }

                    #if _dbg_improvements_only_rule
                    excuse_deterioration = true;
                    #endif

                    a->move(motion);
                    if (audit) sprintf(triedchange, "stochastic linear motion [%f,%f,%f]", motion.x, motion.y, motion.z);

                    tryenerg = cfmol_multibind(a, nearby);
                    

                    #if _dbg_fitness_plummet
                    if (!i) cout << "(linear motion try " << -tryenerg << ") ";
                    #endif

                    if (tryenerg.improved(benerg))
                    {
                        benerg = tryenerg;
                        pib.copy_state(a);
                        test_and_update_absolute_best_poses

                        #if _dbg_linear_motion
                        if (!ares) cout << ">";
                        #endif
                    }
                    else
                    {
                        switch(xyz)
                        {
                            case 0: motion.x *= -2; break;
                            case 1: motion.y *= -2; break;
                            case 2: motion.z *= -2; break;
                            default:
                            ;
                        }

                        #if _dbg_improvements_only_rule
                        excuse_deterioration = true;
                        #endif

                        a->move(motion);
                        if (audit) sprintf(triedchange, "reversed linear motion [%f,%f,%f]", motion.x, motion.y, motion.z);

                        tryenerg = cfmol_multibind(a, nearby);
                        

                        if (tryenerg.improved(benerg))
                        {
                            benerg = tryenerg;
                            pib.copy_state(a);
                            test_and_update_absolute_best_poses

                            #if _dbg_linear_motion
                            if (!ares) cout << "<";
                            #endif
                        }
                        else
                        {
                            pib.restore_state(a);

                            #if _dbg_linear_motion
                            if (!ares) cout << "-" << flush;
                            #endif
                        }
                    }
                }
            }       // If can recenter.
            #if _dbg_linear_motion
            else
            {
                if (!ares) cout << iter << "* ";
            }
            #endif
            #endif

            #if _dbg_fitness_plummet
            if (!i) cout << "linear " << -benerg << " ";
            #endif

            #if _dbg_asunder_atoms
            if (!a->check_Greek_continuity()) throw 0xbadc0de;
            #endif

            /**** Histidine flip ****/
            if (mm[i]->hisflips)
            {
                for (l=0; mm[i]->hisflips[l]; l++)
                {
                    #if _DBG_HISFLIP
                    cout << i << ": Flipping " << mm[i]->name << endl;
                    #endif
                    mm[i]->do_histidine_flip(mm[i]->hisflips[l]);
                    if (audit) sprintf(triedchange, "histidine flip");

                    tryenerg = cfmol_multibind(a, nearby);

                    if (tryenerg.improved(benerg) || tryenerg.attractive > benerg.attractive)
                    {
                        benerg = tryenerg;
                        pib.copy_state(a);
                        test_and_update_absolute_best_poses
                    }
                    else
                    {
                        mm[i]->do_histidine_flip(mm[i]->hisflips[l]);
                    }
                }
            }
            /**** End histidine flip ****/
        
            #if _dbg_asunder_atoms
            if (!a->check_Greek_continuity()) throw 0xbadc0de;
            #endif

            #if allow_axial_tumble
            if ((!iter % 3) && (a->movability & MOV_CAN_AXIAL) && !(a->movability & MOV_FORBIDDEN))
            {
                pib.copy_state(a);
                Point ptrnd(frand(-1,1), frand(-1,1), frand(-1,1));
                if (frand(0,1) < 0.4 && a->best_intera && a->best_other_intera)
                {
                    ptrnd = a->best_other_intera->loc.subtract(a->best_intera->loc);
                }
                if (ptrnd.magnitude())
                {
                    LocatedVector axis = (Vector)ptrnd;
                    axis.origin = a->get_barycenter(true);
                    float theta;

                    theta = frand(-0.5, 0.5)*fiftyseventh*min(10, iter/2);

                    #if _dbg_improvements_only_rule
                    excuse_deterioration = true;
                    #endif

                    a->rotate(&axis, theta);
                    if (audit) sprintf(triedchange, "stochastic rotation %f deg.", theta*fiftyseven);
                    a->enforce_stays(multimol_stays_enforcement);
                    tryenerg = cfmol_multibind(a, nearby);

                    if (tryenerg.improved(benerg))
                    {
                        benerg = tryenerg;
                        pib.copy_state(a);
                        test_and_update_absolute_best_poses
                    }
                    else
                    {
                        #if multimol_stays_allow_revert_worsening
                        pib.restore_state(a);
                        #endif
                    }
                }
            }       // If can axial rotate.
            #endif

            #if _dbg_fitness_plummet
            if (!i) cout << "axial " << -benerg << " ";
            #endif

            #if _dbg_asunder_atoms
            if (!a->check_Greek_continuity()) throw 0xbadc0de;
            #endif

            #if allow_bond_rots
            #if _dbg_mol_flexion
            bool is_flexion_dbg_mol = (ares == 107);
            if (is_flexion_dbg_mol) cout << a->name << " movability " << hex << a->movability << dec << endl << flush;
            #endif
            if (((a->movability & MOV_CAN_FLEX) && !(a->movability & MOV_FORBIDDEN) && a->movability != MOV_PINNED)
                && (!ares || frand(0,1) < sidechain_flexion_frequency)
                )
            {
                pib.copy_state(a);
                #if _dbg_asunder_atoms
                if (!a->check_Greek_continuity()) throw 0xbadc0de;
                #endif

                float self_clash = max(1.25*a->base_internal_clashes, clash_limit_per_aa);
                Bond** bb = a->get_rotatable_bonds(true);
                if (bb)
                {
                    int q, rang=0, qiter;
                    float bbodds = 1;
                    for (qiter=0; qiter<flexion_sub_iterations; qiter++) for (q=0; bb[q]; q++)
                    {
                        if (!bb[q]->atom1 || !bb[q]->atom2) continue;         // Sanity check, otherwise we're sure to get random foolish segfaults.
                        if (bb[q]->atom1->get_Greek() > bb[q]->atom2->get_Greek()) bb[q] = bb[q]->get_reversed();
                        if (!bb[q]->count_moves_with_atom2()) continue;
                        if (bb[q]->atom1->is_backbone && strcmp(bb[q]->atom1->name, "CA")) continue;
                        if (bb[q]->atom2->is_backbone) continue;
                        bbodds *= 0.85;
                        if (frand(0,1) > bbodds) continue;
                        float theta;
                        int heavy_atoms = bb[q]->count_heavy_moves_with_atom2();
                        if (heavy_atoms && (!(a->movability & MOV_CAN_FLEX) || (a->movability & MOV_FORBIDDEN))) continue;

                        #if _dbg_mol_flexion
                        bool is_flexion_dbg_mol_bond = is_flexion_dbg_mol & !strcmp(bb[q]->atom2->name, "OG");
                        #endif

                        benerg.clash += a->total_eclipses();
                        if (do_full_rotation && bb[q]->can_rotate
                            #if fullrot_flex_first_subiter_only
                            && !qiter
                            #endif

                            #if fullrot_forbid_residues
                            && !ares
                            #elif fullrot_flex_residues_only
                            && ares
                            #endif

                            #if fullrot_flex_unfavorable_energy_only
                            && benerg.summed() >= 0
                            #endif
                            )
                        {
                            float best_theta = 0;
                            Pose prior_state;
                            prior_state.copy_state(a);
                            for (theta=_fullrot_steprad; theta < M_PI*2; theta += _fullrot_steprad)
                            {
                                prior_state.restore_state(a);
                                bb[q]->rotate(theta, false);
                                a->enforce_stays(multimol_stays_enforcement);
                                tryenerg = cfmol_multibind(a, nearby);
                                tryenerg.clash += a->total_eclipses();
                                if (audit) sprintf(triedchange, "fullrot flexion %s-%s %f deg.", bb[q]->atom1->name, bb[q]->atom2->name, theta*fiftyseven);

                                #if _dbg_mol_flexion
                                if (is_flexion_dbg_mol_bond) cout << (theta*fiftyseven) << "deg: " << -tryenerg << endl;
                                #endif

                                if (tryenerg.improved(benerg) && a->get_internal_clashes() <= self_clash)
                                {
                                    benerg = tryenerg;
                                    best_theta = theta;
                                    test_and_update_absolute_best_poses
                                }
                            }
                            if (best_theta)
                            {
                                prior_state.restore_state(a);
                                bb[q]->rotate(best_theta, false);
                                a->been_flexed = true;

                                #if _dbg_mol_flexion
                                if (is_flexion_dbg_mol_bond) cout << "Rotating to " << (best_theta*fiftyseven) << "deg." << endl << endl;
                                #endif
                            }
                        }
                        else
                        {
                            theta = frand(-flexion_maxangle, flexion_maxangle);

                            if (!bb[q]->can_rotate)
                            {
                                bb[q]->compute_flip_capability();
                                if (!bb[q]->flip_angle) bb[q]->flip_angle = M_PI;
                            }

                            Ring* isra = bb[q]->atom1->in_same_ring_as(bb[q]->atom2);
                            if (isra)
                            {
                                if (rang) continue;
                                isra->flip_atom(bb[q]->atom1);
                                rang++;
                                flipped_rings = true;
                                if (audit) sprintf(triedchange, "ring flip of %s", bb[q]->atom1->name);
                            }
                            else
                            {
                                bb[q]->rotate(theta, false);
                                if (audit) sprintf(triedchange, "stochastic flexion %s-%s %f deg.", bb[q]->atom1->name, bb[q]->atom2->name, theta*fiftyseven);
                            }

                            a->enforce_stays(multimol_stays_enforcement);
                            tryenerg = cfmol_multibind(a, nearby);
                            tryenerg.clash += a->total_eclipses();

                            #if _dbg_mol_flexion
                            if (is_flexion_dbg_mol_bond) cout << "Trying " << (theta*fiftyseven) << "deg rotation...";
                            #endif

                            if (tryenerg.improved(benerg) && a->get_internal_clashes() <= self_clash)
                            {
                                // Leaving this in case the "nearbys" feature misses any more clashable residues.
                                // If it does, adjust the constants on the cosine in AminoAcid::can_reach().
                                /*if (a->is_residue() == 109)
                                {
                                    cout << endl << "Nearby: ";
                                    int nb;
                                    for (nb=0; nearby[nb]; nb++) cout << nearby[nb]->is_residue() << " ";
                                    cout << endl << "Accepted: " << -tryenerg.attractive << "+" << tryenerg.clash << endl;
                                }*/

                                benerg = tryenerg;
                                pib.copy_state(a);
                                a->been_flexed = true;
                                test_and_update_absolute_best_poses

                                #if _dbg_mol_flexion
                                if (is_flexion_dbg_mol_bond) cout << " energy now " << -tryenerg << ", keeping." << endl << endl;
                                #endif
                            }
                            else
                            {
                                #if multimol_stays_allow_revert_worsening
                                pib.restore_state(a);

                                #if _dbg_mol_flexion
                                if (is_flexion_dbg_mol_bond) cout << " energy from " << -benerg << " to " << -tryenerg << ", reverting." << endl << endl;
                                #endif
                                #endif
                            }
                        }
                    }
                }       // Rotatable bonds.

                #if _dbg_fitness_plummet
                if (!i) cout << "flexion " << -benerg << " ";
                #endif

                mm[i]->lastbind = benerg.summed();
            }   // if MOV_CAN_FLEX
            #endif

            #if _dbg_fitness_plummet
            if (!i) cout << "final " << -benerg << " " << endl << flush;
            #endif

            if (!ares && flipped_rings) a->evolve_structure(100);

            if (!i && !ares)
            {
                float ttl_atom_mtn = a->iterbegan->total_atom_motions() / a->get_heavy_atom_count();
                if (ttl_atom_mtn < iter_lostreturns_threshold) a->iters_without_change++;
                else a->iters_without_change = 0;
                if (iter >= mi && a->iters_without_change >= max_iters_without_ligand_change && mm[0]->lastbind < -5) iter = iters;
            }
            if (a->iterbegan) a->iterbegan->deallocate();
            pib.deallocate();

            if (!(i%8) && progress)
            {
                float f = (float)i / n;
                float fiter = (f + iter) / iters * 100;
                progress(fiter);
            }
        }       // for i

        #if allow_iter_cb
        if (cb) cb(iter+1, mm);
        #endif

        minimum_searching_aniso *= 0.99;

        #if _dbg_fitness_plummet
        if (!i) cout << endl << flush;
        #endif

        if (_kbhit())
        {
            if (_getch() == 'q')
            {
                cout << endl << "Keyboard quit command detected; stopping dock..." << endl;
                end_program = true;
                return;
            }
        }
        if (end_program || end_iterations) return;
    }       // for iter

    for (abc=0; abc<nmm; abc++)
    {
        if (mm[abc]->is_residue()) absolute_best[abc].restore_state_relative(mm[abc], "CA");
        else absolute_best[abc].restore_state(mm[abc]);
    }

    #if _dbg_linear_motion
    cout << endl;
    #endif

    minimum_searching_aniso = 0;
}

float Molecule::get_molecular_wt()
{
    if (!atoms) return 0;

    int i;
    float result=0;
    for (i=0; i<atcount; i++)
    {
        if (!atoms[i]) break;
        result += atoms[i]->get_atomic_weight();
    }

    return result;
}

int Molecule::get_heavy_atom_count() const
{
    if (!atoms) return 0;

    int i, result=0;
    for (i=0; i<atcount; i++)
    {
        if (!atoms[i]) break;
        if (atoms[i]->Z > 1) result++;
    }

    return result;
}

#define dbg_optimal_contact 0
Vector Molecule::motion_to_optimal_contact(Molecule* l)
{
    Pose p;
    p.copy_state(this);

    MovabilityType m = movability, lm = l->movability;
    movability = l->movability = MOV_FLEXONLY;

    Vector total_motion(0,0,0);
    Vector incremental_motion;

    Interaction energy = this->get_intermol_binding(l, true, true);

    incremental_motion = this->get_barycenter().subtract(l->get_barycenter());
    incremental_motion.r = -sgn(energy.summed());

    int i;

    #if dbg_optimal_contact
    cout << "Energy was " << energy << endl;
    #endif

    for (i=0; i<50; i++)
    {
        this->move(incremental_motion, true);
        Interaction new_energy = this->get_intermol_binding(l, true, true);

        if (new_energy.improved(energy))
        {
            energy = new_energy;
            total_motion = total_motion.add(incremental_motion);

            #if dbg_optimal_contact
            cout << "Energy has improved to " << energy << "; keeping." << endl;
            #endif
        }
        else
        {
            incremental_motion.r *= -1;
            this->move(incremental_motion, true);
            incremental_motion.r *= 0.8;

            #if dbg_optimal_contact
            cout << "Energy tried " << energy << endl;
            #endif
        }

        if (fabs(incremental_motion.r) < 0.01) break;
    }

    p.restore_state(this);
    movability = m;
    l->movability = lm;

    return total_motion;
}

void Molecule::wipe_vdw_surface()
{
    if (!vdw_vertex_count) return;
    delete[] vdw_surface;
    vdw_vertex_count = 0;
    vdw_surface = nullptr;
}

const Point* Molecule::obtain_vdW_surface(float d)
{
    if (!atcount || !atoms) return nullptr;
    if (vdw_surface && vdw_vertex_count > 0) return vdw_surface;

    int maxpoints = atcount * d * d / 2 + 16384;
    if (!vdw_surface)
    {
        vdw_surface = new Point[maxpoints];
        vdw_vertex_atom = new Atom*[maxpoints];
    }

    float halfstep = M_PI / d;
    float step = halfstep * 2;

    int i, ivdW = 0;
    Vector v;
    for (i=0; i<atcount && atoms[i]; i++)
    {
        Point aloc = atoms[i]->loc;
        v.r = atoms[i]->vdW_radius;
        float ystep = step / v.r / v.r;
        for (v.theta = -square; v.theta <= square; v.theta += step)
        {
            float xstep = step / v.r / fmax(cos(v.theta), 0.000001);
            float end = M_PI*2-xstep/2;
            for (v.phi = 0; v.phi < end; v.phi += xstep)
            {
                Point pt = aloc.add(v);
                Atom* na = this->get_nearest_atom(pt);
                if (na != atoms[i] && pt.get_3d_distance(na->loc) < na->vdW_radius) continue;
                if (!pt.x && !pt.y && !pt.z) pt = Point(-0.001, 0.001, -0.001);
                vdw_vertex_atom[ivdW] = atoms[i];
                vdw_surface[ivdW++] = pt;
                if (ivdW >= maxpoints)
                {
                    cout << "Too many vdW surface vertices. Please increase limit in code." << endl;
                    throw 0xbadc0de;
                }
            }
        }
    }
    vdw_vertex_count = ivdW;

    return vdw_surface;
}

int Molecule::get_atom_vdW_vertex_count(Atom* a)
{
    if (!atoms) return 0;
    int i=0, result=0;

    if (!vdw_vertex_count || !vdw_surface || !vdw_vertex_atom) obtain_vdW_surface(vdw_surface_density);

    for (i=0; i<vdw_vertex_count; i++)
    {
        if (vdw_vertex_atom[i] == a) result++;
    }

    return result;
}

Atom* numbered[10];
bool ring_warned = false;

bool Molecule::from_smiles(char const * smilesstr, bool use_parser)
{
    if (strchr(smilesstr, '{')) use_parser = false;	    // {AtomName} is a nonstandard feature and cannot be handled by a third party app.

    // obabel is currently exhibiting a chiral enantiomer bug, so we unfortunately cannot use it for chiral molecules.
    // But the in-house molecule assembly code is failing horribly with rings, so we must use obabel for all cyclic molecules.
    if (strchr(smilesstr, '@') && !strchr(smilesstr, '1') ) use_parser = false;

    if (use_parser)
    {
        // Check if OpenBabel is installed.
        FILE* pf = popen(CMD_CHECK_INSTALLED_3P_SMILES_PARSER, "r");
        if (pf)
        {
            char buffer[1024];
            char* got = fgets(buffer, 1022, pf);
            if (strlen(buffer))			// TODO: Change this to employ a regex.
            {
                fclose(pf);
                std::string sdfdat = "";

                // Temporarily reuse buffer as the obabel command.
                sprintf(buffer, CMD_CALL_3P_SMILES_PARSER, smilesstr);
                pf = popen(buffer, "r");

                // Resume using buffer with fgets().
                int lno = 0;
                while (buffer[0] != '$' && !feof(pf))
                {
                    got = fgets(buffer, 1022, pf);
                    if (!got) break;
                    lno++;
                    sdfdat += (std::string)buffer;

                    if (lno == 2) sdfgen_aboutline = buffer;
                }
                fclose(pf);

                int result = from_sdf(sdfdat.c_str());
                return (result > 0);
            }
            buffer[0] = 0;
            fclose(pf);
        }
    }

    if (strchr(smilesstr, '!')) ring_warned = true;

    smlen = strlen(smilesstr);
    paren = new SMILES_Parenthetical[smlen];
    spnum = 0;

    int i;
    for (i=0; i<10; i++) numbered[i] = 0;

    bool retval = from_smiles(smilesstr, nullptr);

    for (i=0; i<spnum; i++)
    {
        retval &= from_smiles(paren[i].smilesstr, paren[i].startsfrom);
    }

    // hydrogenate(true);
    float anomaly = correct_structure();
    // cout << "# Structural anomaly = " << anomaly << endl;
    hydrogenate(false);
    identify_conjugations();
    identify_rings();
    identify_cages();

    return retval;
}

bool Molecule::from_smiles(char const * smilesstr, Atom* ipreva)
{
    if (!smilesstr) return true;
    if (!strlen(smilesstr)) return true;

    Atom* stack[256];
    Atom* sequence[65536];
    bool seqarom[65536];
    int sqidx[10];
    int sp = 0;
    bool bracket=false, prevarom=false;
    Atom* bracketed=0;

    immobile = false;

    int i, j=1, k=0, l, atno=get_atom_count()+1;

    Atom* preva = ipreva;
    float card = ipreva?1:0;
    int len = strlen(smilesstr);
    int lastEZ = 0;

    int numdb = 0, dbi = 0;
    for (i=0; i<len; i++) if (smilesstr[i] == '=') numdb++;

    int EZgiven[numdb+4];
    Atom* EZatom0[numdb+4];
    Atom* EZatom1[numdb+4];

    smlen = strlen(smilesstr);
    paren = new SMILES_Parenthetical[smlen];

    for (i=0; i<len; i++)
    {
        if (smilesstr[i] == '!') continue;

        if (smilesstr[i] == '.')
        {
            card = 0;
            continue;
        }

        if (smilesstr[i] == '-')
        {
            card = 1;
            continue;
        }

        if (smilesstr[i] == ':')
        {
            card = 1.5;
            continue;
        }

        if (smilesstr[i] == '=')
        {
            card = 2;
            EZgiven[dbi] = 0;
            EZatom0[dbi] = preva;
            continue;
        }

        if (smilesstr[i] == '#')
        {
            card = 3;
            continue;
        }

        if (smilesstr[i] == '$')
        {
            card = 4;
            continue;
        }

        if (smilesstr[i] == '(')
        {
            // stack[sp++] = preva;
            paren[spnum].startsfrom = preva;
            paren[spnum].smilesstr = new char[smlen-i+4];
            strcpy(paren[spnum].smilesstr, &smilesstr[++i]);
            // while (smilesstr[i] != ')') i++;

            int level = 1;
            while (level)
            {
                i++;
                if (smilesstr[i] == '(') level++;
                if (smilesstr[i] == ')') level--;
                if (!smilesstr[i]) throw 0xbade9c0d;
            }

            /*for (l=0; paren[spnum].smilesstr[l]; l++)
            {	if (paren[spnum].smilesstr[l] == ')')
            	{	paren[spnum].smilesstr[l+1] = nullptr;
            		break;
            	}
            }*/
            spnum++;
            continue;
        }

        if (smilesstr[i] == ')')
        {
            if (!sp) return false;
            if (ipreva) return true;
            // preva = stack[--sp];
            continue;
        }

        // Nonstandard feature.
        if (smilesstr[i] == '{')
        {
            char* aname = new char[15];
            i++;
            j=0;
            while (smilesstr[i] != '}' && j<15)
            {
                aname[j] = smilesstr[i];
                i++;
                j++;
            }
            aname[j] = 0;
            preva->name = aname;
            //i++;
            continue;
        }

        if (smilesstr[i] == '\\'
                ||
                smilesstr[i] == '/'
           )
        {
            if (lastEZ)
            {
                int EZ = (smilesstr[i] == '/') ? 1 : -1;
                preva->EZ_flip = EZgiven[dbi-1] = sgn(EZ) == sgn(lastEZ);
                // cout << preva->name << " EZ flip: " << preva->EZ_flip << endl;
                lastEZ = 0;
            }
            else
            {
                lastEZ = (smilesstr[i] == '/') ? 1 : -1;
            }
            continue;
        }

        if (smilesstr[i] >= '0' && smilesstr[i] <= '9')
        {
            if (!ring_warned)
            {
                cout << "WARNING: Native support of SMILES is still in development. ";
                cout << "Certain molecules containing rings might not render properly." << endl;
                cout << "If possible, it is recommended to install OpenBabel (e.g. sudo apt-get install openbabel) ";
                cout << "so the integration feature can be used and SMILES strings converted seamlessly to ";
                cout << "the corresponding molecular structures." << endl;
                ring_warned = true;
            }

            j = smilesstr[i] - 48;
            if (!numbered[j])
            {
                numbered[j] = preva;
                sqidx[j] = k-1;
                //cout << "+" << endl;;
                continue;
            }
            else
            {
                //cout << " connecting..." << endl;
                // Rotate bonds to close the loop.
                bool allarom = true;
                Atom* aloop[256];
                for (l=sqidx[j]; l<k; l++)
                {
                    aloop[l-sqidx[j]] = sequence[l];
                    if (!seqarom[l]) allarom = false;
                }
                aloop[l-sqidx[j]] = 0;
                int ringsz = l-sqidx[j];

                add_ring(aloop);
                float anomaly = close_loop(aloop, card);

                if (card) preva->bond_to(numbered[j], card);
                card = 1;

                numbered[j] = 0;

                continue;
            }
        }

        if (smilesstr[i] == '[')
        {
            bracket = true;
            continue;
        }

        if (bracket)
        {
            bool aromatic = false;

            if (smilesstr[i] == '@')
            {
                if (bracketed)
                {
                    bracketed->swap_chirality();
                    continue;
                }
                else throw 0xbade9c0d;
            }

            if (	(smilesstr[i] >= 'A' && smilesstr[i] <= 'Z')
                    ||
                    (smilesstr[i] >= 'a' && smilesstr[i] <= 'z')
               )
            {
                char esym[5];
                int ioff=1;

                esym[0] = smilesstr[i];
                esym[1] = smilesstr[i+1];
                esym[2] = 0;

                if (smilesstr[i+1] < 'a'
                        ||
                        smilesstr[i+1] > 'z'
                        ||
                        !Atom::Z_from_esym(esym)
                   )
                {
                    esym[1] = 0;
                    ioff = 0;
                }

                if (esym[0] >= 'a')
                {
                    aromatic = true;
                    esym[0] &= 0x5f;
                    if (prevarom) card=1.5;
                }

                if (Atom::Z_from_esym(esym))
                {
                    char aname[9];
                    if (!bracketed)
                    {
                        sprintf(aname, "%s%d", esym, (atno++%100));
                        bracketed = add_atom(esym, aname, preva, card);
                        if (aromatic) bracketed->aromatize();
                        if (card == 2)
                        {
                            EZatom1[dbi++] = bracketed;
                        }
                        bracketed->dnh = true;
                        i += ioff;
                        ioff=0;
                    }
                    else
                    {
                        i += ioff;
                        ioff=0;
                        int plex = atoi(smilesstr+i+1);
                        if (!plex) plex=1;

                        for (l=0; l<plex; l++)
                        {
                            sprintf(aname, "%s%d", esym, (atno++%100));
                            Atom* a = add_atom(esym, aname, bracketed, 1);
                            a->dnh = true;
                        }
                    }
                }
                else
                {
                    std::string str = smilesstr, estr = esym;
                    cout << "Bad element symbol " << esym << " (" << i << ")." << endl
                         << str.substr(0, i)
                         << "\x1b[4m" << str.substr(i, estr.length())
                         << "\x1b[0m" << str.substr(i+estr.length())
                         << endl;
                    throw 0xbade9c0d;
                }

                continue;
            }

            if (smilesstr[i] == '+')
            {
                if (bracketed)
                {
                    int plex = atoi(smilesstr+i+1);
                    if (!plex) plex=1;
                    bracketed->increment_charge(plex);
                }
                else throw 0xbade9c0d;
            }

            if (smilesstr[i] == '-')
            {
                if (bracketed)
                {
                    int plex = atoi(smilesstr+i+1);
                    if (!plex) plex=1;
                    bracketed->increment_charge(-plex);
                }
                else throw 0xbade9c0d;
            }

            if (smilesstr[i] == ']')
            {
                bracket = false;
                preva = bracketed;
                bracketed = 0;

                seqarom[k] = aromatic;
                sequence[k++] = preva;
                card = 1;
                prevarom = aromatic;

                continue;
            }

            continue;
        }

        char esym[5], aname[5];
        bool aromatic=false;

        sprintf(esym, "%c", smilesstr[i]);
        if (esym[0] >= 'a' && esym[0] <= 'z')
        {
            esym[0] &= 0x5f;
            aromatic = true;
            if (prevarom) card=1.5;
        }

        sprintf(aname, "%c%d", smilesstr[i], (atno++%100));
        // cout << "Bonding new " << esym << " named " << aname << " to " << (preva?preva->name:"nothing") << " cardinality " << card << endl;
        Atom* a = add_atom(esym, aname, card ? preva : 0, card);
        if (aromatic) a->aromatize();
        if (card == 2)
        {
            EZatom1[dbi++] = a;
        }

        seqarom[k] = aromatic;
        sequence[k++] = a;

        preva = a;
        card = 1;
        prevarom = aromatic;
    }
    atoms[atcount]=0;

    identify_conjugations();
    return true;
}


void Molecule::make_coplanar_ring(Atom** ring_members, int ringid)
{
    if (!ring_members) return;
    int i, j, l, ringsz;

    bool allarom = true;
    // cout << "Making coplanar ring of ";
    for (i=0; ring_members[i]; i++)
    {
    	// cout << ring_members[i]->name << " ";
        ringsz = i+1;
        Bond* ab[16];
        ring_members[i]->fetch_bonds(ab);
        bool haspi = false;
        for (j=0; ab[j]; j++)
            if (ab[j]->cardinality > 1 && ab[j]->cardinality < 2) haspi = true;
        if (!haspi) allarom = false;
    }
    // cout << endl;

    if (ringsz<3) return;
    if (ringsz>6) return;

    Vector normal;
    Point ringcen;

	if (!ring_members || !ring_members[0])
	{
		cout << "Notice: empty ring passed to Molecule::make_coplanar_ring()." << endl;
		return;
	}
    Bond* a0b[16];
    ring_members[0]->fetch_bonds(a0b);
    if (!a0b[0]->atom2 || !a0b[1]->atom2)
    {
        cout << "Attempted to form coplanar ring with starting atom bonded to fewer than two other atoms." << endl;
        throw 0xbad12196;					// If you use your imagination, 12196 spells "ring".
    }

    Point A, B, C;
    A = ring_members[0]->loc;
    B = a0b[0]->atom2->loc;
    C = a0b[1]->atom2->loc;

    normal = compute_normal(&A, &B, &C);
    while (!normal.r)
    {
    	B.x = frand(-0.01, 0.01);
    	B.y = frand(-0.01, 0.01);
    	B.z = frand(-0.01, 0.01);
    	C.x = frand(-0.01, 0.01);
    	C.y = frand(-0.01, 0.01);
    	C.z = frand(-0.01, 0.01);
    	normal = compute_normal(&A, &B, &C);
    }

    // TODO: Un-hardcode the bond lengths below.
    ringcen = A.subtract(&B);
    ringcen.scale(polygon_radius(allarom?1.40:1.54, ringsz));
    ringcen = A.add(ringcen);

    for (l=0; l<ringsz; l++)
    {
        if (l)
        {
            Point lpt = rotate3D(&A, &ringcen, &normal, M_PI*2/ringsz*l);
            ring_members[l]->move(&lpt);
        }

        Bond* b2=0;
        int bgeo = ring_members[l]->get_geometry();
    	// cout << bgeo << ", checking " << ring_members[l]->name << "... " << flush;
        for (i=0; i<bgeo; i++)
        {
        	// cout << i << " ";
            b2 = ring_members[l]->get_bond_by_idx(i);
            if (!b2)
            {
            	// cout << "nullptr bond." << endl;
                continue;
            }
            if (!b2->atom2)
            {
            	// cout << "nullptr atom2." << endl;
                b2=0;
                continue;
            }
            for (j=0; j<ringsz; j++)
            {
            	if (ring_members[j] == b2->atom2)
                {
            		// cout << "atom2 is part of the ring." << endl;
                    b2 = 0;
                    break;
                }
            }
            if (b2) break;
        }
        if (b2 && b2->atom2)
        {
        	// cout << "found " << b2->atom2->name << endl;
            Point ptnew = ring_members[l]->loc.subtract(ringcen);
            ptnew.scale(InteratomicForce::covalent_bond_radius(ring_members[l], b2->atom2, b2->cardinality));
            ptnew = ring_members[l]->loc.add(ptnew);
            b2->atom2->move_assembly(&ptnew, ring_members[l]);
        }
    }
}

float Molecule::fsb_lsb_anomaly(Atom* first, Atom* last, float lcard, float bond_length)
{
    // Last-should-be and first-should-be positions.
    Vector lsbv = first->get_next_free_geometry(lcard);
    Vector fsbv = last->get_next_free_geometry(lcard);
    lsbv.r = fsbv.r = bond_length;
    Point  lsb  = first->loc.add(&lsbv);
    Point  fsb  = last->loc.add(&fsbv);

    return first->loc.get_3d_distance(fsb) + last->loc.get_3d_distance(lsb);
}

float Molecule::close_loop(Atom** path, float lcard)
{
    Bond* rotables[65536];

    if (!path) return 0;

    int i, j, k=0;
    Atom* first, *last;

    first = path[0];
    if (!first) return 0;

    for (i=0; path[i]; i++)
    {
        path[i]->doing_ring_closure = true;
        Bond* b[16];
        path[i]->fetch_bonds(b);
        if (!b[0]) continue;
        int geo = path[i]->get_geometry();

        for (j=0; j<geo; j++)
        {
            if (b[j] && b[j]->atom2
                    &&
                    (	b[j]->can_rotate
                        ||
                        b[j]->can_flip
                    )
                    &&
                    strcmp(b[j]->atom2->name, "N")
               ) rotables[k++] = b[j];
        }

        last = path[i];
    }
    rotables[k] = 0;
    if (_DBGCLSLOOP) cout << "Close Loop: found " << k << " rotables." << endl;

    if (last == first) return 0;
    last->mirror_geo = -1;
    int ringsize = k;

    float bond_length = InteratomicForce::covalent_bond_radius(first, last, lcard);

    if (ringsize < 5)
    {
        // TODO: Make equilateral ring, except accommodating any differences in bond length.
        // But know that cyclopentane and cyclopropane are not coplanar, but rather puckered because of steric hindrance.
    }

    int iter;
    float anomaly = fsb_lsb_anomaly(first, last, lcard, bond_length);
    float iclash = get_internal_clashes();
    float oclash = iclash;
    float bondrot[ringsize];

    for (i=0; i<ringsize; i++) bondrot[i] = M_PI/2*randsgn();

    float allowance = 2.5;

    for (iter=0; iter<250+(20*ringsize); iter++)
    {
        for (i=0; rotables[i]; i++)
        {
            float rr = rotables[i]->atom1->loc.get_3d_distance(rotables[i]->atom2->loc);

            if (fabs(rr-bond_length) > 0.01)
            {
                Point aloc = rotables[i]->atom1->loc;
                Point bloc = rotables[i]->atom1->loc;

                bloc = bloc.subtract(aloc);
                bloc.scale(bond_length);
                bloc = bloc.add(aloc);

                rotables[i]->atom2->move(bloc);
            }

            // issue_5
            if (rotables[i]->rotate(bondrot[i], true))
            {
                float newanom = fsb_lsb_anomaly(first, last, lcard, bond_length);
                float nclash = get_internal_clashes();
                if ((	newanom <= anomaly
                        ||
                        (rotables[i]->can_flip && (rand()%100) < 22)
                    )
                        &&
                        nclash <= allowance * iclash
                   )
                {
                    if (_DBGCLSLOOP) cout << "Anomaly was " << anomaly << " now " << newanom << ", keeping." << endl;
                    anomaly = newanom;
                    iclash = nclash;
                }
                else
                {
                    if (_DBGCLSLOOP) cout << "Anomaly was " << anomaly << " now " << newanom << ", reverting." << endl;
                    rotables[i]->rotate(-bondrot[i]);
                    bondrot[i] *= -0.6;
                }
            }
            else
            {
                bondrot[i] *= -0.5;
            }
        }
        allowance = (allowance-1)*.99+1;
    }

    for (i=0; path[i]; i++)
    {
        path[i]->doing_ring_closure = false;
        path[i]->clear_geometry_cache();
        if (_DBGCLSLOOP) cout << "Reset " << path[i]->name << " geometry." << endl;
    }

    for (i=0; rotables[i]; i++)
        rotables[i]->can_rotate = false;

    return anomaly;
}

Atom** Molecule::get_most_bindable(int max_count)
{
    if (noAtoms(atoms)) return 0;
    if (most_bindable) return most_bindable;

    int i, j=-1, k, l;
    float best[max_count+2];
    most_bindable = new Atom*[max_count+2];

    for (k=0; k<max_count; k++)
    {
        best[k]=0;
        most_bindable[k]=0;
    }

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;

        #if _DBG_MOLBB
        cout << "Probing atom " << atoms[i]->name << endl;
        #endif

        float score = 0;
        atoms[i]->clear_geometry_cache();

        for (k=0; most_bindable[k] && k<max_count; k++)
        {
            if (most_bindable[k] == atoms[i])
            {
                #if _DBG_MOLBB
                cout << "Atom is already in return array." << endl;
                #endif
                goto _resume;
            }
            if (most_bindable[k]->is_bonded_to(atoms[i]) && best[k] >= 10)
            {
                #if _DBG_MOLBB
                cout << "Atom is bonded to " << most_bindable[k]->name << ", already in return array." << endl;
                #endif
                goto _resume;
            }
            if (most_bindable[k]->shares_bonded_with(atoms[i]) && best[k] >= 10)
            {
                #if _DBG_MOLBB
                cout << "Atom shares a bond with " << most_bindable[k]->name << ", already in return array." << endl;
                #endif
                goto _resume;
            }
        }

        if (atoms[i]->get_charge() || atoms[i]->get_acidbase())
            score += 1000;

        if (atoms[i]->is_metal())
            score += 500;

        if (atoms[i]->is_thio())
            score += 20;
        else if (atoms[i]->is_polar())
            score += 100 * fabs(atoms[i]->is_polar());

        if (atoms[i]->is_pi())
            score += 50;

        if (!score) score += 5;		// van der Waals.

        #if _DBG_MOLBB
        cout << "Score is " << score << endl;
        #endif

        for (k=0; k<max_count; k++)
        {
            if (score > best[k])
            {
                #if _DBG_MOLBB
                cout << "Score claims new #" << k << " spot." << endl;
                #endif

                for (l=max_count; l>k; l--)
                {
                    best[l] = best[l-1];
                    most_bindable[l] = most_bindable[l-1];
                }

                best[k] = score;
                most_bindable[k] = atoms[i];
                if (!k) j = i;

                break;
            }
            #if _DBG_MOLBB
            cout << "Score does not exceed previous #" << k << " spot of " << best[k] << endl;
            #endif
        }

    _resume:
        ;
    }
    most_bindable[max_count] = 0;

    if (j < 0) return 0;
    else return most_bindable;
}

#define DBG_BINDABLE 0
Atom** Molecule::get_most_bindable(int max_num, Atom* for_atom)
{
    if (!atoms) return nullptr;
    if (!for_atom) return nullptr;

    #if DBG_BINDABLE
    cout << "Molecule::get_most_bindable( " << max_num << ", " << for_atom->name << " )" << endl;
    #endif

    int mn2 = max_num+2;

    int i, j, k;
    float bb[mn2];
    Atom** bba = new Atom*[mn2];

    for (i=0; i<mn2; i++)
    {
        bb[i] = 0;
        bba[i] = nullptr;
    }

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone) continue;
        InteratomicForce* iff[32];
        InteratomicForce::fetch_applicable(atoms[i], for_atom, iff);
        if (!iff[0]) continue;

        float lbb = 0;
        for (j=0; iff[j]; j++)
        {
            float kj = iff[j]->get_kJmol();

            if (atoms[i]->get_charge()) kj *= fabs(atoms[i]->get_charge());
            if (for_atom->get_charge()) kj *= fabs(for_atom->get_charge());

            lbb += kj;
        }

        #if DBG_BINDABLE
        // Note: A competent programmer would not have to put these debugs all over the place every stinking function.
        cout << "Binding potential for " << atoms[i]->name << " " << lbb << " kJ/mol." << endl;
        #endif

        for (j=0; j<max_num; j++)
        {
            if (lbb > bb[j])
            {
                for (k=max_num-1; k>j; k--)
                {
                    bba[k] = bba[k-1];
                    bb[k] = bb[k-1];
                }

                bba[j] = atoms[i];
                bb[j] = lbb;

                #if DBG_BINDABLE
                cout << "Potential is better than previous result " << j << endl;
                #endif

                break;
            }
        }
    }

    #if DBG_BINDABLE
    cout << "Returning:" << endl << flush;
    for (i=0; i<max_num; i++)
    {
        cout << i << ": " << flush << bba[i] << flush << " ";
        if (bba[i]) cout << bba[i]->name << flush;
        cout << endl << flush;
    }
    cout << endl;
    #endif

    return bba;
}

Point Molecule::get_bounding_box() const
{
    if (noAtoms(atoms)) return 0;

    int i;
    float xmax=0, ymax=0, zmax=0;

    for (i=0; atoms[i]; i++)
    {
        Vector v(atoms[i]->loc);
        float r = atoms[i]->vdW_radius;
        Point pt(v);

        pt.x = fabs(pt.x)+r;
        pt.y = fabs(pt.y)+r;
        pt.z = fabs(pt.z)+r;

        if (pt.x > xmax) xmax = pt.x;
        if (pt.y > ymax) ymax = pt.y;
        if (pt.z > zmax) zmax = pt.z;
    }

    Point pt(xmax, ymax, zmax);
    return pt;
}

float Molecule::get_volume()
{
    if (!atoms) return 0;
    int i, j;
    float result = 0;
    for (i=0; atoms[i]; i++)
    {
        result += atoms[i]->get_sphere().volume();

        for (j=0; atoms[j]; j++)
        {
            if (i==j) continue;
            // TODO: use cap_volume() instead.
            result -= sphere_intersection(atoms[i]->vdW_radius, atoms[j]->vdW_radius, atoms[i]->distance_to(atoms[j]))/2;
        }
    }

    return result;
}

float Molecule::get_surface_area(bool p, bool oaa)
{
    return get_exposed_surface_area(nullptr, p, oaa);
}

float Molecule::get_exposed_surface_area(Molecule** neighbors, bool p, bool oaa)
{
    if (!atoms) return 0;
    int i, j, l, m, n;
    float result = 0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->is_backbone)
        {
            if (oaa) atoms[i]->molsurf_area = 0;
            continue;
        }
        float A, Aeff;
        if (oaa)
        {
            Sphere s;
            s.radius = atoms[i]->vdW_radius;
            A = s.area();
            Aeff = A;

            for (j=0; atoms[j]; j++)
            {
                if (i==j) continue;
                if (!atoms[i]->is_bonded_to(atoms[j])) continue;
                float r = atoms[i]->distance_to(atoms[j]);
                float ca = cap_area(atoms[i]->vdW_radius, atoms[j]->vdW_radius, r);
                if (ca) Aeff -= ca;
                else
                {
                    Aeff -= A * (atoms[j]->vdW_radius / atoms[i]->vdW_radius) * 1.0/(packed_sphere_ligancy/2) / (r*r);
                }
            }

            if (neighbors)
            {
                for (l=0; neighbors[l]; l++)
                {
                    if (neighbors[l] == this) continue;         // already computed own atoms.
                    n = neighbors[l]->get_atom_count();
                    for (m=0; m<n; m++)
                    {
                        float r = atoms[i]->distance_to(neighbors[l]->atoms[n]);
                        float ca = cap_area(atoms[i]->vdW_radius, neighbors[l]->atoms[n]->vdW_radius, r);
                        if (ca) Aeff -= ca;
                        else
                        {
                            Aeff -= A * 1.0 - (neighbors[l]->atoms[n]->vdW_radius / atoms[i]->vdW_radius) 
                                * 1.0/(packed_sphere_ligancy/2) / (r*r);
                        }
                    }
                }
            }

            Aeff = fmax(0, Aeff);
        }
        else Aeff = atoms[i]->molsurf_area;

        atoms[i]->molsurf_area = Aeff;
        result += Aeff * (p ? fabs(atoms[i]->is_polar()) : 1);
    }

    return result;
}

bool Molecule::ring_is_coplanar(int ringid)
{
    if (!rings) return false;
    if (!rings[ringid]) return false;
    return rings[ringid]->is_coplanar();
}

bool Molecule::ring_is_aromatic(int ringid) const
{
    if (!rings) return false;
    if (!rings[ringid]) return false;
    return rings[ringid]->get_type() == AROMATIC;
}

Point Molecule::get_ring_center(int ringid)
{
    if (!rings) return Point(0,0,0);
    if (!rings[ringid]) return Point(0,0,0);
    return rings[ringid]->get_center();
}

Vector Molecule::get_ring_normal(int ringid)
{
    if (!rings) return Vector(0,0,0);
    return rings[ringid]->get_normal();
}

Atom** Molecule::get_ring_atoms(int ringid)
{
    if (!rings) return nullptr;
    return rings[ringid]->get_atoms();
}

int Molecule::get_ring_num_atoms(int ringid)
{
    if (!rings) return 0;
    return rings[ringid]->get_atom_count();
}


void Molecule::recenter_ring(int ringid, Point new_ring_cen)
{
    if (!rings) return;
    Point old_ring_cen = get_ring_center(ringid);
    Vector motion = new_ring_cen.subtract(old_ring_cen);
    int i;
    Atom** ring_atoms = rings[ringid]->get_atoms();
    for (i=0; ring_atoms[i]; i++)
        ring_atoms[i]->move_rel(&motion);

    delete[] ring_atoms;
}

void Molecule::rotate_ring(int ringid, Rotation rot)
{
    if (!rings) return;
    Point origin = get_ring_center(ringid);
    int i;
    Atom** ring_atoms = rings[ringid]->get_atoms();
    for (i=0; ring_atoms[i]; i++)
    {
        Point aloc = ring_atoms[i]->loc;
        aloc = rotate3D(&aloc, &origin, &rot);
        ring_atoms[i]->move(aloc);
    }

    delete[] ring_atoms;
}

int Molecule::get_num_rings() const
{
    if (!rings) return 0;
    int i;
    for (i=0; rings[i]; i++);	// Get count.
    return i;
}

bool Molecule::in_same_ring(Atom* a, Atom* b)
{
    int i;
    Ring** r = a->get_rings();
    if (!r) return false;
    for (i=0; r[i]; i++)
        if (b->is_in_ring(r[i]))
        {
            delete[] r;
            return true;
        }

    delete[] r;
    return false;
}

float Molecule::get_atom_error(int i, LocatedVector* best_lv, bool hemi)
{
    int j;
    float error = 0;
    Point bloc;
    Atom* atom2;
    Bond* b[16];
    int g, bg;
    float b_bond_angle;
    LocatedVector lv;

    // Get the atom's zero-index bonded atom. Call it atom2 (because why not overuse a foolish pun?).
    atoms[i]->fetch_bonds(b);
    if (!b[0]) return error;
    atom2 = b[0]->atom2;
    float card = b[0]->cardinality;
    if (!atom2)
    {
        return error;
    }

    bloc = atom2->loc;
    g = atoms[i]->get_geometry();
    bg = atom2->get_geometry();
    b_bond_angle = atom2->get_geometric_bond_angle();

    // Make an imaginary sphere around atom2, whose radius equals the optimal bond distance.
    lv.origin = bloc;
    if (atoms[i]->Z == 1 || atom2->Z == 1) card = 1;
    float optimal_radius = InteratomicForce::covalent_bond_radius(atoms[i], atom2, card);

    lv = (Vector)atoms[i]->loc.subtract(bloc);
    lv.origin = bloc;

    error += _SANOM_BOND_RADIUS_WEIGHT * (fabs(optimal_radius-lv.r)/optimal_radius);

    error += pow(_SANOM_BOND_ANGLE_WEIGHT*atom2->get_bond_angle_anomaly(lv, atoms[i]), 2);

    float thstep = fiftyseventh*5;
    float besttheta = 0, bestphi = 0, bestscore = 0;
    if (hemi)
    {
        bestscore = -1e9;
        lv.r = InteratomicForce::covalent_bond_radius(atoms[i], atom2, card);
        for (lv.theta = -square; lv.theta <= square; lv.theta += thstep)
        {
            float phstep = M_PI/(20.0*(sin(lv.theta) + 1));
            for (lv.phi = 0; lv.phi < (M_PI*2); lv.phi += phstep)
            {
                // At many points along the sphere, evaluate the goodness-of-fit as a function of:
                // Success in conforming to atom2's geometry;
                // Success in avoiding clashes with atoms not bonded to self or atom2;
                // Success in maintaining optimal binding distances to own bonded atoms.
                // Later, we'll test edge cases where bond strain distorts the usual angles.
                float score = 0;

                score -= _SANOM_BOND_ANGLE_WEIGHT*atom2->get_bond_angle_anomaly(lv, atoms[i]);

                // Avoid clashes with strangers.
                for (j=0; atoms[j]; j++)
                {
                    if (j == i) continue;
                    if (atoms[j]->is_bonded_to(atoms[i])) continue;

                    float r = atoms[j]->loc.get_3d_distance(lv.to_point());
                    score -= _SANOM_CLASHES_WEIGHT/fabs(r+0.000000001);
                }

                // Seek optimal bond radii.
                for (j=1; b[j]; j++)
                {
                    if (!b[j]->atom2) continue;
                    float optimal = InteratomicForce::covalent_bond_radius(atoms[i], b[j]->atom2, b[j]->cardinality);
                    float r = b[j]->atom2->loc.get_3d_distance(lv.to_point());

                    score -= _SANOM_BOND_RAD_WEIGHT * fabs(optimal-r);
                }

                if (score > bestscore)
                {
                    besttheta = lv.theta;
                    bestphi = lv.phi;
                    bestscore = score;
                }
            }
        }

        if (best_lv)
        {
            best_lv->origin = lv.origin;
            best_lv->r = lv.r;
            best_lv->theta = besttheta;
            best_lv->phi = bestphi;
        }
    }

    for (j=0; atoms[j]; j++)
    {
        if (j == i) continue;
        if (atoms[j]->is_bonded_to(atoms[i])) continue;

        float r = atoms[j]->loc.get_3d_distance(lv.to_point());
        error += _SANOM_CLASHES_WEIGHT/fabs(r+0.000000001);
    }

    for (j=1; b[j]; j++)
    {
        if (!b[j]->atom2) continue;
        float optimal = InteratomicForce::covalent_bond_radius(atoms[i], b[j]->atom2, b[j]->cardinality);
        float r = b[j]->atom2->loc.get_3d_distance(lv.to_point());

        error += _SANOM_BOND_RAD_WEIGHT * fabs(optimal-r);
    }

    return fmax(0, error+bestscore);
}


#define _DEV_FIX_MSTRUCT 1
float Molecule::correct_structure(int iters)
{
    if (noAtoms(atoms)) return 0;
    int iter, i, j, k, n;
    Point zero(0,0,0);
    float error = 0;
    Point aloc, bloc;
    Atom* atom2;
    Bond* b[16];
    int g, bg;
    float b_bond_angle;
    LocatedVector lv;

    #if _DEV_FIX_MSTRUCT
    // TODO
    if (n = get_num_rings())            // Assignment, not comparison.
    {
        for (i=0; i<n; i++)
        {
            if (ring_is_aromatic(i) || get_ring_num_atoms(i) < 6)
            {
                make_coplanar_ring(get_ring_atoms(i), i);
            }
        }
    }
    else return 0;			// Non-ring structures work fine. The buggy algorithm is when rings are involved.

    for (iter=0; iter<iters; iter++)
    {
        error = 0;
        for (i=0; atoms[i]; i++)
        {
            // Get the atom's zero-index bonded atom. Call it atom2 (because why not overuse a foolish pun?).
            atoms[i]->fetch_bonds(b);
            if (!b[0]) return error;
            atom2 = b[0]->atom2;
            if (!atom2) return error;

            // TODO
            if (atoms[i]->num_rings() && atoms[i]->is_pi() && in_same_ring(atoms[i], atom2)) continue;

            error += get_atom_error(i, &lv);
            if (!iter) continue;

            // Once a "best fit" point in space is found, move there.
            if (atoms[i]->num_rings() && atoms[i]->is_pi())
            {
                for (j=0; j<n; j++)
                {
                    Atom** ring_atoms_j = get_ring_atoms(j);
                    for (k=0; ring_atoms_j[k]; k++)
                    {
                        if (ring_atoms_j[k] == atoms[i])
                        {
                            // Get distance from atom to ring center.
                            Point rcen = get_ring_center(j);
                            Point aloc = atoms[i]->loc;
                            float rad = rcen.get_3d_distance(aloc);

                            // Center ring at combined distance from atom2.
                            LocatedVector lvr = lv;
                            lvr.r += rad;
                            recenter_ring(j, lvr.to_point());

                            // Rotate ring, and all assemblies, to align atom to atom2.
                            Point atarget = lv.to_point();
                            aloc = atoms[i]->loc;
                            rcen = get_ring_center(j);
                            Rotation rot = align_points_3d(&aloc, &atarget, &rcen);
                            rotate_ring(j, rot);

                            break;
                        }
                    }
                }
            }
            else
            {
                /*Point pt = lv.to_point();
                atoms[i]->move_assembly(&pt, atom2);*/
                atoms[i]->move(lv.to_point());
            }
        }
        // cout << error << " " << atcount << endl;
        // if (error <= 10.0*atcount) break;
    }
    #endif

    error = 0;
    for (i=0; atoms[i]; i++)
    {
        error += get_atom_error(i, &lv);
    }

    return error;
}

bool Molecule::is_thiol()
{
    if (!atoms) return false;
    int i, j;

    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->Z > 10 && atoms[i]->get_family() == CHALCOGEN)
        {
            Atom* H = atoms[i]->is_bonded_to("H");
            if (H) return true;
        }
    }

    return false;
}

void Molecule::identify_cages()
{
    if (!rings) return;

    int i, j, l, n;
    for (i=0; rings[i]; i++)
    {
        Atom** ra = rings[i]->get_atoms();
        if (!ra) continue;
        for (j=0; ra[j]; j++)
        {
            for (l=j+1; ra[l]; l++)
            {
                Ring* other = ra[j]->in_same_ring_as(ra[l], rings[i]);

                #if _dbg_identify_rings
                if (other) cout << ra[j]->name << " in same ring as " << ra[l]->name << ": " << *other << endl;
                else cout << ra[j]->name << " -/- " << ra[l]->name << endl;
                #endif

                if (other)
                {
                    // Any bond between the two atoms cannot flip.
                    Bond* b = ra[j]->get_bond_between(ra[l]);
                    if (b)
                    {
                        b->can_rotate = b->can_flip = false;
                        b = ra[l]->get_bond_between(ra[j]);
                        if (b) b->can_rotate = b->can_flip = false;

                        #if _dbg_identify_rings
                        cout << *b << " cannot flip." << endl;
                        #endif
                    }
                    else
                    {
                        // If no bond between atoms, all bonds in both rings cannot flip.
                        Bond** bb = rings[i]->get_bonds();
                        if (bb) for (n=0; bb[n]; n++)
                        {
                            bb[n]->can_rotate = bb[n]->can_flip = false;
                            bb[n]->caged = true;
                        }
                        delete[] bb;
                        bb = other->get_bonds();
                        if (bb) for (n=0; bb[n]; n++)
                        {
                            bb[n]->can_rotate = bb[n]->can_flip = false;
                            bb[n]->caged = true;
                        }
                        delete[] bb;

                        #if _dbg_identify_rings
                        cout << *rings[i] << " immobilized." << endl;
                        cout << *other << " immobilized." << endl;
                        #endif
                    }
                }
            }
        }
        delete[] ra;
    }
}

float Molecule::get_atom_bond_length_anomaly(Atom* a, Atom* ignore)
{
    if (!atoms) return 0;
    if (!a) return 0;

    Bond* thesebonds[16];
    a->fetch_bonds(thesebonds);
    if (!thesebonds) return 0;
    int i;
    float anomaly = 0;

    int geometry = a->get_geometry();
    for (i=0; i<geometry; i++)
    {
        if (!thesebonds[i]) continue;
        if (thesebonds[i]->atom2)
        {
            if (thesebonds[i]->atom2 == ignore) continue;
            float optimal = InteratomicForce::covalent_bond_radius(a, thesebonds[i]->atom2, thesebonds[i]->cardinality);
            float r = a->distance_to(thesebonds[i]->atom2);
            anomaly += pow(1.0+fabs(r-optimal)/optimal, 4)-1;

            #if _dbg_molstruct_evolution_bond_lengths
            if (fabs(r-optimal) > 0.01) cout << "Atoms " << a->name << " and " << thesebonds[i]->atom2->name
                << ", cardinality = " << thesebonds[i]->cardinality
                << " should be " << optimal << "Å apart, but they're " << r << endl;
            #endif
        }
    }

    return anomaly;
}

float Molecule::evolve_structure(int gens, float mr, int ps)
{
    if (!atoms) return 0;
    int ac = get_atom_count();

    Point parents[2][ac];
    Point population[ps][ac];
    float anomalies[ps];
    int i, j, l, n, gen;
    float r, optimal;

    for (i=0; i<ac; i++)
    {
        parents[0][i] = parents[1][i] = atoms[i]->loc;
    }

    int ibest, i2best;          // Indices of best and second-best anomalies.
    float fbest, f2best;        // Values of best and second-best anomalies.

    // Main loop
    for (gen=1; gen<=gens; gen++)
    {
        ibest = i2best = -1;

        for (i=0; i<ps; i++)
        {
            for (j=0; j<ac; j++)
            {
                float atom_displacement = _evolution_atom_displacement;
                if (atoms[j]->is_pi()) atom_displacement /= _evolution_aromatic_rigidity;

                switch (i)
                {
                    case 0:
                    case 1:
                    population[i][j] = parents[i][j];
                    break;

                    default:
                    int which_parent = rand() & 0x1;
                    population[i][j].x = parents[which_parent][j].x;
                    if (frand(0, 1) < mr) population[i][j].x += frand(-atom_displacement, atom_displacement);
                    population[i][j].y = parents[which_parent][j].y;
                    if (frand(0, 1) < mr) population[i][j].y += frand(-atom_displacement, atom_displacement);
                    population[i][j].z = parents[which_parent][j].z;
                    if (frand(0, 1) < mr) population[i][j].z += frand(-atom_displacement, atom_displacement);

                    if (frand(0,1) < 0.81)
                    {
                        Bond* b0 = atoms[j]->get_bond_by_idx(0);
                        if (b0)
                        {
                            Atom* aparent = b0->atom2;
                            if (aparent)
                            {
                                optimal = InteratomicForce::covalent_bond_radius(atoms[j], aparent, b0->cardinality);
                                Vector v = population[i][j].subtract(aparent->loc);
                                v.r = optimal;
                                population[i][j] = aparent->loc.add(v);
                            }
                        }
                    }
                }

                atoms[j]->move(population[i][j]);
            }

            float anomaly = 0;
            for (j=0; j<ac; j++)
            {
                // anomaly += get_atom_error(j, nullptr, false);
                Bond* b0 = atoms[j]->get_bond_by_idx(0);
                if (!b0) continue;

                Atom* aparent = b0->atom2;
                if (!aparent) continue;

                anomaly += get_atom_bond_length_anomaly(aparent, atoms[j]);

                Vector v = atoms[j]->loc.subtract(aparent->loc);
                anomaly += aparent->get_bond_angle_anomaly(v, atoms[j]);
            }
            anomalies[i] = anomaly;

            if (ibest < 0 || anomaly < fbest)
            {
                i2best = ibest;
                f2best = fbest;
                ibest = i;
                fbest = anomaly;
            }
            else if (i2best < 0 || anomaly < f2best)
            {
                i2best = i;
                f2best = anomaly;
            }
        }

        for (i=0; i<ac; i++)
        {
            parents[0][i] = population[ibest][i];
            parents[1][i] = population[i2best][i];
        }

        #if _dbg_molstruct_evolutions
        cout << ibest << ':' << fbest << ' ' << i2best << ':' << f2best << endl << flush;
        #endif
    }

    for (i=0; i<ac; i++)
    {
        atoms[i]->move(population[ibest][i]);
    }

    return fbest / get_atom_count();
}

bool Molecule::is_chiral()
{
    if (noAtoms(atoms)) return false;

    int i;
    Bond* bndbuf[16];
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->get_geometry() < 4) continue;
        atoms[i]->fetch_bonds(bndbuf);

        if (!bndbuf[0]) continue;
        if (!bndbuf[1]) continue;
        if (!bndbuf[2]) continue;
        if (bndbuf[0]->is_equivalent(bndbuf[1])) continue;
        if (bndbuf[0]->is_equivalent(bndbuf[2])) continue;
        if (bndbuf[1]->is_equivalent(bndbuf[2])) continue;
        if (bndbuf[3])
        {
            if (bndbuf[0]->is_equivalent(bndbuf[3])) continue;
            if (bndbuf[1]->is_equivalent(bndbuf[3])) continue;
            if (bndbuf[2]->is_equivalent(bndbuf[3])) continue;
        }
        return true;
    }

    return false;
}

void Molecule::mirror()
{
    if (noAtoms(atoms)) return;

    if (movability <= MOV_FLEXONLY) return;

    int i;
    Point cen = get_barycenter();
    for (i=0; atoms[i]; i++)
    {
        // if (atoms[i]->residue) return;
        Point loc = atoms[i]->loc;
        loc = loc.subtract(cen);
        loc.x *= -1;
        loc = loc.add(cen);
        atoms[i]->move(&loc);
    }

    if (vdw_surface && vdw_vertex_count)
    {
        for (i=0; i<vdw_vertex_count; i++)
        {
            Point loc = vdw_surface[i];
            loc = loc.subtract(cen);
            loc.x *= -1;
            loc = loc.add(cen);
            vdw_surface[i] = loc;
        }
    }
}

float g_total_mclash(void* mol)
{
    return reinterpret_cast<Molecule*>(mol)->get_total_mclashes();
}

bool mclash_delta(void* mol, float previous_mclashes)
{
    return reinterpret_cast<Molecule*>(mol)->get_total_mclashes() <= previous_mclashes;
}









