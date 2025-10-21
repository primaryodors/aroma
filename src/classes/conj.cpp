
#include "conj.h"


Conjugation::Conjugation()
{
    atoms = new Atom*[1024];
    atoms[0] = nullptr;
}

Conjugation::Conjugation(Atom* from_atom)
{
    atoms = new Atom*[1024];
    atoms[0] = nullptr;
    add_atom(from_atom);
}

Conjugation::Conjugation(Ring *from_ring)
{
    int i, n = from_ring->get_atom_count();
    atoms = new Atom*[n+4];
    nheavy = 0;
    for (i=0; i<n; i++)
    {
        atoms[i] = from_ring->get_atom(i);
        if (atoms[i]->Z > 1) nheavy++;
    }
    atoms[n] = nullptr;
}

Conjugation::~Conjugation()
{
    delete[] atoms;
}

float conj_get_charge_ftn(void* conj)
{
    return reinterpret_cast<Conjugation*>(conj)->get_net_charge();
}

void Conjugation::add_atom(Atom* add)
{
    add_atom(add, add, add);
}

void Conjugation::add_atom(Atom* a, Atom* prev, Atom* orig)
{
    int i, j;

    if (!a->is_pi()) return;

    for (j=0; atoms[j]; j++)
    {
        if (atoms[j] == a) return;
    }

    atoms[j] = a;
    atoms[j+1] = nullptr;
    if (a->Z > 1) nheavy++;
    net_charge_known = false;
    a->conjugation = this;
    conj_get_charge = &conj_get_charge_ftn;

    Bond* b[16];
    for (i=0; i<16; i++) b[i] = nullptr;
    a->fetch_bonds(b);

    for (i=0; b[i]; i++)
    {
        if (!b[i]->atom2) continue;
        if (!b[i]->atom2->is_pi()) continue;
        if (b[i]->atom2 == prev) continue;
        if (b[i]->atom2 == orig) continue;

        add_atom(b[i]->atom2, a, orig);              // RECURSION!
    }
}

float Conjugation::get_net_charge()
{
    if (!net_charge_known)
    {
        int i;
        net_charge = 0;
        for (i=0; atoms[i]; i++)
        {
            if (atoms[i]->Z <= 1) continue;
            net_charge += atoms[i]->get_orig_charge();
        }

        net_charge_known = true;
    }

    return net_charge;
}

bool Conjugation::has_hb_acceptors()
{
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->Z <= 1) continue;
        if (atoms[i]->get_num_lone_pairs()) return true;
    }
    
    return false;
}

bool Conjugation::has_hb_donors()
{
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->Z <= 1) continue;
        Atom* H = atoms[i]->is_bonded_to("H");
        if (H && H->is_polar() >= hydrophilicity_cutoff) return true;
    }
    
    return false;
}

bool Conjugation::contains(Atom *a)
{
    int i;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i] == a) return true;
    }
    return false;
}

bool Conjugation::contains(Conjugation *c)
{
    if (!c->atoms) return false;
    int i, j=0;
    for (i=0; c->atoms[i]; i++)
    {
        if (!contains(c->atoms[i])) return false;
    }
    return true;
}

float Conjugation::get_sum_polarity()
{
    int i;
    float polar = 0;
    for (i=0; atoms[i]; i++)
    {
        if (atoms[i]->Z <= 1) continue;
        polar += fabs(atoms[i]->is_polar());
    }

    return polar;
}

Point Conjugation::get_barycenter()
{
    if (!atoms || !atoms[0]) return Point(0,0,0);

    Point result, pt;
    float divisor = 0;
    int i;

    for (i=0; atoms[i]; i++)
    {
        pt = atoms[i]->loc;
        float atwt = atoms[i]->get_atomic_weight();
        pt.multiply(atwt);
        result = result.add(pt);
        divisor += atwt;
    }

    if (divisor) result.multiply(1.0 / divisor);

    return result;
}

Atom* Conjugation::get_nearest_atom(Point pt)
{
    if (!atoms || !atoms[0]) return nullptr;

    float r, best_r;
    Atom* best_a;
    int i;

    for (i=0; atoms[i]; i++)
    {
        r = atoms[i]->loc.get_3d_distance(pt);
        if (!i || r < best_r)
        {
            best_a = atoms[i];
            best_r = r;
        }
    }

    return best_a;
}

Atom* Conjugation::get_nearest_atom(Conjugation* conj)
{
    /*if (mutual == conj)
    {
        mutual->mutual = nullptr;
        mutual = nullptr;
        return mutual_nearest;
    }*/

    Point bcen1 = get_barycenter(), bcen2 = conj->get_barycenter();
    Atom* a = get_nearest_atom(bcen2), *b = conj->get_nearest_atom(bcen1);

    a = get_nearest_atom(b->loc);
    b = get_nearest_atom(a->loc);
    a = get_nearest_atom(b->loc);
    b = get_nearest_atom(a->loc);
    a = get_nearest_atom(b->loc);
    b = get_nearest_atom(a->loc);

    conj->mutual = this;
    conj->mutual_nearest = b;

    return a;
}

std::string Conjugation::to_std_string()
{
    std::string result;
    int i, n = count_atoms();
    for (i=0; i<n; i++)
    {
        if (i) result += std::string(" ");
        Atom* a = get_atom(i);
        if (a) result += std::string(a->name);
        else result += std::string("(null)");
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const Conjugation& c)
{
    int i, n = c.count_atoms();
    for (i=0; i<n; i++)
    {
        if (i) os << " ";
        Atom* a = c.get_atom(i);
        if (a) os << a->name;
        else os << "(null)";
    }

    return os;
}