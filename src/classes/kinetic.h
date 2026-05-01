
// THIS is the molecular dynamics class.

#ifndef _KINETIC
#define _KINETIC

#include "atom.h"

class MolecularKinetics
{
    protected:
    Vector* momenta = nullptr;
    int momenta_allocated = 0;

    void allocate_momenta();
    bool atom_is_nearby(Atom* test, Atom* rel_to);
    Vector get_energetic_force(Atom* a);

    public:
    MolecularKinetics();
    ~MolecularKinetics();

    static float generate_Boltzmann_velocity(float mass_Daltons);
    void set_Boltzmann_momenta();
    void advance_clock(float femtoseconds);
    void dump();

    AtomCollection atoms;
};



#endif