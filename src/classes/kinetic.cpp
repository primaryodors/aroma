
#include <cmath>
#include "kinetic.h"
#include "intera.h"

void MolecularKinetics::allocate_momenta()
{
    if (momenta_allocated)
    {
        // TODO: Support adding new atoms after initial allocation.
        throw 0x9071e7;
    }
    else
    {
        int n = atoms.size()+1;
        momenta = new Vector[n];
        momenta_allocated = n-1;
    }
}

bool MolecularKinetics::atom_is_nearby(Atom *test, Atom *rel_to)
{
    return test->distance_to(rel_to) < _INTERA_R_CUTOFF;
}

Vector MolecularKinetics::get_energetic_force(Atom *a)
{
    int g, i;
    Vector result(0,0,0);

    g = a->get_geometry();
    for (i=0; i<g; i++)
    {
        Bond *b = a->get_bond_by_idx(i);
        if (b && b->atom2)
        {
            // Find difference between optimal and actual distance
            float displacement = a->distance_to(b->atom2) - InteratomicForce::covalent_bond_radius(a, b->atom2, b->cardinality);

            // Get directional vector
            Vector dir = b->atom2->loc.subtract(a->loc);

            // Get bond strength in kJ/mol.
            float Eopt = InteratomicForce::covalent_bond_energy(a, b->atom2, b->cardinality);
            float Eact = 0.5 * Eopt * pow(displacement, 2);
            float Edelta = Eopt - Eact;                     // should always be positive because we squared the displacement

            // Interaction is in kJ/mol, and 1J = 1 meter in 1 second, so force component strength in m/s = kJ/mol * 1e3 / mass.
            dir.r = Edelta / a->get_atomic_weight();

            result = result.add(dir);
        }
    }

    for (i=0; atoms._atoms[i]; i++)
    {
        if (atoms._atoms[i] == a) continue;
        if (atom_is_nearby(atoms._atoms[i], a))
        {
            // Find displacement and optimal distance
            float optimal_distance = InteratomicForce::optimal_distance(a, atoms._atoms[i]);
            float r = a->distance_to(atoms._atoms[i]);
            float displacement = r - optimal_distance;

            // Find directional vector
            Vector dir = a->loc.subtract(atoms._atoms[i]->loc);

            // Compute energetic anomaly
            float Eopt = InteratomicForce::potential_binding(a, atoms._atoms[i], false);
            float Eact = InteratomicForce::total_binding(a, atoms._atoms[i]).summed();
            float Edelta = fabs(Eopt - Eact) * -sgn(displacement);

            // Interaction is in kJ/mol, and 1J = 1 meter in 1 second, so force component strength in m/s = kJ/mol * 1e-3.
            dir.r = Edelta / a->get_atomic_weight();

            result = result.add(dir);
        }
    }

    return result;
}

MolecularKinetics::MolecularKinetics()
{
}

MolecularKinetics::~MolecularKinetics()
{
}

void MolecularKinetics::set_Boltzmann_momenta()
{
    int i, n = atoms.size();
    if (!n) throw 0xbadc0de;             // have to put atoms in the collection before calling the Boltzmann function
    if (momenta_allocated < n) allocate_momenta();

    for (i=0; i<n; i++)
    {
        double sigma = pow(8.0 * kB * temperature / atoms._atoms[i]->get_atomic_weight() / Dalton / M_PI, 0.5) * 31.5;          // wtf is this 31.5 constant? it seems to make the math work...
        Vector v = Point(frand(-1,1),frand(-1,1),frand(-1,1));
        v.r = generate_gaussian(0, sigma);
        momenta[i] = v;
    }
}

void MolecularKinetics::advance_clock(float femtoseconds)
{
    int i;

    // Calculate interatomic forces.
    Vector forces[atoms.size()+1];
    for (i=0; atoms._atoms[i]; i++)
    {
        forces[i] = get_energetic_force(atoms._atoms[i]);
    }

    // Update locations from existing momenta.
    for (i=0; atoms._atoms[i]; i++)
    {
        // Momenta are in meters per second, and we're using angstroms and femtoseconds,
        // so multiply by 1e10 and divide by 1e15.
        Vector motion = momenta[i];
        motion.r *= femtoseconds / 1e5;

        atoms._atoms[i]->move(atoms._atoms[i]->loc.add(motion));
    }

    // Update momenta according to forces.
    for (i=0; atoms._atoms[i]; i++)
    {
        momenta[i] = momenta[i].add(forces[i]);
    }
}

void MolecularKinetics::dump()
{
    int i, n = atoms.size();

    cout << n << " atoms." << endl;

    for (i=0; i<n; i++)
    {
        Atom* a = atoms._atoms[i];
        if (a->residue) cout << a->residue << ":";
        cout << a->name;
        cout << " location=" << a->loc;
        cout << " momentum=" << (Point)momenta[i];
        cout << endl;
    }
}
