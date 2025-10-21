
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "atom.h"

#ifndef _CONJ
#define _CONJ

class Conjugation
{
    public:
    Conjugation();
    Conjugation(Atom* from_atom);
    Conjugation(Ring* from_ring);
    ~Conjugation();

    void add_atom(Atom* add);
    float get_net_charge();
    float get_sum_polarity();
    Point get_barycenter();
    Atom* get_nearest_atom(Point pt);
    Atom* get_nearest_atom(Conjugation* conj);
    int count_atoms() const { if (!atoms) return 0; int i; for (i=0; atoms[i]; i++); return i; }
    Atom* get_atom(int i) const { if (!atoms) return nullptr; return atoms[i]; }
    bool has_hb_acceptors();
    bool has_hb_donors();
    bool contains(Atom* a);
    bool contains(Conjugation* c);
    std::string to_std_string();

    protected:
    Atom** atoms = nullptr;           // Not using a vector for performance reasons.
    float net_charge = 0;
    bool net_charge_known = false;
    Conjugation* mutual = nullptr;
    Atom* mutual_nearest = nullptr;
    int nheavy = 0;

    void add_atom(Atom* a, Atom* prev, Atom* orig);

public:
    const int& num_heavy_atoms = nheavy;
};

float conj_get_charge_ftn(void* conj);



std::ostream& operator<<(std::ostream& os, const Conjugation& c);

#endif
