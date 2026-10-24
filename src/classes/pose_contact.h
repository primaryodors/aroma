#ifndef AROMA_POSE_CONTACT_H
#define AROMA_POSE_CONTACT_H
#include "protein.h"

class PoseContactSnapshot
{
    struct Contact
    {
        std::string ligand_atom, receptor_atom;
        int residue = 0;
        void capture(Atom* mine, Atom* other)
        {
            ligand_atom = mine ? mine->name : "";
            receptor_atom = other ? other->name : "";
            residue = other ? other->residue : 0;
        }
        void restore(Molecule* ligand, Protein* protein, Atom*& mine, Atom*& other, Molecule*& owner) const
        {
            AminoAcid* aa = protein->get_residue(residue);
            mine = ligand_atom.empty() ? nullptr : ligand->get_atom(ligand_atom.c_str());
            other = aa && !receptor_atom.empty() ? aa->get_atom(receptor_atom.c_str()) : nullptr;
            owner = aa;
            if (!mine || !other) { mine = other = nullptr; owner = nullptr; }
        }
    } primary, secondary;
    float tolerance = 0, optimal = 2, secondary_optimal = 2;
public:
    void capture(Molecule* ligand)
    {
        primary.capture(ligand->stay_close_mine, ligand->stay_close_other);
        secondary.capture(ligand->stay_close2_mine, ligand->stay_close2_other);
        tolerance = ligand->stay_close_tolerance;
        optimal = ligand->stay_close_optimal;
        secondary_optimal = ligand->stay_close2_optimal;
    }
    void restore(Molecule* ligand, Protein* protein) const
    {
        primary.restore(ligand, protein, ligand->stay_close_mine, ligand->stay_close_other, ligand->stay_close_mol);
        secondary.restore(ligand, protein, ligand->stay_close2_mine, ligand->stay_close2_other, ligand->stay_close2_mol);
        ligand->stay_close_tolerance = tolerance;
        ligand->stay_close_optimal = optimal;
        ligand->stay_close2_optimal = secondary_optimal;
    }
};
#endif
