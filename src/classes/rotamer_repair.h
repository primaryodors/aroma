#ifndef AROMA_ROTAMER_REPAIR_H
#define AROMA_ROTAMER_REPAIR_H

#include "protein.h"
#include <string>
#include <vector>

struct ReceptorRotamer
{
    struct Site
    {
        std::string name;
        int element;
        bool backbone;
        Point location;
    };
    char residue_type = 0;
    std::vector<Site> sites;

    void capture(AminoAcid& residue)
    {
        residue_type = residue.get_letter();
        for (int i = 0; i < residue.get_atom_count(); ++i)
        {
            Atom* atom = residue.get_atom(i);
            sites.push_back({atom->name, atom->Z, atom->is_backbone, atom->loc});
        }
    }

    bool repair(AminoAcid& residue, Protein& protein, Molecule& ligand,
                Molecule* metal, Molecule** waters, float clash_limit) const
    {
        if (residue_type != residue.get_letter()
            || sites.size() != static_cast<size_t>(residue.get_atom_count())) return false;

        auto clear_of_ligand = [&] {
            Atom *first = nullptr, *second = nullptr;
            residue.mutual_closest_atoms(&ligand, &first, &second);
            return first && second && first->distance_to(second) >= _INTERA_R_CUTOFF;
        };
        if (!clear_of_ligand()) return false;

        std::vector<Atom*> side_chain;
        std::vector<Point> original, replacement;
        for (const auto& site : sites)
        {
            Atom* atom = residue.get_atom(site.name.c_str());
            if (!atom || atom->Z != site.element || atom->is_backbone != site.backbone) return false;
            // The input rotamer is in the input backbone frame.
            if (site.backbone || site.name == "CB")
            {
                if (atom->loc.get_3d_distance(site.location) > 0.01f) return false;
                continue;
            }
            side_chain.push_back(atom);
            original.push_back(atom->loc);
            replacement.push_back(site.location);
        }
        if (side_chain.empty()) return false;

        struct Clash { Molecule* neighbor; float before; };
        std::vector<Clash> clashes;
        auto record = [&](Molecule* neighbor) {
            if (neighbor && neighbor != &residue)
                clashes.push_back({neighbor, residue.get_intermol_clashes(neighbor)});
        };
        for (int r = protein.get_start_resno(); r <= protein.get_end_resno(); ++r)
            record(protein.get_residue(r));
        record(metal);
        if (waters) for (int i = 0; waters[i]; ++i) record(waters[i]);

        for (size_t i = 0; i < side_chain.size(); ++i) side_chain[i]->move(replacement[i]);
        bool repaired = clear_of_ligand() && residue.get_internal_clashes(false, false) <= clash_limit;
        for (const auto& clash : clashes)
        {
            if (!repaired) break;
            repaired = residue.get_intermol_clashes(clash.neighbor) <= clash.before + 0.1f;
        }
        if (!repaired)
            for (size_t i = 0; i < side_chain.size(); ++i) side_chain[i]->move(original[i]);
        residue.invalidate_memoized_clashes();
        return repaired;
    }
};

#endif
