#ifndef AROMA_RESIDUE_GRID_H
#define AROMA_RESIDUE_GRID_H

#include "spatialgrid.h"

class ResidueGrid
{
public:
    struct Entry
    {
        AminoAcid* residue;
        Point center;
        float reach;
    };
    std::vector<Entry> entries;
    std::unordered_map<VoxelCoord,std::vector<size_t>,VoxelCoordHash> cells;
    float max_reach = 0;
    static constexpr double cell_size = 12.0;

    VoxelCoord coordinate(const Point& point) const
    {
        return {int(floor(point.x/cell_size)),int(floor(point.y/cell_size)),int(floor(point.z/cell_size))};
    }

    void update(AminoAcid** residues)
    {
        size_t count = 0;
        bool changed = false;
        for (; residues && residues[count]; count++)
        {
            Point center = residues[count]->get_CA_location();
            if (count >= entries.size() || entries[count].residue != residues[count]
                || entries[count].center.x != center.x || entries[count].center.y != center.y
                || entries[count].center.z != center.z || entries[count].reach != residues[count]->get_reach())
                changed = true;
        }
        if (!changed && count == entries.size()) return;
        entries.clear();
        cells.clear();
        max_reach = 0;
        for (size_t i=0; i<count; i++)
        {
            Entry entry{residues[i],residues[i]->get_CA_location(),residues[i]->get_reach()};
            entries.push_back(entry);
            cells[coordinate(entry.center)].push_back(i);
            max_reach = fmax(max_reach,entry.reach);
        }
    }

    void select(const Point& center, double radius, std::vector<bool>& selected) const
    {
        radius = fmax(radius,0.0)+0.001;
        VoxelCoord low = coordinate(Point(center.x-radius,center.y-radius,center.z-radius));
        VoxelCoord high = coordinate(Point(center.x+radius,center.y+radius,center.z+radius));
        double volume = double(high.x-low.x+1)*(high.y-low.y+1)*(high.z-low.z+1);
        auto accept = [&](size_t i) {
            if (entries[i].center.get_3d_distance_squared(center) <= radius*radius) selected[i] = true;
        };
        if (volume > entries.size()*4)
        {
            for (size_t i=0; i<entries.size(); i++) accept(i);
            return;
        }
        for (int x=low.x; x<=high.x; x++)
        for (int y=low.y; y<=high.y; y++)
        for (int z=low.z; z<=high.z; z++)
        {
            auto found = cells.find({x,y,z});
            if (found != cells.end()) for (size_t i : found->second) accept(i);
        }
    }

    std::vector<AminoAcid*> query(const Point& center, double radius, const Point* other = nullptr,
                                bool priority = false, const int* additional = nullptr) const
    {
        std::vector<bool> selected(entries.size(),false);
        select(center,radius,selected);
        if (other) select(*other,radius,selected);
        std::vector<AminoAcid*> result;
        for (size_t i=0; i<entries.size(); i++)
        {
            AminoAcid* residue = entries[i].residue;
            bool forced = priority && residue->priority;
            for (int j=0; additional && additional[j]; j++)
                if (additional[j] == residue->get_residue_no()) forced = true;
            if (selected[i] || forced) result.push_back(residue);
        }
        return result;
    }
};

#endif
