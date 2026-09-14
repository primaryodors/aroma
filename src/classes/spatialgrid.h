#ifndef _SPATIAL_GRID_H
#define _SPATIAL_GRID_H

#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
#include "point.h"
#include "atom.h"
#include "aminoacid.h"

struct VoxelCoord
{
    int x, y, z;
    bool operator==(const VoxelCoord& o) const
    {
        return x == o.x && y == o.y && z == o.z;
    }
};

struct VoxelCoordHash
{
    std::size_t operator()(const VoxelCoord& k) const
    {
        return (static_cast<std::size_t>(k.x) * 73856093) ^
               (static_cast<std::size_t>(k.y) * 19349663) ^
               (static_cast<std::size_t>(k.z) * 83492791);
    }
};

class SpatialGrid
{
public:
    float cellSize;
    std::unordered_map<VoxelCoord, std::vector<Atom*>, VoxelCoordHash> atom_cells;

    SpatialGrid(float cell_size = 7.0f) : cellSize(cell_size) {}

    void clear()
    {
        atom_cells.clear();
    }

    VoxelCoord coord_for(const Point& pt) const
    {
        return VoxelCoord
        {
            (int)floor(pt.x / cellSize),
            (int)floor(pt.y / cellSize),
            (int)floor(pt.z / cellSize)
        };
    }

    void insert(Atom* a)
    {
        if (!a) return;
        atom_cells[coord_for(a->loc)].push_back(a);
    }

    void build(AminoAcid** residues)
    {
        clear();
        if (!residues) return;
        for (int i = 0; residues[i]; i++)
        {
            int atcount = residues[i]->get_atom_count();
            for (int j = 0; j < atcount; j++)
            {
                Atom* at = residues[i]->get_atom(j);
                if (at) insert(at);
            }
        }
    }

    Atom* get_nearest_atom(const Point& pt, int sr = 0, int er = 0) const
    {
        if (atom_cells.empty()) return nullptr;
        VoxelCoord center = coord_for(pt);
        Atom* best = nullptr;
        float best_d2 = 1e18f;

        for (int ring = 0; ring <= 8; ring++)
        {
            for (int dx = -ring; dx <= ring; dx++)
            {
                for (int dy = -ring; dy <= ring; dy++)
                {
                    for (int dz = -ring; dz <= ring; dz++)
                    {
                        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) != ring) continue;
                        VoxelCoord vc{center.x + dx, center.y + dy, center.z + dz};
                        auto it = atom_cells.find(vc);
                        if (it == atom_cells.end()) continue;
                        for (Atom* a : it->second)
                        {
                            if (sr && a->residue < sr) continue;
                            if (er && a->residue > er) continue;
                            float d2 = a->loc.get_3d_distance_squared(pt);
                            if (d2 < best_d2)
                            {
                                best_d2 = d2;
                                best = a;
                            }
                        }
                    }
                }
            }
            if (best)
            {
                float min_d_next = ring * cellSize;
                if (min_d_next * min_d_next >= best_d2)
                {
                    break;
                }
            }
        }
        return best;
    }
};

#endif
