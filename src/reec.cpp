
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <string>
#include <regex>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <algorithm>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "classes/search.h"
#include "classes/scoring.h"
#include "classes/cavity.h"
#include "classes/progress.h"

using namespace std;

int main(int argc, char** argv)
{
    cout << "Receptor Energetic Efficacy Calculator" << endl;

    Protein prot("TheProtein");
    Molecule ligand("TheLigand");
    Progressbar pbr;
    Box search_area(0,0,0,0,0,0);
    ResidueAtomPlaceholder bsr[256];
    int nbsr = 0;
    std::string protname, ligname;

    int i, j, l, n;
    FILE* fp;
    char buffer[65536];

    for (i=1; i<argc; i++)
    {
        if (!strcmp(argv[i], "-a") || !strcmp(argv[i], "--area"))
        {
            search_area.x1 = atof(argv[++i]);
            search_area.y1 = atof(argv[++i]);
            search_area.z1 = atof(argv[++i]);
            search_area.x2 = atof(argv[++i]);
            search_area.y2 = atof(argv[++i]);
            search_area.z2 = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bsr"))
        {
            i++;
            while (i<argc && argv[i][0] >= '0' && argv[i][0] <= '9')
            {
                bsr[nbsr++].set(argv[i++]);
            }
            i--;
        }
        else if (file_exists(argv[i]))
        {
            char* fttl = &(strrchr(argv[i], '/')[1]);
            char* ext = strrchr(argv[i], '.');
            if (!strcmp(ext, ".sdf"))
            {
                fp = fopen(argv[i], "rb");
                if (!fp)
                {
                    cerr << "FAILED to open " << argv[i] << " for reading." << endl;
                    return 1;
                }
                n = fread(buffer, sizeof(char), 65534, fp);
                buffer[n] = 0;
                fclose(fp);
                ligand.from_sdf(buffer);

                *ext = 0;
                ligname = fttl;
                *ext = '.';
            }
            else if (!strcmp(ext, ".pdb"))
            {
                fp = fopen(argv[i], "rb");
                if (!fp)
                {
                    cerr << "FAILED to open " << argv[i] << " for reading." << endl;
                    return 1;
                }
                prot.load_pdb(fp);
                fclose(fp);

                *ext = 0;
                protname = fttl;
                *ext = '.';
            }
        }
    }

    if (!prot.get_seq_length())
    {
        cerr << "Error no protein residues." << endl;
        return 1;
    }
    if (!ligand.get_atom_count())
    {
        cerr << "Error no ligand atoms." << endl;
        return 1;
    }

    for (i=0; i<nbsr; i++)
    {
        bsr[i].resolve_resno(&prot);
    }

    if (nbsr && !search_area.volume())
    {
        j=0;
        for (i=0; i<nbsr; i++)
        {
            AminoAcid* aa = prot.get_residue(bsr[i].resno);
            if (!aa) continue;

            Point CA = aa->get_CA_location();
            // cout << aa->get_name() << " CA = " << CA << endl;
            if (!j || CA.x < search_area.x1) search_area.x1 = CA.x;
            if (!j || CA.y < search_area.y1) search_area.y1 = CA.y;
            if (!j || CA.z < search_area.z1) search_area.z1 = CA.z;
            if (!j || CA.x > search_area.x2) search_area.x2 = CA.x;
            if (!j || CA.y > search_area.y2) search_area.y2 = CA.y;
            if (!j || CA.z > search_area.z2) search_area.z2 = CA.z;
            j++;
        }
    }

    cav_xmin = search_area.x1;
    cav_xmax = search_area.x2;
    cav_ymin = search_area.y1;
    cav_ymax = search_area.y2;
    cav_zmin = search_area.z1;
    cav_zmax = search_area.z2;
    pbr.set_color(224, 192, 64);            // curmi uelor
    cout << "Performing cavity search within box "
        << "(" << cav_xmin << "," << cav_ymin << "," << cav_zmin << "), "
        << "(" << cav_xmax << "," << cav_ymax << "," << cav_zmax << ")"
        << "..." << endl;

    Cavity cavities[1029];
    int qfound = Cavity::scan_in_protein(&prot, cavities, 1024, &pbr);
    if (!qfound)
    {
        cerr << "No cavities found." << endl;
        return 1;
    }

    // Biggest cavity
    j = n = 0;
    for (i=0; i<qfound; i++)
    {
        l = cavities[i].count_partials();
        if (l>n)
        {
            n = l;
            j = i;
        }
    }

    cout << "Total cavity volume: " << cavities[j].get_volume() << " Å³" << endl;
    cout << "Ligand volume: " << ligand.get_volume() << " Å³" << endl;

    Point pocketcen = cavities[j].get_center();

    pbr.set_color(160, 224, 64);
    cout << "Tumbling ligand inside cavity..." << endl;
    prot.tumble_ligand_inside_pocket(&ligand, pocketcen, 1, &pbr);

    int bsresno[1024];
    Molecule* bsresm[1024];
    n = cavities[j].resnos_as_array(&prot, bsresno);
    l = 0;
    for (i=0; i<n; i++)
    {
        AminoAcid* aa = prot.get_residue(bsresno[i]);
        if (!aa) continue;
        bsresm[l++] = (Molecule*)aa;
    }
    bsresm[l] = nullptr;

    // TODO: Move individual ligand atoms outward from ligand center to where they're energetically favored.

    pbr.set_color(64, 224, 128);
    cout << "Now the magic happens..." << endl;
    ligand.refine_structure(5000, 0.1, 100, bsresm, &pbr);

    // TODO: Perform bond rotations of heavy atoms bound to hydrogen and only one other heavy atom.

    // TODO: Once the ligand is optimized, individually refine all the BSRs around it.

    // TODO: Then refine the rest of the protein one residue at a time.

    std::string fam = family_from_protid(protname);
    char outfname[1024];

    mode_t permissions = 0755;
    sprintf(outfname, "output/%s", fam.c_str());
    if (!dir_exists(outfname)) if (mkdir(outfname, permissions))
    {
        cerr << "FAILED to create " << outfname << endl;
        return 1;
    }
    sprintf(outfname, "output/%s/%s", fam.c_str(), protname.c_str());
    if (!dir_exists(outfname)) if (mkdir(outfname, permissions))
    {
        cerr << "FAILED to create " << outfname << endl;
        return 1;
    }

    sprintf(outfname, "output/%s/%s/%s~%s.bound.pdb", fam.c_str(), protname.c_str(), protname.c_str(), ligname.c_str());
    fp = fopen(outfname, "wb");
    if (!fp)
    {
        cerr << "FAILED to open " << outfname << " for writing." << endl;
        return 1;
    }
    prot.save_pdb(fp, &ligand);
    fclose(fp);
    cout << "Wrote " << outfname << endl;
}

