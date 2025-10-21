#include "classes/protein.h"

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        cout << "Usage:" << endl << endl;
        cout << "bin/qc path/to/structure.pdb A B" << endl;
        cout << endl << "...replacing A and B with single letters specifying the target strands of the source PDB file." << endl;
        cout << "This program will scan the boundary between specified strands to find the binding energies of contacts." << endl;
        return -1;
    }
    char* fname;
    Protein p("strand1");
    Protein q("strand2");
    char strandid1, strandid2;

    bool dohyd = false, polaronly = false, totalall = false;
    float ttl = 0;

    int i;
    FILE* fp;
    float threshold = -3;
    for (i=1; i<argc; i++)
    {
        if (!strcmp(argv[i], "hydro"))
        {
            dohyd = true;
        }
        else if (!strcmp(argv[i], "nooil"))
        {
            polaronly = true;
        }
        else if (!strcmp(argv[i], "tall"))
        {
            totalall = true;
        }
        else if (atof(argv[i])) threshold = atof(argv[i]);
        else if (strstr(argv[i], ".pdb"))
        {
            fname = argv[i];
            fp = fopen(fname, "rb");
            if (!fp)
            {
                cout << "Failed to open " << fname << " for reading." << endl;
                return -1;
            }
            if (i<argc-1 && strlen(argv[i+1]) == 1) p.load_pdb(fp, 0, argv[++i][0]);
            else return -2;
            strandid1 = argv[i][0];
            cout << "Loaded " << p.get_seq_length() << " residues from strand " << strandid1 << "." << endl;
            rewind(fp);
            if (i<argc-1 && strlen(argv[i+1]) == 1) q.load_pdb(fp, 0, argv[++i][0]);
            else return -2;
            strandid2 = argv[i][0];
            cout << "Loaded " << q.get_seq_length() << " residues from strand " << strandid2 << "." << endl;
            fclose(fp);
        }
    }

    if (dohyd)
    {
        int resno, endres = p.get_end_resno();
        cout << "Hydrogenating...";
        for (resno=1; resno<=endres; resno++)
        {
            AminoAcid* res = p.get_residue(resno);
            if (res)
            {
                res->hydrogenate();
            }
            cout << "." << flush;
        }

        endres = q.get_end_resno();
        for (resno=1; resno<=endres; resno++)
        {
            AminoAcid* res = q.get_residue(resno);
            if (res)
            {
                res->hydrogenate();
            }
            cout << "." << flush;
        }
        cout << endl;
    }

    p.optimize_hydrogens();
    q.optimize_hydrogens();

    int j, m = p.get_end_resno(), n = q.get_end_resno();

    for (i=1; i<=m; i++)
    {
        AminoAcid* aa1 = p.get_residue(i);
        if (!aa1) continue;

        for (j=1; j<=n; j++)
        {
            AminoAcid* aa2 = q.get_residue(j);
            if (!aa2) continue;

            Atom *a, *b;
            aa1->Molecule::mutual_closest_atoms(aa2, &a, &b);
            if (!a || !b) continue;
            float r = a->distance_to(b);
            if (r > _INTERA_R_CUTOFF) continue;

            bool contact_polar = (fabs(a->is_polar()) >= hydrophilicity_cutoff && fabs(a->get_heavy_atom()->is_polar()) >= hydrophilicity_cutoff)
                && (fabs(b->is_polar()) >= hydrophilicity_cutoff && fabs(b->get_heavy_atom()->is_polar()) >= hydrophilicity_cutoff);

            if (polaronly && !totalall)
            {
                if (!contact_polar) continue;
            }

            Interaction e = aa1->get_intermol_binding(aa2, false);
            // if (a->distance_to(b) > (a->vdW_radius+b->vdW_radius)) e.clash = 0;
            float f = e.summed();
            ttl += f;
            if (polaronly)
            {
                if (!contact_polar) continue;
            }

            a = a->get_heavy_atom();
            b = b->get_heavy_atom();
            r = a->distance_to(b);

            if ((threshold < 0 && f <= threshold) || (threshold > 0 && f >= threshold))
            {
                cout << strandid1 << ":" << *aa1;
                if (a) cout << ":" << a->name;
                cout << "-" << strandid2 << ":" << *aa2;
                if (b) cout << ":" << b->name;
                cout << ": " << r << " Å; " << f << " kJ/mol." << endl;
            }
        }
    }

    cout << "Total: " << ttl << " kJ/mol." << endl;
}