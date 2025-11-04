
#include "classes/search.h"
#include "classes/cavity.h"
#include "classes/progress.h"

int main(int argc, char** argv)
{
    Protein p("Test Receptor");
    Protein* protein = &p;
    Molecule m("Test Ligand");
    Molecule* ligand = &m;
    Cavity cvty;
    bool priorities[256];
    std::string outfname = "cavfit.pdb";

    bool save_tmp_pdbs = false;

    int i, j, l, n, nligconf = 8192, iters = 200;
    float minimal_fit_threshold = 0.25, reasonable_fit_threshold = 0.8;
    for (i=0; i<256; i++) priorities[i] = false;

    FILE* fp;
    char buffer[49152];
    int got;
    for (i=1; i<argc; i++)
    {
        char* dot;
        dot = strrchr(argv[i], '.');
        if (dot && !strcasecmp(dot, ".pdb"))
        {
            fp = fopen(argv[i], "rb");
            if (!fp)
            {
                cout << "Failed to open " << argv[i] << " for reading." << endl;
                return -1;
            }
            p.load_pdb(fp);
            fclose(fp);
            cout << "Loaded " << argv[i] << endl;
        }
        else if (dot && !strcasecmp(dot, ".sdf"))
        {
            fp = fopen(argv[i], "rb");
            if (!fp)
            {
                cout << "Failed to open " << argv[i] << " for reading." << endl;
                return -1;
            }
            got = fread(buffer, 1, 49150, fp);
            fclose(fp);
            m.from_sdf(buffer);
            cout << "Ligand molecular weight: " << m.get_molecular_wt() << endl;
        }
        else if (dot && !strcasecmp(dot, ".cvty"))
        {
            if (!file_exists(argv[i]))
            {
                cout << "ERROR file not found: " << argv[i] << endl;
                return -7;
            }
            FILE* fp = fopen(argv[i], "rb");
            if (!fp)
            {
                cout << "FAILED to open " << argv[i] << " for reading." << endl;
                return -7;
            }

            // Load all partials from the cavity file into a single cavity object.
            cout << "Reading " << argv[i] << "..." << endl;
            char buffer[1024];
            while (!feof(fp))
            {
                char* got = fgets(buffer, 1022, fp);
                if (!got) break;
                CPartial cp;
                int cno = cp.from_cvty_line(buffer);
                cvty.add_partial(cp);
            }
            fclose(fp);
            cout << "Read cavities." << endl;
        }
        else if (!strcmp(argv[i], "--bsr"))
        {
            for (i++; argv[i] && argv[i][0] >= '0' && argv[i][0] <= '9'; i++)
            {
                dot = strrchr(argv[i], '.');
                AminoAcid* aa;
                if (dot) aa = p.get_residue_bw(argv[i]);
                else aa = p.get_residue(atoi(argv[i]));
                if (aa) aa->priority = true;
            }
            if (argv[i]) i--;
        }
        else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--out"))
        {
            i++;
            outfname = argv[i];
        }
        else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--confs"))
        {
            i++;
            nligconf = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--reasonable"))
        {
            i++;
            reasonable_fit_threshold = atoi(argv[i]);
        }
        else cout << "Warning: unknown command argument " << argv[i] << endl;
    }

    if (!m.get_atom_count())
    {
        cerr << "ERROR: No ligand." << endl;
        return -2;
    }

    if (!cvty.count_partials())
    {
        cerr << "ERROR: No cavities." << endl;
        return -2;
    }

    cout << "Total cavity volume: " << cvty.get_volume() << " Å³" << endl;
    cout << "Ligand volume: " << m.get_volume() << " Å³" << endl;

    if (!p.get_end_resno())
    {
        cout << "Cavity fit search requires a protein PDB." << endl;
        return 0;
    }


    // TODO: Facility for metal coordination.

    // Generate a large number of randomized ligand conformers.
    Pose  ligconf[nligconf+2];
    bool  cfvalid[nligconf+2];
    float cfscore[nligconf+2];

    Bond** rotb = m.get_rotatable_bonds(true);

    cout << "Generating conformers..." << endl << endl << flush;
    Progressbar pb;
    pb.set_color(192, 64, 224);
    pb.minimum = 0;
    pb.maximum = nligconf-1;
    for (i=0; i<nligconf; i++)
    {
        // Straighten the ligand.
        m.minimize_internal_clashes();

        // Center the conformer in a random partial.
        j = rand() % cvty.count_partials();
        Point cen = cvty.get_partial_by_idx(j)->s.center;
        m.recenter(cen);

        // Crumple some of the conformers.
        if (frand(0,1) < 0.01) m.crumple(frand(0, hexagonal));

        // Perform random rotations.
        Vector axis = Point(1, 0, 0);
        m.rotate(axis, frand(-M_PI, M_PI));
        axis = Point(0, 1, 0);
        m.rotate(axis, frand(-M_PI, M_PI));
        axis = Point(0, 0, 1);
        m.rotate(axis, frand(-M_PI, M_PI));

        // Perform random flexions sometimes.
        if (frand(0,1) < 0.5)
        {
            if (rotb)
            {
                for (j=0; rotb[j]; j++)
                {
                    if (frand(0,1) > 0.03) continue;
                    float theta = pow(frand(-hexagonal, hexagonal), 2);
                    if (frand(0,1) < 0.25) theta += M_PI;
                    rotb[j]->rotate(theta);
                }
            }
        }

        // Perform random ring flips.
        n = m.get_atom_count();
        for (j=0; j<n; j++)
        {
            if (frand(0,1) > 0.03) continue;
            Atom* a = m.get_atom(j);
            if (!a) continue;
            if (!a->num_rings()) continue;
            Ring** rr = a->get_rings();
            if (rr && rr[0]) rr[0]->flip_atom(a);
            delete rr;
        }

        // Store the conformation.
        ligconf[i].copy_state(&m);

        float f = cvty.molecule_inside_pocket(&m);
        cfvalid[i] = (f >= minimal_fit_threshold);

        if (frand(0,1) < 0.1) pb.update(i);
    }
    pb.erase();

    // Narrow down the ligand conformers to those that meet a minimal threshold of molecule_inside_pocket().
    cout << "Filtering by goodness of fit..." << endl << endl << flush;
    pb.set_color(160, 64, 224);
    j = 0;
    for (i=0; i<nligconf; i++)
    {
        if (cfvalid[i])
        {
            cfvalid[j] = true;
            ligconf[j] = ligconf[i];          // We can get away with this because j will never exceed i.
            j++;
        }
        pb.update(i);
    }
    pb.erase();

    int nmatches = j;
    for (; j<nligconf; j++) cfvalid[j] = false;
    cout << endl;

    n = cvty.count_partials();
    CPartial* polar_partials[n+4];
    j = 0;
    for (i=0; i<n; i++)
    {
        CPartial* part = cvty.get_partial_by_idx(i);
        if (!part) break;
        if (part->polar) polar_partials[j++] = part;
    }
    int npolarpart = j;

    n = m.get_atom_count();
    Atom* polar_atoms[n+4];
    j = 0;
    for (i=0; i<n; i++)
    {
        Atom* a = m.get_atom(i);
        if (!a) break;
        if (fabs(a->is_polar()) >= hydrophilicity_cutoff) polar_atoms[j++] = a;
    }
    int npolarat = j;

    // Perform translations, rotations, and flexions to try to optimize each conformer's cavity fit.
    cout << "Optimizing fit..." << endl << endl << flush;
    pb.set_color(128, 64, 224);
    pb.maximum = nmatches;
    for (i=0; i<nmatches; i++)
    {
        ligconf[i].restore_state(&m);
        for (l=0; l<iters; l++)
        {
            Vector mov = Point(frand(-1, 1), frand(-1, 1), frand(-1, 1));

            if (npolarpart && npolarat)
            {
                int ip = rand() % npolarpart, ia = rand() % npolarat;
                Point ptp = polar_partials[ip]->s.center;
                Point pta = polar_atoms[ia]->loc;
                Point ptc = m.get_barycenter();

                float angle = find_3d_angle(ptp, pta, ptc);
                if (angle <= square)
                {
                    mov = ptp.subtract(pta);
                    mov.r /= 2.5;
                    
                    Rotation rot = align_points_3d(pta, ptp, ptc);
                    m.rotate(&rot.v, rot.a, false);
                }
            }

            float before = cvty.molecule_inside_pocket(&m, true);
            m.move(mov);
            float after = cvty.molecule_inside_pocket(&m, true);
            if (after < before)
            {
                mov.r *= -2;
                m.move(mov);
                after = cvty.molecule_inside_pocket(&m, true);
                if (after < before)
                {
                    mov.r *= -0.5;
                    m.move(mov);
                }
            }

            Vector axis;
            std::vector<Atom*> lda;
            if (frand(0,1) < 0.2) lda = m.longest_dimension();
            if (lda.size() > 1) axis = lda[1]->loc.subtract(lda[0]->loc);
            else axis = Point(frand(-1, 1), frand(-1, 1), frand(-1, 1));
            float theta = frand(-hexagonal, hexagonal);

            before = cvty.molecule_inside_pocket(&m, true);
            m.rotate(axis, theta);
            after = cvty.molecule_inside_pocket(&m, true);
            if (after < before)
            {
                m.rotate(axis, -theta*2);
                after = cvty.molecule_inside_pocket(&m, true);
                if (after < before)
                {
                    m.rotate(axis, theta);
                }
            }

            if (rotb)
            {
                for (n=0; rotb[n]; n++);       // count.
                for (j=0; j<n; j++)
                {
                    if (frand(0,1) > 0.02) continue;
                    float theta = pow(frand(-hexagonal, hexagonal), 3);
                    if (frand(0,1) < 0.25) theta += M_PI;

                    before = cvty.molecule_inside_pocket(&m, true) - m.total_eclipses()/10;
                    rotb[j]->rotate(theta);
                    after = cvty.molecule_inside_pocket(&m, true) - m.total_eclipses()/10;

                    if (after < before)
                    {
                        rotb[j]->rotate(-theta*2);
                        after = cvty.molecule_inside_pocket(&m, true) - m.total_eclipses()/10;

                        if (after < before) rotb[j]->rotate(theta);
                    }
                }
            }
        }

        ligconf[i].copy_state(&m);
        if (frand(0,1) < 0.5) pb.update(i);
    }
    pb.erase();

    int mostfit = 0;

    // Narrow the conformers further to those that meet a reasonable threshold.
    cout << "Selecting output conformers..." << endl << endl << flush;
    pb.set_color(96, 64, 224);
    j = 0;
    for (i=0; i<nmatches; i++)
    {
        ligconf[i].restore_state(&m);
        float f = cvty.molecule_inside_pocket(&m, true);
        f -= 0.0001 * (m.get_internal_clashes() + m.total_eclipses());
        if (f >= reasonable_fit_threshold)
        {
            cfvalid[j] = true;
            ligconf[j] = ligconf[i];
            cfscore[j] = f;
            j++;

            int mult = cvty.estimate_multiplicity(ligand);
            if (mult > mostfit) mostfit = mult;
        }
        pb.update(i);
    }
    pb.erase();
    nmatches = j;

    cout << "Cavity can hold up to " << mostfit << " ligands." << endl;

    // Sort the conformers by their cavity fit scores.
    cout << "Sorting..." << endl << endl << flush;
    pb.set_color(64, 64, 224);
    pb.maximum = nmatches;
    for (i=0; i<nmatches; i++)
    {
        for (j=1; j<nmatches; j++)
        {
            if (cfscore[j-1] < cfscore[j])
            {
                float f = cfscore[j-1];
                Pose ps = ligconf[j-1];
                cfscore[j-1] = cfscore[j];
                ligconf[j-1] = ligconf[j];
                cfscore[j] = f;
                ligconf[j] = ps;
            }
        }

        if (frand(0,1) < 0.5) pb.update(i);
    }
    pb.erase();
    cout << endl;

    // for (i=0; i<nmatches; i++) cout << cfscore[i] << endl;
    cout << "Found " << nmatches << " matches." << endl;

    fp = fopen(outfname.c_str(), "w");
    if (!fp)
    {
        cout << "FAILED to open " << outfname << " for writing." << endl;
        return -3;
    }

    n = cvty.count_partials();
    for (i=0; i<n; i++)
    {
        CPartial* cp = cvty.get_partial_by_idx(i);
        strcpy(buffer, "REMARK 821 ");
        cp->write_cvty_line(buffer+11, 0, &p);
        fprintf(fp, "%s", buffer);
    }

    // Write the top 15 conformers as PDB data, giving each the letter code LIG and its own residue number.
    cout << "Writing data..." << endl << flush;
    l = 0;
    for (i=0; i<nmatches && i<15; i++)
    {
        fprintf(fp, "REMARK CFSCORE %f\n", cfscore[i]);
        ligconf[i].restore_state(&m);
        n = m.get_atom_count();
        for (j=0; j<n; j++)
        {
            Atom* a = m.get_atom(j);
            if (!a) continue;
            a->residue = i+1;
            if (!i) strcpy(a->aa3let, "LIG");
        }
        m.save_pdb(fp, l);
        l += n;
    }

    fclose(fp);
    cout << "Saved " << outfname << endl;
    return 0;
}
