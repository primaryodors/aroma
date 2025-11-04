#include "appear.h"

void Appear::update(Protein *prot, int effective_iter)
{
    if (effective_iter != disappear_iter && effective_iter != reappear_iter) return;
    int i, j, m, n = prot->get_end_resno();
    Vector ldisp;
    start.resolve_resno(prot);
    end.resolve_resno(prot);

    #if _dbg_appears
    cout << "Updating " << start.resno << "-" << end.resno << endl << endl;
    #endif

    if (displacement)
    {
        ldisp = prot->get_region_center(start.resno, end.resno).subtract(dispcen);
        ldisp.r = displacement;
    }

    for (i = start.resno; i <= end.resno; i++)
    {
        AminoAcid* aa = prot->get_residue(i);
        if (!aa) continue;

        if (effective_iter == disappear_iter)
        {
            m = aa->get_atom_count();
            for (j=0; j<m; j++)
            {
                Atom* a = aa->get_atom(j);
                a->vanished = true;
            }

            if (aa->coordmtl) aa->coordmtl->vanished = true;
        }
        else
        {
            m = aa->get_atom_count();
            for (j=0; j<m; j++)
            {
                Atom* a = aa->get_atom(j);
                a->vanished = false;
            }

            if (aa->coordmtl) aa->coordmtl->vanished = false;

            if (displacement) 
            {
                MovabilityType mt = aa->movability;
                aa->movability = MOV_ALL;
                aa->aamove(ldisp);
                aa->movability = mt;
            }
        }
    }
}