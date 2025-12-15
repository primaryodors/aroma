
#include "activation.h"

void Activation::load_acvm_file(AcvType acvt)
{
    char infname[256];
    switch (acvt)
    {
        case acv_GPCR:
        strcpy(infname, "data/gpcr.acvm");
        break;

        default:
        cerr << "Unknown activation type." << endl;
        throw 0xbadc0de;
    }

    FILE* fp = fopen(infname, "rb");
    if (!fp)
    {
        cerr << "Failed to open " << infname << " for reading." << endl;
        throw 0xbadf12e;
    }

    char buffer[1024];
    m_acvm = new ActiveMotion[256];
    while (!feof(fp))
    {
        char* ln = fgets(buffer, 1020, fp);
        if (!ln || feof(fp)) break;
        char* hash = strchr(buffer, '#');
        if (hash) *hash = 0;
        char** fields = chop_spaced_words(buffer);
        if (!fields[0] || !strlen(fields[0])) continue;

        if (!strcmp(fields[0], "PIVOT"))                    // https://www.youtube.com/watch?v=n67RYI_0sc0
        {
            m_acvm[nacvm].rap_start.set(fields[1]);
            m_acvm[nacvm].rap_end.set(fields[2]);
            m_acvm[nacvm].rap_fulcrum.set(fields[3]);
            m_acvm[nacvm].rap_index.set(fields[4]);
            m_acvm[nacvm].tgtdist = atof(fields[5]);
            m_acvm[nacvm].rap_target.set(fields[6]);
            nacvm++;
        }
    }
    fclose(fp);
}

void Activation::apply(Protein *p)
{
    int i;
    for (i=0; i<nacvm; i++) m_acvm[i].apply(p);
}

void ActiveMotion::apply(Protein *p)
{
    if (acvmt == acvm_pivot)
    {
        rap_start.resolve_resno(p);
        rap_end.resolve_resno(p);
        rap_fulcrum.resolve_resno(p);
        rap_index.resolve_resno(p);
        rap_target.resolve_resno(p);

        Point pt_fulcrum = rap_fulcrum.loc(), pt_index = rap_index.loc(), pt_target = rap_target.loc();
        Vector tolerance = pt_index.subtract(pt_target);
        tolerance.r = tgtdist;
        pt_target = pt_target.add(tolerance);
        Rotation rot = align_points_3d(pt_index, pt_target, pt_fulcrum);

        p->rotate_piece(rap_start.resno, rap_end.resno, pt_fulcrum, rot.v, rot.a);
    }
}
