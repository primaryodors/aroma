#define main aromadock_cli_main
#include "../aromadock.cpp"
#undef main

#include <cassert>
#include <chrono>
#include <limits>
#include <sys/wait.h>

// Legacy momentum initialization reseeds from wall time. Freeze it in this executable.
extern "C" time_t time(time_t* result) noexcept
{
    time_t fixed = 1700000000;
    if (result) *result = fixed;
    return fixed;
}

struct CallbackMolecule : Molecule
{
    using Molecule::lmx;
    using Molecule::lmy;
    using Molecule::lmz;
};

static void predicate_checks()
{
    int skipped = 0;
    for (float attraction : {-100.0f,0.0f,10.0f,100.0f})
    for (float clash : {0.0f,1.0f,7.0f,8.0f})
    for (float potential : {-10.0f,0.0f,10.0f,100.0f})
    for (float missed : {0.0f,0.5f,1.0f})
    {
        Interaction score;
        score.attractive = attraction;
        score.clash = clash;
        float energy = score.summed();
        bool reference = clash > clash_limit_per_aa || energy > -0.3*potential
            || (clash > 0 && energy > 5) || missed >= 1;
        int calls = 0;
        bool actual = callback_flexion_unfavorable(score,missed,[&]() { calls++; return potential; });
        assert(actual == reference);
        if (!calls) skipped++;
    }
    assert(skipped > 0);
}

struct CallbackResult
{
    std::vector<Point> positions;
    std::vector<Point> initial;
    std::vector<MovabilityType> motions;
    std::vector<float> momentum;
    bool stopped;
    int next_random, calls;
    std::vector<int> random_trace;
    double seconds;
};

static CallbackResult run_callbacks(bool optimized, bool anchored, bool flexible)
{
    Protein receptor;
    CallbackMolecule odorant;
    FILE* input = fopen("pdbs/OR1/OR1A1.active.pdb","r");
    assert(input);
    receptor.load_pdb(input);
    fclose(input);
    srand(1397);
    if (anchored) odorant.from_smiles("CCO",false);
    else
    {
        std::ifstream sdf("sdf/d-limonene.sdf");
        assert(sdf);
        std::string contents((std::istreambuf_iterator<char>(sdf)),std::istreambuf_iterator<char>());
        odorant.from_sdf(contents.c_str());
    }
    Point center = receptor.get_residue(108)->get_CA_location().add(Point(4,0,0));
    odorant.recenter(center);
    protein = &receptor;
    ligand = &odorant;
    defer_callback_potential = optimized;
    flex = flexible;
    iters = 50;
    drift = initial_drift;
    loneliest = ligcen_target = center;
    search_size = Point(6,7.5,5.5);
    softness = 0;
    nsoftrgn = 0;
    nappears = 0;
    waters = nullptr;
    pivotal_hbond_aaa = nullptr;
    pivotal_hbond_la = nullptr;
    g_bbr[0].pri_res = nullptr;
    g_bbr[0].pri_tgt = nullptr;
    end_iterations = false;
    selfish_side_chains = false;
    output_each_iter = false;
    bridges.clear();
    forced_static_resnos.clear();
    worst_mol_clash = 0;
    minimum_searching_aniso = 0.5f;
    for (float& value : total_binding_by_type) value = 0;
    missed_connection = Vector();
    mc_bpotential = 0;

    AminoAcid* residues[SPHREACH_MAX+4];
    int count = receptor.get_residues_can_clash_ligand(residues,&odorant,center,Point(8,8,8),nullptr);
    std::vector<Molecule*> molecules{&odorant};
    for (int i=0; i<count; i++)
    {
        residues[i]->movability = MOV_FLXDESEL;
        molecules.push_back(residues[i]);
    }
    molecules.push_back(nullptr);
    if (anchored) odorant.stay_close_mine = odorant.get_atom(0);
    srand(89273);
    CallbackResult result;
    for (Molecule* molecule : molecules)
    {
        if (!molecule) break;
        for (int i=0; i<molecule->get_atom_count(); i++) result.initial.push_back(molecule->get_atom(i)->loc);
    }
    result.calls = 0;
    auto start = std::chrono::steady_clock::now();
    for (int step=1; step<=4; step++)
    {
        iteration_callback(step,molecules.data());
        result.random_trace.push_back(rand());
        result.calls++;
        if (end_iterations) break;
    }
    result.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    result.next_random = rand();
    result.stopped = end_iterations;
    result.momentum = {odorant.lmx,odorant.lmy,odorant.lmz};
    for (Molecule* molecule : molecules)
    {
        if (!molecule) break;
        result.motions.push_back(molecule->movability);
        for (int i=0; i<molecule->get_atom_count(); i++) result.positions.push_back(molecule->get_atom(i)->loc);
    }
    protein = nullptr;
    ligand = nullptr;
    minimum_searching_aniso = 0;
    return result;
}

template<class T>
static void transfer(int fd, T& value, bool writing)
{
    char* cursor = reinterpret_cast<char*>(&value);
    size_t remaining = sizeof(T);
    while (remaining)
    {
        ssize_t count = writing ? write(fd,cursor,remaining) : read(fd,cursor,remaining);
        assert(count > 0);
        cursor += count;
        remaining -= count;
    }
}

template<class T>
static void transfer(int fd, std::vector<T>& values, bool writing)
{
    size_t count = values.size();
    transfer(fd,count,writing);
    if (!writing) values.resize(count);
    for (T& value : values) transfer(fd,value,writing);
}

static CallbackResult isolated_run(bool optimized, bool anchored, bool flexible)
{
    int descriptors[2];
    assert(pipe(descriptors) == 0);
    pid_t child = fork();
    assert(child >= 0);
    bool writing = child == 0;
    close(descriptors[writing ? 0 : 1]);
    int fd = descriptors[writing ? 1 : 0];
    CallbackResult result;
    if (writing) result = run_callbacks(optimized,anchored,flexible);
    transfer(fd,result.positions,writing);
    transfer(fd,result.initial,writing);
    transfer(fd,result.motions,writing);
    transfer(fd,result.momentum,writing);
    transfer(fd,result.stopped,writing);
    transfer(fd,result.next_random,writing);
    transfer(fd,result.calls,writing);
    transfer(fd,result.seconds,writing);
    transfer(fd,result.random_trace,writing);
    close(fd);
    if (writing) _exit(0);
    int status;
    assert(waitpid(child,&status,0) == child);
    assert(WIFEXITED(status) && WEXITSTATUS(status) == 0);
    return result;
}

static void compare_callbacks(bool anchored, bool flexible)
{
    srand(413);
    CallbackResult reference = isolated_run(false,anchored,flexible);
    srand(413);
    CallbackResult optimized = isolated_run(true,anchored,flexible);
    assert(reference.initial.size() == optimized.initial.size());
    for (size_t i=0; i<reference.initial.size(); i++)
        assert(reference.initial[i].get_3d_distance(optimized.initial[i]) < 0.000001f);
    assert(reference.calls == optimized.calls);
    assert(reference.stopped == optimized.stopped);
    if (reference.next_random != optimized.next_random)
    {
        for (size_t i=0; i<reference.random_trace.size(); i++)
            fprintf(stderr,"Step %zu: %d vs %d\n",i,reference.random_trace[i],optimized.random_trace[i]);
        fprintf(stderr,"Callback random mismatch (anchor %d, flex %d): %d vs %d\n",
            anchored,flexible,reference.next_random,optimized.next_random);
    }
    assert(reference.next_random == optimized.next_random);
    assert(reference.motions == optimized.motions);
    assert(reference.positions.size() == optimized.positions.size());
    for (size_t i=0; i<reference.positions.size(); i++)
    {
        float error = reference.positions[i].get_3d_distance(optimized.positions[i]);
        if (error >= 0.0001f) fprintf(stderr,"Callback geometry mismatch at atom %zu: %g\n",i,error);
        assert(error < 0.0001f);
    }
    for (size_t i=0; i<reference.momentum.size(); i++)
        assert(fabs(reference.momentum[i]-optimized.momentum[i]) < 0.0001f*(1+fabs(reference.momentum[i])));
    printf("Callback comparison (anchor %d, flex %d): %d calls, %.3fs vs %.3fs (%.2fx)\n",
        anchored,flexible,reference.calls,reference.seconds,optimized.seconds,reference.seconds/optimized.seconds);
}

int main()
{
    predicate_checks();
    compare_callbacks(false,true);
    compare_callbacks(true,true);
    compare_callbacks(false,false);
    puts("Callback decisions, random sequence, geometry, and momentum match the reference.");
}
