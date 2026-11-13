#include "../classes/failure_reporting.h"
#include <cassert>

int main()
{
    DockResult result;
    result.kJmol = -10;
    result.proximity = 5;
    result.worst_nrg_aa = 20;
    result.ligand_pocket_occlusion = 0.5;
    assert(output_filter_reasons(result,-10,5,20,true,0.5,1,"kJ/mol").empty());
    result.kJmol = -9;
    result.proximity = 6;
    result.worst_nrg_aa = 21;
    result.ligand_pocket_occlusion = 0.4;
    const auto reasons = output_filter_reasons(result,-10,5,20,true,0.5,1,"kJ/mol");
    for (const char* expected : {"Ligand binding energy -9 exceeds energy limit -10 kJ/mol.",
         "Ligand too far from pocket center: 6 A exceeds limit 5 A.",
         "Residue clashes too great: 21 exceeds limit 20 kJ/mol.",
         "Insufficient occlusion: 0.4 below limit 0.5."})
        assert(reasons.find(expected) != std::string::npos);
    assert(output_filter_reasons(result,-10,5,20,false,0.5,1,"kJ/mol").find("occlusion") == std::string::npos);
    const auto kcal_reasons = output_filter_reasons(result,-10,5,20,true,0.5,0.25,"kcal/mol");
    assert(kcal_reasons.find("Residue clashes too great: 5.25 exceeds limit 5 kcal/mol.") != std::string::npos);

    result.disqualify_reason = "ligand internal clashes too great 12.5 (C1...O3). ";
    append_failure_reason(result.disqualify_reason,reasons);
    const auto original = result.disqualify_reason;
    complete_failure_reason(result,"Different candidate failed cavity containment. ",2,false);
    assert(result.disqualify_reason == original);

    Protein protein;
    assert(protein.add_sequence("AAA"));
    Molecule ligand;
    ligand.from_smiles("CCO",false);
    AminoAcid* residues[] = {nullptr};
    DockResult exported(&protein,&ligand,0,residues);
    exported.disqualified = true;
    exported.disqualify_reason = original;
    exported.include_pdb_data = false;
    std::ostringstream report;
    report << exported;
    assert(report.str().find("Disqualified because: ligand internal clashes too great") != std::string::npos);
    assert(report.str().find("Residue clashes too great: 21") != std::string::npos);
    assert(report.str().find("No candidate passed all output filters.") == std::string::npos);
    assert(result.disqualify_reason.find("internal clashes too great") != std::string::npos);
    append_failure_reason(result.disqualify_reason,reasons);
    assert(result.disqualify_reason == original);

    DockResult fallback;
    complete_failure_reason(fallback,original,2,true);
    assert(fallback.disqualify_reason == original);
    complete_failure_reason(fallback,"Other reason. ",2,true);
    assert(fallback.disqualify_reason == original);
    DockResult no_score;
    complete_failure_reason(no_score,"",0,true);
    assert(no_score.disqualify_reason.find("No scored candidate with saved coordinates") != std::string::npos);
    assert(no_score.disqualify_reason.find("scored candidates: 0") != std::string::npos);
    result.kJmol = nan("");
    assert(output_filter_reasons(result,-10,5,20,false,0.5,1,"kJ/mol").find("Non-finite") != std::string::npos);
    puts("Specific rejection causes, limits, and internal-clash diagnostics survive fallback output.");
}
