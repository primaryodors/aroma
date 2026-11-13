#ifndef AROMA_FAILURE_REPORTING_H
#define AROMA_FAILURE_REPORTING_H

#include "scoring.h"
#include <sstream>

inline void append_failure_reason(std::string& reasons, const std::string& reason)
{
    if (!reason.empty() && reasons.find(reason) == std::string::npos)
    {
        if (!reasons.empty() && reasons.back() != ' ') reasons += ' ';
        reasons += reason;
    }
}

inline std::string output_filter_reasons(const DockResult& result, float energy_limit,
    float proximity_limit, float residue_clash_limit, bool check_occlusion,
    float occlusion_limit, float energy_multiplier, const char* energy_units)
{
    std::ostringstream reasons;
    if (!std::isfinite(result.kJmol))
        reasons << "Non-finite ligand binding energy. ";
    else if (result.kJmol > energy_limit)
        reasons << "Ligand binding energy " << result.kJmol*energy_multiplier
            << " exceeds energy limit " << energy_limit*energy_multiplier << " " << energy_units << ". ";
    if (result.proximity > proximity_limit)
        reasons << "Ligand too far from pocket center: " << result.proximity
            << " A exceeds limit " << proximity_limit << " A. ";
    if (result.worst_nrg_aa > residue_clash_limit)
        reasons << "Residue clashes too great: " << result.worst_nrg_aa*energy_multiplier
            << " exceeds limit " << residue_clash_limit*energy_multiplier << " " << energy_units << ". ";
    if (check_occlusion && result.ligand_pocket_occlusion < occlusion_limit)
        reasons << "Insufficient occlusion: " << result.ligand_pocket_occlusion
            << " below limit " << occlusion_limit << ". ";
    return reasons.str();
}

inline void complete_failure_reason(DockResult& result, const std::string& candidate_reasons,
    int scored_candidates, bool initial_placement)
{
    if (result.disqualify_reason.empty()) append_failure_reason(result.disqualify_reason,candidate_reasons);
    if (result.disqualify_reason.empty())
    {
        std::ostringstream reason;
        if (initial_placement)
            reason << "No scored candidate with saved coordinates was available; initial placement exported";
        else if (result.pose <= 0)
            reason << "Candidate did not receive an output pose rank";
        else if (result.pdbdat.empty())
            reason << "Candidate has no saved PDB data";
        else
            reason << "Candidate path did not complete output selection (pose rank " << result.pose << ")";
        reason << "; scored candidates: " << scored_candidates << ". ";
        result.disqualify_reason = reason.str();
    }
}

#endif
