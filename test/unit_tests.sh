#!/bin/bash
set -u
set -o pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
mkdir -p test/received || exit 1
make tests aromadock || exit 1

failures=0
approval_tests=(point ageo amino aniso bond_rotation conj eclipsing histidine inte mcoord moiety mol_assem pi_stack vdw_vertex probability)
for name in "${approval_tests[@]}"; do
    test_name="${name}_test"
    output="test/received/$test_name.received"
    if ! "test/$test_name" | sed '/^#/d' > "$output"; then
        printf '%s: execution failed\n' "$test_name"
        failures=$((failures+1))
        continue
    fi
    if python3 test/compare_output.py "test/approved/$test_name.approved" "$output"; then
        printf '%s: passed\n' "$test_name"
    else
        printf '%s: output comparison failed\n' "$test_name"
        failures=$((failures+1))
    fi
done

for name in spatial pose_state allocation result_state numeric clash_selection entropy pocket_search smiles_input cavity_residue soft_scoring protonation disulfide atom_target soft_disulfide target_rotation primary_target soft_covalent repack_state cavity_intersection cavity_storage space_containment space_atom; do
    if ! "test/${name}_test"; then
        printf '%s: regression failed\n' "$name"
        failures=$((failures+1))
    fi
done

if ! python3 -m unittest discover -s test -p '*_test.py'; then
    failures=$((failures+1))
fi
if command -v node >/dev/null 2>&1; then
    if ! node test/viewer_client_test.js; then failures=$((failures+1)); fi
fi
printf '%d failing checks\n' "$failures"
if (( failures > 0 )); then exit 1; fi
