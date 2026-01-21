#!/bin/bash

RED="\033[31;1m"
GRN="\033[32m"
NC="\033[0m"

do_test()
{
    local test=$1

    test/$test | sed '/^#/d' > test/received/$test.received
    result=$(diff --unified test/approved/$test.approved test/received/$test.received)
    if [ -z "$result" ]; then
        printf "${GRN}\u2588${NC}"
    else
        printf "\n${RED}${test} FAILED:\n"
        diff --color --unified $REPORT test/approved/$test.approved test/received/$test.received | grep -E '^[+][^+]'
        printf "${NC}"
    fi
}

do_test "point_test"
do_test "ageo_test"
do_test "amino_test"
do_test "aniso_test"
do_test "bond_rotation_test"
do_test "conj_test"
do_test "eclipsing_test"
do_test "histidine_test" # check this, it's probably wrong
do_test "inte_test"
do_test "mcoord_test"
do_test "moiety_test"
do_test "mol_assem_test"
do_test "pi_stack_test" # check this, it's probably wrong
do_test "vdw_vertex_test"
do_test "probability_test"

# TODO: change these to not require a command line arg, else change the test framework to accommodate args.
# do_test "atom_test"
# do_test "backbone_test"
# do_test "bb_test"
# do_test "chirality_test"
# do_test "flexion_test"
# do_test "molecule_test"
# do_test "multimer_test"
# do_test "protein_test"
# do_test "ring_test"
# do_test "solvent_test"


printf "\n\n"