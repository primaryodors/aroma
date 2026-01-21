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
        printf "\n${RED}${test} FAILED:${NC}\n"
        diff --color --unified $REPORT test/approved/$test.approved test/received/$test.received
    fi
}

test="point_test"
do_test $test


printf "\n\n"