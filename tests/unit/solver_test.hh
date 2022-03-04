/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef SOLVER_TEST_H
#define SOLVER_TEST_H

#include "test.hh"
#include "../../src/solver.hh"

class Solver_test : Test
{
private:
    int n;

    void test_solver();

    void run()
    {
        test_solver();
    }

public:
    Solver_test(int n, int tests);
};

#endif
