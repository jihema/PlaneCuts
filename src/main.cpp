//
//  main.cpp
//  AvisFukuda
//
//  Created by Jean-Marie Aubry on 25/03/2017.
//  Copyright © 2017 Jean-Marie Aubry. All rights reserved.
//

#include "SimplexLPSolver.h"
#include "LPSolverTest.h"

int main(int argc, const char * argv[])
{
    using Scalar = double;
    bool okay = true;
    for (int test = 1; test <= 9; ++test)
    {
        okay = LPSolverTest<Scalar>(test).execute() && okay;
        std::cout << '\n';
    }

    if (okay)
    {
        std::cout << "All OK\n";
    }

    return 0;
}
