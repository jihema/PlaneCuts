//
//  main.cpp
//  AvisFukuda
//
//  Created by Jean-Marie Aubry on 25/03/2017.
//  Copyright © 2017 Jean-Marie Aubry. All rights reserved.
//

#include <iostream>

#include "AvisFukudaSolver.h"
#include "SimplexLPSolver.h"
#include "LPSolverTest.h"

int main(int argc, const char * argv[])
{
	for (int test = 1; test <= 5; ++test)
	{
		LPSolverTest(test).execute<SimplexLPSolver<double>>();
	}

	return 0;
}
