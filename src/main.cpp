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
	using Scalar = float;
	for (int test = 1; test <= LPSolverTest<Scalar>::get_num_tests(); ++test)
	{
		LPSolverTest<Scalar>(test).execute();
	}

	return 0;
}
