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

int main(int argc, const char * argv[])
{

	using MatXd = LinearProgrammingSolver<double>::MatX;
	using VecX = LinearProgrammingSolver<double>::VecX;

	static int d = 3; // Dimension of space.
	static int const n0 = 2; // Number of hyperplanes.

	MatXd A_prime(n0, d);
	A_prime << 3, 2, 1, 2, 5, 3;

	VecX b(n0);
	b << 10, 15;

	VecX c(d);
	c << -2, -3, -4;

	/*
	 static int d = 5; // Dimension of space.
	 static int const n0 = 2; // Number of hyperplanes.

	 MatXd A_prime(n0, d);
	 A_prime <<  3, 2, 1, 1, 0, 2, 5, 3, 0, 1;

	 VecX b(n0);
	 b <<  10, 15;

	 VecX c(d);
	 c << -2,-3,-4,0,0;
	 */

	auto slps = SimplexLPSolver<double>(
			LinearProgrammingSolver<double>(A_prime, b, c).get_tableau());

	slps.solve();

	std::cout << "Optimal value = " << slps.get_optimal_value() << '\n';
	std::cout << "Solution:\n" << slps.get_solution() << '\n';

	return 0;
}
