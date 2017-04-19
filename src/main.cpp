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

double sqr(double x)
{
	return x * x;
}

int main(int argc, const char * argv[])
{
	using MatXd = LinearProgrammingSolver<double>::MatX;
	using VecX = LinearProgrammingSolver<double>::VecX;

	MatXd A_prime;
	VecX b;
	VecX c;
	VecX known_solution;
	double known_best_value;

	int d; // Number of variables.
	int n; // Number of constraints.
	int sign;

	for (int test = 1; test <= 5; ++test)
	{

		switch (test) {

		case 1:
			d = 5;
			n = 2;
			sign = 1;
			A_prime.resize(n, d);
			A_prime << .3, .2, .1, .1, 0, 2, 5, 3, 0, 1;
			b.resize(n);
			b << 1.0, 15;
			c.resize(d);
			c << -2, -3, -4, 0, 0;

			known_solution.resize(d);
			known_solution << 0, 0, 5, 5, 0;
			known_best_value = -20;
			break;

		case 2:
			d = 3;
			n = 2;
			sign = 1;

			A_prime.resize(n, d);
			A_prime << 3, 2, 1, 2, 5, 3;
			b.resize(n);
			b << 10, 15;
			c.resize(d);
			c << -2, -3, -4;

			known_solution.resize(d);
			known_solution << 15. / 7., 0, 25. / 7.;
			known_best_value = -130. / 7.;

			break;

		case 3: // Okay.
			d = 6;
			n = 3;
			sign = -1;

			A_prime.resize(n, d);
			A_prime << 2, 1, 1, 1, 0, 0, 1, 2, 3, 0, 1, 0, 2, 2, 1, 0, 0, 1;
			b.resize(n);
			b << 2, 5, 6;
			c.resize(d);
			c << 3, 2, 3, 0, 0, 0;

			known_solution.resize(d);
			known_solution << 1. / 5., 0, 8. / 5., 0, 0, 4;
			known_best_value = 27. / 5.;
			break;

		case 4: // Okay.
			d = 3;
			n = 2;
			sign = -1;
			A_prime.resize(n, d);
			A_prime << 2, 1, -1, 1, 2, 0;
			b.resize(n);
			b << 4, 6;
			c.resize(d);
			c << 1, 1, 0;

			known_solution.resize(d);
			known_solution << 6, 0, 8;
			known_best_value = 6;
			break;

		case 5:
			d = 4;
			n = 2;
			sign = 1;

			A_prime.resize(n, d);
			A_prime << 1, 1, 0, -1, 0, -1, 1, -1;
			b.resize(n);
			b << -1, -3;
			c.resize(d);
			c << 2, 3, -1, 1;

			known_solution.resize(d);
			known_solution << 0, 1, 0, 2;
			known_best_value = 17;
			break;
		}

		auto slps =
				SimplexLPSolver<double>(
						LinearProgrammingSolver<double>(A_prime, b, sign * c).get_tableau());

		slps.solve();

		double constexpr sq_tolerance = std::numeric_limits<double>::epsilon();

		std::cout << "Test " << test << '\n';
		std::cout << "Solution:\t" << slps.get_solution().transpose() << '\n';
		std::cout << "Known solution:\t" << known_solution.transpose() << '\n';
		std::cout << "Optimal value:\t" << sign * slps.get_optimal_value()
				<< '\n';
		std::cout << "Known optimal:\t" << known_best_value << '\n';

		if ((slps.get_solution() - known_solution).squaredNorm() <= sq_tolerance
				&& sqr(sign * slps.get_optimal_value() - known_best_value)
						<= sq_tolerance)
		{
			std::cout << "OK\n";
		} else
		{
			std::cout << "FAILED\n";
		}
		std::cout << '\n';
	}

	return 0;
}
