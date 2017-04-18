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
	A_prime << 3, 2, 1, 2, 5, 0;

	VecX b(n0);
	b << 10, 15;

	VecX c(d);
	c << -2, -3, -4;

	/*
	 static int  d = 6; // Dimension of space.
	 static int const n0 = 3; // Number of hyperplanes.

	 // We consider the convex polyhedron defined by: A' y \leq b.
	 MatXd A_prime;
	 A_prime.resize(n0, d);
	 VecX b;
	 b.resize(n0);

	 A_prime << 1,2,3,4,0,0,
	 0,3,2,1,1,0,
	 0,2,5,3,0,1;
	 b << 0,10,15;

	 VecX c;
	 c.resize(d);
	 c << 0,0,0,0,1,1;
	 */

	auto slps = SimplexLPSolver<double>(
			LinearProgrammingSolver<double>(A_prime, b, c).get_tableau());
	slps.make_canonical();

//	for (int i = 1; i <= n0; ++i)
//	{
//		m_b[i - 1] = i;
//		for (int j = 1; j <= d; ++j)
//		{
//			m_A_prime(i - 1, j - 1) = i * j;
//		}
//	}

//	PlaneCuts<double, d> plane_cuts(A_prime, b);
//
//	AvisFukudaSolver<double, d> afs(plane_cuts);
//
//	std::cout << afs << '\n';
//
//	for (int i = 0; i < 5; ++i)
//	{
//		afs.pivot_criss_cross();
//		std::cout << afs << '\n';
//	}

//	Polytope<double, d> polygon = avis_fukuda_solver.get_polygon();

	return 0;
}
