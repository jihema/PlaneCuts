/*
 * LPSolverTest.h
 *
 *  Created on: Apr 20, 2017
 *      Author: JM
 */

#pragma once

#include "LinearProgrammingSolver.h"

inline double sqr(double x)
{
	return x * x;
}

class LPSolverTest {
public:
	LPSolverTest(int test);

	static int get_num_tests()
	{
		return num_tests;
	}

	template<typename Solver>
	void execute();

private:
	using MatXd = LinearProgrammingSolver<double>::MatX;
	using VecX = LinearProgrammingSolver<double>::VecX;

	void resize(int d, int n);

	int d; // Number of variables.
	int n; // Number of constraints.
	int m_sign;
	int test_id;

	MatXd m_A;
	VecX m_b;
	VecX m_c;
	VecX m_inequalities;
	VecX m_known_solution;

	static const int num_tests = 6;
	static const double Stigler_data[77][9];
	static const double Stigler_nutrients[9];
	static const double Stigler_solution[77];
};

