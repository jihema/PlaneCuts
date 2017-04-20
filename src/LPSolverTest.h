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
	VecX m_known_solution;
	double m_known_best_value;
};

