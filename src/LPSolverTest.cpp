/*
 * LPSolverTest.cpp
 *
 *  Created on: Apr 20, 2017
 *      Author: JM
 */

#include "LPSolverTest.h"
#include "SimplexLPSolver.h"

LPSolverTest::LPSolverTest(int test) :
		test_id(test)
{
	switch (test) {

	case 1:
		resize(5, 2);
		m_sign = 1;
		m_A << .3, .2, .1, .1, 0, 2, 5, 3, 0, 1;
		m_b << 1.0, 15;
		m_c << -2, -3, -4, 0, 0;
		m_known_solution << 0, 0, 5, 5, 0;
		m_known_best_value = -20;
		break;

	case 2:
		resize(3, 2);
		m_sign = 1;
		m_A << 3, 2, 1, 2, 5, 3;
		m_b << 10, 15;
		m_c << -2, -3, -4;
		m_known_solution << 15. / 7., 0, 25. / 7.;
		m_known_best_value = -130. / 7.;
		break;

	case 3:
		resize(6, 3);
		m_sign = -1;
		m_A << 2, 1, 1, 1, 0, 0, 1, 2, 3, 0, 1, 0, 2, 2, 1, 0, 0, 1;
		m_b << 2, 5, 6;
		m_c << 3, 2, 3, 0, 0, 0;
		m_known_solution << 1. / 5., 0, 8. / 5., 0, 0, 4;
		m_known_best_value = 27. / 5.;
		break;

	case 4:
		resize(3, 2);
		m_sign = -1;
		m_A << 2, 1, -1, 1, 2, 0;
		m_b << 4, 6;
		m_c << 1, 1, 0;
		m_known_solution << 6, 0, 8;
		m_known_best_value = 6;
		break;

	case 5:
		resize(4, 2);
		m_sign = 1;
		m_A << 1, 1, 0, -1, 0, -1, 1, -1;
		m_b << -1, -3;
		m_c << 2, 3, -1, 1;
		m_known_solution << 0, 1, 0, 2;
		m_known_best_value = 17;
		break;
	}

}

void LPSolverTest::resize(int d, int n)
{
	m_A.resize(n, d);
	m_b.resize(n);
	m_c.resize(d);
	m_known_solution.resize(d);
}

template<>
void LPSolverTest::execute<SimplexLPSolver<double>>()
{
	auto slps =
			SimplexLPSolver<double>(m_A, m_b, m_sign * m_c);

	slps.solve();

	double constexpr sq_tolerance = std::numeric_limits<double>::epsilon();

	std::cout << "Test " << test_id << '\n';
	std::cout << "Solution:\t" << slps.get_solution().transpose() << '\n';
	std::cout << "Known solution:\t" << m_known_solution.transpose() << '\n';
	std::cout << "Optimal value:\t" << m_sign * slps.get_optimal_value()
			<< '\n';
	std::cout << "Known optimal:\t" << m_known_best_value << '\n';

	if ((slps.get_solution() - m_known_solution).squaredNorm() <= sq_tolerance
			&& sqr(m_sign * slps.get_optimal_value() - m_known_best_value)
					<= sq_tolerance)
	{
		std::cout << "OK\n";
	} else
	{
		std::cout << "FAILED\n";
	}

	std::cout << '\n';
}
