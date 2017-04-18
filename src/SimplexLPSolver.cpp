/*
 * SimplexLPSolver.cpp
 *
 *  Created on: Apr 18, 2017
 *      Author: JM
 */

#include "SimplexLPSolver.h"
#include <iostream>
#include <map>
#include <assert.h>

template<typename Scalar>
SimplexLPSolver<Scalar>::SimplexLPSolver(MatX const& tablo) :
		m_tableau(tablo)
{
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::make_canonical()
{
	std::cout << "Tableau:\n" << m_tableau << '\n';

	search_basic_variables();

	price_out();
	std::cout << "Priced out tableau:\n" << m_tableau << '\n';

	if (is_canonical())
	{
		std::cout << "Canonical: basic variables are\n";
		print_basic_variables(std::cout);
	} else
	{
		std::cout << "Non-canonical: basic variables are\n";
		print_basic_variables(std::cout);

		std::cout << "Extending...\n";

		create_artificial_variables();

		std::cout << "After extension, basic variables are:\n";
		print_basic_variables(std::cout);

		std::cout << "Tableau:\n" << m_tableau << '\n';

		assert(is_canonical());

		price_out();
		std::cout << "Priced out tableau:\n" << m_tableau << '\n';
	}
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::print_basic_variables(std::ostream& os) const
{
	for (auto ij = m_basic_variables.begin(); ij != m_basic_variables.end();
			++ij)
	{
		os << '(' << ij->first << ", " << ij->second << "), ";
	}
	os << '\n';
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::create_artificial_variables()
{
	int const new_variables = m_tableau.rows() - m_basic_variables.size() - 1;

	MatX tableau(m_tableau.rows() + 1, m_tableau.cols() + new_variables + 1);
	tableau(0, 0) = 1;

	auto Ax = tableau.block(1, 1, m_tableau.rows(),
			m_tableau.cols() - 1 + new_variables);
	auto bx = tableau.bottomRightCorner(m_tableau.rows(), 1);
	auto cx = tableau.block(0, 1, 1, m_tableau.cols() - 1 + new_variables);

	Ax.leftCols(m_tableau.cols() - 1) = m_tableau.leftCols(
			m_tableau.cols() - 1);

	// Introduce artificial variables.
	int new_idx = m_tableau.cols() - 1;
	for (int row = 0; row < m_tableau.rows() - 1; ++row)
	{
		if (m_basic_variables.count(row) == 0)
		{
			Ax(row + 1, new_idx) = 1;
			cx(0, new_idx) = -1;
			new_idx++;
		}
	}

	bx = m_tableau.rightCols(1);

	std::swap(m_tableau, tableau);

	search_basic_variables();
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::price_out()
{
	for (auto ij = m_basic_variables.begin(); ij != m_basic_variables.end();
			++ij)
	{
		m_tableau.row(0) -= m_tableau.row(ij->first + 1)
				* m_tableau(0, ij->second + 1);
	}
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::search_basic_variables()
{
	m_basic_variables.clear();

	for (int col = 1; col < m_tableau.cols() - 1; ++col)
	{
		int const one_row = identify_one(
				m_tableau.block(1, col, m_tableau.rows() - 1, 1));
		if (one_row >= 0 && m_basic_variables.count(one_row) == 0) // Not already there
		{
			m_basic_variables[one_row] = col - 1;
		}
	}
}

template<typename Scalar>
int SimplexLPSolver<Scalar>::identify_one(VecX const& x)
{
	int ii = -1;
	for (int i = 0; i < x.size(); ++i)
	{
		if (x[i] == 1 && ii == -1)
		{
			ii = i;
		} else if (x[i] != 0)
		{
			return -1;
		}
	}

	return ii;
}

template class SimplexLPSolver<double> ;
