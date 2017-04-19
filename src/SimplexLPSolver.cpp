/*
 * SimplexLPSolver.cpp
 *
 *  Created on: Apr 18, 2017
 *      Author: JM
 */

#include "SimplexLPSolver.h"

template<typename Scalar>
SimplexLPSolver<Scalar>::SimplexLPSolver(MatX const& tableau) :
		m_tableau(tableau), m_num_extra_variables(0)
{
	search_basic_variables();
//	std::cout << "Tableau:\n" << m_tableau << '\n';
	price_out();
//	std::cout << "Priced out tableau:\n" << m_tableau << '\n';
}

template<typename Scalar>
bool SimplexLPSolver<Scalar>::solve()
{
	Scalar constexpr epsilon = 1.e-12;

	if (!is_canonical())
	{
		make_canonical();
		iterate_pivot();
		if (m_tableau.topRightCorner(1, 1)(0, 0) > epsilon)
		{
			return false;
		}

		MatX tableau(m_tableau.rows() - 1,
				m_tableau.cols() - m_num_extra_variables - 1);

		tableau.leftCols(m_tableau.cols() - m_num_extra_variables - 1) =
				m_tableau.block(1, 1, m_tableau.rows() - 1,
						m_tableau.cols() - m_num_extra_variables - 1);
		tableau.rightCols(1) = m_tableau.bottomRightCorner(m_tableau.rows() - 1,
				1);

		std::swap(m_tableau, tableau);
		search_basic_variables(); // TODO: just update.
	}

	//std::cout << "Phase 2 tableau:\n" << m_tableau << '\n';

	iterate_pivot();

	return true;
}

template<typename Scalar>
typename SimplexLPSolver<Scalar>::VecX SimplexLPSolver<Scalar>::get_solution() const
{
	VecX solution(m_tableau.cols() - 2);
	for (int i = 0; i < m_tableau.cols() - 2; ++i)
	{
		if (m_reverse_basic_variables.count(i))
		{
			solution[i] = m_tableau.rightCols(1)(
					m_reverse_basic_variables.at(i) + 1, 0);
		}
	}

	return solution;
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::iterate_pivot()
{
	for (;;)
	{
		int const col = find_pivot_col();
		if (col == -1)
		{
			break;
		}
		const int row = find_pivot_row(col);
		pivot(row, col);

		//	std::cout << "Pivot: " << row << " " << col << '\n';
		//	std::cout << "After pivot, tableau =\n" << get_tableau() << '\n';
	}
}

template<typename Scalar>
int SimplexLPSolver<Scalar>::find_pivot_col() const
{
	for (int i = 0; i < m_tableau.cols() - 2; ++i)
//	for (int i = m_tableau.cols() - 3; i >= 0; --i)
	{
		//	std::cout << "i = "<< i <<": " << m_tableau(0, i + 1) << '\n';
		if (m_reverse_basic_variables.count(i))
		{
			continue;
		}
		if (m_tableau(0, i + 1) > 0.)
		{
			return i;
		}
	}
	return -1;
}

template<typename Scalar>
int SimplexLPSolver<Scalar>::find_pivot_row(int col) const
{
	// std::cout << "Column " << col << ": " << m_tableau.block(1, col + 1, m_tableau.rows() - 1, 1) << '\n';

	VecX ratios = m_tableau.bottomRightCorner(m_tableau.rows() - 1, 1).array()
			/ m_tableau.block(1, col + 1, m_tableau.rows() - 1, 1).array();

	for (int i = 0; i < ratios.size(); ++i)
	{
		if (ratios[i] <= 0)
		{
			ratios[i] = std::numeric_limits<Scalar>::max();
		}
	}

	// std::cout << "Ratios = " << ratios << '\n';

	int row;
	ratios.minCoeff(&row);
	return row;
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::pivot(int row, int col)
{
	m_tableau.row(row + 1) /= m_tableau(row + 1, col + 1);
	for (int i = 0; i < m_tableau.rows(); ++i)
	{
		if (i != row + 1)
		{
			m_tableau.row(i) -= m_tableau(i, col + 1) * m_tableau.row(row + 1);
		}
	}

	search_basic_variables(); // TODO: just change what's needed.

//	std::cout << "After pivot basic variables = ";
//	print_basic_variables(std::cout);
//	std::cout << '\n';
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::make_canonical()
{
//	std::cout << "Non-canonical: basic variables are\n";
//	print_basic_variables(std::cout);

//	std::cout << "Extending...\n";

	create_artificial_variables();

//	std::cout << "After extension, basic variables are:\n";
//	print_basic_variables(std::cout);
//
//	std::cout << "Tableau:\n" << m_tableau << '\n';

	assert(is_canonical());

	price_out();
//	std::cout << "Priced out tableau:\n" << m_tableau << '\n';
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
	m_num_extra_variables = m_tableau.rows() - m_basic_variables.size() - 1;

	MatX tableau(m_tableau.rows() + 1, m_tableau.cols() + m_num_extra_variables + 1);
	tableau(0, 0) = 1;

	auto Ax = tableau.block(1, 1, m_tableau.rows(),
			m_tableau.cols() - 1 + m_num_extra_variables);
	auto bx = tableau.bottomRightCorner(m_tableau.rows(), 1);
	auto cx = tableau.block(0, 1, 1, m_tableau.cols() - 1 + m_num_extra_variables);

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
	m_reverse_basic_variables.clear();

	for (int col = 1; col < m_tableau.cols() - 1; ++col)
	{
		int const one_row = identify_one(
				m_tableau.block(1, col, m_tableau.rows() - 1, 1));
		if (one_row >= 0 && m_basic_variables.count(one_row) == 0) // Not already there
		{
			m_basic_variables[one_row] = col - 1;
			m_reverse_basic_variables[col - 1] = one_row;
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
