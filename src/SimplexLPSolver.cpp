/*
 * SimplexLPSolver.cpp
 *
 *  Created on: Apr 18, 2017
 *      Author: JM
 */

#include "SimplexLPSolver.h"
#include <iostream>

template<>
float const SimplexLPSolver<float>::s_epsilon = 1.e-5;
template<>
double const SimplexLPSolver<double>::s_epsilon = 1.e-12;

template<typename Scalar>
SimplexLPSolver<Scalar>::SimplexLPSolver(MatX const& A, VecX const& b,
        VecX const& c, VecX const& inequalities, int const num_free_variables) :
		SimplexLPSolver<Scalar>(make_tableau(A, b, c, inequalities))
{
	// Record number of slack variables so we can omit them from the solution.
	m_num_slack_variables = m_tableau.cols() - A.cols() - 2;

	set_aside_free_variables(num_free_variables);
}

template<typename Scalar>
SimplexLPSolver<Scalar>::SimplexLPSolver(MatX const& tableau) :
		m_tableau(tableau), //
		m_new_tableau(m_tableau.rightCols(m_tableau.cols() - 1)), m_num_slack_variables(
		        0), //
		m_num_extra_variables(0), //
		m_num_free_variables(0)
{
	std::cout << "Made tableau\n" << m_tableau << '\n';
	std::cout << "------------\n" << m_new_tableau << '\n';

	make_b_non_negative();

	new_search_basic_variables();

	price_out();
	new_price_out();

	std::cout << "After pricing out\n" << m_tableau << '\n';
	std::cout << "------------\n" << m_new_tableau << '\n';
	assert(m_new_tableau == m_tableau.rightCols(m_tableau.cols() - 1));
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::set_aside_free_variables(
        int const num_free_variable)
{
	m_num_free_variables = num_free_variable;
	assert(m_num_free_variables <= m_tableau.cols() - 2);

	int variable_idx = 0;
	for (; variable_idx < m_num_free_variables; ++variable_idx)
	{
		// Find a pivot for this variable.
		int constraint_idx = variable_idx;
		for (;
		        constraint_idx < m_tableau.rows() - 1
		                && m_tableau(constraint_idx + 1, variable_idx + 1) == 0;
		        ++constraint_idx)
		{
		}

		if (constraint_idx < m_tableau.rows() - 1) // Pivot found.
		{
			// Swap so the pivot is in diagonal position.
			if (constraint_idx > variable_idx)
			{
				VecX const tmp = m_tableau.row(variable_idx + 1);
				m_tableau.row(variable_idx + 1) = m_tableau.row(
				        constraint_idx + 1);
				m_tableau.row(constraint_idx + 1) = tmp;
			}

			// Eliminate the pivot variable from other rows.
			m_tableau.row(variable_idx + 1) /= m_tableau(variable_idx + 1,
			        variable_idx + 1);
			for (int row = 0; row < m_tableau.rows(); ++row)
			{
				if (row != variable_idx + 1)
				{
					m_tableau.row(row) -= m_tableau(row, variable_idx + 1)
					        * m_tableau.row(variable_idx + 1);
				}
			}
		} else // Pivot not found, we just stop here.
		{
			break;
		}
	}
	int const effective_free_variables = variable_idx;

	// variable_idx now contains the number of effective free variables
	// i.e. the ones that will need to be computed in the end.
	// Other free variables can be simply set to zero.
	m_free_variable_equations = m_tableau.block(1, 1 + m_num_free_variables,
	        effective_free_variables,
	        m_tableau.cols() - 1 - m_num_free_variables);

	// Create a new tableau for the non-free variables (standard form) linear programming problem.
	MatX new_tableau(m_tableau.rows() - effective_free_variables,
	        m_tableau.cols() - m_num_free_variables);
	new_tableau(0, 0) = 1;
	const int new_cols = m_tableau.cols() - m_num_free_variables - 1;
	new_tableau.block(0, 1, 1, new_cols) = m_tableau.block(0,
	        m_num_free_variables + 1, 1, new_cols);
	new_tableau.block(1, 1, m_tableau.rows() - effective_free_variables - 1,
	        new_cols) = m_tableau.block(effective_free_variables + 1,
	        m_num_free_variables + 1,
	        m_tableau.rows() - effective_free_variables - 1, new_cols);
	std::swap(m_tableau, new_tableau);

	// Change row signs again if necessary, as b may have become negative during pivoting.
	make_b_non_negative();

	// Prepare for canonization.
	search_basic_variables();
	price_out();
}

template<typename Scalar>
typename SimplexLPSolver<Scalar>::MatX SimplexLPSolver<Scalar>::make_tableau(
        MatX const& A, VecX const& b, VecX const& c, VecX const& inequalities)
{
	assert(A.rows() == b.size());
	assert(A.cols() == c.size());

	// Add slack variables if necessary.
	int const num_slacks = (inequalities.array() != 0).count();
	MatX extra_slack(A.rows(), num_slacks);
	for (int s = 0, si = 0; s < inequalities.size(); ++s)
	{
		if (inequalities[s] != 0)
		{
			extra_slack(s, si++) = inequalities[s] > 0 ? 1 : -1;
		}
	}

	MatX tableau(A.rows() + 1, A.cols() + num_slacks + 2);
	tableau(0, 0) = 1;
	tableau.block(1, 1, A.rows(), A.cols()) = A;
	tableau.block(1, A.cols() + 1, A.rows(), num_slacks) = extra_slack;
	tableau.block(0, 1, 1, A.cols()) = -c.transpose();
	tableau.block(1, A.cols() + num_slacks + 1, A.rows(), 1) = b;

	return tableau;
}

template<typename Scalar>
bool SimplexLPSolver<Scalar>::solve()
{
	if (!is_canonical())
	{
		make_canonical();

		iterate_pivot(); // Phase 1.
		//	new_iterate_pivot();

		std::cout << "After iterate\n" << m_tableau << '\n';
		std::cout << "------------\n" << m_new_tableau << '\n';
		assert(m_new_tableau == m_tableau.rightCols(m_tableau.cols() - 1));

		if (m_new_tableau.topRightCorner(1, 1)(0, 0) > s_epsilon)
		{
			return false;
		}

		// Build augmented tableau.
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

	iterate_pivot(); // Phase 2.

	return true;
}

template<typename Scalar>
typename SimplexLPSolver<Scalar>::VecX SimplexLPSolver<Scalar>::get_solution() const
{
	int const num_non_free_non_slack_variables = m_tableau.cols()
	        - m_num_slack_variables - 2;

	VecX solution(num_non_free_non_slack_variables + m_num_free_variables);

	// Copy non-free variables from the tableau.
	for (int i = 0; i < num_non_free_non_slack_variables; ++i)
	{
		if (m_reverse_basic_variables.count(i))
		{
			solution[i + m_num_free_variables] = m_tableau.rightCols(1)(
			        m_reverse_basic_variables.at(i) + 1, 0);
		}
	}

	// Compute free variables.
	if (m_free_variable_equations.rows() > 0)
	{
		int const num_non_free_variables = m_tableau.cols() - 2;
		VecX solution_with_slack(num_non_free_variables + 1);
		for (int i = 0; i < num_non_free_variables; ++i)
		{
			if (m_reverse_basic_variables.count(i))
			{
				solution_with_slack[i] = m_tableau.rightCols(1)(
				        m_reverse_basic_variables.at(i) + 1, 0);
			}
		}

		solution_with_slack[num_non_free_variables] = -1;

		solution.head(m_free_variable_equations.rows()) =
		        -m_free_variable_equations * solution_with_slack;
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
		new_pivot(row, col);
	}
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::new_iterate_pivot()
{
	for (;;)
	{
		int const col = new_find_pivot_col();
		if (col == -1)
		{
			break;
		}
		const int row = new_find_pivot_row(col);
		new_pivot(row, col);
	}
}

template<typename Scalar>
int SimplexLPSolver<Scalar>::find_pivot_col() const
{
	for (int i = 0; i < m_tableau.cols() - 2; ++i)
	{
		if (m_reverse_basic_variables.count(i))
		{
			continue;
		}
		if (m_tableau(0, i + 1) > s_epsilon) // Not sure why we need an epsilon here, but if we don't a cycle may occur.
		{
			return i;
		}
	}
	return -1;
}

template<typename Scalar>
int SimplexLPSolver<Scalar>::new_find_pivot_col() const
{
	for (int i = 0; i < m_new_tableau.cols() - 1; ++i)
	{
		if (m_reverse_basic_variables.count(i))
		{
			continue;
		}
		if (m_tableau(0, i) > s_epsilon) // Not sure why we need an epsilon here, but if we don't a cycle may occur.
		{
			return i;
		}
	}
	return -1;
}

template<typename Scalar>
int SimplexLPSolver<Scalar>::find_pivot_row(int col) const
{
	VecX ratios = m_tableau.bottomRightCorner(m_tableau.rows() - 1, 1).array()
	        / m_tableau.block(1, col + 1, m_tableau.rows() - 1, 1).array();

	for (int i = 0; i < ratios.size(); ++i)
	{
		if (ratios[i] <= 0)
		{
			ratios[i] = std::numeric_limits<Scalar>::max();
		}
	}

	int row;
	ratios.minCoeff(&row);
	return row;
}

template<typename Scalar>
int SimplexLPSolver<Scalar>::new_find_pivot_row(int col) const
{
	VecX ratios =
	        m_new_tableau.bottomRightCorner(m_tableau.rows() - 1, 1).array()
	                / m_new_tableau.block(1, col, m_tableau.rows() - 1, 1).array();

	for (int i = 0; i < ratios.size(); ++i)
	{
		if (ratios[i] <= 0)
		{
			ratios[i] = std::numeric_limits<Scalar>::max();
		}
	}

	int row;
	ratios.minCoeff(&row);
	return row;
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::pivot(int row, int col)
{
	std::cout << "Pivot: " << row << " " << col << '\n';

	m_tableau.row(row + 1) /= m_tableau(row + 1, col + 1);
	for (int i = 0; i < m_tableau.rows(); ++i)
	{
		if (i != row + 1)
		{
			m_tableau.row(i) -= m_tableau(i, col + 1) * m_tableau.row(row + 1);
		}
	}

	// Update basic variables.
//	m_reverse_basic_variables.erase(m_basic_variables[row]);
//	m_basic_variables[row] = col;
//	m_reverse_basic_variables[col] = row;
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::new_pivot(int row, int col)
{
	std::cout << "New pivot: " << row << " " << col << '\n';

	m_new_tableau.row(row + 1) /= m_new_tableau(row + 1, col);
	for (int i = 0; i < m_new_tableau.rows(); ++i)
	{
		if (i != row + 1)
		{
			m_new_tableau.row(i) -= m_new_tableau(i, col)
			        * m_new_tableau.row(row + 1);
		}
	}

	// Update basic variables.
	m_reverse_basic_variables.erase(m_basic_variables[row]);
	m_basic_variables[row] = col;
	m_reverse_basic_variables[col] = row;
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::make_canonical()
{
	create_artificial_variables();

	assert(is_canonical());

	price_out();
	new_price_out();
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::print_basic_variables(std::ostream& os) const
{
	for (auto col = m_basic_variables.begin(); col != m_basic_variables.end();
	        ++col)
	{
		if (*col >= 0)
		{
			os << '(' << col - m_basic_variables.begin() << ", " << *col
			        << "), ";
		}
	}
	os << '\n';
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::create_artificial_variables()
{
	m_num_extra_variables = m_tableau.rows() - m_reverse_basic_variables.size()
	        - 1;

	{
		MatX tableau(m_tableau.rows() + 1,
		        m_tableau.cols() + m_num_extra_variables + 1);
		tableau(0, 0) = 1;

		auto Ax = tableau.block(1, 1, m_tableau.rows(),
		        m_tableau.cols() - 1 + m_num_extra_variables);
		auto bx = tableau.bottomRightCorner(m_tableau.rows(), 1);
		auto cx = tableau.block(0, 1, 1,
		        m_tableau.cols() - 1 + m_num_extra_variables);

		Ax.leftCols(m_tableau.cols() - 1) = m_tableau.leftCols(
		        m_tableau.cols() - 1);

		// Introduce artificial variables.
		int new_idx = m_tableau.cols() - 1;
		for (int row = 0; row < m_tableau.rows() - 1; ++row)
		{
			if (m_basic_variables[row] < 0) // Not already a basic variable.
			{
				Ax(row + 1, new_idx) = 1;
				cx(0, new_idx) = -1;
				new_idx++;
			}
		}

		bx = m_tableau.rightCols(1);

		std::swap(m_tableau, tableau);
	}

	{
		MatX new_tableau(m_new_tableau.rows() + 1,
		        m_new_tableau.cols() + m_num_extra_variables + 1);

		//	new_tableau = m_tableau.rightCols(m_tableau.cols() - 1);

		new_tableau(1, 0) = 1;
		new_tableau.block(1, 1, m_new_tableau.rows(), m_new_tableau.cols() - 1) =
		        m_new_tableau.leftCols(m_new_tableau.cols() - 1);
		new_tableau.block(1, m_new_tableau.cols() + m_num_extra_variables,
		        m_new_tableau.rows(), 1) = m_new_tableau.rightCols(1);

		int new_idx = m_new_tableau.cols() - 1;
		for (int row = 0; row < m_new_tableau.rows() - 1; ++row)
		{
			if (m_basic_variables[row] < 0) // Not already a basic variable.
			{
				new_tableau(row + 2, new_idx + 1) = 1;
				new_tableau(0, new_idx + 1) = -1;
				new_idx++;
			}
		}

		std::swap(m_new_tableau, new_tableau);
	}

//	std::cout << "After artificial variables\n" << m_tableau << '\n';
//	std::cout << "------------\n" << m_new_tableau << '\n';
//	assert(m_new_tableau == m_tableau.rightCols(m_tableau.cols() - 1));

	search_basic_variables();
	new_search_basic_variables();
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::price_out()
{
	for (auto ij = m_basic_variables.begin(); ij != m_basic_variables.end();
	        ++ij)
	{
		if (*ij >= 0)
		{
			m_tableau.row(0) -= m_tableau.row(
			        ij - m_basic_variables.begin() + 1) * m_tableau(0, *ij + 1);
		}
	}
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::new_price_out()
{
	for (auto ij = m_basic_variables.begin(); ij != m_basic_variables.end();
	        ++ij)
	{
		if (*ij >= 0)
		{
			m_new_tableau.row(0) -= m_new_tableau.row(
			        ij - m_basic_variables.begin() + 1) * m_new_tableau(0, *ij);
		}
	}
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::search_basic_variables()
{
	m_basic_variables.clear();
	m_basic_variables.resize(m_tableau.cols() - 2, -1);
	m_reverse_basic_variables.clear();

	for (int col = 1; col < m_tableau.cols() - 1; ++col)
	{
		int const one_row = identify_one(
		        m_tableau.block(1, col, m_tableau.rows() - 1, 1));
		if (one_row >= 0 // Found one
		&& m_tableau.rightCols(1)(one_row + 1, 0) >= 0 // Acceptable
		&& m_basic_variables[one_row] < 0  // Not already there
		        )
		{
			m_basic_variables[one_row] = col - 1;
			m_reverse_basic_variables[col - 1] = one_row;
			m_tableau.row(one_row + 1) /= m_tableau(one_row + 1, col);
		}
	}
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::new_search_basic_variables()
{
	m_basic_variables.clear();
	m_basic_variables.resize(m_new_tableau.cols() - 1, -1);
	m_reverse_basic_variables.clear();

	for (int col = 0; col < m_new_tableau.cols(); ++col)
	{
		int const one_row = identify_one(
		        m_new_tableau.block(1, col, m_tableau.rows() - 1, 1));
		if (one_row >= 0 // Found one
		&& m_tableau.rightCols(1)(one_row + 1, 0) >= 0 // Acceptable
		&& m_basic_variables[one_row] < 0  // Not already there
		        )
		{
			m_basic_variables[one_row] = col - 1;
			m_reverse_basic_variables[col - 1] = one_row;
			m_new_tableau.row(one_row + 1) /= m_new_tableau(one_row + 1, col);
		}
	}
}

template<typename Scalar>
int SimplexLPSolver<Scalar>::identify_one(VecX const& x)
{
	int ii = -1;
	for (int i = 0; i < x.size(); ++i)
	{
		if (x[i] > 0 && ii == -1)
		{
			ii = i;
		} else if (x[i] != 0)
		{
			return -1;
		}
	}

	return ii;
}

template class SimplexLPSolver<float> ;
template class SimplexLPSolver<double> ;
