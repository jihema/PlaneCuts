/*
 * SimplexLPSolver.h
 *
 *  Created on: Apr 17, 2017
 *      Author: JM
 */

#pragma once

#include "EigenIncludes.h"

template<typename Scalar>
class SimplexLPSolver {
public:
	using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
	using MatX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

	/**
	 * \brief Constructor from standard form.
	 *
	 * \param A and \param b define the constraints: A x <?> b where the meaning of <?> is given
	 * by \param inequalities: +1 for <= constraint, -1 for >= constraint, 0 for == constraint.
	 * \param c defines the objective function x -> c.x to be minimized.
	 */
	SimplexLPSolver(MatX const& A, VecX const& b, VecX const& c,
			VecX const& inequalities = VecX());
	SimplexLPSolver(MatX const& tableau);

	bool solve();
	void print_basic_variables(std::ostream& os) const;
	VecX get_solution() const;
	const MatX& get_tableau() const
	{
		return m_tableau;
	}
	Scalar get_optimal_value() const
	{
		return m_tableau(0, m_tableau.cols() - 1);
	}

private:
	static Scalar const s_epsilon;

	/**
	 * \brief Builds the initial tableau, adding slack variables if necessary.
	 */
	static MatX make_tableau(MatX const& A, VecX const& b, VecX const& c,
			VecX const& inequalities);
	/**
	 *  \brief If x has exactly one positive coefficient and the rest is zero,
	 *  return its index; otherwise returns -1.
	 */
	static int identify_one(VecX const& x);

	void iterate_pivot();
	void make_canonical();
	/**
	 * \brief Finds the pivot column according to Bland's rule.
	 */
	int find_pivot_col() const;
	/**
	 * \brief Finds the pivot row for the given column according to Bland's rule.
	 */
	int find_pivot_row(int col) const;
	/**
	 * \brief Do the pivot, setting A(row, col) to 1.
	 */
	void pivot(int row, int col);
	/**
	 * \brief Creates artificial variables as needed during phase 1 to put the
	 * problem in canonical form.
	 */
	void create_artificial_variables();
	void price_out();
	/**
	 * \brief Checks if the problem is in canonical form.
	 *
	 * Basic variables bookkeeping must be up to date, or call search_basic_variables().
	 */
	bool is_canonical() const
	{
		return m_reverse_basic_variables.size() == m_tableau.rows() - 1;
	}
	/**
	 * \brief Identify existing basic variables.
	 */
	void search_basic_variables();

	MatX m_tableau;
	int m_num_slack_variables;
	std::vector<int> m_basic_variables; ///< Maps basic variable indices to their column index in A.
	std::map<int, int> m_reverse_basic_variables; ///< Reverse map of the basic variables.
	int m_num_extra_variables;
};
