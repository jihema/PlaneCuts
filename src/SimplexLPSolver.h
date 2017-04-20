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

	SimplexLPSolver(MatX const& A, VecX const& b, VecX const& c);
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
	void iterate_pivot();
	void make_canonical();
	int find_pivot_col() const;
	int find_pivot_row(int col) const;
	void pivot(int row, int col);
	void create_artificial_variables();
	void price_out();
	bool is_canonical() const
	{
		return m_basic_variables.size() == m_tableau.rows() - 1;
	}
	void search_basic_variables();

	static int identify_one(VecX const& x);
	static MatX make_tableau(MatX const& A, VecX const& b, VecX const& c);

	MatX m_tableau;
	std::map<int, int> m_basic_variables; ///< Maps basic variable indices to their column index in A. TODO: replace this with a vector.
	std::map<int, int> m_reverse_basic_variables;
	int m_num_extra_variables;
};
