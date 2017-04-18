/*
 * SimplexLPSolver.h
 *
 *  Created on: Apr 17, 2017
 *      Author: JM
 */

#pragma once

#include "LinearProgrammingSolver.h"

template<typename Scalar>
class SimplexLPSolver {
public:
	using MatX= typename LinearProgrammingSolver<Scalar>::MatX;
	using VecX = typename LinearProgrammingSolver<Scalar>::VecX;

	SimplexLPSolver(MatX const& tablo);

	void make_canonical();

	void print_basic_variables(std::ostream& os) const;

protected:

	void create_artificial_variables();

	void price_out();

	bool is_canonical() const
	{
		return m_basic_variables.size() == m_tableau.rows() - 1;
	}

	void search_basic_variables();

	static int identify_one(VecX const& x);

	std::map<int, int> m_basic_variables; ///< Maps basic variable indices to their column index in A.
	MatX m_tableau;
};
