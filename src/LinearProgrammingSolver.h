/*
 * LinearProgrammingSolver.h
 *
 *  Created on: Apr 17, 2017
 *      Author: JM
 */

#ifndef SRC_LINEARPROGRAMMINGSOLVER_H_
#define SRC_LINEARPROGRAMMINGSOLVER_H_

#include "EigenIncludes.h"

/**
 * \brief Virtual base class for linear programming solvers.
 */
template<typename Scalar>
class LinearProgrammingSolver {
public:
	using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
	using MatX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

	/**
	 * \brief Constructor in standard form.
	 *
	 * The feasible region is: $A x = b, x \geq 0$.
	 * The objective function to be maximized is: $c \cdot x$.
	 */
	LinearProgrammingSolver(MatX const& A, VecX const& b, VecX const& c) :
			m_A(A), m_b(b), m_c(c)
	{
		assert(A.rows() == b.size());
		assert(A.cols() == c.size());
	}

	virtual ~LinearProgrammingSolver()
	{
	}




	MatX get_tableau() const
	{
		MatX tableau(m_A.rows() + 1, m_A.cols() + 2);
		tableau(0, 0) = 1;
		tableau.block(1, 1, m_A.rows(), m_A.cols()) = m_A;
		tableau.block(0, 1, 1, m_A.cols()) = -m_c.transpose();
		tableau.block(1, m_A.cols() + 1, m_A.rows(), 1) = m_b;

		return tableau;
	}

protected:
	MatX const m_A; ///< (n x d)-matrix of constraint directions.
	VecX const m_b; ///< n-vector of constraint values.
	VecX const m_c; ///< n-vector of optimization direction.

};

#endif /* SRC_LINEARPROGRAMMINGSOLVER_H_ */
