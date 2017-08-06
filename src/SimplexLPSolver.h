/*
 * SimplexLPSolver.h
 *
 *  Created on: Apr 17, 2017
 *      Author: JM
 */

#pragma once

#include <iostream>
#include "EigenIncludes.h"

template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> x)
{
    std::copy(x.begin(), x.end(), std::ostream_iterator<T>(os, ", "));
    return os;
}

template<typename Scalar>
class SimplexLPSolver {
public:
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MatX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    /**
     * \brief Constructor for the linear programming problem.
     *
     * The objective function is x -> c.x, to be minimized under constraints A x <?> b,
     * where the meaning of <?> is given by the inequalities vector:
     * +1 for <= constraint, -1 for >= constraint, 0 for == constraint.
     */
    SimplexLPSolver(MatX const& A, VecX const& b, VecX const& c,
            VecX const& inequalities = VecX(),
            int const num_free_variables = 0);

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

    int num_variables() const
    {
        return m_tableau.cols() - 1;
    }
    int num_constraints() const
    {
        return m_tableau.rows() - 1;
    }

    /**
     * \brief Builds the initial tableau.
     *
     * For each non-zero entry in inequalities, adds a slack variable to turn that
     * inequality constraint into an equality.
     */
    static MatX make_tableau(MatX const& A, VecX const& b, VecX const& c,
            VecX const& inequalities);

    /**
     * \brief If x has exactly one positive coefficient and the rest is zero,
     *  return its index; otherwise returns -1.
     */
    static int identify_one(VecX const& x);

    /**
     * \brief Remove free variables from linear programming problem, save their
     * binding equations.
     */
    void set_aside_free_variables();

    /**
     * Make sure that b >= 0, this is assumed by search_basic_variables().
     */
    void make_b_non_negative();
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
     * \brief Creates extra variables as needed during phase 1 to put the
     * problem in canonical form.
     */
    void create_extra_variables();
    /**
     * \brief Eliminate extra variables, which are not needed anymore..
     */
    void eliminate_extra_variables();
    void price_out();
    /**
     * \brief Identify existing basic variables.
     */
    void search_basic_variables();
    bool is_basic_variable(int i) const
    {
        return m_reverse_basic_variables.count(i) > 0;
    }
    int num_basic_variables() const
    {
        return m_reverse_basic_variables.size();
    }
    /**
     * \brief Checks if the problem is in canonical form.
     *
     * Basic variables bookkeeping must be up to date, or call search_basic_variables().
     */
    bool is_canonical() const
    {
        return num_basic_variables() == num_constraints();
    }

    MatX m_tableau;
    int m_num_slack_variables;
    std::vector<int> m_basic_variables; ///< Maps basic variable indices to their column index in A.
    std::map<int, int> m_reverse_basic_variables; ///< Reverse map of the basic variables.
    int m_num_extra_variables;
    int m_num_free_variables;
    MatX m_free_variable_equations;
};
