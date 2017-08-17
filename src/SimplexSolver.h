/*
 * SimplexSolver.h
 *
 *  Created on: Aug 12, 2017
 *      Author: JM
 */

#pragma once

#include <iostream>

#include "utilities.h"

template<typename Scalar>
class SimplexSolver {
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
    SimplexSolver(MatX const& A, VecX const& b, VecX const& c,
            VecX const& inequalities = VecX(),
            long const num_free_variables = 0);

    virtual ~SimplexSolver()
    {
    }

    bool solve();

    VecX get_solution() const;

    Scalar get_optimal_value() const
    {
        return m_tableau.bottomRightCorner(1, 1)(0, 0);
    }

    void reverse_solve();

    long getReverseSolveCounter() const
    {
        return m_reverse_solve_counter;
    }

    void setReverseSolveCounter(long reverseSolveCounter)
    {
        m_reverse_solve_counter = reverseSolveCounter;
    }

private:
    static Scalar const s_epsilon;

    /**
     * \brief Builds the initial tableau.
     *
     * For each non-zero entry in inequalities, adds a slack variable to turn that
     * inequality constraint into an equality.
     */
    static MatX make_tableau(MatX const& A, VecX const& b, VecX const& c,
            VecX const& inequalities);

    /**
     * \brief Remove free variables from linear programming problem, save their
     * binding equations.
     */
    void set_aside_free_variables();

    /**
     * Make sure that b >= 0, this is assumed by search_basic_variables().
     */
    void make_b_non_negative();

    /**
     * \brief If x has exactly one positive coefficient and the rest is zero,
     *  return its index; otherwise returns -1.
     */
    static long identify_one(VecX const& x);

    /**
     * Price out all basic variables.
     */
    void price_out();

    /**
     * Transforms the problem into canonical form, by adding extra variables.
     */
    void make_canonical();
    /**
     * \brief Creates extra variables as needed during phase 1 to put the
     * problem in canonical form.
     */
    void create_extra_variables();
    /**
     * \brief Eliminate extra variables, when they are not needed anymore.
     */
    void eliminate_extra_variables();
    /**
     * \brief Checks if the problem is in canonical form.
     *
     * Basic variables bookkeeping must be up to date, or call search_basic_variables().
     */
    bool is_canonical() const
    {
        return num_basic_variables() == num_constraints();
    }

    /**
     * \brief Identify existing basic variables.
     */
    void search_basic_variables();
    bool is_basic_variable(long i) const
    {
        return m_reverse_basic_variables.count(i) > 0;
    }

    /**
     * Iterates pivot algorithm until optimal solution is found.
     */
    void iterate_pivot();
    /**
     * \brief Finds the pivot column according to Bland's rule.
     */
    long find_pivot_col() const;
    /**
     * \brief Finds the pivot row for the given column according to Bland's rule.
     */
    long find_pivot_row(long col) const;
    /**
     * \brief Do the pivot, setting A(row, col) to 1.
     */
    void pivot(long row, long col);

    long num_basic_variables() const
    {
        return m_reverse_basic_variables.size();
    }
    long num_variables() const
    {
        return m_tableau.cols() - 1;
    }
    long num_constraints() const
    {
        return m_tableau.rows() - m_num_objectives;
    }

    MatX m_tableau;
    long m_num_slack_variables;
    long m_num_extra_variables;
    long m_num_free_variables;
    long m_num_objectives;
    long m_reverse_solve_counter;
    MatX m_free_variable_equations;
    std::vector<long> m_basic_variables; ///< Maps basic variable indices to their column index in A.
    std::map<long, long> m_reverse_basic_variables; ///< Reverse map of the basic variables.
};

