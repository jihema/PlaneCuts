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
     * \brief Builds the initial tableau, adding slack variables if necessary.
     */
    void make_tableau(MatX const& A, VecX const& b, VecX const& c,
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
    void set_aside_free_variables(int const num_free_variables);
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

    static int constexpr first_constraint_row = 1;
    inline int objective_row_idx() const
    {
        return 0;
    }
    inline static int constraint_row_idx(int constraint)
    {
        return constraint + first_constraint_row;
    }
    inline typename MatX::RowXpr constraint_row(int constraint)
    {
        return m_tableau.row(constraint_row_idx(constraint));
    }
    inline typename MatX::ConstRowXpr constraint_row(int constraint) const
    {
        return m_tableau.row(constraint_row_idx(constraint));
    }
    inline typename MatX::RowXpr objective_row()
    {
        return m_tableau.row(objective_row_idx());
    }
    inline typename MatX::ConstRowXpr objective_row() const
    {
        return m_tableau.row(objective_row_idx());
    }
    inline auto constraints_rhs()
    {
        return m_tableau.block(first_constraint_row, num_variables(),
                num_constraints(), 1);
    }
    inline auto constraints_rhs() const
    {
        return m_tableau.block(first_constraint_row, num_variables(),
                num_constraints(), 1);
    }

    MatX m_tableau;
    int m_num_slack_variables;
    std::vector<int> m_basic_variables; ///< Maps basic variable indices to their column index in A.
    std::map<int, int> m_reverse_basic_variables; ///< Reverse map of the basic variables.
    int m_num_extra_variables;
    int m_num_free_variables;
    MatX m_free_variable_equations;
};
