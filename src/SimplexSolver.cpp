/*
 * SimplexSolver.cpp
 *
 *  Created on: Aug 12, 2017
 *      Author: JM
 */

#include "SimplexSolver.h"

template<>
float const SimplexSolver<float>::s_epsilon = 1.e-5;
template<>
double const SimplexSolver<double>::s_epsilon = 1.e-12;

template<typename Scalar>
SimplexSolver<Scalar>::SimplexSolver(MatX const& A, VecX const& b,
        VecX const& c, VecX const& inequalities, long const num_free_variables) :
        m_tableau(make_tableau(A, b, c, inequalities)), //
        m_num_slack_variables(num_variables() - A.cols()), //
        m_num_extra_variables(0), //
        m_num_free_variables(num_free_variables), //
        m_num_objectives(1), //
        m_reverse_solve_counter(0)
{
    set_aside_free_variables();
    make_b_non_negative();
    search_basic_variables();
    price_out();
}

template<typename Scalar>
bool SimplexSolver<Scalar>::solve()
{
    if (!is_canonical()) // Phase 1.
    {
        make_canonical();
        iterate_pivot();

        // TODO: understand why the following test fails (test 9 & 10) and why it doesn't matter.
//        if (fabs(m_tableau.bottomRightCorner(1, 1)(0, 0)) > s_epsilon) // Successful phase 1 should leave objective function value to zero.
//        {
//            std::cout << "Phase 1 failed: "
//                    << m_tableau.bottomRightCorner(1, 1)(0, 0) << '\n';
//            return false;
//        }

        eliminate_extra_variables();
    }

    iterate_pivot(); // Phase 2.

    return true;
}

template<typename Scalar>
void SimplexSolver<Scalar>::reverse_solve()
{
    assert(is_canonical());

    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
            ", ", ", ", "", "", " << ", ";");
    std::cout << get_solution().format(CommaInitFmt) << '\n';
    m_reverse_solve_counter++;

    for (long variable = 0; variable < num_variables(); ++variable)
    {
        if (m_reverse_basic_variables.count(variable) > 0)
        {
            continue;
        }

        for (long constraint = 0; constraint < num_constraints(); ++constraint)
        {
            if (m_tableau(constraint, variable) == 0)
            {
                continue;
            }

            long const basic_variable = m_basic_variables[constraint];
            pivot(constraint, variable);
//            std::cout << "After pivot " << constraint << " " << variable << '\n' << m_tableau << '\n' << '\n';

            VecX const sol = m_tableau.topRightCorner(num_constraints(),1);
            if ((sol.array() >= 0).all())
            {
                long const back_variable = find_pivot_col();
                if (back_variable >= 0)
                {
                    long const back_constraint = find_pivot_row(back_variable);
                    if (back_variable == basic_variable
                            && back_constraint == constraint)
                    {
                        reverse_solve();
                    }
                }
            }

            pivot(constraint, basic_variable);
        }
    }
}

template<typename Scalar>
void SimplexSolver<Scalar>::make_canonical()
{
    create_extra_variables();
    assert(is_canonical());
    price_out();
}

template<typename Scalar>
void SimplexSolver<Scalar>::create_extra_variables()
{
    m_num_extra_variables = num_constraints() - num_basic_variables();

    MatX tableau(m_tableau.rows() + 1,
            m_tableau.cols() + m_num_extra_variables);

    tableau.block(0, 0, m_tableau.rows(), num_variables()) = m_tableau.leftCols(
            num_variables());

    tableau.block(0, m_tableau.cols() + m_num_extra_variables - 1,
            m_tableau.rows(), 1) = m_tableau.rightCols(1);

    long idx = num_variables();
    for (long constraint = 0; constraint < num_constraints(); ++constraint)
    {
        if (m_basic_variables[constraint] < 0) // Not already a basic variable.
        {
            tableau(constraint, idx) = 1;
            tableau.bottomRows(1)(0, idx) = -1;
            idx++;
        }
    }

    std::swap(m_tableau, tableau);
    m_num_objectives++;

    search_basic_variables();
}

template<typename Scalar>
void SimplexSolver<Scalar>::eliminate_extra_variables()
{
    MatX tableau(m_tableau.rows() - 1,
            m_tableau.cols() - m_num_extra_variables);

    tableau.leftCols(num_variables() - m_num_extra_variables) =
            m_tableau.topLeftCorner(m_tableau.rows() - 1,
                    num_variables() - m_num_extra_variables);

    tableau.rightCols(1) = m_tableau.topRightCorner(m_tableau.rows() - 1, 1);

    std::swap(m_tableau, tableau);
    m_num_objectives--;

    search_basic_variables();
}

template<typename Scalar>
typename SimplexSolver<Scalar>::MatX SimplexSolver<Scalar>::make_tableau(
        MatX const& A, VecX const& b, VecX const& c, VecX const& inequalities)
{
    long const num_constraints = b.size();
    long const num_variables = c.size(); // Without slack variables.
    assert(A.rows() == num_constraints);
    assert(A.cols() == num_variables);

    // Add slack variables if necessary to turn inequality constraints into equality constraints.
    long const num_slacks = (inequalities.array() != 0).count();
    MatX slack_constraints(num_constraints, num_slacks);
    for (long s = 0, si = 0; s < inequalities.size(); ++s)
    {
        if (inequalities[s] != 0)
        {
            slack_constraints(s, si++) = inequalities[s] > 0 ? 1 : -1;
        }
    }

    // Build the tableau with objective function on the first row,
    // original problem + slack variables constraints on subsequent rows.
    MatX tableau(num_constraints + 1, num_variables + num_slacks + 1);
    tableau.block(0, 0, num_constraints, num_variables) = A;
    tableau.block(0, num_variables, num_constraints, num_slacks) =
            slack_constraints;
    tableau.block(0, num_variables + num_slacks, A.rows(), 1) = b;
    tableau.block(num_constraints, 0, 1, num_variables) = -c.transpose();

    return tableau;
}

template<typename Scalar>
void SimplexSolver<Scalar>::set_aside_free_variables()
{
    assert(m_num_free_variables <= num_variables());

    int variable;
    for (variable = 0; variable < m_num_free_variables; ++variable)
    {
        // Find a pivot for this variable.
        int constraint;
        for (constraint = variable;
                constraint < num_constraints()
                        && m_tableau(constraint, variable) == 0; ++constraint)
        {
        }

        if (constraint < num_constraints()) // Pivot found.
        {
            // Swap so the pivot is in diagonal position.
            if (constraint > variable)
            {
                VecX const tmp = m_tableau.row(variable);
                m_tableau.row(variable) = m_tableau.row(constraint);
                m_tableau.row(constraint) = tmp;
                constraint = variable;
            }

            // Eliminate the pivot variable from other rows.
            m_tableau.row(constraint) /= m_tableau(constraint, variable);
            for (int row = 0; row < m_tableau.rows(); ++row)
            {
                if (row != variable)
                {
                    m_tableau.row(row) -= m_tableau(row, variable)
                            * m_tableau.row(variable);
                }
            }
        } else // Pivot not found, we just stop here.
        {
            std::cout << __PRETTY_FUNCTION__ << "Warning: pivot not found.\n";
            break;
        }
    }

    // variable now contains the number of effective free variables
    // i.e. the ones that will need to be computed in the end.
    // Other free variables can be simply set to zero.
    int const effective_free_variables = variable;
    m_free_variable_equations = m_tableau.block(0, m_num_free_variables,
            effective_free_variables, m_tableau.cols() - m_num_free_variables);

    // Create a new tableau for the non-free variables (standard form) linear programming problem.
    long const cols = m_tableau.cols() - m_num_free_variables;
    long const rows = m_tableau.rows() - effective_free_variables;
    MatX tableau(rows, cols);
    tableau.block(rows - 1, 0, 1, cols) = m_tableau.block(num_constraints(),
            m_num_free_variables, 1, cols);
    tableau.block(0, 0, rows - 1, cols) = m_tableau.block(
            effective_free_variables, m_num_free_variables, rows - 1, cols);
    std::swap(m_tableau, tableau);
}

template<typename Scalar>
void SimplexSolver<Scalar>::make_b_non_negative()
{
    for (int constraint = 0; constraint < num_constraints(); ++constraint)
    {
        if (m_tableau.rightCols(1)(constraint, 0) < 0)
        {
            m_tableau.row(constraint) *= -1;
        }
    }
}

template<typename Scalar>
void SimplexSolver<Scalar>::search_basic_variables()
{
    m_basic_variables.clear();
    m_basic_variables.resize(num_constraints(), -1);
    m_reverse_basic_variables.clear();

    for (int variable = 0; variable < num_variables(); ++variable)
    {
        long const one_row = identify_one(
                m_tableau.block(0, variable, num_constraints(), 1));
        if (one_row >= 0 // Found one.
        && m_tableau.rightCols(1)(one_row, 0) >= 0 // Acceptable.
        && m_basic_variables[one_row] < 0 // Not already a basic variable.
                )
        {
            m_basic_variables[one_row] = variable;
            m_reverse_basic_variables[variable] = one_row;
            m_tableau.row(one_row) /= m_tableau(one_row, variable);
        }
    }
}

template<typename Scalar>
long SimplexSolver<Scalar>::identify_one(VecX const& x)
{
    long ii = -1;
    for (long i = 0; i < x.size(); ++i)
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

template<typename Scalar>
void SimplexSolver<Scalar>::price_out()
{
    for (auto variable = m_basic_variables.begin();
            variable != m_basic_variables.end(); ++variable)
    {
        if (*variable >= 0)
        {
            m_tableau.bottomRows(1) -= m_tableau.row(
                    variable - m_basic_variables.begin())
                    * m_tableau.bottomRows(1)(0, *variable);
        }
    }
}

template<typename Scalar>
void SimplexSolver<Scalar>::iterate_pivot()
{
    long col = -1;
    while ((col = find_pivot_col()) != -1)
    {
        long const row = find_pivot_row(col);
        pivot(row, col);
    }
}

template<typename Scalar>
long SimplexSolver<Scalar>::find_pivot_col() const
{
    for (long variable = 0; variable < num_variables(); ++variable)
    {
        if (!is_basic_variable(variable)
                && m_tableau.bottomRows(1)(0, variable) > s_epsilon) // Not sure why we need an epsilon here, but if we don't a cycle may occur.
        {
            return variable;
        }
    }

    // No pivot found.
    return -1;
}

template<typename Scalar>
long SimplexSolver<Scalar>::find_pivot_row(long variable) const
{
    VecX ratios = m_tableau.topRightCorner(num_constraints(), 1).array()
            / m_tableau.block(0, variable, num_constraints(), 1).array();

    for (long constraint = 0; constraint < ratios.size(); ++constraint)
    {
        if (m_tableau(constraint, variable) <= 0)
        {
            ratios[constraint] = std::numeric_limits<Scalar>::max();
        }
    }

    long row;
    ratios.minCoeff(&row); // NB For Bland's rule, we rely here on minCoeff returning the first one in case of a tie.

    return row;
}

template<typename Scalar>
void SimplexSolver<Scalar>::pivot(long constraint, long variable)
{
    m_tableau.row(constraint) /= m_tableau(constraint, variable);
    for (long i = 0; i < m_tableau.rows(); ++i)
    {
        if (i != constraint)
        {
            m_tableau.row(i) -= m_tableau(i, variable)
                    * m_tableau.row(constraint);
        }
    }

    // Update basic variables.
    m_reverse_basic_variables.erase(m_basic_variables[constraint]);
    m_basic_variables[constraint] = variable;
    m_reverse_basic_variables[variable] = constraint;
}

template<typename Scalar>
typename SimplexSolver<Scalar>::VecX SimplexSolver<Scalar>::get_solution() const
{
    long const num_non_free_non_slack_variables = num_variables()
            - m_num_slack_variables;

    VecX solution(num_non_free_non_slack_variables + m_num_free_variables);

    // Copy non-free variables from the tableau.
    for (long i = 0; i < num_non_free_non_slack_variables; ++i)
    {
        if (is_basic_variable(i))
        {
            solution[i + m_num_free_variables] = m_tableau.rightCols(1)(
                    m_reverse_basic_variables.at(i), 0);
        }
    }

    // Compute free variables.
    if (m_free_variable_equations.rows() > 0)
    {
        VecX solution_with_slack(m_tableau.cols());
        for (long i = 0; i < num_variables(); ++i)
        {
            if (is_basic_variable(i))
            {
                solution_with_slack[i] = m_tableau.rightCols(1)(
                        m_reverse_basic_variables.at(i), 0);
            }
        }
        solution_with_slack[num_variables()] = -1;

        solution.head(m_free_variable_equations.rows()) =
                -m_free_variable_equations * solution_with_slack;
    }

    return solution;
}

template class SimplexSolver<float> ;
template class SimplexSolver<double> ;

