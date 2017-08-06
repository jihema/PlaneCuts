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
        m_tableau(make_tableau(A, b, c, inequalities)), //
        m_num_slack_variables(num_variables() - A.cols()), //
        m_num_extra_variables(0), //
        m_num_free_variables(num_free_variables)
{
    set_aside_free_variables();
    make_b_non_negative();
    search_basic_variables();
    price_out();
}

template<typename Scalar>
typename SimplexLPSolver<Scalar>::MatX SimplexLPSolver<Scalar>::make_tableau(
        MatX const& A, VecX const& b, VecX const& c, VecX const& inequalities)
{
    int const num_constraints = b.size();
    int const num_variables = c.size(); // Without slack variables.
    assert(A.rows() == num_constraints);
    assert(A.cols() == num_variables);

    // Add slack variables if necessary to turn inequality constraints into equality constraints.
    int const num_slacks = (inequalities.array() != 0).count();
    MatX slack_constraints(num_constraints, num_slacks);
    for (int s = 0, si = 0; s < inequalities.size(); ++s)
    {
        if (inequalities[s] != 0)
        {
            slack_constraints(s, si++) = inequalities[s] > 0 ? 1 : -1;
        }
    }

    // Build the tableau with objective function on the first row,
    // original problem + slack variables constraints on subsequent rows.
    MatX tableau(num_constraints + 1, num_variables + num_slacks + 1);
    tableau.block(0, 0, 1, num_variables) = -c.transpose();
    tableau.block(1, 0, num_constraints, num_variables) = A;
    tableau.block(1, num_variables, num_constraints, num_slacks) =
            slack_constraints;
    tableau.block(1, num_variables + num_slacks, A.rows(), 1) = b;

    return tableau;
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::make_b_non_negative()
{
    for (int row = 1; row < m_tableau.rows(); ++row)
    {
        if (m_tableau.rightCols(1)(row, 0) < 0)
        {
            m_tableau.row(row) *= -1;
        }
    }
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::search_basic_variables()
{
    m_basic_variables.clear();
    m_basic_variables.resize(num_variables(), -1);
    m_reverse_basic_variables.clear();

    for (int variable = 0; variable < num_variables(); ++variable)
    {
        int const one_row = identify_one(
                m_tableau.block(1, variable, num_constraints(), 1));
        if (one_row >= 0 // Found one.
        && m_tableau.rightCols(1)(one_row + 1, 0) >= 0 // Acceptable.
        && m_basic_variables[one_row] < 0 // Not already a basic variable.
                )
        {
            m_basic_variables[one_row] = variable;
            m_reverse_basic_variables[variable] = one_row;
            m_tableau.row(one_row + 1) /= m_tableau(one_row + 1, variable);
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

template<typename Scalar>
void SimplexLPSolver<Scalar>::price_out()
{
    for (auto variable = m_basic_variables.begin();
            variable != m_basic_variables.end(); ++variable)
    {
        if (*variable >= 0)
        {
            m_tableau.row(0) -= m_tableau.row(
                    variable - m_basic_variables.begin() + 1)
                    * m_tableau(0, *variable);
        }
    }
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::set_aside_free_variables()
{
    assert(m_num_free_variables <= num_variables());

    int variable;
    for (variable = 0; variable < m_num_free_variables; ++variable)
    {
        // Find a pivot for this variable.
        int constraint;
        for (constraint = variable;
                constraint < num_constraints()
                        && m_tableau(constraint + 1, variable) == 0;
                ++constraint)
        {
        }

        if (constraint < num_constraints()) // Pivot found.
        {
            // Swap so the pivot is in diagonal position.
            if (constraint > variable)
            {
                VecX const tmp = m_tableau.row(variable + 1);
                m_tableau.row(variable + 1) = m_tableau.row(constraint + 1);
                m_tableau.row(constraint + 1) = tmp;
                constraint = variable;
            }

            // Eliminate the pivot variable from other rows.
            m_tableau.row(constraint + 1) /= m_tableau(constraint + 1,
                    variable);
            for (int row = 0; row < m_tableau.rows(); ++row)
            {
                if (row != variable + 1)
                {
                    m_tableau.row(row) -= m_tableau(row, variable)
                            * m_tableau.row(variable + 1);
                }
            }
        } else // Pivot not found, we just stop here.
        {
            std::cout << "Pivot not found\n";
            break;
        }
    }

    // variable now contains the number of effective free variables
    // i.e. the ones that will need to be computed in the end.
    // Other free variables can be simply set to zero.
    int const effective_free_variables = variable;
    m_free_variable_equations = m_tableau.block(1, m_num_free_variables,
            effective_free_variables, m_tableau.cols() - m_num_free_variables);

    // Create a new tableau for the non-free variables (standard form) linear programming problem.
    int const cols = m_tableau.cols() - m_num_free_variables;
    int const rows = m_tableau.rows() - effective_free_variables;
    MatX tableau(rows, cols);
    tableau.block(0, 0, 1, cols) = m_tableau.block(0, m_num_free_variables, 1,
            cols);
    tableau.block(1, 0, rows - 1, cols) = m_tableau.block(
            effective_free_variables + 1, m_num_free_variables, rows - 1, cols);
    std::swap(m_tableau, tableau);
}

template<typename Scalar>
bool SimplexLPSolver<Scalar>::solve()
{
    if (!is_canonical()) // Phase 1.
    {
        make_canonical();
        iterate_pivot();
        if (fabs(m_tableau.topRightCorner(1, 1)(0, 0)) > s_epsilon) // Successful phase 1 should leave objective function value to zero.
        {
            return false;
        }
        eliminate_extra_variables();
    }

    iterate_pivot(); // Phase 2.

    return true;
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::make_canonical()
{
//    std::cout << "Tableau before: \n" << m_tableau << '\n';
//    std::cout << "Basic variables: " << m_basic_variables << '\n';
    create_extra_variables();
//    std::cout << "Tableau after: \n" << m_tableau << '\n';
//    std::cout << "Basic variables: " << m_basic_variables << '\n';
    assert(is_canonical());
    price_out();
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::create_extra_variables()
{
    m_num_extra_variables = num_constraints() - num_basic_variables();

    {
        MatX tableau(m_tableau.rows() + 1,
                m_tableau.cols() + m_num_extra_variables + 1);

        tableau(1, 0) = 1; // TODO: for historical reasons, this variable is created to the left of the existing ones, should be to the right like the other extra variables below.
        tableau.block(1, 1, m_tableau.rows(), num_variables()) =
                m_tableau.leftCols(num_variables());
        tableau.block(1, m_tableau.cols() + m_num_extra_variables,
                m_tableau.rows(), 1) = m_tableau.rightCols(1);

        int idx = num_variables();
        for (int constraint = 0; constraint < num_constraints(); ++constraint)
        {
            if (m_basic_variables[constraint] < 0) // Not already a basic variable.
            {
                tableau(constraint + 2, idx + 1) = 1;
                tableau(0, idx + 1) = -1;
                idx++;
            }
        }

        std::swap(m_tableau, tableau);
    }

    search_basic_variables();
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::eliminate_extra_variables()
{
    MatX tableau(num_constraints(), num_variables() - m_num_extra_variables);
    tableau.leftCols(num_variables() - m_num_extra_variables) = m_tableau.block(
            1, 1, num_constraints(), num_variables() - m_num_extra_variables);
    tableau.rightCols(1) = m_tableau.bottomRightCorner(num_constraints(), 1);
    std::swap(m_tableau, tableau);

    search_basic_variables();
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::iterate_pivot()
{
    int col = -1;
    while ((col = find_pivot_col()) != -1)
    {
        pivot(find_pivot_row(col), col);
    }
}

template<typename Scalar>
int SimplexLPSolver<Scalar>::find_pivot_col() const
{
    for (int variable = 0; variable < num_variables(); ++variable)
    {
        if (!is_basic_variable(variable) && m_tableau(0, variable) > s_epsilon) // Not sure why we need an epsilon here, but if we don't a cycle may occur.
        {
            return variable;
        }
    }

    // No pivot found.
    return -1;
}

template<typename Scalar>
int SimplexLPSolver<Scalar>::find_pivot_row(int variable) const
{
    VecX ratios = m_tableau.bottomRightCorner(num_constraints(), 1).array()
            / m_tableau.block(1, variable, num_constraints(), 1).array();

    for (int constraint = 0; constraint < ratios.size(); ++constraint)
    {
        if (ratios[constraint] <= 0)
        {
            ratios[constraint] = std::numeric_limits<Scalar>::max();
        }
    }

    int row;
    ratios.minCoeff(&row); // NB For Bland's rule, we rely here on minCoeff returning the first one in case of a tie.
    return row;
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::pivot(int contraint, int variable)
{
    m_tableau.row(contraint + 1) /= m_tableau(contraint + 1, variable);
    for (int i = 0; i < m_tableau.rows(); ++i)
    {
        if (i != contraint + 1)
        {
            m_tableau.row(i) -= m_tableau(i, variable)
                    * m_tableau.row(contraint + 1);
        }
    }

    // Update basic variables.
    m_reverse_basic_variables.erase(m_basic_variables[contraint]);
    m_basic_variables[contraint] = variable;
    m_reverse_basic_variables[variable] = contraint;
}

template<typename Scalar>
typename SimplexLPSolver<Scalar>::VecX SimplexLPSolver<Scalar>::get_solution() const
{
    int const num_non_free_non_slack_variables = num_variables()
            - m_num_slack_variables;

    VecX solution(num_non_free_non_slack_variables + m_num_free_variables);

    // Copy non-free variables from the tableau.
    for (int i = 0; i < num_non_free_non_slack_variables; ++i)
    {
        if (is_basic_variable(i))
        {
            solution[i + m_num_free_variables] = m_tableau.rightCols(1)(
                    m_reverse_basic_variables.at(i) + 1, 0);
        }
    }

    // Compute free variables.
    if (m_free_variable_equations.rows() > 0)
    {
        VecX solution_with_slack(m_tableau.cols());
        for (int i = 0; i < num_variables(); ++i)
        {
            if (is_basic_variable(i))
            {
                solution_with_slack[i] = m_tableau.rightCols(1)(
                        m_reverse_basic_variables.at(i) + 1, 0);
            }
        }
        solution_with_slack[num_variables()] = -1;

        solution.head(m_free_variable_equations.rows()) =
                -m_free_variable_equations * solution_with_slack;
    }

    return solution;
}

template<typename Scalar>
void SimplexLPSolver<Scalar>::print_basic_variables(std::ostream& os) const
{
    for (auto variable = m_basic_variables.begin();
            variable != m_basic_variables.end(); ++variable)
    {
        if (*variable >= 0)
        {
            os << '(' << variable - m_basic_variables.begin() << ", "
                    << *variable << "), ";
        }
    }
    os << '\n';
}

template class SimplexLPSolver<float> ;
template class SimplexLPSolver<double> ;
