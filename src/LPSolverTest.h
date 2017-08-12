/*
 * LPSolverTest.h
 *
 *  Created on: Apr 20, 2017
 *      Author: JM
 */

#pragma once

#include "utilities.h"

template<typename Scalar>
inline Scalar sqr(Scalar x)
{
    return x * x;
}

template<typename Scalar>
class LPSolverTest {
public:
    LPSolverTest(int test);

    static int get_num_tests()
    {
        return num_tests;
    }

    bool execute();

private:
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MatX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    void resize(int num_constraints, int num_variables);

    int m_sign; ///< Set to +1 for minimization, -1 for maximization of the objective function.
    int m_test_id;

    MatX m_A;
    VecX m_b;
    VecX m_c;
    VecX m_inequalities;
    int m_num_free_variables;
    VecX m_known_solution;
    bool m_verbose;

    static const int num_tests = 9;
    static const Scalar Stigler_data[77][9];
    static const Scalar Stigler_nutrients[9];
    static const Scalar Stigler_solution[77];
};

