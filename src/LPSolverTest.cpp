/*
 * LPSolverTest.cpp
 *
 *  Created on: Apr 20, 2017
 *      Author: JM
 */

#include "LPSolverTest.h"

#include "SimplexSolver.h"

template<typename Scalar>
LPSolverTest<Scalar>::LPSolverTest(int test_id) :
        m_test_id(test_id), //
        m_num_free_variables(0), //
        m_verbose(false)
{
    switch (m_test_id) {

    case 1:
        resize(2, 3); // 2 constraints, 3 variables.
        m_sign = 1; // Objective function is to be minimized.
        m_A << 3, 2, 1, // Constraint lhs coefficients
        2, 5, 3;
        m_b << 10, 15; // Constraint rhs coefficients
        m_c << -2, -3, -4; // Linear objective function coefficients
        m_inequalities << 1, 1; // Both constraints are <= inequalities.
        m_known_solution << 0, 0, 5; // For validation.
        break;

    case 2:
        resize(2, 3);
        m_sign = 1;
        m_A << 3, 2, 1, //
        2, 5, 3;
        m_b << 10, 15;
        m_c << -2, -3, -4;
        // If m_inequalities is not specified, default is 0 which means that all constraints are equalities.
        m_known_solution << 15. / 7., 0, 25. / 7.;
        m_verbose = false;
        break;

    case 3:
        resize(3, 6);
        m_sign = -1; // Objective function is to be maximized.
        m_A << 2, 1, 1, 1, 0, 0, //
        1, 2, 3, 0, 1, 0, //
        2, 2, 1, 0, 0, 1;
        m_b << 2, 5, 6;
        m_c << 3, 2, 3, 0, 0, 0;
        // Again, all constraints are equalities.
        m_known_solution << 1. / 5., 0, 8. / 5., 0, 0, 4;
        break;

    case 4:
        resize(2, 3);
        m_sign = -1;
        m_A << 2, 1, -1, 1, 2, 0;
        m_b << 4, 6;
        m_c << 1, 1, 0;
        m_known_solution << 6, 0, 8;
        break;

    case 5:
        resize(2, 4);
        m_sign = 1;
        m_A << 1, 1, 0, -1, 0, -1, 1, -1;
        m_b << -1, -3;
        m_c << 2, 3, -1, 1;
        m_known_solution << 0, 1, 0, 2;
        break;

    case 6:  // Google 'Stigler diet" for the history of this test.
        resize(9, 77);
        m_sign = 1;
        m_A = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic,
                Eigen::RowMajor>::Map(Stigler_data[0], 77, 9).transpose();
        m_b = VecX::Map(Stigler_nutrients, 9);
        m_c.setOnes();
        m_inequalities.setConstant(-1); // All constraints are >= (minimum amount of nutrient required).
        m_known_solution = VecX::Map(Stigler_solution, 77);
        break;

    case 7:
        resize(3, 5);
        m_sign = 1;
        m_A << -1, 1, 2, 1, 2, -1, 2, 3, 1, 1, -1, 1, 1, 2, 1;
        m_b << 7, 6, 4;
        m_c << -2, 4, 7, 1, 5;
        m_num_free_variables = 1; // The first variable is free, i.e. the implicit >= 0 constraint is lifted.
        m_known_solution << -1, 0, 1, 0, 2;
        break;

    case 8: // Free cube.
        resize(6, 3);
        m_A << 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1;
        m_b << 1, 1, 1, -1, -1, -1;
        m_c << 1, 1, -2;
        m_sign = -1;
        m_num_free_variables = 3;
        m_inequalities << 1, 1, 1, -1, -1, -1;
        m_known_solution << 1, 1, -1;
        //m_verbose = true;
        break;

    case 9:
        resize(7, 3);
        m_A << 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1;
        m_b << 1, 1, 1, -1, -1, -1, 1;
        m_c << 1, .9, 1;
        m_sign = -1;
        m_num_free_variables = 3;
        m_inequalities << 1, 1, 1, -1, -1, -1, 1;
        m_known_solution << 1, -1, 1;
        //m_verbose = true;
        break;

    case 10:
        resize(5, 2);
        m_A << 1, 0, 0, 1, 1, 0, 0, 1, 2, 1;
        m_b << 1, 1, -1, -1, 1;
        m_inequalities << 1, 1, -1, -1, 1;
        m_c << 1, 1;
        m_sign = -1;
        m_num_free_variables = 2;
        m_known_solution << 0, 1;
        //m_verbose = true;
        break;

    case 11:
        resize(4, 3);
        m_A << 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1;
        m_b << 1, 1, 1, 2.5;
        m_inequalities << 1, 1, 1, 1;
        m_c << 1, .9, 1;
        m_sign = -1;
        m_known_solution << 1, 0.5, 1;
        m_verbose = false;
        break;

    case 12:
        resize(4, 3);
        m_A << 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1;
        m_b << 1, 1, 1, 2.5;
        m_inequalities << 1, 1, 1, -1;
        m_c << 1, 1.1, 1;
        m_sign = +1;
        m_known_solution << 1, .5, 1;
        m_verbose = false;
        break;

    case 13: // Free cube with cut corners.
        resize(14, 3);
        m_A << 1, 0, 0, 0, 1, 0, 0, 0, 1 //
        , 1, 0, 0, 0, 1, 0, 0, 0, 1 //
        , 1, 1, 1, 1, 1, 1 //
        , 1, 1, -1, 1, 1, -1 //
        , -1, 1, 1, -1, 1, 1 //
        , 1, -1, 1, 1, -1, 1;
        m_b << 1, 1, 1 //
        , -1, -1, -1 //
        , 2.75, 0.25 //
        , 2.75, 0.25 //
        , 2.75, 0.25 //
        , 2.75, 0.25;
        m_c << 1, 0, 0;
        m_sign = -1;
        m_num_free_variables = 3;
        m_inequalities << 1, 1, 1, -1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1;
        m_known_solution << 1, .75, 1;
        m_verbose = false;
        break;

    case 14: // Free square with cut corners.
        resize(8, 2);
        m_A << 1, 0, 0, 1 //
        , 1, 0, 0, 1 //
        , 1, 1, 1, 1 //
        , 1, -1, 1, -1;
        m_b << 1, 1 //
        , -1, -1 //
        , 1.5, -1.5 //
        , 1.5, -1.5;
        m_c << 1, 1;
        m_sign = 1;
        m_num_free_variables = 2;
        m_inequalities << 1, 1, -1, -1, 1, -1, 1, -1;
        m_known_solution << -1, -.5;
        m_verbose = false;
        break;

    }
}

template<typename Scalar>
void LPSolverTest<Scalar>::resize(int num_constraints, int num_variables)
{
    m_A.resize(num_constraints, num_variables);
    m_b.resize(num_constraints);
    m_c.resize(num_variables);
    m_inequalities.resize(num_constraints);
    m_known_solution.resize(num_variables);
}

template<typename Scalar>
bool LPSolverTest<Scalar>::execute()
{
    std::cout << "Test number " << m_test_id << ": ";

    SimplexSolver<Scalar> simplex_solver = SimplexSolver<Scalar>(m_A, m_b,
            m_sign * m_c, m_inequalities, m_num_free_variables);

    if (!simplex_solver.solve())
    {
        std::cout << "SimplexSolver solve failed!\n";
        return false;
    }

//    // Uncomment this to test reverse solve as well.
//    std::cout << '\n';
//    simplex_solver.reverse_solve();
//    std::cout << "Number of vertices found: "
//            << simplex_solver.getReverseSolveCounter() << '\n';

    typename SimplexSolver<Scalar>::VecX const solution =
            simplex_solver.get_solution();

    Scalar constexpr sq_tolerance = std::numeric_limits<Scalar>::epsilon();

    if (m_verbose)
    {
        std::cout << '\n';
        std::cout << "Solution:\t" << solution.transpose() << '\n';
        std::cout << "Known solution:\t" << m_known_solution.transpose()
                << '\n';
        std::cout << "Optimal value:\t"
                << m_sign * simplex_solver.get_optimal_value() << '\n';
        std::cout << "Known optimal:\t" << m_known_solution.dot(m_c) << '\n';
    }

    if ((solution - m_known_solution).squaredNorm() <= sq_tolerance
            && sqr(
                    m_sign * simplex_solver.get_optimal_value()
                            - m_known_solution.dot(m_c)) <= sq_tolerance)
    {
        std::cout << "OK\n";
        return true;
    } else
    {
        std::cout << "FAILED\n";
        return false;
    }
}

template<typename Scalar>
Scalar const LPSolverTest<Scalar>::Stigler_data[77][9] = { //
        { 44.7, 1411, 2, 365, 0, 55.4, 33.3, 441, 0 }, //
                { 11.6, 418, 0.7, 54, 0, 3.2, 1.9, 68, 0 }, //
                { 11.8, 377, 14.4, 175, 0, 14.4, 8.8, 114, 0 }, //
                { 11.4, 252, 0.1, 56, 0, 13.5, 2.3, 68, 0 }, //
                { 36., 897, 1.7, 99, 30.9, 17.4, 7.9, 106, 0 }, //
                { 28.6, 680, 0.8, 80, 0, 10.6, 1.6, 110, 0 }, //
                { 21.2, 460, 0.6, 41, 0, 2, 4.8, 60, 0 }, //
                { 25.3, 907, 5.1, 341, 0, 37.1, 8.9, 64, 0 }, //
                { 15., 488, 2.5, 115, 0, 13.8, 8.5, 126, 0 }, //
                { 12.2, 484, 2.7, 125, 0, 13.9, 6.4, 160, 0 }, //
                { 12.4, 439, 1.1, 82, 0, 9.9, 3, 66, 0 }, //
                { 8., 130, 0.4, 31, 18.9, 2.8, 3, 17, 0 }, //
                { 12.5, 288, 0.5, 50, 0, 0, 0, 0, 0 }, //
                { 6.1, 310, 10.5, 18, 16.8, 4, 16, 7, 177 }, //
                { 8.4, 422, 15.1, 9, 26, 3, 23.5, 11, 60 }, //
                { 10.8, 9, 0.2, 3, 44.2, 0, 0.2, 2, 0 }, //
                { 20.6, 17, 0.6, 6, 55.8, 0.2, 0, 0, 0 }, //
                { 2.9, 238, 1., 52, 18.6, 2.8, 6.5, 1, 0 }, //
                { 7.4, 448, 16.4, 19, 28.1, 0.8, 10.3, 4, 0 }, //
                { 3.5, 49, 1.7, 3, 16.9, 0.6, 2.5, 0, 17 }, //
                { 15.7, 661, 1., 48, 0, 9.6, 8.1, 471, 0 }, //
                { 8.6, 18, 0.2, 8, 2.7, 0.4, 0.5, 0, 0 }, //
                { 20.1, 0, 0, 0, 0, 0, 0, 0, 0 }, //
                { 41.7, 0, 0, 0, 0.2, 0, 0.5, 5, 0 }, //
                { 2.9, 166, 0.1, 34, 0.2, 2.1, 2.9, 69, 0 }, //
                { 2.2, 214, 0.1, 32, 0.4, 2.5, 2.4, 87, 0 }, //
                { 3.4, 213, 0.1, 33, 0, 0, 2, 0, 0 }, //
                { 3.6, 309, 0.2, 46, 0.4, 1, 4, 120, 0 }, //
                { 8.5, 404, 0.2, 62, 0, 0.9, 0, 0, 0 }, //
                { 2.2, 333, 0.2, 139, 169.2, 6.4, 50.8, 316, 525 }, //
                { 3.1, 245, 0.1, 20, 0, 2.8, 3.9, 86, 0 }, //
                { 3.3, 140, 0.1, 15, 0, 1.7, 2.7, 54, 0 }, //
                { 3.5, 196, 0.2, 30, 0, 17.4, 2.7, 60, 0 }, //
                { 4.4, 249, 0.3, 37, 0, 18.2, 3.6, 79, 0 }, //
                { 10.4, 152, 0.2, 23, 0, 1.8, 1.8, 71, 0 }, //
                { 6.7, 212, 0.2, 31, 0, 9.9, 3.3, 50, 0 }, //
                { 18.8, 164, 0.1, 26, 0, 1.4, 1.8, 0, 0 }, //
                { 1.8, 184, 0.1, 30, 0.1, 0.9, 1.8, 68, 46 }, //
                { 1.7, 156, 0.1, 24, 0, 1.4, 2.4, 57, 0 }, //
                { 5.8, 705, 6.8, 45, 3.5, 1, 4.9, 209, 0 }, //
                { 5.8, 27, 0.5, 36, 7.3, 3.6, 2.7, 5, 544 }, //
                { 4.9, 60, 0.4, 30, 17.4, 2.5, 3.5, 28, 498 }, //
                { 1., 21, 0.5, 14, 0, 0.5, 0, 4, 952 }, //
                { 2.2, 40, 1.1, 18, 11.1, 3.6, 1.3, 10, 1998 }, //
                { 2.4, 138, 3.7, 80, 69, 4.3, 5.8, 37, 862 }, //
                { 2.6, 125, 4., 36, 7.2, 9, 4.5, 26, 5369 }, //
                { 2.7, 73, 2.8, 43, 188.5, 6.1, 4.3, 89, 608 }, //
                { 0.9, 51, 3., 23, 0.9, 1.4, 1.4, 9, 313 }, //
                { 0.4, 27, 1.1, 22, 112.4, 1.8, 3.4, 11, 449 }, //
                { 5.8, 166, 3.8, 59, 16.6, 4.7, 5.9, 21, 1184 }, //
                { 14.3, 336, 1.8, 118, 6.7, 29.4, 7.1, 198, 2522 }, //
                { 1.1, 106, 0, 138, 918.4, 5.7, 13.8, 33, 2755 }, //
                { 9.6, 138, 2.7, 54, 290.7, 8.4, 5.4, 83, 1912 }, //
                { 3.7, 20, 0.4, 10, 21.5, 0.5, 1, 31, 196 }, //
                { 3., 8, 0.3, 8, 0.8, 0.8, 0.8, 5, 81 }, //
                { 2.4, 16, 0.4, 8, 2, 2.8, 0.8, 7, 399 }, //
                { 0.4, 33, 0.3, 12, 16.3, 1.4, 2.1, 17, 272 }, //
                { 1., 54, 2, 65, 53.9, 1.6, 4.3, 32, 431 }, //
                { 7.5, 364, 4, 134, 3.5, 8.3, 7.7, 56, 0 }, //
                { 5.2, 136, 0.2, 16, 12, 1.6, 2.7, 42, 218 }, //
                { 2.3, 136, 0.6, 45, 34.9, 4.9, 2.5, 37, 370 }, //
                { 1.3, 63, 0.7, 38, 53.2, 3.4, 2.5, 36, 1253 }, //
                { 1.6, 71, 0.6, 43, 57.9, 3.5, 2.4, 67, 862 }, //
                { 8.5, 87, 1.7, 173, 86.8, 1.2, 4.3, 55, 57 }, //
                { 12.8, 99, 2.5, 154, 85.7, 3.9, 4.3, 65, 257 }, //
                { 13.5, 104, 2.5, 136, 4.5, 6.3, 1.4, 24, 136 }, //
                { 20., 1367, 4.2, 345, 2.9, 28.7, 18.4, 162, 0 }, //
                { 17.4, 1055, 3.7, 459, 5.1, 26.9, 38.2, 93, 0 }, //
                { 26.9, 1691, 11.4, 792, 0, 38.4, 24.6, 217, 0 }, //
                { 0, 0, 0, 0, 0, 4, 5.1, 50, 0 }, //
                { 0, 0, 0, 0, 0, 0, 2.3, 42, 0 }, //
                { 8.7, 237, 3, 72, 0, 2, 11.9, 40, 0 }, //
                { 8., 77, 1.3, 39, 0, 0.9, 3.4, 14, 0 }, //
                { 34.9, 0, 0, 0, 0, 0, 0, 0, 0 }, //
                { 14.7, 0, 0.5, 74, 0, 0, 0, 5, 0 }, //
                { 9., 0, 10.3, 244, 0, 1.9, 7.5, 146, 0 }, //
                { 6.4, 11, 0.4, 7, 0.2, 0.2, 0.4, 3, 0 } //
        };

template<typename Scalar>
Scalar const LPSolverTest<Scalar>::Stigler_nutrients[9] = { 3, 70, 0.8, 12, 5,
        1.8, 2.7, 18, 75 };

template<typename Scalar>
Scalar const LPSolverTest<Scalar>::Stigler_solution[77] = {
        0.029519061676488264, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0.0018925572907052778, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0.01121443524614487, 0., 0., 0., 0., 0.,
        0.005007660466725201, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.06102856352669324, 0., 0., 0., 0., 0., 0., 0., 0. };

template class LPSolverTest<float> ;
template class LPSolverTest<double> ;
