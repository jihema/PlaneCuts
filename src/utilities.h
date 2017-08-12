/*
 * utilities.h
 *
 *  Created on: Apr 17, 2017
 *      Author: JM
 */

#pragma once

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <eigen3/Eigen/Eigen>

template<typename T>
inline std::ostream& operator<<(std::ostream& os, std::vector<T> x)
{
    std::copy(x.begin(), x.end(), std::ostream_iterator<T>(os, ", "));
    return os;
}
