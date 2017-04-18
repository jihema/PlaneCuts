//
//  PlaneCuts.h
//  AvisFukuda
//
//  Created by Jean-Marie Aubry on 25/03/2017.
//  Copyright © 2017 Jean-Marie Aubry. All rights reserved.
//

#ifndef PlaneCuts_h
#define PlaneCuts_h

#include <iostream>
#include <set>

#include "EigenIncludes.h"

template<typename Scalar, int d>
class PlaneCuts {
public:
	PlaneCuts(Eigen::Matrix<double, Eigen::Dynamic, d> A_prime,
			Eigen::Matrix<double, Eigen::Dynamic, 1> b)
	{
		assert(A_prime.rows() == b.size());

		m_n0 = b.size();

		m_m = m_n0 + 1;
		m_n = m_n0 + d + 2;
		m_f = m_n - 2;
		m_g = m_n-1;

		m_A.resize(m_m, m_n);
		m_A.block(0, m_n0, 1, d + 1).setConstant(1.);
		m_A.block(1, 0, m_n0, m_n0).setIdentity();
		m_A.block(1, m_n0, m_n0, d) = A_prime;
		m_A.block(1, m_n0 + d + 1, m_n0, 1) = -b;

		for (int i = 0; i < m_n0; ++i) {
			m_B.push_back(i);
		}
		m_B.push_back(m_f);

		for (int i = 0; i < d; ++i) {
			m_N.push_back(m_n0 + i);
		}
		m_N.push_back(m_g);

//		std::cout << "A = \n" << m_A << '\n';
	}

	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> m_A;
	int m_n0, m_m, m_n, m_f, m_g;
	std::vector<int> m_B, m_N;
};

#endif /* PlaneCuts_h */
