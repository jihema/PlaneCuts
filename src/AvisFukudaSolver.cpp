//
//  AvisFukudaSolver.cpp
//  AvisFukuda
//
//  Created by Jean-Marie Aubry on 25/03/2017.
//  Copyright © 2017 Jean-Marie Aubry. All rights reserved.
//

#include "AvisFukudaSolver.h"

template<typename Scalar, int d>
void AvisFukudaSolver<Scalar, d>::set_B(std::vector<int> const& B)
{
	// Store basis indices.
	m_B = B;

	auto const AB = extract_B_matrix();
//	std::cout << "AB = \n" << AB <<'\n';

// Basis solver.
	Eigen::FullPivLU<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> B_solver(
			AB);

	m_A_bar.resize(m_B.size(), m_A.cols() - m_B.size());

	// Store complementary indices and solve for the corresponding columns of m_A_bar.
	m_N.clear();
	for (int j = 0; j < m_A.cols(); ++j)
	{
		if (std::find(m_B.begin(), m_B.end(), j) == m_B.end())
		{
			m_A_bar.col(m_N.size()) = B_solver.solve(m_A.col(j));
			m_N.push_back(j);
		}
	}
}

template<typename Scalar, int d>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> AvisFukudaSolver<Scalar, d>::extract_B_matrix() const
{
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> B_matrix(m_B.size(),
			m_B.size());

	int j = 0;
	for (auto i = m_B.begin(); i != m_B.end(); ++i)
	{
		B_matrix.col(j++) = m_A.col(*i);
	}

	return B_matrix;
}

template<typename Scalar, int d>
void AvisFukudaSolver<Scalar, d>::print(std::ostream& os) const
{
	os << "B = ";
	std::copy(m_B.begin(), m_B.end(), std::ostream_iterator<int>(os, " "));
	std::cout << " (f = " << m_f << ")\n";

	os << "N = ";
	std::copy(m_N.begin(), m_N.end(), std::ostream_iterator<int>(os, " "));
	std::cout << " (g = " << m_g << ")\n";

	os << "A_bar = \n" << m_A_bar << '\n';

	std::cout << "Primal feasible: " << is_primal_feasible() << '\n';
	std::cout << "Dual feasible: " << is_dual_feasible() << '\n';

	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> x_N(m_N.size());
	x_N[m_N.size() - 1] = 1;
	os << "Basic solution:\n" << (m_A_bar * x_N).head(d) << '\n';
}

template<typename Scalar, int d>
bool AvisFukudaSolver<Scalar, d>::is_primal_feasible() const
{
	std::vector<int>::const_iterator const found_g = std::find(m_N.begin(),
			m_N.end(), m_g);
	if (found_g == m_N.end())
	{
		return false; // In fact an exception.
	}
	int const index_g = found_g - m_N.begin();

	for (int i = 0; i < m_B.size(); ++i)
	{
		if (m_B[i] != m_f && m_A_bar(i, index_g) < 0)
		{
			return false;
		}
	}
	return true;
}

template<typename Scalar, int d>
bool AvisFukudaSolver<Scalar, d>::is_dual_feasible() const
{
	std::vector<int>::const_iterator const found_f = std::find(m_B.begin(),
			m_B.end(), m_f);
	if (found_f == m_B.end())
	{
		return false; // In fact an exception.
	}
	int const index_f = found_f - m_B.begin();

	for (int j = 0; j < m_N.size(); ++j)
	{
		if (m_N[j] != m_g && m_A_bar(index_f, j) > 0)
		{
			return false;
		}
	}
	return true;
}

template<typename Scalar, int d>
void AvisFukudaSolver<Scalar, d>::pivot_Bland()
{
	std::vector<int>::const_iterator const found_f = std::find(m_B.begin(),
			m_B.end(), m_f);
	if (found_f == m_B.end())
	{
		return;
	}
}

template<typename Scalar, int d>
void AvisFukudaSolver<Scalar, d>::pivot_criss_cross()
{
	if (is_optimal())
	{
		return;
	}

	std::vector<int>::const_iterator const found_g = std::find(m_N.begin(),
			m_N.end(), m_g);
	if (found_g == m_N.end())
	{
		std::cout << "Criss cross rule failed\n";
		return; // In fact an exception.
	}
	int const index_g = found_g - m_N.begin();

	for (int index_r = 0; index_r < m_B.size(); ++index_r)
	{
		if (m_B[index_r] != m_f && m_A_bar(index_r, index_g) < 0) // Found a primal unfeasible
		{
			std::cout << "Found primal unfeasible r = " << m_B[index_r] << '\n';
			std::cout << m_A_bar.row(index_r)<<'\n';
			for (int index_s = 0; index_s < m_N.size(); ++index_s)
			{
				if (m_A_bar(index_r, index_s) > 0)
				{
					std::cout << "Pivot on primal unfeasible " << index_r << " "
							<< index_s << '\n';
					pivot(m_B[index_r], m_N[index_s]);
					return;
				}
			}
			std::cout << "Nope\n";
		}
	}

	///

	std::vector<int>::const_iterator const found_f = std::find(m_B.begin(),
			m_B.end(), m_f);
	if (found_f == m_B.end())
	{
		std::cout << "Criss cross rule failed\n";
		return;
	}
	int const index_f = found_f - m_B.begin();

	for (int index_s = 0; index_s < m_N.size(); ++index_s)
	{
		if (m_N[index_s] != m_g && m_A_bar(index_f, index_s) > 0)
		{
			std::cout << "Found dual unfeasible s = " << m_N[index_s] << '\n';
			std::cout << m_A_bar.col(index_s) <<'\n';
			for (int index_r = 0; index_r < m_B.size(); ++index_r)
			{
				if (m_A_bar(index_r, index_s) < 0)
				{
					std::cout << "Pivot on dual unfeasible " << index_r << " "
							<< index_s << '\n';
					pivot(m_B[index_r], m_N[index_s]);
					return;
				}
			}
			std::cout << "Nope\n";
		}
	}

	std::cout << "Criss cross rule failed\n";

}

template<typename Scalar, int d>
void AvisFukudaSolver<Scalar, d>::pivot(int r, int s)
{
	std::cout << "Pivoting r = " << r << " and s = " << s << '\n';

	std::vector<int>::iterator const found_r = std::find(m_B.begin(), m_B.end(),
			r);
	std::vector<int>::iterator const found_s = std::find(m_N.begin(), m_N.end(),
			s);
	if (found_r == m_B.end() || found_s == m_N.end()) // r must be in B and s must be in N.
	{
		return;
	}
	int const r_index = found_r - m_B.begin();
	int const s_index = found_s - m_N.begin();

	// Update A_bar.
	auto A_bar = m_A_bar;
	A_bar(r_index, s_index) = 1. / m_A_bar(r_index, s_index);

	for (int i_B = 0; i_B < m_B.size(); ++i_B)
	{
		if (i_B == r_index)
		{
			continue;
		}

		A_bar(i_B, s_index) *= A_bar(r_index, s_index);

		for (int j_N = 0; j_N < m_N.size(); ++j_N)
		{
			if (j_N == s_index)
			{
				continue;
			}

			A_bar(i_B, j_N) -= m_A_bar(i_B, s_index) * m_A_bar(r_index, j_N)
					* A_bar(r_index, s_index);
		}

	}
	for (int j_N = 0; j_N < m_N.size(); ++j_N)
	{
		if (j_N == s_index)
		{
			continue;
		}
		A_bar(r_index, j_N) *= -A_bar(r_index, s_index);
	}

	// Finalize swap.
	m_A_bar = A_bar;
	*found_s = r;
	*found_r = s;
}

template class AvisFukudaSolver<double, 2> ;
template class AvisFukudaSolver<double, 3> ;
