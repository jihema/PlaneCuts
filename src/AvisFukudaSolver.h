//
//  AvisFukudaSolver.h
//  AvisFukuda
//
//  Created by Jean-Marie Aubry on 25/03/2017.
//  Copyright © 2017 Jean-Marie Aubry. All rights reserved.
//

#ifndef AvisFukudaSolver_h
#define AvisFukudaSolver_h

#include "Polytope.h"
#include "PlaneCuts.h"
#include <vector>
#include <ostream>

template<typename Scalar, int d>
class AvisFukudaSolver {
public:
	AvisFukudaSolver(PlaneCuts<Scalar, d> const& plane_cuts) :
			m_A(plane_cuts.m_A)
	{
		set_B(plane_cuts.m_B);
		m_f = m_B.size() + 1;
		m_g = m_f + 1;
	}

	Polytope<Scalar, d> get_polygon()
	{
		Polytope<Scalar, d> polytope;

		// Some work goes here.

		return polytope;
	}


	void pivot_Bland();
	void pivot_criss_cross();
	void pivot(int r, int s);
	void print(std::ostream& os) const;

	bool is_primal_feasible() const;
	bool is_dual_feasible() const;
	bool is_optimal() const{return is_primal_feasible() && is_dual_feasible();}

private:
	void find_optimum();
	void set_B(std::vector<int> const& B); ///< Sets B, its complementary N and computes the corresponding A_bar.
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> extract_B_matrix() const; /// Returns a matrix containing the B columns of A;

public:
	std::vector<int> m_B;
	std::vector<int> m_N;
	int m_f, m_g;
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> const& m_A;
	Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> m_A_bar; ///< Dimension b x (n-b).
};

template<typename Scalar, int d>
std::ostream& operator<<(std::ostream& os,
		AvisFukudaSolver<Scalar, d> const& afs)
{
	afs.print(os);
	return os;
}

#endif /* AvisFukudaSolver_h */
