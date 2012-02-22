//============================================================================
// Adopt - Environment for nonlinear optimization
// Copyright (C) 2010 Eduardo Oda
//
// This file is part of Adopt.
//
// Adopt is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Adopt is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Adopt.  If not, see <http://www.gnu.org/licenses/>.
//
//
// Name        : optcontrol.hpp
// Author      : Eduardo Oda
// Version     : 0.1
// Description : Consider the optimum control problem:
//
// 		        _T
// 		       /
// 		min   / c(x,u) dt
// 		    _/
// 		     0
//
//		s.t. x'=f(x,u)
//			 x_L < x < x_U
// 			 u_L < u < u_U
// 			 x(0)=x0
// 			 x(T)=xF
//
//		x \in R^M      u \in R^C
//
// 		We use a approximation of the cost function and a numeric integration
//	of the dynamical restriction to obtain the nonlinear problem:
//
//		  x0   x1            ...             xN=xF
//		    u0   u1          ...          u6
//		  |----|----|----|----|----|----|----|------------>t
//		  0    1    2    3    4    5    6    7=T     (N=7)
//        /----/
//		    h = T/N
//
//            __N-1
//		min  \    c(x_k,u_k)*h
//           /__
//            k=0
//	    s.t.
//			 x_0=x0
//			 x_{k+1} = x_k + h*PHI(x_k,u_k,t_k), k=0..N-1
//			 x_N=xF
//			 x_L < x_k < x_U
//			 u_L < u_k < u_U
//
//		x=(h, x_0, x_1, ..., x_N, x_N, u_0, ..., u_{N-1})
//
//		Where PHI denotes the Classic Runge-Kutta method.
//
// Prereqs     : CppAD
// 				 Boost uBlas
//============================================================================

#ifndef OPTCONTROLPROBLEM_HPP_
#define OPTCONTROLPROBLEM_HPP_

#include "adopt.hpp"
#include "boost/numeric/ublas/vector_proxy.hpp"

class OptControlProblem : public Problem{
private:
	// Helper methods
	void get_x_k(vector< AD<double> >& x, size_t k, vector< AD<double> >& x_k);
	void set_x_k(vector<double>& x, size_t k, vector<double>& x_k);
	void get_u_k(vector< AD<double> >& x, size_t k, vector< AD<double> >& u_k);
	void set_u_k(vector<double>& x, size_t k, vector<double>& u_k);
	AD<double> get_h(vector< AD<double> >& x);
	void set_h(vector<double>& x, double h);

	// Some integration methods
	void EULER(vector< AD<double> >& x, size_t k, vector< AD<double> >& phi_k);
	void HEUN(vector< AD<double> >& x, size_t k, vector< AD<double> >& phi_k);
	void RK44(vector< AD<double> >& x, size_t k, vector< AD<double> >& phi_k);

	// These can/must be provided by the user
	virtual void field(AD<double> t, vector< AD<double> >& x, vector< AD<double> >& u, vector< AD<double> >& y) =0;
	virtual AD<double> cost(vector< AD<double> >& x, vector< AD<double> >& u) =0;
	virtual void get_state_L(vector<double>& x_L) =0;
	virtual void get_state_U(vector<double>& x_U) =0;
	virtual void get_control_L(vector<double>& u_L) =0;
	virtual void get_control_U(vector<double>& u_U) =0;
	virtual void get_initial_state(vector<double>& x0) =0;
	virtual void get_final_state(vector<double>& xF) =0;
	virtual void get_initial_trajectory(vector<double>& x, vector<double>& u);
	virtual size_t get_N() =0;
	virtual size_t get_M() =0;
	virtual size_t get_C() =0;

public:
	// Methods from Problem class
	bool eval_f(vector< AD<double> >& x, AD<double>& y);
	bool eval_g(vector< AD<double> >& x, vector< AD<double> >& z);
	bool get_g_L(vector<double>& g_L);
	bool get_g_U(vector<double>& g_U);
	bool get_x_L(vector<double>& x_L);
	bool get_x_U(vector<double>& x_U);
	virtual bool get_initial_point(vector<double>& x0);
	size_t get_n();
	size_t get_m();
};

#endif /* OPTCONTROLPROBLEM_HPP_ */
