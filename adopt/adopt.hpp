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
// Name        : adop.hpp
// Author      : Eduardo Oda
// Version     : 0.1
// Description : A nonlinear problem:
//					min	 f(x)
//					s.t. g_L <= g(x) <= g_U
//						 x_L <=   x  <= x_U
//					x \in R^n
//					f: R^n --> R
//					g: R^n --> R^m
//				 can be solved by implementation of the classes: Problem and
//				 Solver.
//
// 				 Your implementations of Problem can optionally redefine the
// 				 get_initial_multipliers method to set initial values to
// 				 Lagrangian multipliers.
//
// 				 Derivatives are evaluated through Automatic Derivation (CppAD).
//
// 				 There are implementations of Solver to the following solvers:
// 					Ipopt
//
// Prereqs     : CppAD
// 				 Boost uBlas
//============================================================================

#ifndef ODA_HPP_
#define ODA_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cppad/cppad.hpp>

#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

using std::cout;
using std::endl;
using std::set;
using CppAD::AD;
using CppAD::ADFun;
using boost::numeric::ublas::vector;

class Problem{
private:
	ADFun<double> F,G;
	bool get_initial_point_AD(vector< AD<double> >& x);

public:
	double time_hess_L;
	double time_F;
	double time_grad_F;
	double time_G;
	double time_jac_G;
	vector< set<size_t> > F_grad_sparsity_pattern;
	size_t F_grad_nne;
	vector< set<size_t> > F_hess_sparsity_pattern;
	vector< set<size_t> > G_jac_sparsity_pattern;
	size_t G_jac_nne;
	vector< vector< set<size_t> > > Gi_hess_sparsity_pattern;
	vector< set<size_t> > G_hess_sparsity_pattern;
	vector< set<size_t> > L_hess_sparsity_pattern;
	size_t L_hess_nne;
	bool init();
	bool tape_f();
	bool tape_f(vector< AD<double> > X);
	bool tape_g();
	bool tape_g(vector< AD<double> > X);

	bool f(vector<double>& x, double& y);
	bool grad_f(vector<double>& x, vector<double>& grad);
	bool grad_f_sparse(vector<double>& x, vector<double>& grad);
	bool hess_f(vector<double>& x, vector<double>& hess);
	bool hess_f_sparse(vector<double>& x, vector<double>& hess);
	bool g(vector<double>& x, vector<double>& z);
	bool jac_g(vector<double>& x, vector<double>& jac);
	bool jac_g_sparse(vector<double>& x, vector<double>& jac);
	bool hess_g(vector<double>& x, vector<double>& l, vector<double>& hess);
	bool hess_g_sparse(vector<double>& x, vector<double>& l, vector<double>& hess);
	// L = r.f + < l , g > : Lagrangian
	bool hess_L(vector<double>& x, double r, vector<double>& l, vector<double>& hess);
	bool hess_L_sparse(vector<double>& x, double r, vector<double>& l, vector<double>& hess);

	virtual bool eval_f(vector< AD<double> >& x, AD<double>& y) =0;
	virtual bool eval_g(vector< AD<double> >& x, vector< AD<double> >& z) =0;
	virtual bool get_initial_point(vector< double >& x) =0;
	virtual bool get_initial_multipliers(vector< double >& lambda);
	virtual bool get_g_L(vector<double>& g_L) =0;
	virtual bool get_g_U(vector<double>& g_U) =0;
	virtual bool get_x_L(vector<double>& x_L) =0;
	virtual bool get_x_U(vector<double>& x_U) =0;
	virtual size_t get_n() =0;
	virtual size_t get_m() =0;
	void process_mem_usage(double& vm_usage, double& resident_set);
};

class Solver{
protected:
	Problem *problem;
public:
	Solver(Problem &p){
		problem = &p;
		problem->init();
	}
	virtual bool solve() =0;
	virtual double get_optval() =0;
	virtual void get_sol(vector< double >& x) =0;
	virtual void get_multipliers(vector< double >& l) =0;
};

#endif /* ODA_HPP_ */
