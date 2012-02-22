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
// Name        : IpoptSolver.hpp
// Author      : Eduardo Oda
// Version     : 0.1
// Description : A simple Ipopt interface
// Prereqs     : CppAD
// 				 Boost uBlas
//============================================================================

#ifndef IPOPTSOLVER_HPP_
#define IPOPTSOLVER_HPP_

#include "adopt/adopt.hpp"
#include "IpTNLP.hpp"

using namespace Ipopt;

class IpoptSolver:public Solver, public TNLP{
public:
	// TODO: do we really need this constructor?
	IpoptSolver(Problem &p):Solver(p){}
	bool solve();
	double get_optval();
	void get_sol(vector< double >& x);
	void get_multipliers(vector< double >& l);

	// Ipopt Interface specific methods
	bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style);
	bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u);
	bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda);
	bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);
	bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);
	bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
	bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values);
	bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values);
	void finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq);
};

#endif /* IPOPTSOLVER_HPP_ */
