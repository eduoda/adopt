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
// Name        : IpoptSolver.cpp
// Author      : Eduardo Oda
// Version     : 0.1
//============================================================================

#include "IpoptSolver.hpp"
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
using namespace Ipopt;

bool IpoptSolver::solve(){
	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	// Initialize the IpoptApplication and process the options
	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded) {
		printf("\n\n*** Error during initialization!\n");
		return (int) status;
	}
	status = app->OptimizeTNLP(this);
	if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
		// Retrieve some statistics about the solve
		Index iter_count = app->Statistics()->IterationCount();
		printf("\n\n*** The problem solved in %d iterations!\n", iter_count);

		Number final_obj = app->Statistics()->FinalObjective();
		printf("\n\n*** The final value of the objective function is %e.\n", final_obj);

	}

	return status;
}

double IpoptSolver::get_optval(){
	return 0.;
}

void IpoptSolver::get_sol(vector< double >& x){

}

void IpoptSolver::get_multipliers(vector< double >& l){

}


// Ipopt specific methods
bool IpoptSolver::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style){
	n = problem->get_n();
	m = problem->get_m();

	nnz_jac_g = 0;
	for(int i=0; i<m; i++){
		std::set<size_t> s = problem->G_jac_sparsity_pattern[i];
		nnz_jac_g+=s.size();
	}

	nnz_h_lag=0;
	for(int i=0; i<n; i++){
		std::set<size_t> s = problem->L_hess_sparsity_pattern[i];
		nnz_h_lag+=s.size();
	}

	// TODO: what are the styles?  FORTRAN_STYLE
	index_style = TNLP::C_STYLE;

	return true;
}

bool IpoptSolver::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u){
	vector<double> x(n);
	problem->get_x_L(x);
	std::copy(x.begin(),x.end(),x_l);
	x.clear();
	problem->get_x_U(x);
	std::copy(x.begin(),x.end(),x_u);
	x.clear();

	vector<double> z(m);
	problem->get_g_L(z);
	std::copy(z.begin(),z.end(),g_l);
	z.clear();
	problem->get_g_U(z);
	std::copy(z.begin(),z.end(),g_u);
	z.clear();

	return true;
}

bool IpoptSolver::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda){
	// Here, we assume we only have starting values for x, if you code
	// your own NLP, you can provide starting values for the others if
	// you wish.
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);

	vector<double> x_v(problem->get_n());
	problem->get_initial_point(x_v);
	std::copy(x_v.begin(),x_v.end(),x);

	return true;
}

bool IpoptSolver::eval_f(Index n, const Number* x, bool new_x, Number& obj_value){
	vector<double> x_v(problem->get_n());
	for(int i=0; i<n; i++)
		x_v[i]=x[i];
	double y;
	problem->f(x_v,y);
	obj_value = y;

	return true;
}

bool IpoptSolver::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){
	vector<double> x_v(problem->get_n());
	for(int i=0; i<n; i++)
		x_v[i]=x[i];

	vector<double> grad(problem->get_n());
	problem->grad_f(x_v,grad);
	for(int i=0; i<n; i++)
		grad_f[i]=grad[i];

	return true;
}

bool IpoptSolver::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
	vector<double> x_v(n);
	for(int i=0; i<n; i++)
		x_v[i]=x[i];
	vector<double> z_v(m);
	problem->g(x_v,z_v);
	for(int i=0; i<m; i++)
		g[i]=z_v[i];
	return true;
}

bool IpoptSolver::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values){
	if (values == NULL) {
		int k = 0;
		for(int i=0; i<m; i++){
			std::set<size_t> s = problem->G_jac_sparsity_pattern[i];
			for(std::set<size_t>::iterator it=s.begin(); it!=s.end(); it++) {
				iRow[k] = i;
				jCol[k] = (int)*it;
				k++;
			}
		}
	} else {
		vector<double> x_v(n);
		for(int i=0; i<n; i++)
			x_v[i]=x[i];

		vector<double> jac(nele_jac);
		problem->jac_g_sparse(x_v,jac);
		for(int i=0; i<nele_jac; i++)
			values[i]=jac[i];

//		vector<double> jac(n*m);
//		problem->jac_g(x_v,jac);
//		int k = 0;
//		for(int i=0; i<m; i++){
//			std::set<size_t> s = problem->G_jac_sparsity_pattern[i];
//			for(std::set<size_t>::iterator it=s.begin(); it!=s.end(); it++) {
//				values[k] = jac[i*n + (int)*it];
//				k++;
//			}
//		}

	}
	return true;
}

bool IpoptSolver::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values){
	if (values == NULL) {
		int k = 0;
		for(int i=0; i<n; i++){
			std::set<size_t> s = problem->L_hess_sparsity_pattern[i];
			for(std::set<size_t>::iterator it=s.begin(); it!=s.end(); it++) {
				iRow[k] = i;
				jCol[k] = (int)*it;
				k++;
			}
		}
	} else {
		vector<double> x_v(n);
		for(int i=0; i<n; i++)
			x_v[i]=x[i];
		vector<double> l_v(m);
		for(int i=0; i<m; i++)
			l_v[i]=lambda[i];

		vector<double> hess(nele_hess);
		problem->hess_L_sparse(x_v,obj_factor,l_v,hess);
		for(int i=0; i<nele_hess; i++)
			values[i]=hess[i];

//		vector<double> hess(n*n);
//		problem->hess_L(x_v,obj_factor,l_v,hess);
//		int k = 0;
//		for(int i=0; i<n; i++){
//			std::set<size_t> s = problem->L_hess_sparsity_pattern[i];
//			for(std::set<size_t>::iterator it=s.begin(); it!=s.end(); it++) {
//				values[k] = hess[i*n + (int)*it];
//				k++;
//			}
//		}

	}
	return true;
}

void IpoptSolver::finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq){
	cout << endl;
	cout << "X=[" << endl;
	for(int i=0; i<n; i++){
		cout << x[i] << ";";
	}
	cout << "];" << endl;

	cout << "L=[" << endl;
	for(int i=0; i<m; i++){
		cout << lambda[i] << "; ";
	}
	cout << "];" << endl;
}
