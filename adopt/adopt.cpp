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
// Name        : adop.cpp
// Author      : Eduardo Oda
// Version     : 0.1
//============================================================================

#include "adopt.hpp"
#include "time.h"

bool Problem::init(){
	bool ret;
	time_hess_L=0;
	time_F=0;
	time_grad_F=0;
	time_G=0;
	time_jac_G=0;

	ret = tape_f() && tape_g();

	F_grad_nne = F_grad_sparsity_pattern[0].size();

	G_jac_nne=0;
	for(size_t i=0; i<get_m(); i++){
		G_jac_nne+=G_jac_sparsity_pattern[i].size();
	}

	L_hess_nne=0;
	for(size_t i=0; i<get_n(); i++){
		L_hess_nne+=L_hess_sparsity_pattern[i].size();
	}

	return ret;
}

bool Problem::tape_f(){
	vector< AD<double> > X(get_n());
	get_initial_point_AD(X);
	return tape_f(X);
}

bool Problem::tape_f(vector< AD<double> > X){
	vector< AD<double> > Y(1);
	CppAD::Independent(X);
	eval_f(X, Y[0]);
	F.Dependent(X, Y);
	F.optimize();

	vector< std::set<size_t> > Id_n(get_n());
	for(size_t i = 0; i < get_n(); i++)
		Id_n[i].insert(i);
	vector< std::set<size_t> > Id_1(1);
	Id_1[0].insert(0);

	F.ForSparseJac(get_n(), Id_n);
	F_hess_sparsity_pattern = F.RevSparseHes(get_n(), Id_1);
	F_grad_sparsity_pattern = F.RevSparseJac(1, Id_1);

	for(size_t i=0; i<get_n(); i++){
		set<size_t> *s = &F_hess_sparsity_pattern[i];
		for(set<size_t>::iterator it=s->begin(); it!=s->end(); it++){
			if(*it>i){
				s->erase(it,s->end());
				break;
			}
		}
	}

	return true;
}

bool Problem::tape_g(){
	vector< AD<double> > X(get_n());
	get_initial_point_AD(X);
	return tape_g(X);
}

bool Problem::tape_g(vector< AD<double> > X){
	vector< AD<double> > Z(get_m());
	CppAD::Independent(X);
	eval_g(X, Z);
	G.Dependent(X, Z);
	G.optimize();

	vector< std::set<size_t> > Id_n(get_n());
	for(size_t i = 0; i < get_n(); i++)
		Id_n[i].insert(i);
	vector< std::set<size_t> > Id_m(get_m());
	for(size_t i = 0; i < get_m(); i++)
		Id_m[i].insert(i);
	vector< std::set<size_t> > Id_1(1);
	Id_1[0].insert(0);

	if( get_n() <= get_m() ) {
		// are you crazy? more restrictions then variables?!?
		// check if your problem satisfies the kkt conditions
		G_jac_sparsity_pattern = G.ForSparseJac(get_n(), Id_n);
	} else {
		G_jac_sparsity_pattern = G.RevSparseJac(get_m(), Id_m);
		G.ForSparseJac(get_n(), Id_n);
	}
	vector < vector< std::set<size_t> > > e_i(get_m());
	Gi_hess_sparsity_pattern.resize(get_m());
	for(size_t i = 0; i < get_m(); i++){
		e_i[i].resize(1);
		e_i[i][0].insert(i);
		Gi_hess_sparsity_pattern[i] = G.RevSparseHes(get_n(), e_i[i]);
		e_i[i][0].clear();
	}

	// we also calculate the Hess L sparsity pattern merging G and F sparsity patterns
	L_hess_sparsity_pattern.resize(get_n());
	G_hess_sparsity_pattern.resize(get_n());
	for(size_t i = 0; i < get_n(); i++){
		set<size_t> fs = F_hess_sparsity_pattern[i];
		for(size_t j = 0; j < get_m(); j++){
			set<size_t> gis = Gi_hess_sparsity_pattern[j][i];
			for(set<size_t>::iterator it=gis.begin(); it!=gis.end() && *it<=i; it++){
				L_hess_sparsity_pattern[i].insert(*it);
			}
		}
		G_hess_sparsity_pattern[i]=L_hess_sparsity_pattern[i];
		for(set<size_t>::iterator it=fs.begin(); it!=fs.end() && *it<=i; it++){
			L_hess_sparsity_pattern[i].insert(*it);
		}
	}
//	for(size_t i = 0; i < get_n(); i++){
//		set<size_t> s = F_hess_sparsity_pattern[i];
//		for(set<size_t>::iterator it=s.begin(); it!=s.end(); it++){
//			cout << *it << ", ";
//		}
//		cout << endl;
//	}

	return true;
}


bool Problem::f(vector<double>& x, double& y) {
	clock_t start, end;
	start = clock();

	y = (F.Forward(0,x))[0];

	end = clock();
	time_F+=(double)(end-start)/CLOCKS_PER_SEC;

	return true;
}

bool Problem::grad_f(vector<double>& x, vector<double>& grad) {
	clock_t start, end;
	start = clock();

	grad = F.SparseJacobian(x,F_grad_sparsity_pattern);

	end = clock();
	time_grad_F+=(double)(end-start)/CLOCKS_PER_SEC;

	return true;
}

bool Problem::grad_f_sparse(vector<double>& x, vector<double>& grad) {
	clock_t start, end;
	start = clock();

	F.Forward(0,x);
	std::set<size_t> s = F_grad_sparsity_pattern[0];
	grad.resize(s.size());
	size_t k=0;
	vector<double> dx(get_n(), 0.0);
	for(std::set<size_t>::iterator it=s.begin(); it!=s.end(); it++){
		dx[*it]=1.0;
		grad[k]=(F.Forward(1, dx))[0];
		dx[*it]=0.0;
		k++;
	}

	end = clock();
	time_grad_F+=(double)(end-start)/CLOCKS_PER_SEC;

	return true;
}

bool Problem::hess_f(vector<double>& x, vector<double>& hess) {
	vector<double> l(get_m(),0.0);
	hess_L(x,1.0,l,hess);
	return true;
}

bool Problem::hess_f_sparse(vector<double>& x, vector<double>& hess) {
	vector<double> l(get_m(),0.0);
	hess_L_sparse(x,1.0,l,hess);
//	vector<double> w(1,1.);
//	hess = F.SparseHessian(x,w,F_hess_sparsity_pattern);
//	hess = F.Hessian(x,w);
	return true;
}

bool Problem::g(vector<double>& x, vector<double>& z) {
	clock_t start, end;
	start = clock();

	z = G.Forward(0,x);

	end = clock();
	time_G+=(double)(end-start)/CLOCKS_PER_SEC;

	return true;
}

bool Problem::jac_g(vector<double>& x, vector<double>& jac) {
	clock_t start, end;
	start = clock();

	jac = G.SparseJacobian(x,G_jac_sparsity_pattern);

	end = clock();
	time_jac_G+=(double)(end-start)/CLOCKS_PER_SEC;

	return true;
}

bool Problem::jac_g_sparse(vector<double>& x, vector<double>& jac) {
	clock_t start, end;
	start = clock();

	G.Forward(0,x);
	vector<double> dw(get_m(), 0.0);
	vector<double> dy(get_m());
	size_t k=0;
	for(size_t i=0; i<get_m(); i++){
		std::set<size_t> s = G_jac_sparsity_pattern[i];
		dw[i]=1.0;
		dy=G.Reverse(1, dw);
		for(std::set<size_t>::iterator it=s.begin(); it!=s.end(); it++){
			jac[k]=dy[*it];
			k++;
		}
		dw[i]=0.0;
	}

	end = clock();
	time_jac_G+=(double)(end-start)/CLOCKS_PER_SEC;

	return true;
}

bool Problem::hess_g(vector<double>& x, vector<double>& l, vector<double>& hess){
	hess_L(x,0.0,l,hess);
	return true;
}

bool Problem::hess_g_sparse(vector<double>& x, vector<double>& l, vector<double>& hess){
	hess_L_sparse(x,0.0,l,hess);
	return true;
}

bool Problem::hess_L(vector<double>& x, double r, vector<double>& l, vector<double>& hess) {
	clock_t start, end;
	start = clock();

//	vector<double> w(1,r);
//	hess = F.SparseHessian(x,w,F_hess_sparsity_pattern);
//	hess = hess + G.SparseHessian(x,l);
//	hess = F.Hessian(x,w);
//	hess = hess + G.Hessian(x,l);

	vector<double> hess_sparse(L_hess_nne);
	hess_L_sparse(x,r,l,hess_sparse);
	size_t k=0;
	for(size_t i=0; i<get_n(); i++){
		for(size_t j=0; j<get_n(); j++)
			hess[i*get_n()+j]=0.0;
		set<size_t> s = L_hess_sparsity_pattern[i];
		for(set<size_t>::iterator it=s.begin(); it!=s.end(); it++){
			hess[i*get_n()+*it]=hess_sparse[k];
			hess[*it*get_n()+i]=hess_sparse[k];
			k++;
		}
	}

	end = clock();
	time_hess_L+=(double)(end-start)/CLOCKS_PER_SEC;
	return true;
}

bool Problem::hess_L_sparse(vector<double>& x, double r, vector<double>& l, vector<double>& hess) {
	clock_t start, end;
	start = clock();

	F.Forward(0,x);
	G.Forward(0,x);

	vector<double> ddG(2*get_n(),0.0);
	vector<double> dx(get_n(), 0.0);
	size_t k=0;
	for(size_t i=0; i<get_n(); i++){
		vector<double> ddF(2*get_n(),0.0);
		std::set<size_t> s = L_hess_sparsity_pattern[i];
		if(s.size()==0)
			continue;
		dx[i] = 1.0;
		G.Forward(1,dx);
		F.Forward(1,dx);
		ddG = G.Reverse(2,l);
		if(r!=0.0){
			vector<double> r_v(1,r);
			ddF = F.Reverse(2,r_v);
		}
		for(std::set<size_t>::iterator it=s.begin(); it!=s.end(); it++){
			hess[k]=ddF[*it*2 + 1] + ddG[*it*2 + 1];
			k++;
		}
		dx[i] = 0.0;
	}
	end = clock();
	time_hess_L+=(double)(end-start)/CLOCKS_PER_SEC;

	return true;
}

bool Problem::get_initial_point_AD(vector< AD<double> >& X){
	vector<double> x_v(get_n());
	get_initial_point(x_v);
	for(size_t i=0; i<get_n(); i++)
		X[i]=x_v[i];
	return true;
}
bool Problem::get_initial_multipliers(vector< double >& lambda){
	return false;
}


//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void Problem::process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}
