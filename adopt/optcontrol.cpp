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
// Name        : optcontrol.cpp
// Author      : Eduardo Oda
// Version     : 0.1
//============================================================================

#include "optcontrol.hpp"

bool OptControlProblem::eval_f(vector< AD<double> >& x, AD<double>& y) {
	y=0;
	vector< AD<double> > x_k(get_M());
	vector< AD<double> > u_k(get_C());
	for(size_t k=0; k<get_N(); k++){
		get_x_k(x,k,x_k);
		get_u_k(x,k,u_k);
		y+=cost(x_k,u_k)*get_h(x);
	}
	return true;
}

bool OptControlProblem::eval_g(vector< AD<double> >& x, vector< AD<double> >& z) {
	vector<double> x0(get_M());
	vector<double> xF(get_M());
	vector< AD<double> > x_k(get_M());
	vector< AD<double> > x_k_1(get_M());
	vector< AD<double> > phi_k(get_M());
	get_initial_state(x0);
	get_final_state(xF);
	get_x_k(x,0,x_k);
	subrange(z,0,get_M()) = x0 - x_k;
	for(size_t k=0;k<get_N();k++){
		get_x_k(x,k,x_k);
		get_x_k(x,k+1,x_k_1);
//		EULER(x, k, phi_k);
		RK44(x, k, phi_k);
		subrange(z,(k+1)*get_M(),(k+1)*get_M()+get_M()) = x_k_1 - x_k - get_h(x)*phi_k;
	}
	get_x_k(x,get_N(),x_k);
	subrange(z,(get_N()+1)*get_M(),(get_N()+2)*get_M()) = xF - x_k;

	return true;
}

bool OptControlProblem::get_g_L(vector<double>& g_L) {
	for(size_t i=0;i<get_m();i++)
		g_L[i]=0.;
	return true;
}

bool OptControlProblem::get_g_U(vector<double>& g_U) {
	for(size_t i=0;i<get_m();i++)
		g_U[i]=0.;
	return true;
}

bool OptControlProblem::get_x_L(vector<double>& x_L) {
	set_h(x_L,0.0);
//	x_L[0]=0.0;
	vector<double> v1(get_M());
	vector<double> v2(get_C());
	get_state_L(v1);
	get_control_L(v2);
	for(size_t k=0;k<get_N();k++){
		set_x_k(x_L,k,v1);
		set_u_k(x_L,k,v2);
	}
	set_x_k(x_L,get_N(),v1);
	return true;
}

bool OptControlProblem::get_x_U(vector<double>& x_U) {
	set_h(x_U,1.0);
//	x_U[0]=1.0;
	vector<double> v1(get_M());
	vector<double> v2(get_C());
	get_state_U(v1);
	get_control_U(v2);
	for(size_t k=0;k<get_N();k++){
		set_x_k(x_U,k,v1);
		set_u_k(x_U,k,v2);
	}
	set_x_k(x_U,get_N(),v1);
	return true;
}

size_t OptControlProblem::get_n() {
	return get_N()*(get_M()+get_C())+get_M()+1;
}

size_t OptControlProblem::get_m() {
	return (get_N()+2)*get_M();
}

void OptControlProblem::EULER(vector< AD<double> >& x, size_t k, vector< AD<double> >& phi_k){
	AD<double> t_k;
	vector< AD<double> > x_k(get_M());
	vector< AD<double> > u_k(get_C());
	t_k=get_h(x);
	get_x_k(x,k,x_k);
	get_u_k(x,k,u_k);
	field(t_k,x_k,u_k,phi_k);
}

void OptControlProblem::HEUN(vector< AD<double> >& x, size_t k, vector< AD<double> >& phi_k){
	AD<double> t_k;
	vector< AD<double> > x_k(get_M());
	vector< AD<double> > u_k(get_C());
	vector< AD<double> > u_k_1(get_C());
	t_k=get_h(x);
	get_x_k(x,k,x_k);
	get_u_k(x,k,u_k);
	if(k<get_N()-1){
		get_u_k(x,k+1,u_k_1);
	} else {
		u_k_1=u_k;
	}
	vector< AD<double> > K1(get_M());
	vector< AD<double> > K2(get_M());
	vector< AD<double> > aux(get_M());
	field(t_k, x_k, u_k, K1);
	aux = x_k + get_h(x)*K1;
	field(t_k+get_h(x), aux, u_k_1, K2);
	phi_k = 0.5*(K1+K2);
}

void OptControlProblem::RK44(vector< AD<double> >& x, size_t k, vector< AD<double> >& phi_k){
	AD<double> t_k;
	vector< AD<double> > x_k(get_M());
	vector< AD<double> > u_k(get_C());
	vector< AD<double> > u_k_1(get_C());
	t_k=get_h(x);
	get_x_k(x,k,x_k);
	get_u_k(x,k,u_k);
	if(k<get_N()-1){
		get_u_k(x,k+1,u_k_1);
	} else {
		u_k_1=u_k;
	}

	vector< AD<double> > K1(get_M());
	vector< AD<double> > K2(get_M());
	vector< AD<double> > K3(get_M());
	vector< AD<double> > K4(get_M());

	vector< AD<double> > aux(get_M());
	field(t_k, x_k, u_k, K1);
	aux = x_k + 0.5*get_h(x)*K1;
	field(t_k+0.5*get_h(x), aux, u_k, K2);
	aux = x_k + 0.5*get_h(x)*K2;
	field(t_k+0.5*get_h(x), aux, u_k, K3);
	aux = x_k + get_h(x)*K3;
	field(t_k+get_h(x), aux, u_k_1, K4);
	phi_k = 0.5*(K1+K2);
}

void OptControlProblem::get_x_k(vector< AD<double> >& x, size_t k, vector< AD<double> >& x_k){
//	using namespace boost::numeric::ublas;
//	size_t a = k*(get_M()+get_C())+1;
//	size_t b = a+get_M();
	size_t a = k*get_M();
	size_t b = a+get_M();

	size_t i,j=0;
	for(i=a; i<b; i++){
		x_k[j]=x[i];
		j++;
	}
//	x_k = subrange(x,a,b);
}

void OptControlProblem::set_x_k(vector<double>& x, size_t k, vector<double>& x_k){
//	using namespace boost::numeric::ublas;
//	size_t a = k*(get_M()+get_C())+1;
//	size_t b = a+get_M();
	size_t a = k*get_M();
	size_t b = a+get_M();
//	subrange(x,a,b)=x_k;
	size_t i,j=0;
	for(i=a; i<b; i++){
		x[i]=x_k[j];
		j++;
	}
}

//(x01, x02, x03, x11, x12, x13, x21, x22, x23, u01, u02, u11, u12, h)
//(0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13)
//(h, x01, x02, x03, u01, u02, x11, x12, x13, u11, u12, x21, x22, x23, u21, u22, ...)
//(0, 1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15,  ...)
void OptControlProblem::get_u_k(vector< AD<double> >& x, size_t k, vector< AD<double> >& u_k){
	using namespace boost::numeric::ublas;
//	size_t a = k*(get_M()+get_C())+get_M()+1;
//	size_t b = a+get_C();
	size_t a = (get_N()+1)*get_M()+k*get_C();
	size_t b = a+get_C();
//	u_k = subrange(x,a,b);
	size_t i,j=0;
	for(i=a; i<b; i++){
		u_k[j]=x[i];
		j++;
	}
}

void OptControlProblem::set_u_k(vector<double>& x, size_t k, vector<double>& u_k){
	using namespace boost::numeric::ublas;
//	size_t a = k*(get_M()+get_C())+get_M()+1;
//	size_t b = a+get_C();
	size_t a = (get_N()+1)*get_M()+k*get_C();
	size_t b = a+get_C();
//	subrange(x,a,b)=u_k;
	size_t i,j=0;
	for(i=a; i<b; i++){
		x[i]=u_k[j];
		j++;
	}

}

AD<double> OptControlProblem::get_h(vector< AD<double> >& x){
//	return x[0];
	return x[(get_N()+1)*get_M() + get_N()*get_C()];
}

void OptControlProblem::set_h(vector<double>& x, double h){
	x[(get_N()+1)*get_M() + get_N()*get_C()]=h;
//	x[0]=h;
}

void OptControlProblem::get_initial_trajectory(vector<double>& x, vector<double>& u){
	for(size_t i=0; i<get_M()*(get_N()+1); i++){
		x[i]=0.;
	}
	for(size_t i=0; i<get_C()*get_N(); i++){
		u[i]=0.;
	}
}

bool OptControlProblem::get_initial_point(vector<double>& x) {
	vector<double> x_I(get_M()*(get_N()+1));
	vector<double> u_I(get_C()*get_N());
	get_initial_trajectory(x_I, u_I);
	vector<double> x_k(get_M());
	vector<double> u_k(get_C());
	for(size_t k=0; k<get_N(); k++){
		x_k = subrange(x_I,k*get_M(),(k+1)*get_M());
		u_k = subrange(u_I,k*get_C(),(k+1)*get_C());
		set_x_k(x,k,x_k);
		set_u_k(x,k,u_k);
	}
	x_k = subrange(x_I,get_N()*get_M(),(get_N()+1)*get_M());
	set_x_k(x,get_N(),x_k);
	set_h(x, 0.2);
//	x[0]=0.2;
	return true;
}
