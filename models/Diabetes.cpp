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
// Name        : Diabetes.cpp
// Author      : Eduardo Oda
// Version     : 0.1
//============================================================================

#include "Diabetes.hpp"

void Diabetes::field(AD<double> t, vector< AD<double> >& x, vector< AD<double> >& u, vector< AD<double> >& y){
	double M1 = 0.0009;
	double M2 = 0.0031;
	double M3 = 0.0415;
	y[0] = -M1*x[0]-M2*x[1];
	y[1] = -M3*x[1]+u[0];
}
AD<double> Diabetes::cost(vector< AD<double> >& x, vector< AD<double> >& u){
	return (exp((x[0]-180)/100)+exp((-x[0]+70))/100);
}
void Diabetes::get_state_L(vector<double>& x_L){
	for(size_t i=0; i<get_M(); i++){
		x_L[i]=-1.0e19;
	}
}
void Diabetes::get_state_U(vector<double>& x_U){
	for(size_t i=0; i<get_M(); i++){
		x_U[i]=1.0e19;
	}
}
void Diabetes::get_control_L(vector<double>& u_L){
	for(size_t i=0; i<get_C(); i++){
		u_L[i]=0;
	}
}
void Diabetes::get_control_U(vector<double>& u_U){
	for(size_t i=0; i<get_C(); i++){
		u_U[i]=100;
	}
}
void Diabetes::get_initial_state(vector<double>& x0){
	x0[0]=230;
	x0[1]=0;
}
void Diabetes::get_final_state(vector<double>& xF){
	xF[0]=100;
	xF[1]=0;
}
size_t Diabetes::get_N(){
	return 240; // discretization steps
}
size_t Diabetes::get_M(){
	return 2; // dimension of states
}
size_t Diabetes::get_C(){
	return 1; // dimension of control
}
