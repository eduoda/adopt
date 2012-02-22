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
// Name        : simple_car.cpp
// Author      : Eduardo Oda
// Version     : 0.1
//============================================================================

#include "simple_car.hpp"

void Car::field(AD<double> t, vector< AD<double> >& x, vector< AD<double> >& u, vector< AD<double> >& y){
	y[0]=x[1];
	y[1]=u[0];
}
AD<double> Car::cost(vector< AD<double> >& x, vector< AD<double> >& u){
	return 1;
}
void Car::get_state_L(vector<double>& x_L){
	for(size_t i=0; i<get_M(); i++){
		x_L[i]=-1.0e19;
	}
}
void Car::get_state_U(vector<double>& x_U){
	for(size_t i=0; i<get_M(); i++){
		x_U[i]=1.0e19;
	}
}
void Car::get_control_L(vector<double>& u_L){
	for(size_t i=0; i<get_C(); i++){
		u_L[i]=-1;
	}
}
void Car::get_control_U(vector<double>& u_U){
	for(size_t i=0; i<get_C(); i++){
		u_U[i]=1;
	}
}
void Car::get_initial_state(vector<double>& x0){
	x0[0]=-1;
	x0[1]=0;
}
void Car::get_final_state(vector<double>& xF){
	xF[0]=0;
	xF[1]=0;
}
size_t Car::get_N(){
	return 100; // discretization steps
}
size_t Car::get_M(){
	return 2; // dimension of states
}
size_t Car::get_C(){
	return 1; // dimension of control
}
