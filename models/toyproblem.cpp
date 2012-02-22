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
// Name        : toyproblem.cpp
// Author      : Eduardo Oda
// Version     : 0.1
//============================================================================

#include "toyproblem.hpp"

bool ToyProblem::eval_f(vector< AD<double> >& x, AD<double>& y) {
//	y = -(x[1] - 2.0) * (x[1] - 2.0);
	y = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
	return true;
}
bool ToyProblem::eval_g(vector< AD<double> >& x, vector< AD<double> >& z) {
//	z[0] = -(x[0]*x[0] + x[1] - 1);
	z[0] = x[0] * x[1] * x[2] * x[3];
	z[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
	return true;
}
bool ToyProblem::get_initial_point(vector< double >& x) {
//	x[0] = 0.5;
//	x[1] = 1.5;
	x[0] = 1.0;
	x[1] = 5.0;
	x[2] = 5.0;
	x[3] = 1.0;
	return true;
}
bool ToyProblem::get_g_L(vector< double >& g_L) {
//	g_L[0] = 0;
	g_L[0] = 25;
	g_L[1] = 40;
	return true;
}
bool ToyProblem::get_g_U(vector< double >& g_U) {
//	g_U[0] = 0;
	g_U[0] = 1.0e19;
	g_U[1] = 40;
	return true;
}
bool ToyProblem::get_x_L(vector< double >& x_L) {
//	x_L[0] = -1.;
//	x_L[1] = -1.0e19;
	x_L[0] = 1.;
	x_L[1] = 1.;
	x_L[2] = 1.;
	x_L[3] = 1.;
	//x_L[1] = -1.0e19;
	return true;
}
bool ToyProblem::get_x_U(vector< double >& x_U) {
//	x_U[0] = 1.;
//	x_U[1] = 1.0e19;
	x_U[0] = 5.;
	x_U[1] = 5.;
	x_U[2] = 5.;
	x_U[3] = 5.;
//		x_U[1] = +1.0e19;
	return true;
}
size_t ToyProblem::get_n() {
	return 4;
}
size_t ToyProblem::get_m() {
	return 2;
}
