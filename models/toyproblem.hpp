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
// Name        : toyproblem.hpp
// Author      : Eduardo Oda
// Version     : 0.1
// Description : This file implements the problem:
//     min   x1*x4*(x1 + x2 + x3)  +  x3
//     s.t.  x1*x2*x3*x4                   >=  25
//           x1**2 + x2**2 + x3**2 + x4**2  =  40
//           1 <=  x1,x2,x3,x4  <= 5
//
//     Starting point:
//        x = (1, 5, 5, 1)
//
//     Optimal solution:
//        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
//============================================================================

#ifndef TOYPROBLEM_HPP_
#define TOYPROBLEM_HPP_

#include "adopt/adopt.hpp"

class ToyProblem : public Problem{
private:
	bool eval_f(vector< AD<double> >& x, AD<double>& y);
	bool eval_g(vector< AD<double> >& x, vector< AD<double> >& z);
public:
	bool get_initial_point(vector< double >& x);
	bool get_g_L(vector< double >& g_L);
	bool get_g_U(vector< double >& g_U);
	bool get_x_L(vector< double >& x_L);
	bool get_x_U(vector< double >& x_U);
	size_t get_n();
	size_t get_m();
};

#endif /* TOYPROBLEM_HPP_ */
