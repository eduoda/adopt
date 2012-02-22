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
// Name        : Diabetes.hpp
// Author      : Eduardo Oda
// Version     : 0.1
// Description : Optimal control of insulin.
//============================================================================

#ifndef DIABETES_HPP_
#define DIABETES_HPP_

#include "adopt/optcontrol.hpp"

class Diabetes : public OptControlProblem{
public:
	void field(AD<double> t, vector< AD<double> >& x, vector< AD<double> >& u, vector< AD<double> >& y);
	AD<double> cost(vector< AD<double> >& x, vector< AD<double> >& u);
	void get_state_L(vector<double>& x_L);
	void get_state_U(vector<double>& x_U);
	void get_control_L(vector<double>& u_L);
	void get_control_U(vector<double>& u_U);
	void get_initial_state(vector<double>& x0);
	void get_final_state(vector<double>& xF);
	size_t get_N();
	size_t get_M();
	size_t get_C();
};

#endif /* DIABETES_HPP_ */
