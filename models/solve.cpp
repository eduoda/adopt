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
// Name        : solver.hpp
// Author      : Eduardo Oda
// Version     : 0.1
// Description :
// Prereqs     :
//============================================================================

//#include "AUV.hpp"
//#include "Diabetes.hpp"
#include "toyproblem.hpp"
//#include "simple_car.hpp"
//#include "adopt/optcontrol.hpp"
#include "adopt/solvers/IpoptSolver.hpp"

int main(void) {
	Problem *a = new ToyProblem();
//	Problem *a = new AUV();
//	Problem *a = new Car();
//	Problem *a = new Diabetes();
	Solver *s = new IpoptSolver(*a);
	s->solve();
	cout << "Time required for execution: " << endl
		<< "F:      " << a->time_F << " seconds" << endl
		<< "G:      " << a->time_G << " seconds" << endl
		<< "grad F: " << a->time_grad_F << " seconds" << endl
		<< "jac G:  " << a->time_jac_G << " seconds" << endl
		<< "hess L: " << a->time_hess_L << " seconds" << endl
		<< "TOTAL: " <<
			a->time_F+
			a->time_G+
			a->time_grad_F+
			a->time_jac_G+
			a->time_hess_L
			<< " seconds" << endl;

	return 0;
}
