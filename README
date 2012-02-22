Adopt - Environment for nonlinear optimization
Copyright (C) 2010 Eduardo Oda

Adopt is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Adopt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Adopt.  If not, see <http://www.gnu.org/licenses/>.


INTRODUCTION
============
Given a nonlinear problem:
	min	 f(x)
	s.t. g_L <= g(x) <= g_U
		 x_L <=   x  <= x_U
	x \in R^n
	f: R^n --> R
	g: R^n --> R^m
it can be solved by implementing two classes: Problem and Solver.

Your implementation of Problem can optionally redefine the
get_initial_multipliers method to set initial values to
Lagrangian multipliers.

You can also use an implementation of the class OptControlProblem to solve an
optimum control problem:

	        _T
	       /
	min   / c(x,u) dt
	    _/
	     0

	s.t.	x'=f(x,u)
		x_L < x < x_U
		u_L < u < u_U
		x(0)=x0
		x(T)=xF

	x \in R^M      u \in R^C

We use a approximation of the cost function and a numeric integration of the
dynamical restriction to obtain the nonlinear problem:

	  x0   x1            ...             xN=xF
	    u0   u1          ...          u6
	  |----|----|----|----|----|----|----|------------>t
	  0    1    2    3    4    5    6    7=T     (N=7)
	  /----/
	    h = T/N

	      __N-1
	min  \    c(x_k,u_k)*h
	     /__
	      k=0
	s.t.
		x_0=x0
		x_{k+1} = x_k + h*PHI(x_k,u_k,t_k), k=0..N-1
		x_N=xF
		x_L < x_k < x_U
		u_L < u_k < u_U

	x=(h, x_0, x_1, ..., x_N, x_N, u_0, ..., u_{N-1})

Where PHI denotes the Classic Runge-Kutta method.

Derivatives are evaluated through Automatic Derivation (CppAD).

There are implementations of Solver to the following solvers:
	Ipopt

Prereqs:
	CppAD
	Boost uBlas
	Ipopt

	
INSTALLATION - GENERAL INSTRUCTIONS
===================================
OBS: See above an example of an installation.

After installing the pre-required software, checkout the latest version from
repository:

	$ svn co http://www.labmap.ime.usp.br/svn/adopt/trunk/ adopt
	$ cd adopt

If you are a Adopt developer you may need to run: 
	$ ./clear.sh
	$ ./autogen.sh

Configure and compile the Adopt:
	$ ./configure                           	\
		--with-ipopt=[IPOPT INSTALL DIR]	\
		--with-cppad=[CPPAD INSTALL DIR]	\
		--prefix=/home/oda64/opt/adopt		\
		CFLAGS="-DNDEBUG -O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math"
	$ make

Run the included examples:
	$ make run

Install adopt;
	$ make install
	

INSTALLATION EXAMPLE
====================
Here we give detailed instructions to have a fully working user installation of
Adopt with Ipopt. We assume that you will install everything in ~/adopt

1. Creating diretories:
	$ cd
	$ mkdir -p			\
		adopt/adopt/build	\
		adopt/ipopt/build	\
		adopt/cppad/build

2. Installing Ipopt:
	$ cd ~/adopt/ipopt
	$ svn co https://projects.coin-or.org/svn/Ipopt/stable/3.8 src
	$ cd src/ThirdParty
	$ cd Blas
	$ ./get.Blas
	$ cd ../Lapack
	$ ./get.Lapack
	$ cd ../ASL
	$ ./get.ASL
	$ cd ../Mumps
	$ ./get.Mumps
	$ cd ../Metis
	$ ./get.Metis

Ipopt requires some HSL Subroutines. They can be downloaded from:
http://hsl.rl.ac.uk/archive/hslarchive.html

Save MA27 as ma27ad.f and MC19 as mc19ad.f, both in ~/adopt/ipopt/src/ThirParty/HSL/ 

Compile and install:
	$ cd ~/adopt/ipopt/build
	$ ../src/configure
	$ make
	$ make test
	$ make install

3. Installing CppAD:
CppAD requires Boost. Usually modern linux distribuitions come with Boost out
of the box. Searsh for boost directory in your includes directories (Slackware
installs Boost in /usr/include).
	$ cd ~/adopt/cppad
	$ svn co https://projects.coin-or.org/svn/CppAD/trunk src
	$ ./configure				\
		--prefix=~/adopt/cppad/build	\
		BOOST_DIR=/usr/include/		\
		IPOPT_DIR=~/adopt/ipopt/build/	\
		CXX_FLAGS="-O3 -DNDEBUG -fPIC"
	$ make
	$ make install

4. Installing Adopt:
	$ svn co http://www.labmap.ime.usp.br/svn/adopt/trunk/ src
	$ cd src
	$ ./configure					\
		--with-ipopt=~/adopt/ipopt/build/	\
		--with-cppad=~/adopt/cppad/build/	\
		--prefix=~/adopt/adopt/build		\
		CXXFLAGS="-DNDEBUG -O3 -fPIC"
	$ make
	$ make run
	$ make install


PROBLEM IMPLEMENTATION
======================
The best way to implement your problem is to take one example as base.

Copy solve.cpp, toyproblem.cpp and toyproblem.hpp to some directory. Edit 
the file toyproblem.cpp do fit your problem.

All functions, except get_n and get_m, return true if the evaluation was ok and
false otherwise. At this time this return value is ignored, so it's safe to
always return true.

	bool eval_f(vector< AD<double> >& x, AD<double>& y);
	Objective function  y=f(x)
	
	bool eval_g(vector< AD<double> >& x, vector< AD<double> >& z)
	Restriction function z=g(x)
	
	bool get_initial_point(vector< double >& x);
	Initial guess
	
	bool get_g_L(vector< double >& g_L);
	Restriction lower bound g_L <= g(x)
	
	bool get_g_U(vector< double >& g_U);
	Restriction upper bound g(x) <= g_U
	
	bool get_x_L(vector< double >& x_L);
	Variable lower bound x_L <= x
	
	bool get_x_U(vector< double >& x_U);
	Variable upper bound x <= x_U
	
	size_t get_n();
	Dimension of x
	
	size_t get_m();
	Dimension of z=g(x) 


OPTIMAL CONTROL PROBLEM IMPLEMENTATION
======================================
You can also take an example as base.

Copy solve.cpp, simplecar.cpp and simplecar.hpp to some directory. Edit 
the file simplecar.cpp do fit your problem.

        void field(AD<double> t, vector< AD<double> >& x, vector< AD<double> >& u, vector< AD<double> >& y);
        The field x'=field(t,x,u)=y
        
        AD<double> cost(vector< AD<double> >& x, vector< AD<double> >& u);
        Cost function
        
        void get_state_L(vector<double>& x_L);
        State lower bound
        
        void get_state_U(vector<double>& x_U);
        State upper bound
        
        void get_control_L(vector<double>& u_L);
        Control lower bound
        
        void get_control_U(vector<double>& u_U);
        Control upper bound
        
        void get_initial_state(vector<double>& x0);
        State at t=0
        
        void get_final_state(vector<double>& xF);
        State to achive
        
        size_t get_N();
        Discretizations steps
        
        size_t get_M();
        Dimension of states
        
        size_t get_C();
	Dimension of control 

You also need to edit solve.cpp commenting the lines:
	#include "toyproblem.hpp"
	Problem *a = new ToyProblem();

And uncommenting the lines:
	#include "simple_car.hpp"
	#include "adopt/optcontrol.hpp"
	Problem *a = new Car();

COMPILING YOUR PROBLEM
======================
To compile yout problem implementation run the following commands:
	$ g++ -c solve.cpp toyproblem.cpp		\
		-I~/adopt/adopt/build/include/		\
		-I~/adopt/cppad/build/include/		\
		-I~/adopt/ipopt/build/include/coin/	\
		-O3 -DNDEBUG -fPIC
	$ g++ -o solve solve.o toyproblem.o			\
		-L~/adopt/ipopt/build/lib/			\
		-L~/adopt/adopt/build/lib/			\
		-lipopt						\
		-ladopt -ladopt-optcontrol -ladopt-ipopt	\
		-lgfortran -lpthread				\
		-O3 -DNDEBUG -fPIC
		
Set the LD_LIBRARY_PATH and run the program:
	$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/adopt/ipopt/build/lib:~/adopt/adopt/build/lib/
	$ ./solve
	
