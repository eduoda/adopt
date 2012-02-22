AUTOMAKE_OPTIONS = subdir-objects
EXTRA_DIST = autogen.sh clear.sh plot.sce

adopt_libincludedir = $(libdir)/adopt/include
adopt_includedir = $(includedir)

nodist_adopt_libinclude_HEADERS = config.h
nobase_adopt_include_HEADERS = adopt/adopt.hpp adopt/optcontrol.hpp adopt/solvers/IpoptSolver.hpp

lib_LTLIBRARIES = libadopt.la libadopt-optcontrol.la libadopt-ipopt.la

libadopt_la_SOURCES = adopt/adopt.cpp adopt/adopt.hpp 
libadopt_la_LDFLAGS = -version-info 1:0:0

libadopt_optcontrol_la_SOURCES = adopt/optcontrol.cpp adopt/optcontrol.hpp
libadopt_optcontrol_la_LDFLAGS = -version-info 1:0:0

libadopt_ipopt_la_SOURCES = adopt/solvers/IpoptSolver.cpp adopt/solvers/IpoptSolver.hpp
libadopt_ipopt_la_LDFLAGS =  -version-info 1:0:0

AM_CXXFLAGS = -pipe
#	-Wall -pedantic -Wcast-qual -Wwrite-strings -Wconversion -Wno-unknown-pragmas \
# -fomit-frame-pointer -g

LDADD = -lgfortranbegin -lgfortran -lgcc_s -lipopt -lpthread -lm -ldl

AM_LDFLAGS = -L${IPOPT_PREFIX}/lib

AM_CPPFLAGS = -I${CPPAD_PREFIX}/include -I${IPOPT_PREFIX}/include/coin

# Tell make to ignore clean.sh
.PHONY : clean

#
EXTRA_PROGRAMS = solve
solve_SOURCES = models/solve.cpp \
	models/simple_car.cpp models/simple_car.hpp \
	models/AUV.cpp models/AUV.hpp   \
	models/toyproblem.cpp models/toyproblem.hpp \
	models/Diabetes.cpp models/Diabetes.hpp
solve_LDADD = libadopt.la libadopt-optcontrol.la libadopt-ipopt.la $(LDADD)
run: solve
	solve

