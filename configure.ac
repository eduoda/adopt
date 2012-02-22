AC_INIT(adopt, 0.1, oda@ime.usp.br)
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_CONFIG_HEADER(config.h)
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL
LT_INIT
AC_ARG_WITH(
	ipopt,
	[  --with-ipopt=<path>     Ipopt installation path: <path>/include
                                                   <path>/lib ],
	[IPOPT_PREFIX=$with_ipopt],
	AC_MSG_ERROR([You must call configure with the --with-ipopt option.])
)
AC_SUBST(IPOPT_PREFIX)
AC_ARG_WITH(
	cppad,
	[  --with-cppad=<path>     CppAD installation path: <path>/include],
	[CPPAD_PREFIX=$with_cppad],
	AC_MSG_ERROR([You must call configure with the --with-cppad option.])
)
AC_SUBST(CPPAD_PREFIX)

AC_OUTPUT(Makefile)
