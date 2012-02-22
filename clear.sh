#!/bin/sh

make clean
rm solve
rm -vrf aclocal.m4 autom4te.cache							\
	libtool													\
	stamp-h*												\
	configure config.guess config.h.in config.sub ltmain.sh	\
	Makefile.in install-sh missing depcomp					\
	Makefile config.log config.status config.h*				\
	src/Makefile.in src/Makefile							\
	models/Makefile.in models/Makefile

