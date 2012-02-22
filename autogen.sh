#!/bin/sh

aclocal
libtoolize --force --copy
autoheader
touch stamp-h
autoconf
automake --add-missing --copy
rm -r autom4te.cache