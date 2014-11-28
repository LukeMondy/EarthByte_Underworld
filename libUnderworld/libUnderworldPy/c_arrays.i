/* -*- C -*-  (not really, but good for syntax highlighting) */

%module c_arrays

/* this guy is an interface to malloc */
%include "carrays.i"
%include "cdata.i"
%array_class(double, DoubleArray)
%array_class(float, FloatArray)
%array_class(int, IntArray)
%array_class(unsigned, UnsignedArray)
