/**
 * Wrapper for R.dll, Converted to D from R.h by htod
 **/
 
module r.r;

//D 2.0 std imports
private import std.loader;
private import std.stdio;
private import std.conv;

import qtl.core.libs.libload;

import std.c.stdlib;
import std.c.stdio;
import std.c.math;

extern (C):
  //These are the functions we map to D from Rmath.h
  double function(double, double, double, int) dnorm;
  double function(double, double, double, int, int) qf;

//Load the functions when the module is loaded
static this(){
  HXModule lib = load_library("R");
  load_function(dnorm)(lib,"Rf_dnorm4");
  load_function(qf)(lib,"Rf_qf");
  writefln("mapped R.dll");
}

unittest {
  writefln("Unit test %s", __FILE__);
  //No chance of finding a 100 in a Normal distribution (m=0,std=1)
  assert(dnorm(100.0,0.0,1.0,0)==0.0);
}
