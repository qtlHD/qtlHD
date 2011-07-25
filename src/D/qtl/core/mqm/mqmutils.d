/**
 * Utilities for the multiple QTL modelling and scanning routine
 **/
 
module qtl.core.mqm.mqmutils;

import std.container;
import qtl.core.primitives;
import std.stdio;
import std.algorithm;
import std.math;
import std.string;
import std.conv;
import qtl.core.genotype;
import qtl.plugins.input.read_csv;

version (Windows) {
  import qtl.core.libs.libload;
  import std.loader;
}else{
  pragma(lib, "r");
}

const HAVE_F77_UNDERSCORE = 1;
const IEEE_754 = 1;
const SUPPORT_UTF8 = 1;
const SUPPORT_MBCS = 1;
const ENABLE_NLS = 1;

version (Windows) {
  extern(C){
    double function(double, double, double, int) dnorm;
    double function(double, double, double, int, int) qf;
  }
  
  static this(){
    HXModule lib = load_library("R");
    load_function(dnorm)(lib,"Rf_dnorm4");
    load_function(qf)(lib,"Rf_qf");
    writeln("Loaded R functionality");
  }
  
}else{
  extern(C){
    double dnorm(double, double, double, int);
    double qf(double, double, double, int, int);
  }
}

double LogNormal(double residual, double variance){
  return dnorm(residual,0,sqrt(variance),0);
}

double InverseF(int df1, int df2, double alfa){
  return qf(1-alfa,df1,df2,0,0);
}

unittest{
  writeln("Unit test " ~ __FILE__);
  assert(to!string(InverseF(100,98,0.05))[0..8]    == "0.717548");
  assert(to!string(InverseF(200,50,0.05))[0..8]    == "0.706904");
  assert(to!string(InverseF(1000,900,0.05))[0..8]  == "0.898742");
  assert(to!string(InverseF(60,2,0.05))[0..8]      == "0.317419");
  assert(to!string(dnorm(3.0,5.0,2.0,1))[0..8]     == "-2.11209");
  assert(to!string(LogNormal(3.1,5.0))[0..8]       == "0.068244");
  assert(to!string(LogNormal(13.0,15.7))[0..8]     == "0.000462");
  assert(to!string(LogNormal(23.0,2.0))[0..8]      == "1.03502e");
  assert(to!string(LogNormal(3.3,45.0))[0..8]      == "0.052693");
  assert(to!string(LogNormal(3.0,5.6))[0..8]       == "0.075479");
}
