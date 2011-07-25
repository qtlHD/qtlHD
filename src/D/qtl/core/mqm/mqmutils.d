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


unittest{
  writeln("Unit test " ~ __FILE__);
  assert(to!string(dnorm(3.0,5.0,2.0,1))[0..8] == "-2.11209");
}
