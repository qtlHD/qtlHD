/*
 * Rdistributions: wrappers for distribution functions (dnorm, qnorm, etc), from R
 *
 *     rdmd -L/usr/lib/libR.dylib --main -unittest math/distributions/Rdistributions.d 
 */

module math.distributions.Rdistributions;

import std.algorithm;
import std.math;
import std.range;
import std.random;
import std.stdio;
import std.string;
import std.conv;

import std.c.stdlib;

extern(C) {
  alias float f_float;
  alias double f_double;
  alias int f_int;
  alias char f_char;
}

version(Windows){
  private import std.loader;
  private import arch.windows;
  private import qtl.plugins.renv.libload;

  extern (C) double function(double, double, double, int) Rf_dnorm4;
  extern (C) double function(double, double, double, int, int) Rf_pnorm5;
  extern (C) double function(double, double, double, int, int) Rf_qnorm5;

  static this(){
    HXModule lib = load_library("R");
    load_function(Rf_dnorm4)(lib, "Rf_dnorm4");
    load_function(Rf_pnorm5)(lib, "Rf_pnorm5");
    load_function(Rf_qnorm5)(lib, "Rf_qnorm5");
    writeln("Loaded R functionality");
  }

}else{
  pragma(lib, "R");

  extern (C) double Rf_dnorm4(double, double, double, int);
  extern (C) double Rf_pnorm5(double, double, double, int, int);
  extern (C) double Rf_qnorm5(double, double, double, int, int);

  // also do  exp, chisq, f, gamma, pois
}


double dnorm(double x, double mu, double sigma, bool give_log)
{
  return Rf_dnorm4(x, mu, sigma, cast(int)give_log);
}

double pnorm(double x, double mu, double sigma, bool lower_tail, bool give_log)
{
  return Rf_pnorm5(x, mu, sigma, cast(int)lower_tail, cast(int)give_log);
}

double qnorm(double p, double mu, double sigma, bool lower_tail, bool log_p)
{
  return Rf_qnorm5(p, mu, sigma, cast(int)lower_tail, cast(int)log_p);
}

unittest {
  writeln("Unit test " ~ __FILE__);

  double[] x = [-0.6584, 6.3812, 6.1650, 5.1128, 3.9026];
  double[] R_dn = [0.022453690350634571832, 0.119608077529658976546, 0.123322523794599181457,
                   0.132886791910129337113, 0.124374806177956639952];
  double[] R_pn = [0.029638494723723464441, 0.677385432655885244557, 0.651115313701991826889,
                   0.514996696038364620840, 0.357257159939336377263];
  
  foreach(i, xv; x) {
    assert(abs(dnorm(xv, 5.0, 3.0, false) - R_dn[i]) < 1e-12);
    assert(abs(pnorm(xv, 5.0, 3.0, true, false) - R_pn[i]) < 1e-12);
    assert(abs(qnorm(R_pn[i], 5.0, 3.0, true, false) - xv) < 1e-12);

    assert(abs(dnorm(xv, 5.0, 3.0, true) - log(R_dn[i])) < 1e-12);
    assert(abs(pnorm(xv, 5.0, 3.0, true, true) - log(R_pn[i])) < 1e-12);
    assert(abs(qnorm(log(R_pn[i]), 5.0, 3.0, true, true) - xv) < 1e-12);

    double y1 = dnorm(xv, 5.0, 3.0, false);
    double y2 = pnorm(xv, 5.0, 3.0, true, false);
    double z = R_pn[i];
    double y3 = qnorm(z, 5.0, 3.0, true, false);

    writeln("This is okay:");
    writefln("%9.5f %9.5f %9.5f %9.5f", xv, y1, y2, y3);

    writeln("But this gives a seg fault:");
    writefln("%9.5f %9.5f %9.5f %9.5f", xv, y1, y2, qnorm(z, 5.0, 3.0, true, false));
  }
}