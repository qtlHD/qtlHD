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
import std.c.stdio;
import std.stdio;
import std.string;
import std.conv;

import std.c.stdlib;

version(Windows){
  private import std.loader;
  private import arch.windows;
  private import qtl.plugins.renv.libload;

  extern (C) double function(double, double, double, int) Rf_dnorm4;
  extern (C) double function(double, double, double, int, int) Rf_pnorm5;
  extern (C) double function(double, double, double, int, int) Rf_qnorm5;

  extern (C) double function(double, double, int) Rf_dpois;
  extern (C) double function(double, double, int, int) Rf_ppois;
  extern (C) double function(double, double, int, int) Rf_qpois;

  static this(){
    HXModule lib = load_library("R");

    load_function(Rf_dnorm4)(lib, "Rf_dnorm4");
    load_function(Rf_pnorm5)(lib, "Rf_pnorm5");
    load_function(Rf_qnorm5)(lib, "Rf_qnorm5");

    load_function(Rf_dpois)(lib, "Rf_dpois");
    load_function(Rf_ppois)(lib, "Rf_ppois");
    load_function(Rf_qpois)(lib, "Rf_qpois");

    writeln("Loaded R functionality");
  }

}else{
  extern (C) double Rf_dnorm4(double, double, double, int);
  extern (C) double Rf_pnorm5(double, double, double, int, int);
  extern (C) double Rf_qnorm5(double, double, double, int, int);

  extern (C) double Rf_dpois(double, double, int);
  extern (C) double Rf_ppois(double, double, int, int);
  extern (C) double Rf_qpois(double, double, int, int);
}


double dnorm(double x, double mu, double sigma, int give_log)
{
  return Rf_dnorm4(x, mu, sigma, give_log);
}

double pnorm(double x, double mu, double sigma, int lower_tail, int give_log)
{
  return Rf_pnorm5(x, mu, sigma, lower_tail, give_log);
}

double qnorm(double p, double mu, double sigma, int lower_tail, int log_p)
{
  return Rf_qnorm5(p, mu, sigma, lower_tail, log_p);
}




double dpois(double x, double lambda, int give_log)
{
  return Rf_dpois(x, lambda, give_log);
}

double ppois(double x, double lambda, int lower_tail, int give_log)
{
  return Rf_ppois(x, lambda, lower_tail, give_log);
}

double qpois(double p, double lambda, int lower_tail, int log_p)
{
  return Rf_qpois(p, lambda, lower_tail, log_p);
}



unittest {
  writeln("Unit test " ~ __FILE__);

  double y=0.0;
  foreach(i; 0..10) {
    writefln("%9.5f %9.5f", y, dnorm(y, 5.0, 2.0, 0));
    y += 1.0;
  }

}