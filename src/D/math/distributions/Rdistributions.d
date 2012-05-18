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

  extern (C) double function(double, double, double, int, int) Rf_qf;
  extern (C) double function(double, double, double, int, int) Rf_pgamma;

  static this(){
    HXModule lib = load_library("R");

    load_function(Rf_dnorm4)(lib, "Rf_dnorm4");
    load_function(Rf_pnorm5)(lib, "Rf_pnorm5");
    load_function(Rf_qnorm5)(lib, "Rf_qnorm5");

    load_function(Rf_dpois)(lib, "Rf_dpois");
    load_function(Rf_ppois)(lib, "Rf_ppois");
    load_function(Rf_qpois)(lib, "Rf_qpois");

    load_function(Rf_qf)(lib, "Rf_qf");
    load_function(Rf_pgamma)(lib, "Rf_pgamma");

    writeln("Loaded R functionality");
  }

}else{
  extern (C) double Rf_dnorm4(double, double, double, int);
  extern (C) double Rf_pnorm5(double, double, double, int, int);
  extern (C) double Rf_qnorm5(double, double, double, int, int);

  extern (C) double Rf_dpois(double, double, int);
  extern (C) double Rf_ppois(double, double, int, int);
  extern (C) double Rf_qpois(double, double, int, int);

  extern (C) double Rf_qf(double, double, double, int, int);
  extern (C) double Rf_pgamma(double, double, double, int, int);
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

double qf(double p, double df1, double df2, int lower_tail, int log_p)
{
  return Rf_qf(p, df1, df2, lower_tail, log_p);
}

double pgamma(double x, double alph, double scale, int lower_tail, int give_log)
{
  return Rf_pgamma(x, alph, scale, lower_tail, give_log);
}



unittest {
  writeln("Unit test " ~ __FILE__);

  double y=0.0;
  foreach(i; 0..10) {
    double z1 = dpois(y, 5.0, 0);
    double z2 = ppois(y, 5.0, 1, 0);
    double logz2 = ppois(y, 5.0, 1, 1);
    //    double z3 = qpois(y/10.0, 5.0, 1, 0);
    writefln("%9.5f %9.5f %9.5f %9.5f", y, z1, z2, logz2);
    y += 1.0;
  }

  writeln();
  for(double p=0.1; p < 0.95; p += 0.1) {
    y = qf(p, 10.0, 1.0, 1, 0);
    writefln("%9.5f %9.5f", p, y);
  }

  writeln();
  for(double p=0.0; p < 9.5; p += 1.0) {
    y = pgamma(5.0, p+1.0, 1.0, 0, 0);
    writefln("%9.5f %9.5f", p, y);
  }

}