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

  extern (C) double function(double, double, int) Rf_dchisq;
  extern (C) double function(double, double, int, int) Rf_pchisq;
  extern (C) double function(double, double, int, int) Rf_qchiq;

  extern (C) double function(double, double, int) Rf_dexp;
  extern (C) double function(double, double, int, int) Rf_pexp;
  extern (C) double function(double, double, int, int) Rf_qexp;

  extern (C) double function(double, double, double, int) Rf_df;
  extern (C) double function(double, double, double, int, int) Rf_pf;
  extern (C) double function(double, double, double, int, int) Rf_qf;

  extern (C) double function(double, double, double, int) Rf_dgamma;
  extern (C) double function(double, double, double, int, int) Rf_pgamma;
  extern (C) double function(double, double, double, int, int) Rf_qgamma;

  extern (C) double function(double, double, int) Rf_dpois;
  extern (C) double function(double, double, int, int) Rf_ppois;
  extern (C) double function(double, double, int, int) Rf_qpois;

  static this(){
    HXModule lib = load_library("R");

    load_function(Rf_dnorm4)(lib, "Rf_dnorm4");
    load_function(Rf_pnorm5)(lib, "Rf_pnorm5");
    load_function(Rf_qnorm5)(lib, "Rf_qnorm5");

    load_function(Rf_dnorm)(lib, "Rf_dnorm");
    load_function(Rf_pnorm)(lib, "Rf_pnorm");
    load_function(Rf_qnorm)(lib, "Rf_qnorm");

    load_function(Rf_dchisq)(lib, "Rf_dchisq");
    load_function(Rf_pchisq)(lib, "Rf_pchisq");
    load_function(Rf_qchisq)(lib, "Rf_qchisq");

    load_function(Rf_dexp)(lib, "Rf_dexp");
    load_function(Rf_pexp)(lib, "Rf_pexp");
    load_function(Rf_qexp)(lib, "Rf_qexp");

    load_function(Rf_df)(lib, "Rf_df");
    load_function(Rf_pf)(lib, "Rf_pf");
    load_function(Rf_qf)(lib, "Rf_qf");

    load_function(Rf_dgamma)(lib, "Rf_dgamma");
    load_function(Rf_pgamma)(lib, "Rf_pgamma");
    load_function(Rf_qgamma)(lib, "Rf_qgamma");

    load_function(Rf_dpois)(lib, "Rf_dpois");
    load_function(Rf_ppois)(lib, "Rf_ppois");
    load_function(Rf_qpois)(lib, "Rf_qpois");

    writeln("Loaded R functionality");
  }

}else{
  extern (C) double Rf_dnorm4(double, double, double, int);
  extern (C) double Rf_pnorm5(double, double, double, int, int);
  extern (C) double Rf_qnorm5(double, double, double, int, int);

  extern (C) double Rf_dchisq(double, double, int);
  extern (C) double Rf_pchisq(double, double, int, int);
  extern (C) double Rf_qchisq(double, double, int, int);

  extern (C) double Rf_dexp(double, double, int);
  extern (C) double Rf_pexp(double, double, int, int);
  extern (C) double Rf_qexp(double, double, int, int);

  extern (C) double Rf_df(double, double, double, int);
  extern (C) double Rf_pf(double, double, double, int, int);
  extern (C) double Rf_qf(double, double, double, int, int);

  extern (C) double Rf_dgamma(double, double, double, int);
  extern (C) double Rf_pgamma(double, double, double, int, int);
  extern (C) double Rf_qgamma(double, double, double, int, int);

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




double dchisq(double x, double df, int give_log)
{
  return Rf_dchisq(x, df, give_log);
}

double pchisq(double x, double df, int lower_tail, int give_log)
{
  return Rf_pchisq(x, df, lower_tail, give_log);
}

double qchisq(double p, double df, int lower_tail, int log_p)
{
  return Rf_qchisq(p, df, lower_tail, log_p);
}




double dexp(double x, double rate, int give_log)
{
  return Rf_dexp(x, 1.0/rate, give_log);
}

double pexp(double x, double rate, int lower_tail, int give_log)
{
  return Rf_pexp(x, 1.0/rate, lower_tail, give_log);
}

double qexp(double p, double rate, int lower_tail, int log_p)
{
  return Rf_qexp(p, 1.0/rate, lower_tail, log_p);
}




double df(double x, double df1, double df2, int give_log)
{
  return Rf_df(x, df1, df2, give_log);
}

double pf(double x, double df1, double df2, int lower_tail, int give_log)
{
  return Rf_pf(x, df1, df2, lower_tail, give_log);
}

double qf(double p, double df1, double df2, int lower_tail, int log_p)
{
  return Rf_qf(p, df1, df2, lower_tail, log_p);
}




double dgamma(double x, double shape, double rate, int give_log)
{
  return Rf_dgamma(x, shape, rate, give_log);
}

double pgamma(double x, double shape, double rate, int lower_tail, int give_log)
{
  return Rf_pgamma(x, shape, rate, lower_tail, give_log);
}

double qgamma(double p, double shape, double rate, int lower_tail, int log_p)
{
  return Rf_qgamma(p, shape, rate, lower_tail, log_p);
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

  double[] x = [-0.6584, 6.3812, 6.1650, 5.1128, 3.9026];
  double[] R_dn = [0.022453690350634571832, 0.119608077529658976546, 0.123322523794599181457,
                   0.132886791910129337113, 0.124374806177956639952];
  double[] R_pn = [0.029638494723723464441, 0.677385432655885244557, 0.651115313701991826889,
                   0.514996696038364620840, 0.357257159939336377263];

  foreach(i, xv; x) {
    assert(abs(dnorm(xv, 5.0, 3.0, 0) - R_dn[i]) < 1e-12);
    assert(abs(pnorm(xv, 5.0, 3.0, 1, 0) - R_pn[i]) < 1e-12);
    assert(abs(qnorm(R_pn[i], 5.0, 3.0, 1, 0) - xv) < 1e-12);

    assert(abs(dnorm(xv, 5.0, 3.0, 1) - log(R_dn[i])) < 1e-12);
    assert(abs(pnorm(xv, 5.0, 3.0, 1, 1) - log(R_pn[i])) < 1e-12);
    assert(abs(qnorm(log(R_pn[i]), 5.0, 3.0, 1, 1) - xv) < 1e-12);
  }

  // take absolute values
  foreach(ref xv; x) xv = abs(xv);

  double[] R_dcs = [0.051115639783442655408, 0.088202799976698234574, 0.093320206692885265820,
                       0.119274225280151260908, 0.145673665710276567520];
  double[] R_pcs = [0.014826460654626052657, 0.729127114691895505949, 0.709506763671853946107,
                       0.597730271220054287795, 0.436477204281711328449];

  foreach(i, xv; x) {
    assert(abs(dchisq(xv, 5.0, 0) - R_dcs[i]) < 1e-8);
    writeln(to!string(i) ~ " " ~ to!string(xv) ~ "     " ~
            to!string(pchisq(xv, 5.0, 1, 0)) ~ " " ~ to!string(R_pcs[i]));
    //    assert(abs(pchisq(xv, 5.0, 1, 0) - R_pcs[i]) < 1e-8);
    //    writeln(qchisq(R_pcs[i], 5.0, 1, 0), " ", xv);
    //    assert(abs(qchisq(R_pcs[i], 5.0, 1, 0) - xv) < 1e-8);

    //    assert(abs(dchisq(xv, 5.0, 1) - log(R_dcs[i])) < 1e-8);
    //    assert(abs(pchisq(xv, 5.0, 1, 1) - log(R_pcs[i])) < 1e-8);
    //    assert(abs(qchisq(log(R_pcs[i]), 5.0, 1, 1) - xv) < 1e-8);
  }

  double[] xi = [2.0, 3.0, 4.0, 5.0, 9.0];

  foreach(xiv; xi) {
    writeln(to!string(xiv) ~ " " ~ to!string(dpois(xiv, 5.0, 0)) ~ "     " ~
            to!string(xiv) ~ " " ~ to!string(ppois(xiv, 5.0, 1, 0)));
    //    writeln(qpois(0.5, 5.0, 1, 0));
  }

  foreach(i, xv; x) {
    assert(abs(pexp(xv, 0.2, 1, 0) - (1.0-exp(-xv * 0.2))) < 1e-12);
  }

}