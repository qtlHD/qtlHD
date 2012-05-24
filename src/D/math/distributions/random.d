/*
 * random: functions to generate random draws from various distributions
 */

module math.distributions.random;

import math.distributions.distributions;
import std.stdio;
import std.math;
import std.mathspecial;
import std.algorithm;
import std.random;
import std.string;
import std.conv;


// For each of these, I'm just using the "inverse method"
//   (draw from uniform(0,1) and then find corresponding quantile)
double rnorm(in double mu, in double sigma, ref Random gen)
{
  double r = uniform(0.0, 1.0, gen);
  return( qnorm(r, mu, sigma, true) );
}

double rgamma(in double shape, in double scale, ref Random gen)
{
  double r = uniform(0.0, 1.0, gen);
  return( qgamma(r, shape, scale, true) );
}

uint rpois(in double lambda, ref Random gen)
{
  double r = uniform(0.0, 1.0, gen);
  return( qpois(r, lambda, true) );
}

unittest {
  writeln("Unit test " ~ __FILE__);

  Random gen;
  gen.seed(unpredictableSeed);

  writeln("10 draws from normal(5, 3):");
  foreach(i; 0..10) {
    writef("%9.5f ", rnorm(5.0, 3.0, gen));
  }
  writeln;

  writeln("10 draws from gamma(3,3):");
  foreach(i; 0..10) {
    writef("%9.5f ", rgamma(3.0, 3.0, gen));
  }
  writeln;

  writeln("10 draws from Poisson(4):");
  foreach(i; 0..10) {
    writef("%9d ", rpois(4.0, gen));
  }
  writeln;
}
