/**
 * hmm_util
 */

module qtl.core.hmm.hmm_util;

import std.stdio;
import std.math;

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
double addlog(double a, double b)
{
  enum TOL = 200.0;

  if(b > a + TOL) return(b);
  else if(a > b + TOL) return(a);
  else return(a + log1p(exp(b-a)));
}

import std.conv;

unittest {
  writeln("Unit test " ~ __FILE__);
  writeln("    unit test addlog");
  double a=50, b=60, d=2;
  auto value = abs(addlog(a,b) - log(exp(a)+exp(b)));
  assert(abs(addlog(a,b) - log(exp(a)+exp(b))) < 1e-15, to!string(value));
  assert(abs(addlog(b,a) - log(exp(a)+exp(b))) < 1e-15);
  assert(abs(addlog(a,d) - log(exp(a)+exp(d))) < 1e-15);
  assert(abs(addlog(d,a) - log(exp(a)+exp(d))) < 1e-15);
  assert(abs(addlog(b,d) - log(exp(b)+exp(d))) < 1e-15);
  assert(abs(addlog(d,b) - log(exp(b)+exp(d))) < 1e-15);
  assert(addlog(a,a+300) == a+300);
  assert(addlog(a,a-300) == a);
  assert(addlog(a+300,a) == a+300);
  assert(addlog(a-300,a) == a);
}
