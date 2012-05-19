/*
 Code to simulate meiosis and basic crosses
*/

module qtl.core.simulate.meiosis;

import math.distributions.random;
import std.algorithm;
import std.stdio;
import std.string;
import std.conv;
import std.random;
import std.math;

// simulate crossover locations on a chromosome under no interference
// chrlen = chromosome length in cM
double[] meiosisNI(double chrlen, ref Random gen)
in {
  assert(chrlen > 0, "chrlen must be > 0");
}
body {
  int n_xo = rpois(chrlen/100.0, gen);

  auto result = new double[](n_xo);
  foreach(i; 0..n_xo)
    result[i] = uniform(0.0, chrlen, gen);

  sort(result);
  return(result);
}



// m = interference parameter
// p = proportion of crossovers coming from no interference model
double[] meiosis(double chrlen, int m, double p, ref Random gen)
in {
  assert(chrlen > 0, "chrlen must be > 0");
  assert(m >= 0, "m must be >= 0");
  assert(p >= 0 && p <= 1, "p must be in [0,1]");
}
body {
  // no interference case is simple
  if(m==0 || p==0) return(meiosisNI(chrlen, gen));

  // number of crossovers and intermediates (on 4-strand bundle)
  int n_points_interf = rpois(chrlen*(m+1)/50.0*(1.0-p), gen );
  // number of crossovers from no interference mechanism (on final product)
  int n_points_ni = rpois(chrlen/100.0*p, gen);

  double[] result;
  result.reserve(n_points_interf/(m+1) + n_points_ni);

  // locations of crossovers and intermediates
  auto pts_interf = new double[](n_points_interf);
  foreach(i; 0..n_points_interf)
    pts_interf[i] = uniform(0.0, chrlen, gen);
  sort(pts_interf);

  // every (m+1)st point is chiasma; thin to give crossovers
  int first = cast(int)uniform(0.0, m+1, gen);
  for(auto i=first; i<n_points_interf; i+= m+1)
    if(uniform(0.0, 1.0, gen) < 0.5) result ~= pts_interf[i];

  auto pts_ni = new double[](n_points_ni);
  foreach(i; 0..n_points_ni)
    result ~= uniform(0.0, chrlen, gen);

  sort(result);
  return(result);
}

unittest {
  writeln("Unit test " ~ __FILE__);

  Random gen;
  gen.seed(unpredictableSeed);

  double[] product;

  foreach(m; 0..10) {
    product = meiosis(200.0, m, 0.05, gen);
    writeln("No. crossovers = ", product.length);
    foreach(p; product)
      writef("%6.3f ", p);
    writeln();
  }
}