/*
 Code to simulate meiosis and basic crosses
*/

module qtl.core.simulate.meiosis;

import std.algorithm, std.stdio, std.string, std.conv;
import qtl.core.simulate.rlib;

// simulate crossover locations on a chromosome under no interference
// chrlen = chromosome length in cM
double[] meiosisNI(double chrlen)
in {
  assert(chrlen > 0, "chrlen must be > 0");
}
body {
  int n_xo = rpois(chrlen/100.0);

  auto result = new double[](n_xo);
  foreach(i; 0..n_xo)
    result[i] = unif_rand()*chrlen;
    
  sort(result);
  return(result);
}



// m = interference parameter
// p = proportion of crossovers coming from no interference model


double[] meiosis(double chrlen, int m, double p)
in {
  assert(chrlen > 0, "chrlen must be > 0");
  assert(m >= 0, "m must be >= 0");
  assert(p >= 0 && p <= 1, "p must be in [0,1]");
}
body {
  // no interference case is simple
  if(m==0 || p==0) return(meiosisNI(chrlen));

  // number of crossovers and intermediates (on 4-strand bundle)
  int n_points_interf = rpois(chrlen*(m+1)/50.0*(1.0-p) );
  // number of crossovers from no interference mechanism (on final product)
  int n_points_ni = rpois(chrlen/100.0*p);
  
  double[] result;
  result.reserve(n_points_interf/(m+1) + n_points_ni);

  auto pts_interf = new double[](n_points_interf);
  foreach(i; 0..n_points_interf)
    pts_interf[i] = unif_rand()*chrlen;
  sort(pts_interf);

  auto pts_ni = new double[](n_points_ni);
  foreach(i; 0..n_points_ni)
    result ~= unif_rand()*chrlen;

  return(result);
}