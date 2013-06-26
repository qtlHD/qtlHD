/**
 * Statistics for the multiple QTL modelling and scanning routine
 **/

module qtl.core.mqm.statistics;

import std.math, std.stdio, std.conv;
import qtl.core.mqm.matrix;

/*
 * Structure holding soem basic statistics data
 */
struct Stats{
  real     variance = 0.0;
  double[] fit;
  double[] residual;
}

/**
 * Function to calculate the stats structure ( variance, fit and residual )
 **/
Stats calculateStatistics(size_t nv, size_t ns, in double[][] xt, in double[] xtwy, in double[] y, in double[] w){
  Stats s = Stats(0.0, newvector!double(ns, 0.0), newvector!double(ns, 0.0));

  for(size_t i=0; i < ns; i++){
    s.fit[i] = 0.0;
    s.residual[i] = 0.0;
    for(size_t j=0; j < nv; j++){
      s.fit[i] += xt[j][i] * xtwy[j];
    }
    s.residual[i] = y[i] - s.fit[i];
    s.variance += w[i] * pow(s.residual[i], 2.0);
  }
  s.variance /= ns;
  return s;
}

