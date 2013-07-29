/**
 * Multiple linear regression for the multiple QTL modelling and scanning routine
 **/

module qtl.core.mqm.regression;

import std.math, std.stdio, std.conv;
import qtl.core.mqm.matrix, qtl.core.mqm.statistics, qtl.core.mqm.mqmutils;
import qtl.core.mqm.LUdecomposition;

/*
 * Structure holding our MQM model
 */
struct Model{
  double   logL = 0.0;
  double[] Fy;
  double[] indL;
  double[] params;
  Stats    stats;
}

/*
 * Give a Likelihood of odds score between two models
 */
double toLOD(Model model, Model nullmodel){ 
  return(abs((2.0 * model.logL) - (2.0 * nullmodel.logL)) / 4.60517); 
}

/*
 * Calculate the model loglikelihood
 */
Model calculateLogLikelihood(uint nsamples, in double[] residual, in double[] w, real variance){
  Model model = Model(0.0, newvector!double(nsamples, 0.0), newvector!double(nsamples, 0.0));

  for(size_t i = 0; i < nsamples; i++){
    model.Fy[i]    = LogNormal(residual[i], variance);
    model.indL[i] += w[i] * model.Fy[i];
    model.logL    += log(model.indL[i]);
  }
  return model;
}

/*
 * Weighted multivariate regression
 */
Model multiVariateRegression(in double[][] x, in double[] w, in double[] y, in int[] nullmodel = []){
  size_t nvariables = x[0].length;
  size_t nsamples = x.length;
  double[][] Xt = translate!double(x);
  double[] XtWY = calculateParameters(nvariables,nsamples,Xt, w, y);
  if(nullmodel.length != 0){ // If we have the a NULL model, drop the estimated params to 0.0
    for(size_t i = 1; i < nvariables; i++){
       // The nullmodel has always 1 parameter less Y = M + F1..Fn + Error
       // (The first parameter is the estimated mean)
      if(nullmodel[(i-1)] == 1) XtWY[i] = 0.0;
    }
  }
  Stats s = calculateStatistics(nvariables, nsamples, Xt, XtWY, y, w);
  Model f = calculateLogLikelihood(nsamples, s.residual, w, s.variance);
  f.params = copyvector(XtWY);
  f.stats = s;
  return f;
}

/*
 * Calculate model likelihood by EM
 */
Model likelihoodbyEM(in double[][] x, double[] w, in double[] y, bool verbose = false){
  uint nvariables = cast(uint)x[0].length;
  uint nsamples = cast(uint)x.length;
  uint maxemcycles = 1000;
  uint emcycle = 0;
  double delta = 1.0f;
  double logL = 0.0f;
  double logLprev = 0.0f;
  
  Model f;
  if(verbose) writefln("Starting EM:");
  while((emcycle < maxemcycles) && (delta > 1.0e-10)){
    f = multiVariateRegression(x, w, y);

    for(size_t s = 0; s < nsamples; s++){
      if(w[s] != 0) w[s] = (w[s] + f.Fy[s])/w[s];
    }
    delta = fabs(f.logL - logLprev);
    logLprev=f.logL;
    emcycle++;
  }
  
  if(verbose) writefln("[REGRESSION] EM took %d/%d cyclies", emcycle, maxemcycles);
  return f;
}

