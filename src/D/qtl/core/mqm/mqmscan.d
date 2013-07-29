/**
 * Multiple QTL modelling and scanning routine
 **/

module qtl.core.mqm.mqmscan;
 
import std.container;
import qtl.core.primitives;
import std.stdio;
import std.algorithm;
import std.math;
import std.string;
import qtl.plugins.csv.read_csv;
import qtl.core.mqm.regression;
import qtl.core.mqm.mqmutils; //Load the R-bindings

Model backwardselect(){

}

// scanmqm, many positions - returns matrix of residual sums of squares [position][phenotype]
double[][] scanmqm_rj(in Probability[][][] genoprobs, in Phenotype[][] pheno, in double[][] addcovar, 
                   in double[][] intcovar, in double[] weights, in double tol=1e-8){
  if(genoprobs.length == 0)
    throw new Exception("genoprobs is empty");
  if(genoprobs[0].length != pheno.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and pheno");
  if(weights.length > 0 && genoprobs[0].length != weights.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and weights");
  if(addcovar.length > 0 && addcovar[0].length > 0 && genoprobs[0].length != addcovar.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and addcovar");
  if(intcovar.length > 0 && intcovar[0].length > 0 && genoprobs[0].length != intcovar.length)
    throw new Exception("Mismatch in no. individuals in genoprobs and intcovar");

  double[][] rss = new double[][](genoprobs.length, pheno[0].length);
  return rss;
}

unittest{

}

