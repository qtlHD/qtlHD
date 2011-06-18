/**
 * Module for Haley-Knott scanone (single QTL mapping by HK regression).
 *
 * Run unit tests with:
 *
 *   rm read_csv ; dmd -unittest qtl/plugins/input/read_csv.d qtl/core/*.d ; time ./read_csv
 * 
 */

module qtl.core.scaneone_hk;

import std.container;
import qtl.core.primitives;

/*
 * Performs genome scan using the Haley-Knott regression method
 * (regressing phenotypes on conditional genotype probabilities; the
 * multipoint genotype probabilities have already been calculated in
 * calc.genoprob)
 * 
 * The old R/qtl interface read:
 *
 *   n_ind        Number of individuals
 *   n_pos        Number of marker positions
 *   n_gen        Number of different genotypes
 *   Genoprob     Array of conditional genotype probabilities
 *   Addcov       Matrix of additive covariates: Addcov[cov][ind]
 *   n_addcov     Number of columns of Addcov
 *   Intcov       Number of interactive covariates: Intcov[cov][ind]
 *   n_intcov     Number of columns of Intcov
 *   pheno        Phenotype data, as a vector
 *   nphe         Number of phenotypes
 *   weights      Vector of positive weights, of length n_ind
 *   ind_noqtl    Indicators (0/1) of which individuals should be excluded 
 *                from QTL effects.  
 *   Result       Result matrix of size [n_pos x (nphe)] containing the
 *                LOD scores for each phenotype
 *
 * The new (tentative) interface is:
 *
 * In:
 *
 *   markers          List of markers and positions, including refs to
 *   -> genotype[][]     Matrix of genotypes (one or more)
 *   phenotypes       List of phenotypes
 *   individuals      List of individuals with refs to
 *   -> weights
 *   -> exclude
 *   covariates       List of covariates 
 *
 * Out:
 *  
 *   mappedqtls       The LOD scores for each phenotype
 *
 **********************************************************************/

import std.stdio;
import std.conv;
import std.string;
import std.path;

import qtl.core.primitives;
// import qtl.core.chromosome;
// import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.plugins.input.read_csv;

class M {
}

class Mref(T) {
  Genotype!T[][] genotypes;
}

alias Mref!F2 F2Mref;
alias SList!F2Mref MarkerRefs;

MappedQTLs scanone_hk(in MarkerRefs markers) /*
                      in PhenotypeContainer phenotypes, 
                      in Individual[] individuals,
                      in CovariateContainer covariates) */
{
  return null;
}

unittest {
  writeln("Unit test " ~ __FILE__);
  Genotype!F2[][] genotypes;
  auto mref = new Mref!F2;
  mref.genotypes = genotypes;

  alias F2 XType;
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","listeria.csv");
  writeln("  - reading CSV " ~ fn);
  auto data = new ReadSimpleCSV!XType(fn);
  assert(data.markers.length == 133, to!string(data.markers.length));
  auto covariates = new Covariates!double;
  MarkerRefs markerrefs;
  auto result = scanone_hk(markerrefs); // ,data.phenotypes,data.individuals,covariates);
}

