/**
 * Data manipulation
 */

module qtl.core.util.data_manip;

import std.conv, std.stdio, std.string, std.exception, std.path, std.file;
import qtl.core.primitives;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.plugins.qtab.read_qtab;

// return boolean vector indicating whether any
bool[] is_any_phenotype_missing(T)(Phenotype!T[][] pheno)
{
  auto ret = new bool[](pheno.length);
  foreach(i; 0..pheno.length) {
    ret[i] = false;
    foreach(j; 0..pheno[i].length) {
      if(isNA(pheno[i][j])) ret[i] = true;
      break;
    }
  }
  return ret;
}

// omit individuals from phenotype data
Phenotype!T[][] omit_ind_from_phenotypes(T)(Phenotype!T[][] pheno, bool[] to_omit)
{
  if(pheno.length != to_omit.length)
    throw new Exception("no. individuals in pheno (" ~ to!string(pheno.length) ~ 
                        ") doesn't match length of to_omit (" ~ to!string(to_omit.length) ~ ")");

  Phenotype!T[][] ret;

  foreach(i; 0..to_omit.length) {
    if(!to_omit[i])
      ret ~= pheno[i];
  }

  return ret;
}

// omit individuals from genotype matrix
GenotypeCombinator[][] omit_ind_from_genotypes(GenotypeCombinator[][] geno, bool[] to_omit)
{
  if(geno.length != to_omit.length)
    throw new Exception("no. individuals in geno (" ~ to!string(geno.length) ~ 
                        ") doesn't match length of to_omit (" ~ to!string(to_omit.length) ~ ")");

  GenotypeCombinator[][] ret;

  foreach(i; 0..to_omit.length) {
    if(!to_omit[i])
      ret ~= geno[i];
  }

  return ret;
}

unittest {
  writeln("Unit test " ~ __FILE__);

  // Phenotype reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data"));
  auto pheno_fn = to!string(buildPath(dir,"regression","test_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];

  // find individuals with missing phenotypes
  auto has_missing = is_any_phenotype_missing(pheno);

  // count them
  auto n_missing = count(has_missing, true);

  // check results
  assert(n_missing == 4);
  auto has_missing_shouldbe = new bool[](120);
  foreach(ref val; has_missing_shouldbe) val=false;
  has_missing_shouldbe[29] = true;
  has_missing_shouldbe[71] = true;
  has_missing_shouldbe[75] = true;
  has_missing_shouldbe[76] = true;
  assert(has_missing == has_missing_shouldbe);

  // omit individuals with missing phenotypes
  auto pheno_rev = omit_ind_from_phenotypes(pheno, has_missing);
  assert(pheno_rev.length == 116);
  auto has_missing_rev = is_any_phenotype_missing(pheno_rev);
  assert(count(has_missing_rev, true) == 0);
  
  // First read symbol information (the GenotypeCombinators)
  auto symbol_fn = to!string(buildPath(dir,"regression","test_symbol.qtab"));
  writeln("reading ",symbol_fn);
  auto f = File(symbol_fn,"r");
  auto symbols = read_genotype_symbol_qtab(f);
  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"regression","test_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto f1 = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(f1, symbols);
  auto individuals = ret[0];
  auto genotype_matrix = ret[1];

  assert(genotype_matrix.length == 120);
  // omit individuals with missing phenotypes
  auto genotype_matrix_rev = omit_ind_from_genotypes(genotype_matrix, has_missing);
  assert(genotype_matrix_rev.length == 116);
}